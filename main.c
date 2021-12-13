#include <stdio.h>
#include <string.h>
#include <stdlib.h>

#include "getopt.h"

#include "dna.h"

#include "bsalign.h"
#include "bspoa.h"
#include <stdlib.h>
#include <regex.h>

#include "kvec.h"
#include "kthread.h"
#include "kstring.h"
#include "bamlite.h"
#include "seqio.h"
#include "khash.h"

KHASH_SET_INIT_STR(strset)

SEQIO_INIT(gzFile, gzread)

void kt_for(int n_threads, void (*func)(void *, long, int), void *data, long n);
void kt_pipeline(int n_threads, void *(*func)(void *, int, void *),
                 void *shared_data, int n_steps);

typedef struct __pipeline_t
{
    size_t min_subread_len, max_subread_len;
    int min_fulllen_count;
    int verbose;
    int split_subread;
    int nthreads;
    size_t chunk_size;
    FILE *fp_out;
    kseqs_zmw_t zmw_seqs_in;
    khash_t(strset) *hole_set;
    BSPOA **gg;     /*ccs comput ,per thread*/
} pipeline_t;

typedef struct
{
    kstring_t movie_name, hole, seqs;
    kvec_t(int) lens;
    kstring_t ccsseq;
} zmw_t;

static inline void zmw_initialize(zmw_t *zmw_p, int nseqs)
{
    kv_init(zmw_p->lens);
    kv_resize(int, zmw_p->lens, nseqs);
    zmw_p->ccsseq.m = zmw_p->ccsseq.l = 0, zmw_p->ccsseq.s = 0;
    zmw_p->hole.m = zmw_p->hole.l = 0, zmw_p->hole.s = 0;
    zmw_p->movie_name.m = zmw_p->movie_name.l = 0, zmw_p->movie_name.s = 0;
    zmw_p->seqs.m = zmw_p->seqs.l = 0, zmw_p->seqs.s = 0;
}

static inline zmw_t *zmw_init(int nseqs)
{
    zmw_t *zmw_p = (zmw_t *)calloc(1, sizeof(zmw_t));
    zmw_initialize(zmw_p, nseqs);
    return zmw_p;
}

static inline void zmw_destroy(zmw_t *zmw_p)
{
    if(!zmw_p)
        return;
    char *s = NULL;
    kv_destroy(zmw_p->lens);
    s = ks_release(&zmw_p->ccsseq);
    free(s);
    s = ks_release(&zmw_p->seqs);
    free(s);
    s = ks_release(&zmw_p->hole);
    free(s);
    s = ks_release(&zmw_p->movie_name);
    free(s);
}

typedef kvec_t(zmw_t) vec_zmw_t;
typedef struct
{
    vec_zmw_t zmws;     /*zmw seqs*/
    BSPOA **gg;          /*ccs comput ,per thread*/
    int verbose;
} step_t;

static inline step_t   *step_init(int chunk_size)
{
    step_t *step_p = calloc(1, sizeof(step_t));
    step_p->gg = NULL;
    kv_init(step_p->zmws);
    kv_resize(zmw_t, step_p->zmws, chunk_size);
    return step_p;
}

static inline void step_destroy(step_t *s)
{
    size_t i = 0;
    for(i = 0; i < kv_size(s->zmws); ++i)
    {
        zmw_t *zmw_p = &kv_A(s->zmws, i);
        zmw_destroy(zmw_p);
    }
    kv_destroy(s->zmws);
    if(s)
    {
        free(s);
        s = 0;
    }
}

static void ccs_for(void *_data, long idx, int tid) // kt_for() callback
{
    size_t l;
    step_t *step = (step_t *)_data;
    zmw_t *zmw_p = &kv_A(step->zmws, idx);
    BSPOA *g = step->gg[tid];
    if(kv_size(zmw_p->lens) < 3)
    {
        return;
    }

    int core_len = ks_len(&zmw_p->seqs) - kv_A(zmw_p->lens, 0) - kv_A(zmw_p->lens, kv_size(zmw_p->lens) - 1);
    int averge_core_len = core_len / (kv_size(zmw_p->lens) - 2);
    typedef struct
    {
        int offs;
        int len;
        char reverse;
    } segment_t;
    kvec_t(segment_t) segments;
    kv_init(segments);
    kv_resize(segment_t, segments, kv_size(zmw_p->lens));
    segment_t seg;
    int core_offset = 0;
    int zmw_len;
    for(l = 1; l < kv_size(zmw_p->lens) - 1; ++l)
    {
        zmw_len = kv_A(zmw_p->lens, l);
        if(zmw_len < averge_core_len / 2)
        {
            core_offset += zmw_len;
            continue;
        }
        else
        {
            seg.offs = core_offset + kv_A(zmw_p->lens, 0);
            seg.len = (zmw_len > averge_core_len * 3 / 2) ? averge_core_len : zmw_len;
            seg.reverse = (core_offset + seg.len / 2) / averge_core_len % 2;
            kv_push(segment_t, segments, seg);
        }
        core_offset += zmw_len;
    }

    //the first record
    zmw_len = kv_A(zmw_p->lens, 0);
    if(zmw_len > averge_core_len / 2)
    {
        seg.reverse = 1;
        seg.len = zmw_len > averge_core_len ? averge_core_len : zmw_len;
        seg.offs = 0;
        kv_push(segment_t, segments, seg);
    }


    if(step->verbose > 1)
        fprintf(stderr, "poa begin %s\n", ks_str(&zmw_p->hole));

    beg_bspoa(g);
    /*push seq */
    segment_t *seg_p = NULL;
    for(l = 0; l < kv_size(segments); ++l)
    {
        seg_p = &kv_A(segments, l);
        if(seg_p->reverse)
        {
            seq_reverse_comp(seg_p->len, (unsigned char *)ks_str(&zmw_p->seqs) + seg_p->offs);
        }

        if(step->verbose)
            fprintf(stderr, ">%s_%lu/%lu strand=%d len=%d \n%.*s\n",
                    ks_str(&zmw_p->hole), l, kv_size(segments), seg_p->reverse, seg_p->len, seg_p->len, ks_str(&zmw_p->seqs) + seg_p->offs);
        push_bspoa(g, ks_str(&zmw_p->seqs) + seg_p->offs, seg_p->len);
    }
    end_bspoa(g);

    kv_destroy(segments);

    /*save ccs seq*/
    char *res = malloc(g->cns->size + 1);
    for(l = 0; l < g->cns->size; ++l)
    {
        res[l] = bit_base_table[g->cns->buffer[l]];
    }
    res[l] = 0;
    zmw_p->ccsseq.l = g->cns->size;
    zmw_p->ccsseq.m = g->cns->size + 1;
    zmw_p->ccsseq.s = res;

    if(step->verbose > 1)
        fprintf(stderr, "poa end %s\n", ks_str(&zmw_p->hole));
}

static void ccs_for2(void *_data, long idx, int tid) // kt_for() callback
{
    size_t l;
    step_t *step = (step_t *)_data;
    zmw_t *zmw_p = &kv_A(step->zmws, idx);
    BSPOA *g = step->gg[tid];
    if(kv_size(zmw_p->lens) < 3)
    {
        return;
    }

    int core_len = ks_len(&zmw_p->seqs) - kv_A(zmw_p->lens, 0) - kv_A(zmw_p->lens, kv_size(zmw_p->lens) - 1);
    int averge_core_len = core_len / (kv_size(zmw_p->lens) - 2);
    typedef struct
    {
        int offs;
        int len;
        char reverse;
    } segment_t;
    kvec_t(segment_t) segments;
    kv_init(segments);
    kv_resize(segment_t, segments, kv_size(zmw_p->lens));
    segment_t seg;
    int core_offset = 0;
    int zmw_len;
    for(l = 1; l < kv_size(zmw_p->lens); ++l)
    {
        zmw_len = kv_A(zmw_p->lens, l);
        if(zmw_len < averge_core_len / 2)
        {
            core_offset += zmw_len;
            continue;
        }
        else
        {
            seg.offs = core_offset + kv_A(zmw_p->lens, 0);
            seg.len = (zmw_len > averge_core_len * 3 / 2) ? averge_core_len : zmw_len;
            seg.reverse = (core_offset + seg.len / 2) / averge_core_len % 2;
            kv_push(segment_t, segments, seg);
        }
        core_offset += zmw_len;
    }
    //the first record
    zmw_len = kv_A(zmw_p->lens, 0);
    if(zmw_len > averge_core_len / 2)
    {
        seg.reverse = 1;
        seg.len = zmw_len > averge_core_len ? averge_core_len : zmw_len;
        seg.offs = 0;
        kv_push(segment_t, segments, seg);
    }

    u4i nseq = kv_size(segments);
    String **seq = calloc(kv_size(segments), sizeof(String *));
    segment_t *seg_p = NULL;
    for(l = 0; l < kv_size(segments); ++l)
    {
        seg_p = &kv_A(segments, l);
        if(seg_p->reverse)
        {
            seq_reverse_comp(seg_p->len, (unsigned char *)ks_str(&zmw_p->seqs) + seg_p->offs);
        }
        seq[l] = init_string(seg_p->len);
        append_string(seq[l], ks_str(&zmw_p->seqs) + seg_p->offs, seg_p->len);
    }

    zmw_p->ccsseq.l = zmw_p->ccsseq.m = 0, zmw_p->ccsseq.s = NULL;
    ks_resize(&zmw_p->ccsseq, 10000);



    //
    //fprintf(stderr, "split beg\n");
    //
    char *str;
    int initlen, addlen, minlen, strcap, strsize, pushbeg[nseq], pushlen[nseq], flag;
    double rowrate, colrate;
    u1i *col;
    u4i i, j, k, window, minwin, nogwin, mrow;
    initlen = 1000;
    addlen = 1000;
    minlen = 200;
    strcap = 0;
    strsize = 0;
    flag = 1;
    rowrate = 0.8;
    colrate = 0.8;
    window = 10;
    minwin = 5;
    u4i rowcnt[nseq], colcnt[window], rowmin, colmin, rowidx, colidx, nogcol[window];
    for(i = 0; i < nseq; i++)
    {
        pushbeg[i] = 0;
        pushlen[i] = initlen;
        if(strcap < seq[i]->size)
        {
            strcap = seq[i]->size;
        }
        if(pushbeg[i] + pushlen[i] + minlen > seq[i]->size)
        {
            flag = 0;
        }
    }
    if(flag == 0)
    {
        for(i = 0; i < nseq; i++)
        {
            pushlen[i] = seq[i]->size - pushbeg[i];
        }
    }
    strcap = strcap * 2;
    str = (char *)calloc(strcap, sizeof(char));
    while(flag == 1)
    {
        beg_bspoa(g);
        for(i = 0; i < nseq; i++)
        {
            /*fprintf(stderr, "seq%d %d %d ", i, pushbeg[i], pushlen[i]);
            for(j = pushbeg[i]; j < pushbeg[i] + pushlen[i]; j++)
            {
                fprintf(stderr, "%c", seq[i]->string[j]);
            }
            fprintf(stderr, "\n");*/
            push_bspoa(g, seq[i]->string + pushbeg[i], pushlen[i]);
        }
        end_bspoa(g);
        tidy_msa_bspoa(g);
        mrow = nseq + 1 + 3;
        for(i = g->msaidxs->size - window; i >= 1; i--)
        {
            col = g->msacols->buffer + g->msaidxs->buffer[i] * mrow;
            if(col[nseq + 1] >= 4)
            {
                continue;
            }
            nogwin = 1;
            for(j = i + 1; j < i + window; j++)
            {
                col = g->msacols->buffer + g->msaidxs->buffer[j] * mrow;
                if(col[nseq + 1] < 4)
                {
                    nogwin++;
                }
            }
            if(nogwin < minwin)
            {
                continue;
            }
            for(rowidx = 0; rowidx < nseq; rowidx++)
            {
                rowcnt[rowidx] = 0;
            }
            for(colidx = 0; colidx < window; colidx++)
            {
                colcnt[colidx] = 0;
            }
            for(j = i; j < i + window; j++)
            {
                colidx = j - i;
                col = g->msacols->buffer + g->msaidxs->buffer[j] * mrow;
                if(col[nseq + 1] < 4)
                {
                    for(k = 1; k <= nseq; k++)
                    {
                        rowidx = k - 1;
                        if(col[k] == col[nseq + 1])
                        {
                            rowcnt[rowidx]++;
                            colcnt[colidx]++;
                        }
                    }
                    nogcol[colidx] = 1;
                }
                else
                {
                    nogcol[colidx] = 0;
                }
            }
            /*for(j = i; j < i + window; j++)
                        {
                                col = g->msacols->buffer + g->msaidxs->buffer[j] * mrow;
                                for(k = 1; k <= nseq; k++)
                                {
                                        fprintf(stderr, "%d", col[k]);
                                }
                                fprintf(stderr, "%d\n", col[nseq + 1]);
                        }*/
            rowmin = window;
            for(rowidx = 0; rowidx < nseq; rowidx++)
            {
                //fprintf(stderr, "rowidx=%d, rowcnt[%d]=%d\n", rowidx, rowidx, rowcnt[rowidx]);
                if(rowmin > rowcnt[rowidx])
                {
                    rowmin = rowcnt[rowidx];
                }
            }
            colmin = nseq;
            for(colidx = 0; colidx < window; colidx++)
            {
                //fprintf(stderr, "colidx=%d, colcnt[%d]=%d\n", colidx, colidx, colcnt[colidx]);
                if(nogcol[colidx] == 1 && colmin > colcnt[colidx])
                {
                    colmin = colcnt[colidx];
                }
            }
            //fprintf(stderr, "rowmin=%d/nogwin=%d, colmin=%d/nseq=%d\n", rowmin, nogwin, colmin, nseq);
            if((float)rowmin / (float)nogwin >= rowrate && (float)colmin / (float)nseq >= colrate)
            {
                //fprintf(stderr, "yes\n");
                break;
            }
        }
        if(i >= 1)
        {
            for(j = 0; j < i; j++)
            {
                col = g->msacols->buffer + g->msaidxs->buffer[j] * mrow;
                for(k = 1; k <= nseq; k++)
                {
                    if(col[k] < 4)
                    {
                        pushbeg[k - 1]++;
                    }
                }
                if(col[nseq + 1] < 4)
                {
                    str[strsize] = bit_base_table[col[nseq + 1]];
                    strsize++;
                }
            }
            for(j = 0; j < nseq; j++)
            {
                pushlen[j] = initlen;
                if(pushbeg[j] + pushlen[j] + minlen > seq[j]->size)
                {
                    flag = 0;
                }
            }
            if(flag == 0)
            {
                for(j = 0; j < nseq; j++)
                {
                    pushlen[j] = seq[j]->size - pushbeg[j];
                }
            }
        }
        else
        {
            for(j = 0; j < nseq; j++)
            {
                pushlen[j] = pushlen[j] + addlen;
                if(pushbeg[j] + pushlen[j] + minlen > seq[j]->size)
                {
                    flag = 0;
                }
            }
            if(flag == 0)
            {
                for(j = 0; j < nseq; j++)
                {
                    pushlen[j] = seq[j]->size - pushbeg[j];
                }
            }
        }
    }
    beg_bspoa(g);
    for(i = 0; i < nseq; i++)
    {
        /*fprintf(stderr, "seq%d %d %d ", i, pushbeg[i], pushlen[i]);
        for(j = pushbeg[i]; j < pushbeg[i] + pushlen[i]; j++)
        {
            fprintf(stderr, "%c", seq[i]->string[j]);
        }
        fprintf(stderr, "\n");*/
        push_bspoa(g, seq[i]->string + pushbeg[i], pushlen[i]);
    }
    end_bspoa(g);
    tidy_msa_bspoa(g);
    mrow = nseq + 1 + 3;
    for(i = 0; i < g->msaidxs->size; i++)
    {
        col = g->msacols->buffer + g->msaidxs->buffer[i] * mrow;
        if(col[nseq + 1] < 4)
        {
            str[strsize] = bit_base_table[col[nseq + 1]];
            strsize++;
        }
    }
    str[strsize] = 0;
    //
    //fprintf(stderr, "split end\n");
    //
    kputs(str, &zmw_p->ccsseq);
    free(str);
    //destroy
    for(i = 0; i < nseq; i++)
    {
        if(seq[i])
            free_string(seq[i]);
    }
    free(seq);
    kv_destroy(segments);
    //fprintf(stderr,"end %lu\n",zmw_p->hole);
}

static void *worker_pipeline(void *shared, int step, void *in) // kt_pipeline() callback
{
    pipeline_t *p = (pipeline_t *)shared;
    if(step == 0)    // step 0: read zmw into the buffer
    {
        step_t *s = step_init(p->chunk_size);
        s->verbose = p->verbose;
        s->gg = p->gg;
        kseqs_zmw_t *zmw_seqs = &p->zmw_seqs_in;
        int l;
        while((l = kseq_zmw_read(zmw_seqs)) >= 0)
        {
            if(l < p->min_fulllen_count + 2)
                continue;

            size_t averge_len = (ks_len(&zmw_seqs->seqs) / kv_size(zmw_seqs->lens));
            if(averge_len > p->max_subread_len || averge_len < p->min_subread_len)
            {
                continue;
            }

            if(p->hole_set)
            {
                khash_t(strset)* hole_set = p->hole_set;
                if(kh_get(strset, hole_set, ks_str(&zmw_seqs->hole)) != kh_end(hole_set))
                {
                    continue;
                }
            }

            zmw_t zmw;
            zmw_initialize(&zmw, 8);
            kputsn(ks_str(&zmw_seqs->hole), ks_len(&zmw_seqs->hole), &zmw.hole);
            kputsn(ks_str(&zmw_seqs->movie_name), ks_len(&zmw_seqs->movie_name), &zmw.movie_name);
            kputsn(ks_str(&zmw_seqs->seqs), ks_len(&zmw_seqs->seqs), &zmw.seqs);
            kv_copy(int, zmw.lens, zmw_seqs->lens);
            kv_push(zmw_t, s->zmws, zmw);
            if(kv_size(s->zmws) >= p->chunk_size)
            {
                if(p->chunk_size < 16384)
                {
                    p->chunk_size *= 4;
                }
                break;
            }
        }

        if(kv_size(s->zmws))
            return s;
        else
            step_destroy(s);
    }
    else if(step == 1)      // step 1: ccs
    {
        step_t *s = (step_t *)in;
        if(s)
        {
            if(p->split_subread)
                kt_for(p->nthreads, ccs_for2, in, kv_size(s->zmws));
            else
                kt_for(p->nthreads, ccs_for, in, kv_size(s->zmws));
        }
        return in;
    }
    else if(step == 2)      // step 2: write ccs seq to output
    {
        step_t *s = (step_t *)in;
        if(s)
        {
            size_t i = 0;
            for(i = 0; i < kv_size(s->zmws); ++i)
            {
                zmw_t *zmw_p = &kv_A(s->zmws, i);
                if(ks_len(&zmw_p->ccsseq))
                    fprintf(p->fp_out, ">%s/%s/ccs\n%s\n", ks_str(&zmw_p->movie_name), ks_str(&zmw_p->hole), ks_str(&zmw_p->ccsseq));
            }
            step_destroy(s);
        }
    }
    return 0;
}


static int usage()
{
    fprintf(stdout,
            "Program: ccsx\n"
            "Version: 1.0.0\n"
            "Usage  : ccsx  [options] <INPUT> <OUTPUT>\n"
            "Generate circular consensus sequences (ccs) from subreads.\n"
            "\n"
            "Options:\n"
            "-h             Output this help \n"
            "-v             debug \n"
            "-m     <int>   Minimum length of subreads to use for generating CCS. [10] \n"
            "-M     <int>   Maximum length of subreads to use for generating CCS. [50000] \n"
            "-c     <int>   Minimum number of subreads required to generate CCS. [3] \n"
            "-A             For fasta/fastq input,gzip allowed  \n"
            "-P             primitive bsalign,subread shred by default \n"
            "-X		<str>   Exclude ZMWs from output file,a comma-separated list of ID \n"
            "-j     <int>   Number of threads to use. [2] \n"
            "\n"
            "Arguments:\n"
            "input          Input file.\n"
            "output         Output file.\n"
            "\n"
           );
    return 1;
}


int main(int argc, char **const argv)
{
    int c, verbose = 0, min_subread_len = 10, max_subread_len = 50000,
           min_fulllen_count = 3, nthreads = 1, isbam = 1, split_subread = 1;

    khash_t(strset)* hole_set = NULL;
    kstring_t strholes = {0};
    while((c = getopt(argc, argv, "hm:M:c:j:X:PAv")) != -1)
    {
        switch(c)
        {
            case 'm':
                min_subread_len = atoi(optarg);
                break;
            case 'M':
                max_subread_len = atoi(optarg);
                break;
            case 'P':
                split_subread = 0;
                break;
            case 'A':
                isbam = 0;/*fasta/q*/
                break;
            case 'X':
            {
                kputs(optarg, &strholes);
                int n, i, *fields;
                fields = ksplit(&strholes, ',', &n);
                hole_set = kh_init(strset); // allocate a hash table
                for(i = 0; i < n; ++i)
                {
                    int ret;
                    kh_put(strset, hole_set, ks_str(&strholes) + fields[i], &ret);
                }
                free(fields);
            }
            break;
            case 'c':
                min_fulllen_count = atoi(optarg);
                if(min_fulllen_count < 3)
                {
                    fprintf(stderr, "Error! min fulllen count=[%d] (>=3) !\n", min_fulllen_count);
                    return -1;
                }
                break;
            case 'v':
                verbose++;
                break;
            case 'j':
                nthreads = atoi(optarg);
                break;
            default:
                return usage();
        }
    }

    FILE *fp_out = NULL;
    gzFile fp_in;
    if(argc - optind == 0)
    {
        fp_in = gzdopen(fileno(stdin), "rb") ;
        fp_out = stdout;
    }
    else if(argc - optind == 1)
    {
        fp_in = (strcmp(argv[optind], "-") == 0) ? gzdopen(fileno(stdin), "rb") : gzopen(argv[optind], "rb");
        fp_out = stdout;
    }
    else if(argc - optind == 2)
    {
        fp_in = (strcmp(argv[optind], "-") == 0) ? gzdopen(fileno(stdin), "rb") : gzopen(argv[optind], "rb");
        fp_out = (strcmp(argv[optind + 1], "-") == 0) ? stdout : fopen(argv[optind + 1], "w+");
    }
    else
    {
        return usage();
    }

    /*read infile*/
    if(!fp_in)
    {
        fprintf(stderr, "Error: Failed to open infile!\n");
        return 1;
    }

    if(!fp_out)
    {
        fprintf(stderr, "Cannot open file for write!\n");
        return 1;
    }

    pipeline_t pl;
    pl.max_subread_len = max_subread_len;
    pl.min_subread_len = min_subread_len;
    pl.min_fulllen_count = min_fulllen_count;
    pl.nthreads = nthreads;
    pl.chunk_size = 1024;
    pl.split_subread = split_subread;
    pl.fp_out = fp_out;
    pl.verbose = verbose;
    pl.hole_set = hole_set;
    int i = 0;
    pl.gg = calloc(nthreads, sizeof(BSPOA *));
    for(i = 0; i < nthreads; ++i)
    {
        BSPOAPar par = DEFAULT_BSPOA_PAR;
        par.M = 2;
        par.X = -6;
        par.O = -3;
        par.E = -2;
        par.Q = 0;
        par.P = 0;
        pl.gg[i] = init_bspoa(par);
    }

    kseqs_zmw_initialize(&pl.zmw_seqs_in, fp_in, isbam);

    kt_pipeline(2, worker_pipeline, &pl, 3);
    for(i = 0; i < nthreads; ++i)
    {
        free_bspoa(pl.gg[i]);
    }
    free(pl.gg);
    zmw_seqs_release(&pl.zmw_seqs_in);
    if(pl.hole_set)
    {
        kh_destroy(strset, hole_set);
    }
    free(ks_release(&strholes));
    fclose(fp_out);
    gzclose(fp_in);

    return 0;
}
