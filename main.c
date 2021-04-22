#include <stdio.h>
#include <string.h>
#include <stdlib.h>

#include "getopt.h"

#include "dna.h"
#include "filereader.h"
#include "bsalign.h"
#include "bspoa.h"
#include <stdlib.h>
#include <regex.h>

#include "kvec.h"
#include "kthread.h"
#include "kstring.h"
#include "bamlite.h"

/*https://pacbiofileformats.readthedocs.io/en/9.0/BAM.html*/
#define LCONTEXT_NULL (0)
#define LCONTEXT_ADAPTER_BEFORE (1)
#define LCONTEXT_ADAPTER_AFTER (2)
#define LCONTEXT_BARCODE_AFTER (4)
#define LCONTEXT_BARCODE_BEFORE (8)
#define LCONTEXT_FORWARD_PASS (16)
#define LCONTEXT_REVERSE_PASS (32)
#define LCONTEXT_ADAPTER_BEFORE_BAD (64)
#define LCONTEXT_ADAPTER_AFTER_BAD (128)


#define MAX_MOVIE_NAME_LEN (1024)
int8_t seq_comp_table[16] = { 0, 8, 4, 12, 2, 10, 6, 14, 1, 9, 5, 13, 3, 11, 7, 15 };
const char seq_nt16_str[] = "=ACMGRSVTWYHKDBN";

static char *get_read_seq(uint8_t *seq, int len, int breverse)
{
    int i, j;
    char *s = calloc(1, len + 1);
    if (s) {
        if (breverse) {
            for (i = 0, j = len - 1; i < len; i++, j--)
                s[i] = seq_nt16_str[seq_comp_table[bam1_seqi(seq, j)]];
        } else {
            for (i = 0; i < len; i++)
                s[i] = seq_nt16_str[bam1_seqi(seq, i)];
        }
    }
    return s;
}

void kt_for(int n_threads, void (*func)(void *, long, int), void *data, long n);
void kt_pipeline(int n_threads, void *(*func)(void *, int, void *),
                 void *shared_data, int n_steps);

typedef struct {
    int min_subread_len, max_subread_len;
    int min_fulllen_count;
    int verbose;
    int nthreads;
    size_t chunk_size;
    FILE *fp_out;
    bamFile fp_in;
    BSPOA **gg;     /*ccs comput ,per thread*/
    BSPOA **lgg;     /*ccs comput ,per thread*/
} pipeline_t;

#define SEQ_PARTIAL  (1)
#define SEQ_SHORT (2)
#define SEQ_LONG (4)

typedef struct {
    uint32_t len;
    uint8_t flags;
    uint8_t *buf;
} seq_t;

typedef kvec_t(seq_t) vec_seq_t;
typedef struct {
    vec_seq_t seqs;
    uint64_t hole;
    int npass;
    kstring_t ccsseq;
} zmw_t;

static inline void zmw_initialize(zmw_t *zmw_p, int nseqs)
{
    zmw_p->hole = 0;
    zmw_p->npass = 0;
    kv_init(zmw_p->seqs);
    kv_resize(seq_t, zmw_p->seqs, nseqs);
    zmw_p->ccsseq.m = zmw_p->ccsseq.l = 0, zmw_p->ccsseq.s = 0;
}

static inline void zmw_destroy(zmw_t *zmw_p)
{
    size_t i = 0;
    char *s = NULL;
    seq_t *seq;
    for (i = 0; i < kv_size(zmw_p->seqs); ++i) {
        seq = &kv_A(zmw_p->seqs, i);
        free(seq->buf);
    }

    kv_destroy(zmw_p->seqs);
    s = ks_release(&zmw_p->ccsseq);
    free(s);
}

typedef kvec_t(zmw_t) vec_zmw_t;
typedef struct {
    vec_zmw_t zmws;     /*zmw seqs*/
    BSPOA **gg;          /*ccs comput ,per thread*/
    BSPOA **lgg;          /*ccs comput ,per thread*/
} step_t;

static inline void step_initialize(step_t *step_p, int chunk_size)
{
    step_p->gg = step_p->lgg = NULL;

    kv_init(step_p->zmws);
    kv_resize(zmw_t, step_p->zmws, chunk_size);
}

static inline void step_destroy(step_t *s)
{
        size_t i = 0;
        for (i = 0; i < kv_size(s->zmws); ++i) {
            zmw_t *zmw_p = &kv_A(s->zmws, i);
            zmw_destroy(zmw_p);
        }
        kv_destroy(s->zmws);
        free(s);
}

static void ccs_for(void *_data, long i, int tid) // kt_for() callback
{
    step_t *step = (step_t *)_data;
    zmw_t *zmw_p = &kv_A(step->zmws, i);

    BSPOA *g = step->gg[tid];
    BSPOA *lg = step->lgg[tid];

    if (kv_size(zmw_p->seqs) < 3) {
        return;
    }

    beg_bspoa(g);

    /*push seq */
    size_t k;
    seq_t *str;
    for (k = 0; k < kv_size(zmw_p->seqs); ++k) {
        str = &kv_A(zmw_p->seqs, k);
        char *seq = get_read_seq(str->buf, str->len, k % 2);
        if (!str->flags)
            push_bspoa(g, seq, str->len);
        free(seq);
    }

    end_bspoa(g);
    remsa_bspoa(g, lg);

    /*save ccs seq*/
    char *res = malloc(g->cns->size + 1);
    for (k = 0; k < g->cns->size; ++k) {
        res[k] = bit_base_table[g->cns->buffer[k]];
    }
    res[k] = 0;
    zmw_p->ccsseq.l = g->cns->size;
    zmw_p->ccsseq.m = g->cns->size + 1;
    zmw_p->ccsseq.s = res;
}

static void ccs_for2(void *_data, long idx, int tid) // kt_for() callback
{
    step_t *step = (step_t *)_data;
    zmw_t *zmw_p = &kv_A(step->zmws, idx);

    BSPOA *g = step->gg[tid];
    BSPOA *lg = step->lgg[tid];

    if (kv_size(zmw_p->seqs) < 3) {
        return;
    }

    u4i nseq = 0;
    String **seq = calloc(kv_size(zmw_p->seqs) ,sizeof(String *));

    /*push seq */
    size_t kidx;
    for (kidx = 0; kidx < kv_size(zmw_p->seqs); ++kidx) {
        seq_t *str = &kv_A(zmw_p->seqs, kidx);
        char *read = get_read_seq(str->buf, str->len, kidx % 2);
        if (!str->flags) {
            seq[nseq ] = init_string(10240);
            append_string(seq[nseq], read, str->len);
	    nseq++;
        }
        free(read);
    }

    if(nseq<3) {
        u4i i=0;
    	for (i = 0; i < kv_size(zmw_p->seqs); i++) {
        	if(seq[i]) free_string(seq[i]);
    	}
    	free(seq);
        return;
    }
    //fprintf(stderr,"begin %lu\n",zmw_p->hole);	
    u4i splitlen = 1000, minlen = 200, addlen = 1000, kmersize = 8;
    float rowrate = 0.775, colrate = 0.8;
    zmw_p->ccsseq.l = zmw_p->ccsseq.m = 0, zmw_p->ccsseq.s = NULL;
    ks_resize(&zmw_p->ccsseq, 10000);

    //
    u4i done, i, beg[nseq], len[nseq];
    char *str;
    u4i j, k, off[nseq], off_nseq, ch_i, ch_nseq, k_off[nseq][kmersize], k_eq[nseq][kmersize], rowcnt, colcnt, rowtot, coltot, find;
    //
    done = 0;
    for (i = 0; i < nseq; i++) {
        beg[i] = 0;
        len[i] = splitlen;
        done = seq[i]->size - beg[i] < len[i] + minlen ? 1 : done;
    }
    for (i = 0; i < nseq; i++) {
        len[i] = done == 1 ? seq[i]->size - beg[i] : len[i];
    }
    while (done == 0) {
        beg_bspoa(g);
        for (i = 0; i < nseq; i++) {
            push_bspoa(g, seq[i]->string + beg[i], len[i]);
        }
        end_bspoa(g);
        remsa_bspoa(g, lg);
        //print_msa_mline_bspoa(g, stdout);
        str = str_msa_bspoa(g, g->strs);
        //
        j = g->msaidxs->size - 1;
        for (i = 0; i < nseq; i++) {
            off[i] = 0;
        }
        off_nseq = 0;
        //
        for (k = 0; k < kmersize; k++) {
            ch_nseq = 45;
            for (; j > 0 && ch_nseq == 45; j--) {
                ch_nseq = str[(nseq * (g->msaidxs->size + 1)) + j];
                off_nseq = ch_nseq != 45 ? off_nseq + 1 : off_nseq;
                for (i = 0; i < nseq; i++) {
                    ch_i = str[(i * (g->msaidxs->size + 1)) + j];
                    off[i] = ch_i != 45 ? off[i] + 1 : off[i];
                    k_off[i][k] = off[i];
                    k_eq[i][k] = ch_i == ch_nseq ? 1 : 0;
                }
            }
        }
        find = 1;
        for (i = 0; i < nseq && find == 1; i++) {
            rowcnt = 0;
            for (k = 0; k < kmersize; k++) {
                rowcnt = rowcnt + k_eq[i][k];
            }
            rowtot = kmersize;
            for (k = 0; k <= kmersize - 2; k++) {
                rowtot = k_off[i][k + 1] - k_off[i][k] >= 2 ? rowtot + k_off[i][k + 1] - k_off[i][k] - 1 : rowtot;
            }
            find = (float)rowcnt / (float)rowtot >= rowrate ? 1 : 0;
        }
        for (k = 0; k < kmersize && find == 1; k++) {
            colcnt = 0;
            for (i = 0; i < nseq; i++) {
                colcnt = colcnt + k_eq[i][k];
            }
            coltot = nseq;
            find = (float)colcnt / (float)coltot >= colrate ? 1 : 0;
        }
        //
        while (j > 0 && find == 0) {
            for (i = 0; i < nseq; i++) {
                for (k = 0; k <= kmersize - 2; k++) {
                    k_off[i][k] = k_off[i][k + 1];
                    k_eq[i][k] = k_eq[i][k + 1];
                }
            }
            ch_nseq = 45;
            for (; j > 0 && ch_nseq == 45; j--) {
                ch_nseq = str[(nseq * (g->msaidxs->size + 1)) + j];
                off_nseq = ch_nseq != 45 ? off_nseq + 1 : off_nseq;
                for (i = 0; i < nseq; i++) {
                    ch_i = str[(i * (g->msaidxs->size + 1)) + j];
                    off[i] = ch_i != 45 ? off[i] + 1 : off[i];
                    k_off[i][k] = off[i];
                    k_eq[i][k] = ch_i == ch_nseq ? 1 : 0;
                }
            }
            find = 1;
            for (i = 0; i < nseq && find == 1; i++) {
                rowcnt = 0;
                for (k = 0; k < kmersize; k++) {
                    rowcnt = rowcnt + k_eq[i][k];
                }
                rowtot = kmersize;
                for (k = 0; k <= kmersize - 2; k++) {
                    rowtot = k_off[i][k + 1] - k_off[i][k] >= 2 ? rowtot + k_off[i][k + 1] - k_off[i][k] - 1 : rowtot;
                }
                find = (float)rowcnt / (float)rowtot >= rowrate ? 1 : 0;
            }
            for (k = 0; k < kmersize && find == 1; k++) {
                colcnt = 0;
                for (i = 0; i < nseq; i++) {
                    colcnt = colcnt + k_eq[i][k];
                }
                coltot = nseq;
                find = (float)colcnt / (float)coltot >= colrate ? 1 : 0;
            }
        }
        //
        for (i = 0; i < nseq; i++) {
            beg[i] = find == 1 ? beg[i] + len[i] - off[i] : beg[i];
            len[i] = find == 0 ? len[i] + addlen : splitlen;
            done = seq[i]->size - beg[i] < len[i] + minlen ? 1 : done;
        }
        for (i = 0; i < nseq; i++) {
            len[i] = done == 1 ? seq[i]->size - beg[i] : len[i];
        }
        //
        for (i = 0; i < g->cns->size - off_nseq && find == 1; i++) {
            str[i] = bit_base_table[g->cns->buffer[i]];
        }
        str[i] = 0;
        str[0] = find == 0 ? 0 : str[0];
        kputs(str, &zmw_p->ccsseq);
    }
    beg_bspoa(g);
    for (i = 0; i < nseq; i++) {
        push_bspoa(g, seq[i]->string + beg[i], len[i]);
    }
    end_bspoa(g);
    remsa_bspoa(g, lg);
    //print_msa_mline_bspoa(g, stdout);
    str = str_msa_bspoa(g, g->strs);
    for (i = 0; i < g->cns->size; i++) {
        str[i] = bit_base_table[g->cns->buffer[i]];
    }
    str[i] = 0;
    kputs(str, &zmw_p->ccsseq);

    for (i = 0; i < kv_size(zmw_p->seqs); i++) {
        if(seq[i]) free_string(seq[i]);
    }
    free(seq);    
    //fprintf(stderr,"end %lu\n",zmw_p->hole);
}


static void *worker_pipeline(void *shared, int step, void *in) // kt_pipeline() callback
{
    pipeline_t *p = (pipeline_t *)shared;
    if (step == 0) { // step 0: read zmw into the buffer
        step_t *s;
        s = calloc(1, sizeof(step_t));
        step_initialize(s, p->chunk_size);
        s->gg = p->gg;
        s->lgg = p->lgg;

        bam1_t *b = NULL;
        int rc = 0;
        zmw_t zmw;
        int iseof;
        int pre_zm = 0;
        while (1) {
            if (!b)
                b = bam_init1();
            if (b == NULL) {
                fprintf(stderr, "Error: Failed to init BAM block!\n");
                rc = -1;
                break;
            }

            if ((rc = bam_read1(p->fp_in, b)) < -1) {
                fprintf(stderr, "Error: truncated file \n");
                break;
            } else if (rc == -1) {
                iseof = 1;
                if (!pre_zm)
                    break;
            } else {
                iseof = 0;
            }

            int seq_len = 0, cx = 0, zm = 0;
            if (!iseof) {
                seq_len = b->core.l_qseq;
                uint8_t *tag = NULL;
                if ((tag = bam_aux_get(b, "cx")))
                    cx = bam_aux2i(tag);
                if ((tag = bam_aux_get(b, "zm")))
                    zm = bam_aux2i(tag);
            }

            if (iseof || !pre_zm || pre_zm != zm) {
                /*handle the old*/
                if (pre_zm) {
                    if (zmw.npass >= p->min_fulllen_count && zmw.npass && kv_size(zmw.seqs) >= 3)    {
                        kv_push(zmw_t, s->zmws, zmw);
                    }
                    if (iseof || kv_size(s->zmws) >= p->chunk_size)
                        break;
                }

                /*restart the new*/
                if (!iseof) {
                    pre_zm = zm;

                    zmw_initialize(&zmw, 4);/*4 average*/
                    zmw.npass = 0;
                    zmw.hole = pre_zm;
                }
            }

            if (!iseof) {
                /*update the old*/
                size_t buflen = (b->core.l_qseq + 1) >> 1;
                seq_t seq = {b->core.l_qseq, 0, malloc(buflen)};
                if (seq.buf) memcpy(seq.buf, bam1_seq(b), buflen);

                if (cx & LCONTEXT_ADAPTER_AFTER && cx & LCONTEXT_ADAPTER_BEFORE) {
                    zmw.npass++;
                } else {
                    seq.flags |= SEQ_PARTIAL;
                }
                if (seq_len < p->min_subread_len) {
                    seq.flags |= SEQ_SHORT;
                } else if (seq_len > p->max_subread_len) {
                    seq.flags |= SEQ_LONG;
                }

                kv_push(seq_t, zmw.seqs, seq);
                if (p->verbose) {
                    size_t i, len = b->core.l_qseq;
                    char *tmp = malloc(len + 1);
                    for (i = 0; i < len; i++)
                        tmp[i] = seq_nt16_str[bam1_seqi(bam1_seq(b), i)];
                    tmp[len] = '\0';
                    fprintf(stderr, ">%s\n%s\n", bam1_qname(b), tmp);
                    free(tmp);
                }
            }

        }

        if (b) {
            bam_destroy1(b);
            b = NULL;
        }
        if (kv_size(s->zmws))
            return s;
        else
            step_destroy(s);
    } else if (step == 1) { // step 1: ccs
        step_t *s = (step_t *)in;
        kt_for(p->nthreads, ccs_for, in, kv_size(s->zmws));
        return in;
    } else if (step == 2) { // step 2: write ccs seq to output
        step_t *s = (step_t *)in;
        size_t i = 0;
        for (i = 0; i < kv_size(s->zmws); ++i) {
            zmw_t *zmw_p = &kv_A(s->zmws, i);
            fprintf(p->fp_out, ">%ld\n%s\n", zmw_p->hole, ks_str(&zmw_p->ccsseq));
            zmw_destroy(zmw_p);
        }
        kv_destroy(s->zmws);
        free(s);
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
            "-j     <int>   Number of threads to use, 0 means autodetection. [2] \n"
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
           min_fulllen_count = 3, nthreads = 1;
    while ((c = getopt(argc, argv, "hm:M:c:j:v")) != -1) {
        switch (c) {
            case 'm':
                min_subread_len = atoi(optarg);
                break;
            case 'M':
                max_subread_len = atoi(optarg);
                break;
            case 'c':
                min_fulllen_count = atoi(optarg);
                if (min_fulllen_count < 1) {
                    fprintf(stderr, "Error! min fulllen count=[%d] (>=1) !\n", min_fulllen_count);
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

    char *infile;
    FILE *fp_out = NULL;

    if (argc - optind == 1) {
        infile = argv[optind];
        fp_out = open_file_for_write("-", NULL, 1);;
    } else if (argc - optind == 2) {
        infile = argv[optind];
        fp_out = open_file_for_write(argv[optind + 1], NULL, 1);
    } else {
        return usage();
    }

    /*read bam file*/
    bamFile fp_in = bam_open(infile, "rb");
    if (!fp_in) {
        fprintf(stderr, "Error: Failed to open '%s'!\n", argv[1]);
        return 1;
    }

    bam_header_t  *hdr = bam_header_read(fp_in);
    if (!hdr) {
        fprintf(stderr, "Error: Failed to read BAM head!\n");
        bam_close(fp_in);
        return 1;
    }
    bam_header_destroy(hdr);

    pipeline_t pl;
    pl.max_subread_len = max_subread_len;
    pl.min_subread_len = min_subread_len;
    pl.min_fulllen_count = min_fulllen_count;
    pl.nthreads = nthreads;
    pl.chunk_size = 10240;
    pl.fp_in = fp_in;
    pl.fp_out = fp_out;
    pl.verbose = verbose;
    int i = 0;
    pl.gg = calloc(nthreads, sizeof(BSPOA *));
    pl.lgg = calloc(nthreads, sizeof(BSPOA *));
    for (i = 0; i < nthreads; ++i) {
        BSPOAPar par = DEFAULT_BSPOA_PAR;
        par.M = 2;
        par.X = -6;
        par.O = -3;
        par.E = -2;
        par.Q = 0;
        par.P = 0;
        par.bwtrigger = 10;
        pl.gg[i] = init_bspoa(par);

        BSPOAPar rpar = DEFAULT_BSPOA_PAR;
        rpar.M = 1;
        rpar.X = -2;
        rpar.O = 0;
        rpar.E = -1;
        rpar.Q = 0;
        rpar.P = 0;
        pl.lgg[i] = init_bspoa(rpar);
    }
    kt_pipeline(nthreads, worker_pipeline, &pl, 3);
    for (i = 0; i < nthreads; ++i) {
        free_bspoa(pl.gg[i]);
        free_bspoa(pl.lgg[i]);
    }
    free(pl.gg);free(pl.lgg);

    fclose(fp_out);
    bam_close(fp_in);
    return 0;
}

