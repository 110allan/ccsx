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

typedef struct __pipeline_t {
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

typedef struct {
    kstring_t movie_name, hole;
    kstring_t seqs;
    kvec_t(u4i) lens;
    kvec_t(u4i) offs;
    kstring_t ccsseq;
} zmw_t;

static inline void zmw_initialize(zmw_t *zmw_p, int nseqs)
{
    kv_init(zmw_p->lens);
    kv_resize(u4i, zmw_p->lens, nseqs);
    kv_init(zmw_p->offs);
    kv_resize(u4i, zmw_p->offs, nseqs);
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
    if (!zmw_p)
        return;
    char *s = NULL;
    kv_destroy(zmw_p->lens);
    kv_destroy(zmw_p->offs);
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
typedef struct {
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
    size_t i = 0,i_size=kv_size(s->zmws);
    for (i = 0; i <i_size ; ++i) {
        zmw_t *zmw_p = &kv_A(s->zmws, i);
        zmw_destroy(zmw_p);
    }
    kv_destroy(s->zmws);
    if (s) {
        free(s);
        s = 0;
    }
}

typedef kvec_t(int) vec_int;
typedef struct {
    vec_int group_ids;
    size_t sum_len;
} group_result_t;

typedef kvec_t(group_result_t) vec_group_result_t;

static inline int len_in_group(group_result_t *g, u4i len, int tolerance_pct)
{
    size_t tmp = len * kv_size(g->group_ids);
    size_t diff_len = tmp > g->sum_len ? (tmp - g->sum_len) : (g->sum_len - tmp);
    return diff_len * 100 < tolerance_pct * g->sum_len;
}

static inline int group_in_group(group_result_t *g, group_result_t *g_qry, int tolerance_pct)
{
    size_t mean_group_g = g->sum_len * kv_size(g_qry->group_ids);
    size_t mean_group_g_qry = g_qry->sum_len * kv_size(g->group_ids);
    size_t diff_len = (mean_group_g > mean_group_g_qry) ? (mean_group_g - mean_group_g_qry) : (mean_group_g_qry - mean_group_g);
    return diff_len * 100 < mean_group_g * tolerance_pct;
}

static vec_group_result_t init_group_lens(u4i *array, int n, int tolerance_pct)
{
    group_result_t groups[n];
    int i, j;
    for (j = 0; j < n; ++j) {
        kv_init(groups[j].group_ids);
        groups[j].sum_len = 0;
    }
    for (i = 0; i < n; ++i) {
        //check groups
        for (j = 0; j < i; ++j) {
            if (!groups[j].sum_len) {
                continue;
            }
            if (len_in_group(groups + j, array[i], tolerance_pct)) {
                //add group j with a[i],idx=i
                kv_push(int, groups[j].group_ids, i);
                groups[j].sum_len += array[i];
                break;
            }
        }
        if (j < i) { //find the right group
            continue;
        }
        // cannot find group
        //create group j with a[i],idx=i
        kv_push(int, groups[j].group_ids, i);
        groups[j].sum_len = array[i];
    }

    //merge groups
    int flag = 1;
    while (flag) {
        flag = 0;
        for (j = 0; j < n; ++j) {
            if (kv_size(groups[j].group_ids) == 0) {
                continue;
            }
            //check groups
            for (int k = 0; k < j; ++k) {
                if (kv_size(groups[k].group_ids) && group_in_group(groups + k, groups + j, tolerance_pct)) {
                    //extend group k(mean=mean_group_k) with group j(mean=mean_group_j)
                    size_t extend_size = kv_size(groups[k].group_ids) + kv_size(groups[j].group_ids);
                    kv_resize(int, groups[k].group_ids, extend_size);
                    memcpy(groups[k].group_ids.a + groups[k].group_ids.n,
                           groups[j].group_ids.a, sizeof(int) * groups[j].group_ids.n);
                    kv_size(groups[k].group_ids) = extend_size;
                    groups[k].sum_len += groups[j].sum_len;
                    kv_destroy(groups[j].group_ids);
                    kv_init(groups[j].group_ids);
                    groups[j].sum_len = 0;
                    flag = 1;
                    break;
                }
            }
        }
    }
    vec_group_result_t ret_groups;
    kv_init(ret_groups);
    for (j = 0; j < n; ++j) {
        if (kv_size(groups[j].group_ids) == 0) {
            kv_destroy(groups[j].group_ids);
            continue;
        }
        kv_push(group_result_t, ret_groups, groups[j]);
    }

    //sort ret_groups
    if (kv_size(ret_groups) > 1) {
        bubble_sort_array(ret_groups.a, kv_size(ret_groups), group_result_t,
                          (kv_size(b.group_ids) > kv_size(a.group_ids)));
    }
    return ret_groups;
}

static inline void destroy_group_lens(vec_group_result_t *groups)
{
    for (size_t i = 0,i_size=kv_size(*groups); i < i_size; ++i) {
        kv_destroy(kv_A(*groups, i).group_ids);
    }
    kv_destroy(*groups);
}

static inline void recap_base_bit_u1v(u1v **seq, char *buffer, size_t len, int reverse)
{
    size_t i;
    if (!(*seq))
        *seq = init_u1v(len);

    clear_and_encap_u1v(*seq,  len);
    if (reverse) {
        for (i = 0; i < len; ++i) {
            (*seq)->buffer[i] = 3 - base_bit_table[(int)buffer[len - i - 1]];
        }

    } else {

        for (i = 0; i < len; ++i) {
            (*seq)->buffer[i] = base_bit_table[(u1i)buffer[i]];
        }
    }
    (*seq)->size = len;
}

static inline u1v *init_base_bit_u1v(char *buffer, size_t len, int reverse)
{
    u1v *seq = init_u1v(len);
    recap_base_bit_u1v(&seq, buffer, len, reverse);
    return seq;
}

static inline void destroy_base_bit_u1v(u1v *seq)
{
    free_u1v(seq);
}

static inline int strand_match(u1v *qseq, u1v *tseq, int similarity_pct, seqalign_result_t *rs_p)
{
    seqalign_result_t rs;
    u4v *cigars = init_u4v(64);;
    b1v *memp = adv_init_b1v(1024, 0, WORDSIZE, 0);
    size_t qlen = qseq->size;
    size_t tlen = tseq->size;

    //first align
    rs = kmer_striped_seqedit_pairwise(13, qseq->buffer, qlen, tseq->buffer, tlen, memp, cigars, 0);
#if 0
    {
        char *alnstr[3];
        alnstr[0] = malloc(rs.aln + 1);
        alnstr[1] = malloc(rs.aln + 1);
        alnstr[2] = malloc(rs.aln + 1);
        seqalign_cigar2alnstr(qseq->buffer, tseq->buffer, &rs, cigars, alnstr, rs.aln);
        fprintf(stdout, "%s\t%d\t+\t%d\t%d\t%s\t%d\t+\t%d\t%d\t", "QRY", Int(qseq->size), rs.qb, rs.qe, "REF", Int(tseq->size), rs.tb, rs.te);
        fprintf(stdout, "%d\t%.3f\t%d\t%d\t%d\t%d\n", rs.score, 1.0 * rs.mat / rs.aln, rs.mat, rs.mis, rs.ins, rs.del);
        fprintf(stdout, "%s\n%s\n%s\n", alnstr[0], alnstr[2], alnstr[1]);
        free(alnstr[0]);
        free(alnstr[1]);
        free(alnstr[2]);
    }
#endif
    if ((rs.aln * 2 > (qlen > tlen ? (int)tlen : (int)qlen)) && (rs.mat * 100 >= rs.aln * similarity_pct)) {
        if (rs_p) *rs_p = rs;
        free_u4v(cigars);
        free_b1v(memp);
        return 1;
    }

    free_u4v(cigars);
    free_b1v(memp);
    return 0;
}

typedef struct {
    u4i offs;
    u4i len;
    u1i reverse;
    u4i pos;
} segment_t;
typedef kvec_t(segment_t) vec_segment_t;

static inline int get_template_grp(zmw_t *zmw_p, vec_group_result_t *groups_p)
{
    u4i template_grp = 0;
    if ((kv_size(kv_A(*groups_p, template_grp).group_ids) < 2)) {
        return 0;
    }

    u1v *border_seq = NULL, *main_seq = NULL;
    u4i candidate_i, candidate_len;
    // adjust template
    //check groups with 80% of the first group
    for (u4i candidate_grp = 1; candidate_grp < kv_size(*groups_p); ++candidate_grp) {
        if ((kv_size(kv_A(*groups_p, candidate_grp).group_ids) < 2)
                || (kv_size(kv_A(*groups_p, candidate_grp).group_ids) * 5 < 4 * kv_size(kv_A(*groups_p, 0).group_ids))) {
            continue;
        }

        candidate_i = kv_A(*groups_p, candidate_grp).group_ids.a[kv_size(kv_A(*groups_p, candidate_grp).group_ids) / 2];
        candidate_len = kv_A(zmw_p->lens, candidate_i);
        if ((candidate_len <= kv_A(zmw_p->lens, kv_A(*groups_p, template_grp).group_ids.a[kv_size(kv_A(*groups_p, template_grp).group_ids) / 2]))
                || candidate_len <= 2000) {
            continue;
        }

        recap_base_bit_u1v(&border_seq, ks_str(&zmw_p->seqs) + kv_A(zmw_p->offs, candidate_i), 1000, 1);
        recap_base_bit_u1v(&main_seq, ks_str(&zmw_p->seqs) + kv_A(zmw_p->offs, candidate_i) + 1000, candidate_len - 1000, 0);
        if (strand_match(border_seq, main_seq, 70, NULL)) {
            //fprintf(stderr,"head match while change to %u hole=%s\n",candidate_grp,ks_str(&zmw_p->hole));
            continue;
        }
        recap_base_bit_u1v(&border_seq, ks_str(&zmw_p->seqs) + kv_A(zmw_p->offs, candidate_i) + candidate_len - 1000, 1000, 1);
        recap_base_bit_u1v(&main_seq, ks_str(&zmw_p->seqs) + kv_A(zmw_p->offs, candidate_i), candidate_len - 1000, 0);
        if (strand_match(border_seq, main_seq, 70, NULL)) {
            //fprintf(stderr,"tail match while change to %u hole=%s\n",candidate_grp,ks_str(&zmw_p->hole));
            continue;
        }
        template_grp = candidate_grp;
        //fprintf(stderr,"change to %u hole=%s\n",template_grp,ks_str(&zmw_p->hole));
    }
    if (border_seq) destroy_base_bit_u1v(border_seq);
    if (main_seq) destroy_base_bit_u1v(main_seq);
    return template_grp;
}

static inline vec_segment_t ccs_prepare(zmw_t *zmw_p)
{
    //init
    vec_segment_t segments;
    kv_init(segments);
    kv_resize(segment_t, segments, kv_size(zmw_p->lens));
    const int tolerance_pct = 10;
    //get groups
    vec_group_result_t groups = init_group_lens(zmw_p->lens.a, kv_size(zmw_p->lens), tolerance_pct);
    u4i map_group[kv_size(zmw_p->lens)];
    u4i i, j,i_size,j_size;
    for (i = 0,i_size=kv_size(groups); i <i_size ; ++i) {
        for (j = 0,j_size=kv_size(kv_A(groups, i).group_ids); j < j_size; ++j) {
            map_group[kv_A(kv_A(groups, i).group_ids, j) ] = i;
        }
    }

    seqalign_result_t rs;
    //
    const u4i template_grp = get_template_grp(zmw_p, &groups);
    const u4i template_i = kv_A(groups, template_grp).group_ids.a[kv_size(kv_A(groups, template_grp).group_ids) / 2];
    const u4i template_offs = kv_A(zmw_p->offs, template_i);
    const u4i template_len = kv_A(zmw_p->lens, template_i);


    u1v *tseq = NULL, *t2seq = NULL, *qseq = NULL;

    char reverse = 0, strand_adjust = 0;
    segment_t seg = {template_offs, template_len, reverse, 0};
    kv_push(segment_t, segments, seg);
    for (int k = template_i - 1; k >= 0; --k) {
        reverse = (reverse == 0) ? 1 : 0;
        seg.offs = kv_A(zmw_p->offs, k);
        seg.len = kv_A(zmw_p->lens, k);
        seg.reverse = reverse;
        if (map_group[k] != template_grp) { //not in the template_grp means abnormal length
            strand_adjust = 1;
            if (seg.len < template_len) continue;
        } else if (!strand_adjust) {
            kv_push(segment_t, segments, seg);
            continue;
        }
        if (!tseq) {
            tseq = init_base_bit_u1v(ks_str(&zmw_p->seqs) + template_offs, template_len, 0);
            t2seq = init_base_bit_u1v(ks_str(&zmw_p->seqs) + template_offs, template_len, 1);
        }

        qseq = init_base_bit_u1v(ks_str(&zmw_p->seqs) + seg.offs, seg.len, 0);
        if (strand_match(qseq, tseq, 75, &rs)) {
            reverse = 0;
            seg.offs += rs.qb, seg.len = rs.qe - rs.qb, seg.reverse = reverse;
            if (len_in_group(&kv_A(groups, template_grp), seg.len, tolerance_pct))
                kv_push(segment_t, segments, seg);
            strand_adjust = (map_group[k] != template_grp);
        } else if (strand_match(qseq, t2seq, 75, &rs)) {
            reverse = 1;
            seg.offs += rs.qb, seg.len = rs.qe - rs.qb, seg.reverse = reverse;
            if (len_in_group(&kv_A(groups, template_grp), seg.len, tolerance_pct))
                kv_push(segment_t, segments, seg);
            strand_adjust = (map_group[k] != template_grp);
        } else {
            strand_adjust = 1;//cannot be aligned
        }
        destroy_base_bit_u1v(qseq);
    }

    reverse = 0, strand_adjust = 0;
    for (i = template_i + 1; i < kv_size(zmw_p->lens); ++i) {
        reverse = (reverse == 0) ? 1 : 0;
        seg.offs = kv_A(zmw_p->offs, i);
        seg.len = kv_A(zmw_p->lens, i);
        seg.reverse = reverse;
        if (map_group[i] != template_grp) { //not in template_grp  means abnormal length
            strand_adjust = 1;
            if (seg.len < template_len) continue;
        } else if (!strand_adjust) {
            kv_push(segment_t, segments, seg);
            continue;
        }

        if (!tseq) {
            tseq = init_base_bit_u1v(ks_str(&zmw_p->seqs) + template_offs, template_len, 0);
            t2seq = init_base_bit_u1v(ks_str(&zmw_p->seqs) + template_offs, template_len, 1);
        }
        qseq = init_base_bit_u1v(ks_str(&zmw_p->seqs) + seg.offs, seg.len, 0);
        if (strand_match(qseq, tseq, 75, &rs)) {
            reverse = 0;
            seg.offs += rs.qb, seg.len = rs.qe - rs.qb, seg.reverse = reverse;
            if (len_in_group(&kv_A(groups, template_grp), seg.len, tolerance_pct))
                kv_push(segment_t, segments, seg);
            strand_adjust = (map_group[i] != template_grp);
        } else if (strand_match(qseq, t2seq, 75, &rs)) {
            reverse = 1;
            seg.offs += rs.qb, seg.len = rs.qe - rs.qb, seg.reverse = reverse;
            if (len_in_group(&kv_A(groups, template_grp), seg.len, tolerance_pct))
                kv_push(segment_t, segments, seg);
            strand_adjust = (map_group[i] != template_grp);
        } else {
            strand_adjust = 1;//cannot be aligned
        }
        destroy_base_bit_u1v(qseq);
    }

    destroy_group_lens(&groups);
    if (tseq) {
        destroy_base_bit_u1v(tseq);
        destroy_base_bit_u1v(t2seq);
    }
    return segments;
}

static void ccs_for(void *_data, long idx, int tid) // kt_for() callback
{
    step_t *step = (step_t *)_data;
    zmw_t *zmw_p = &kv_A(step->zmws, idx);
    BSPOA *g = step->gg[tid];
    if (kv_size(zmw_p->lens) < 3) {
        return;
    }

    //done
    vec_segment_t segments = ccs_prepare(zmw_p);
    if (step->verbose > 1)
        fprintf(stderr, "poa begin %s\n", ks_str(&zmw_p->hole));

    // adjust strand
    segment_t *seg_p = NULL;
    for (u4i l = 0,l_size=kv_size(segments); l <l_size ; ++l) {
        seg_p = &kv_A(segments, l);
        if (seg_p->reverse) {
            seq_reverse_comp(seg_p->len, (unsigned char *)ks_str(&zmw_p->seqs) + seg_p->offs);
        }

        if (step->verbose)
            fprintf(stderr, ">%s_%d/%lu strand=%d len=%d \n%.*s\n", ks_str(&zmw_p->hole), l, kv_size(segments),
                    seg_p->reverse, seg_p->len, seg_p->len, ks_str(&zmw_p->seqs) + seg_p->offs);
    }

    zmw_p->ccsseq.l = zmw_p->ccsseq.m = 0, zmw_p->ccsseq.s = NULL;
    ks_resize(&zmw_p->ccsseq, 10240);


    beg_bspoa(g);
    /*push seq */
    for (u4i l = 0,l_size=kv_size(segments); l <l_size; ++l) {
        seg_p = &kv_A(segments, l);
        push_bspoa(g, ks_str(&zmw_p->seqs) + seg_p->offs, seg_p->len);
    }
    end_bspoa(g);

    /*save ccs seq*/
    char *res = malloc(g->cns->size + 1);
    for (u4i l = 0,l_size=g->cns->size; l <l_size ; ++l) {
        res[l] = bit_base_table[g->cns->buffer[l]];
    }
    res[g->cns->size] = 0;
    zmw_p->ccsseq.l = g->cns->size;
    zmw_p->ccsseq.m = g->cns->size + 1;
    zmw_p->ccsseq.s = res;

    kv_destroy(segments);

    if (step->verbose > 1)
        fprintf(stderr, "poa end %s\n", ks_str(&zmw_p->hole));
}

static void ccs_for2(void *_data, long idx, int tid) // kt_for() callback
{
    step_t *step = (step_t *)_data;
    zmw_t *zmw_p = &kv_A(step->zmws, idx);
    BSPOA *g = step->gg[tid];
    if (kv_size(zmw_p->lens) < 3) {
        return;
    }

    //done
    vec_segment_t segments = ccs_prepare(zmw_p);
    if (step->verbose > 1)
        fprintf(stderr, "poa begin %s\n", ks_str(&zmw_p->hole));


    // adjust strand
    segment_t *seg_p = NULL;
    for (u4i l = 0; l < kv_size(segments); ++l) {
        seg_p = &kv_A(segments, l);
        if (seg_p->reverse) {
            seq_reverse_comp(seg_p->len, (unsigned char *)ks_str(&zmw_p->seqs) + seg_p->offs);
        }

        if (step->verbose)
            fprintf(stderr, ">%s_%d/%lu strand=%d len=%d \n%.*s\n", ks_str(&zmw_p->hole), l, kv_size(segments),
                    seg_p->reverse, seg_p->len, seg_p->len, ks_str(&zmw_p->seqs) + seg_p->offs);
    }

    zmw_p->ccsseq.l = zmw_p->ccsseq.m = 0, zmw_p->ccsseq.s = NULL;
    ks_resize(&zmw_p->ccsseq, 10240);

    const u4i window = 10, addlen = 2000, minlen = 1000, initlen = 2000,  minwin = 5;
    u4i rowrate = 80, colrate = 80;
    const u4i nseq = kv_size(segments);
    const u4i mrow = nseq + 1 + 3;
    u4i flag = 1;
    if (nseq < 10) colrate = 60;
    while (flag) {
        u4i i,i_size;
        //get ccs of cur window,window size can be extended
        for (u4i window_size = initlen; ; window_size += addlen) {
            //ccs
            beg_bspoa(g);
            for (i = 0,i_size=kv_size(segments); i <i_size ; ++i) {
                seg_p = &kv_A(segments, i);
                if (seg_p->pos + window_size + minlen >= seg_p->len) {
                    break;
                }
            }
            if ((i < i_size) || (i_size < 3)) { //too short or too few
                flag = 0;
                for (i = 0; i < i_size; ++i) {
                    seg_p = &kv_A(segments, i);
                    push_bspoa(g, ks_str(&zmw_p->seqs) + seg_p->offs + seg_p->pos, seg_p->len - seg_p->pos);
                }
            } else {
                for (i = 0; i < i_size; ++i) {
                    seg_p = &kv_A(segments, i);
                    push_bspoa(g, ks_str(&zmw_p->seqs) + seg_p->offs + seg_p->pos, window_size);
                }
            }
            end_bspoa(g);
            tidy_msa_bspoa(g);

            if (!flag) {
                i = g->msaidxs->size;
                break;
            }

            //find breakpoint  around the window boundary
            for (i = g->msaidxs->size - window; i >= 1; i--) {
                u4i j, k, nogwin = 0;
                u1i rowcnt[nseq];
                memset(rowcnt, 0, nseq * sizeof(u1i));
                for (j = i; j < i + window; ++j) {
                    u1i *col = g->msacols->buffer + g->msaidxs->buffer[j] * mrow;
                    if (col[nseq + 1] >= 4) {
                        if (nogwin) continue;
                        else break;
                    }
                    ++nogwin;
                    u4i colcnt = 0;
                    for (k = 0; k < nseq; k++) {
                        if (col[k + 1] == col[nseq + 1]) {
                            colcnt++;
                            rowcnt[k]++;
                        }
                    }
                    if (colcnt * 100  < colrate * nseq)
                        break;

                }
                if ((j < i + window) || (nogwin < minwin)) {
                    continue;
                }

                for (k = 0; k < nseq; ++k) {
                    if (rowcnt[k] * 100  < rowrate * nogwin)
                        break;
                }
                if (k >= nseq) //find cluster
                    break;
            }
            if (i >= 1) { //now find the breakpoint i
                break;
            }
        }

        //collect result from breakpoint
        if (step->verbose > 2)
            fprintf(stdout, "breakpoint=%u maplen=%llu nseq=%u hole=%s\n", i, g->msaidxs->size,  nseq, ks_str(&zmw_p->hole));

        for (u4i j = 0; j < i; j++) {
            u1i *col = g->msacols->buffer + g->msaidxs->buffer[j] * mrow;
            // if not finished,adjust the pos
            if (flag) {
				const u4i k_size=kv_size(segments);
                for (u4i k = 0; k < k_size; ++k) {
                    seg_p = &kv_A(segments, k);
                    if (col[k + 1] < 4) {
                        ++seg_p->pos;
                    }
                }
            }

            if (col[nseq + 1] < 4) {
                kputc(bit_base_table[col[nseq + 1]], &zmw_p->ccsseq);
            }
        }

        //done
    }

    kv_destroy(segments);

    if (step->verbose > 1)
        fprintf(stderr, "poa end %s\n", ks_str(&zmw_p->hole));
}

static void *worker_pipeline(void *shared, int step, void *in) // kt_pipeline() callback
{
    pipeline_t *p = (pipeline_t *)shared;
    if (step == 0) { // step 0: read zmw into the buffer
        step_t *s = step_init(p->chunk_size);
        s->verbose = p->verbose;
        s->gg = p->gg;
        kseqs_zmw_t *zmw_seqs = &p->zmw_seqs_in;
        int l;
        while ((l = kseq_zmw_read(zmw_seqs)) >= 0) {
            if (l < p->min_fulllen_count + 2)
                continue;

            size_t total_len = ks_len(&zmw_seqs->seqs);
            if (total_len > p->max_subread_len || total_len < p->min_subread_len) {
                continue;
            }

            if (p->hole_set) {
                khash_t(strset)* hole_set = p->hole_set;
                if (kh_get(strset, hole_set, ks_str(&zmw_seqs->hole)) != kh_end(hole_set)) {
                    continue;
                }
            }

            zmw_t zmw;
            zmw_initialize(&zmw, 8);
            kputsn(ks_str(&zmw_seqs->hole), ks_len(&zmw_seqs->hole), &zmw.hole);
            kputsn(ks_str(&zmw_seqs->movie_name), ks_len(&zmw_seqs->movie_name), &zmw.movie_name);
            kputsn(ks_str(&zmw_seqs->seqs), ks_len(&zmw_seqs->seqs), &zmw.seqs);
            kv_copy(u4i, zmw.lens, zmw_seqs->lens);
            int i, offs = 0, len,i_size=kv_size(zmw_seqs->lens);
            for (i = 0; i < i_size; ++i, offs += len) {
                len = kv_A(zmw_seqs->lens, i);
                kv_push(u4i, zmw.offs, offs);
            }
            kv_push(zmw_t, s->zmws, zmw);
            if (kv_size(s->zmws) >= p->chunk_size) {
                if (p->chunk_size < 16384) {
                    p->chunk_size *= 4;
                }
                break;
            }
        }

        if (kv_size(s->zmws))
            return s;
        else
            step_destroy(s);
    } else if (step == 1) { // step 1: ccs
        step_t *s = (step_t *)in;
        if (s) {
            if (p->split_subread)
                kt_for(p->nthreads, ccs_for2, in, kv_size(s->zmws));
            else
                kt_for(p->nthreads, ccs_for, in, kv_size(s->zmws));
        }
        return in;
    } else if (step == 2) { // step 2: write ccs seq to output
        step_t *s = (step_t *)in;
        if (s) {
            size_t i = 0,i_size=kv_size(s->zmws);
            for (i = 0; i < i_size; ++i) {
                zmw_t *zmw_p = &kv_A(s->zmws, i);
                if (ks_len(&zmw_p->ccsseq))
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
            "-m     <int>   Minimum total length of subreads in a hole to use for generating CCS. [5000] \n"
            "-M     <int>   Maximum total length of subreads in a hole to use for generating CCS. [500000] \n"
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
    int c, verbose = 0, min_subread_len = 5000, max_subread_len = 500000,
           min_fulllen_count = 3, nthreads = 1, isbam = 1, split_subread = 1;

    khash_t(strset)* hole_set = NULL;
    kstring_t strholes = {0};
    while ((c = getopt(argc, argv, "hm:M:c:j:X:PAv")) != -1) {
        switch (c) {
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
            case 'X': {
                kputs(optarg, &strholes);
                int n, i, *fields;
                fields = ksplit(&strholes, ',', &n);
                hole_set = kh_init(strset); // allocate a hash table
                for (i = 0; i < n; ++i) {
                    int ret;
                    kh_put(strset, hole_set, ks_str(&strholes) + fields[i], &ret);
                }
                free(fields);
            }
            break;
            case 'c':
                min_fulllen_count = atoi(optarg);
                if (min_fulllen_count < 3) {
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
    if (argc - optind == 0) {
        fp_in = gzdopen(fileno(stdin), "rb") ;
        fp_out = stdout;
    } else if (argc - optind == 1) {
        fp_in = (strcmp(argv[optind], "-") == 0) ? gzdopen(fileno(stdin), "rb") : gzopen(argv[optind], "rb");
        fp_out = stdout;
    } else if (argc - optind == 2) {
        fp_in = (strcmp(argv[optind], "-") == 0) ? gzdopen(fileno(stdin), "rb") : gzopen(argv[optind], "rb");
        fp_out = (strcmp(argv[optind + 1], "-") == 0) ? stdout : fopen(argv[optind + 1], "w+");
    } else {
        return usage();
    }

    /*read infile*/
    if (!fp_in) {
        fprintf(stderr, "Error: Failed to open infile!\n");
        return 1;
    }

    if (!fp_out) {
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
    for (i = 0; i < nthreads; ++i) {
        BSPOAPar par = DEFAULT_BSPOA_PAR;
        par.M = 2;
        par.X = -6;
        par.O = -3;
        par.E = -2;
        par.Q = 0;
        par.P = 0;
        par.editbw = 32;
        par.bandwidth = 128;
        //par.refmode=1;
        pl.gg[i] = init_bspoa(par);
    }

    kseqs_zmw_initialize(&pl.zmw_seqs_in, fp_in, isbam);

    kt_pipeline(2, worker_pipeline, &pl, 3);
    for (i = 0; i < nthreads; ++i) {
        free_bspoa(pl.gg[i]);
    }
    free(pl.gg);
    zmw_seqs_release(&pl.zmw_seqs_in);
    if (pl.hole_set) {
        kh_destroy(strset, hole_set);
    }
    free(ks_release(&strholes));
    fclose(fp_out);
    gzclose(fp_in);

    return 0;
}
