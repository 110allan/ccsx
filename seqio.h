#ifndef _SEQ_IO_H
#define _SEQ_IO_H

#include "kseq.h"
#include "bamlite.h"

#define __SEQIO_TYPE(type_t) \
        typedef struct __kseq_extend_t {\
            kseq_t *seq;\
            int is_bam;\
            bam1_t *b;\
        } kseq_extend_t;\
        \
        typedef struct __kseqs_zmw_t {\
            kvec_t(int) lens;\
        kseq_extend_t seq_extend;\
            kstring_t movie_name, hole, seqs;\
            kstring_t last_movie_name, last_hole, last_seqs;\
        } kseqs_zmw_t;

#define __SEQIO_BASIC(SCOPE,type_t) \
        SCOPE inline void kseq_extend_initialize(kseq_extend_t *s, type_t fd, int isbam)\
        {\
            s->seq = kseq_init(fd);\
            s->is_bam = isbam;\
            if (isbam) {\
                bam_header_t *h = bam_header_read(fd);\
                if (!h) {\
                    fprintf(stderr, "[bam_header_read] invalid BAM header.\n");\
                    kseq_destroy(s->seq);\
                }\
                bam_header_destroy(h);\
                s->b = bam_init1();\
            }\
        }\
        \
        SCOPE inline void kseq_extend_release(kseq_extend_t *ks_extend)\
        {\
            if (!ks_extend)\
                return;\
            kseq_destroy(ks_extend->seq);\
            if (ks_extend->is_bam) {\
                bam_destroy1(ks_extend->b);\
            }\
        }\
        \
        SCOPE inline void kseqs_zmw_initialize(kseqs_zmw_t *s, type_t fd, int isbam)\
        {\
            kseq_extend_initialize(&s->seq_extend, fd, isbam);\
            kv_init(s->lens);\
            kv_resize(int, s->lens, 8);\
            \
            s->movie_name.l = s->hole.l = s->seqs.l = 0;\
            s->movie_name.m = s->hole.m = s->seqs.m = 0;\
            s->movie_name.s = s->hole.s = s->seqs.s = 0;\
            ks_resize(&s->movie_name, 32);\
            ks_resize(&s->hole, 16);\
            ks_resize(&s->seqs, 16384);\
            \
            s->last_movie_name.l = s->last_hole.l = s->last_seqs.l = 0;\
            s->last_movie_name.s = s->last_hole.s = s->last_seqs.s = 0;\
            s->last_movie_name.m = s->last_hole.m = s->last_seqs.m = 0;\
            ks_resize(&s->last_movie_name, 32);\
            ks_resize(&s->last_hole, 16);\
            ks_resize(&s->last_seqs, 16384);\
        }\
        \
        SCOPE inline void zmw_seqs_release(kseqs_zmw_t *zmw_seqs)\
        {\
            if (!zmw_seqs)\
                return;\
            kseq_extend_release(&zmw_seqs->seq_extend);\
            kv_destroy(zmw_seqs->lens);\
            char *s = NULL;\
            s = ks_release(&zmw_seqs->movie_name);\
            free(s);\
            s = ks_release(&zmw_seqs->hole);\
            free(s);\
            s = ks_release(&zmw_seqs->seqs);\
            free(s);\
            s = ks_release(&zmw_seqs->last_movie_name);\
            free(s);\
            s = ks_release(&zmw_seqs->last_hole);\
            free(s);\
            s = ks_release(&zmw_seqs->last_seqs);\
            free(s);\
        }



#define __SEQIO_READ(SCOPE) \
        SCOPE const char seq_nt16_str[] = "=ACMGRSVTWYHKDBN";\
        SCOPE int kseq_extend_read(kseq_extend_t *ks_extend)\
        {\
            if (!ks_extend->is_bam) {\
                return kseq_read(ks_extend->seq);\
            }\
            kstream_t *ks = ks_extend->seq->f;\
            bam1_t *b = ks_extend->b;\
            int bytes_read = bam_read1(ks->f, b);\
            if (bytes_read <= 0)\
                return -1;\
            uint8_t *s = bam1_seq(b);\
            uint8_t *q = bam1_qual(b);\
            int len_qseq = b->core.l_qseq;\
            ks_extend->seq->comment.l = ks_extend->seq->name.l = 0;\
            kputs((const char *)bam1_qname(b), &ks_extend->seq->name);\
            ks_resize(&ks_extend->seq->seq, len_qseq + 1);\
            ks_resize(&ks_extend->seq->qual, len_qseq + 1);\
            int i;\
            for (i = 0; i != len_qseq; ++i){\
                ks_extend->seq->seq.s[i] = seq_nt16_str[(int)bam1_seqi(s, i)];\
                ks_extend->seq->qual.s[i] = q[i] + 33 < 126 ? q[i] + 33 : 126;\
            }\
            ks_extend->seq->seq.s[len_qseq] = ks_extend->seq->qual.s[len_qseq] = 0;\
            ks_extend->seq->seq.l = ks_extend->seq->qual.l = len_qseq;\
            return len_qseq;\
        }\
        \
        SCOPE unsigned char seq_comp_table[256] = {\
                                                   0,    1,   2,   3,   4,   5,   6,   7,   8,   9,  10,  11,  12,  13,  14,  15, \
                                                   16,  17,  18,  19,   20,  21,  22,  23,  24,  25,  26,  27,  28,  29,  30,  31, \
                                                   32,  33,  34,  35,   36,  37,  38,  39,  40,  41,  42,  43,  44,  45,  46,  47, \
                                                   48,  49,  50,  51,   52,  53,  54,  55,  56,  57,  58,  59,  60,  61,  62,  63, \
                                                   64, 'T', 'V', 'G', 'H', 'E', 'F', 'C', 'D', 'I', 'J', 'M', 'L', 'K', 'N', 'O', \
                                                   'P', 'Q', 'Y', 'S', 'A', 'A', 'B', 'W', 'X', 'R', 'Z',   91,  92,  93,  94,  95, \
                                                   96, 't', 'v', 'g', 'h', 'e', 'f', 'c', 'd', 'i', 'j', 'm', 'l', 'k', 'n', 'o', \
                                                   'p', 'q', 'y', 's', 'a', 'a', 'b', 'w', 'x', 'r', 'z', 123, 124, 125, 126, 127, \
                                                   128, 129, 130, 131, 132, 133, 134, 135, 136, 137, 138, 139, 140, 141, 142, 143, \
                                                   144, 145, 146, 147, 148, 149, 150, 151, 152, 153, 154, 155, 156, 157, 158, 159, \
                                                   160, 161, 162, 163, 164, 165, 166, 167, 168, 169, 170, 171, 172, 173, 174, 175, \
                                                   176, 177, 178, 179, 180, 181, 182, 183, 184, 185, 186, 187, 188, 189, 190, 191, \
                                                   192, 193, 194, 195, 196, 197, 198, 199, 200, 201, 202, 203, 204, 205, 206, 207, \
                                                   208, 209, 210, 211, 212, 213, 214, 215, 216, 217, 218, 219, 220, 221, 222, 223, \
                                                   224, 225, 226, 227, 228, 229, 230, 231, 232, 233, 234, 235, 236, 237, 238, 239, \
                                                   240, 241, 242, 243, 244, 245, 246, 247, 248, 249, 250, 251, 252, 253, 254, 255\
                                                  };\
        SCOPE inline void seq_reverse_comp(int l, unsigned char *seq)\
        {\
            int i, t;\
            for (i = 0; i < l >> 1; ++i) {\
                t = seq[l - 1 - i];\
                seq[l - i - 1] = seq_comp_table[seq[i]];\
                seq[i] = seq_comp_table[t];\
            }\
            if (l & 1)\
                seq[l >> 1] = seq_comp_table[seq[l >> 1]];\
        }\


#define __SEQIO_ZMW_READ(SCOPE) \
        SCOPE int kseq_zmw_read(kseqs_zmw_t *zmw_p)\
        {\
            int l;\
            /*clear current zmw*/\
            kv_size(zmw_p->lens) = 0;\
            zmw_p->movie_name.l = zmw_p->hole.l = zmw_p->seqs.l = 0;\
            if (zmw_p->last_movie_name.l != 0) {\
                kputsn((const char *)(zmw_p->last_hole.s), ks_len(&zmw_p->last_hole), &zmw_p->hole);\
                kputsn((const char *)(zmw_p->last_movie_name.s), ks_len(&zmw_p->last_movie_name), &zmw_p->movie_name);\
                kputsn((const char *)(zmw_p->last_seqs.s), ks_len(&zmw_p->last_seqs), &zmw_p->seqs);\
                kv_push(int, zmw_p->lens, ks_len(&zmw_p->last_seqs));\
            }\
            while ((l = kseq_extend_read(&zmw_p->seq_extend)) >= 0) {\
                int n;\
                kseq_t *seq = zmw_p->seq_extend.seq;\
                int *fields = ksplit(&seq->name, '/', &n);\
                if (n != 3) {\
                    free(fields);\
                    fprintf(stderr, "invalid zmw name :%s\n", seq->name.s);\
                    return -1;\
                }\
                if (zmw_p->last_movie_name.l == 0) {\
                    kputs((const char *)ks_str(&seq->name) + fields[0], &zmw_p->movie_name);\
                    kputs((const char *)ks_str(&seq->name) + fields[1], &zmw_p->hole);\
                    kputsn((const char *)ks_str(&seq->seq), ks_len(&seq->seq), &zmw_p->seqs);\
                    kv_push(int, zmw_p->lens, ks_len(&seq->seq));\
                    \
                    zmw_p->last_hole.l = zmw_p->last_movie_name.l = zmw_p->last_seqs.l = 0;\
                    kputsn((const char *)(zmw_p->hole.s), ks_len(&zmw_p->hole), &zmw_p->last_hole);\
                    kputsn((const char *)(zmw_p->movie_name.s), ks_len(&zmw_p->movie_name), &zmw_p->last_movie_name);\
                    kputsn((const char *)ks_str(&seq->seq), ks_len(&zmw_p->seqs), &zmw_p->last_seqs);\
                } else if (strcmp(zmw_p->last_hole.s, seq->name.s + fields[1]) || strcmp(zmw_p->movie_name.s, seq->name.s + fields[0])) {\
                    zmw_p->last_movie_name.l = zmw_p->last_hole.l = zmw_p->last_seqs.l = 0;\
                    kputs((const char *)ks_str(&seq->name) + fields[0], &zmw_p->last_movie_name);\
                    kputs((const char *)ks_str(&seq->name) + fields[1], &zmw_p->last_hole);\
                    kputsn((const char *)ks_str(&seq->seq), seq->seq.l, &zmw_p->last_seqs);\
                    free(fields);\
                    return kv_size(zmw_p->lens);\
                } else {\
                    kputsn((const char *)(seq->seq.s), ks_len(&seq->seq), &zmw_p->seqs);\
                    kv_push(int, zmw_p->lens, ks_len(&seq->seq));\
                }\
                free(fields);\
            }\
            zmw_p->last_movie_name.l = zmw_p->last_hole.l = zmw_p->last_seqs.l = 0;\
            if (kv_size(zmw_p->lens) > 0) {\
                return kv_size(zmw_p->lens);\
            }\
            return -1;\
        }



#define SEQIO_INIT2(SCOPE, type_t,__read) \
        KSEQ_INIT(type_t, __read)\
        __SEQIO_TYPE(type_t) \
        __SEQIO_BASIC(SCOPE, type_t) \
        __SEQIO_READ(SCOPE)\
        __SEQIO_ZMW_READ(SCOPE)

#define SEQIO_INIT(type_t,__read) SEQIO_INIT2(static,type_t, __read)


/*
SEQIO_INIT(gzFile, gzread)
int main(int argc, char *argv[])
{
    gzFile fp;
    zmw_seqs_t *zmw;
    int l;
    if (argc == 1) {
        fprintf(stderr, "Usage: %s <in.fasta>\n", argv[0]);
        return 1;
    }
    fp = gzopen(argv[1], "r");
    zmw = zmw_seqs_init(fp, 1);
    while ((l = zmw_read(zmw)) >= 0) {
        printf(">%s_%s\n", ks_str(&zmw->movie_name),ks_str(&zmw->hole));
        printf("%s\n", ks_str(&zmw->seqs));
    }
    printf("return value: %d\n", l);
    zmw_seqs_destroy(zmw);
    gzclose(fp);
    return 0;
}

*/

#endif
