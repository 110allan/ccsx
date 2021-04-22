#include <stdlib.h>
#include <ctype.h>
#include <string.h>
#include <stdio.h>
#include <errno.h>
#include "bamlite.h"

#ifdef USE_MALLOC_WRAPPERS
#  include "malloc_wrap.h"
#endif

/*********************
 * from bam_endian.c *
 *********************/

static inline int bam_is_big_endian()
{
	long one= 1;
	return !(*((char *)(&one)));
}
static inline uint16_t bam_swap_endian_2(uint16_t v)
{
	return (uint16_t)(((v & 0x00FF00FFU) << 8) | ((v & 0xFF00FF00U) >> 8));
}
static inline void *bam_swap_endian_2p(void *x)
{
	*(uint16_t*)x = bam_swap_endian_2(*(uint16_t*)x);
	return x;
}
static inline uint32_t bam_swap_endian_4(uint32_t v)
{
	v = ((v & 0x0000FFFFU) << 16) | (v >> 16);
	return ((v & 0x00FF00FFU) << 8) | ((v & 0xFF00FF00U) >> 8);
}
static inline void *bam_swap_endian_4p(void *x)
{
	*(uint32_t*)x = bam_swap_endian_4(*(uint32_t*)x);
	return x;
}
static inline uint64_t bam_swap_endian_8(uint64_t v)
{
	v = ((v & 0x00000000FFFFFFFFLLU) << 32) | (v >> 32);
	v = ((v & 0x0000FFFF0000FFFFLLU) << 16) | ((v & 0xFFFF0000FFFF0000LLU) >> 16);
	return ((v & 0x00FF00FF00FF00FFLLU) << 8) | ((v & 0xFF00FF00FF00FF00LLU) >> 8);
}
static inline void *bam_swap_endian_8p(void *x)
{
	*(uint64_t*)x = bam_swap_endian_8(*(uint64_t*)x);
	return x;
}

/**************
 * from bam.c *
 **************/

int bam_is_be;

bam_header_t *bam_header_init()
{
	bam_is_be = bam_is_big_endian();
	return (bam_header_t*)calloc(1, sizeof(bam_header_t));
}

void bam_header_destroy(bam_header_t *header)
{
	int32_t i;
	if (header == 0) return;
	if (header->target_name) {
		for (i = 0; i < header->n_targets; ++i)
			if (header->target_name[i]) free(header->target_name[i]);
		if (header->target_len) free(header->target_len);
		free(header->target_name);
	}
	if (header->text) free(header->text);
	free(header);
}

bam_header_t *bam_header_read(bamFile fp)
{
	bam_header_t *header;
	char buf[4];
	int magic_len;
	int32_t i = 1, name_len;
	// read "BAM1"
	magic_len = bam_read(fp, buf, 4);
	if (magic_len != 4 || strncmp(buf, "BAM\001", 4) != 0) {
		fprintf(stderr, "[bam_header_read] invalid BAM binary header (this is not a BAM file).\n");
		return NULL;
	}
	header = bam_header_init();
	// read plain text and the number of reference sequences
	if (bam_read(fp, &header->l_text, 4) != 4) goto fail; 
	if (bam_is_be) bam_swap_endian_4p(&header->l_text);
	header->text = (char*)calloc(header->l_text + 1, 1);
	if (bam_read(fp, header->text, header->l_text) != header->l_text) goto fail;
	if (bam_read(fp, &header->n_targets, 4) != 4) goto fail;
	if (bam_is_be) bam_swap_endian_4p(&header->n_targets);
	// read reference sequence names and lengths
	header->target_name = (char**)calloc(header->n_targets, sizeof(char*));
	header->target_len = (uint32_t*)calloc(header->n_targets, 4);
	for (i = 0; i != header->n_targets; ++i) {
		if (bam_read(fp, &name_len, 4) != 4) goto fail;
		if (bam_is_be) bam_swap_endian_4p(&name_len);
		header->target_name[i] = (char*)calloc(name_len, 1);
		if (bam_read(fp, header->target_name[i], name_len) != name_len) {
			goto fail;
		}
		if (bam_read(fp, &header->target_len[i], 4) != 4) goto fail;
		if (bam_is_be) bam_swap_endian_4p(&header->target_len[i]);
	}
	return header;
 fail:
	bam_header_destroy(header);
	return NULL;
}

static void swap_endian_data(const bam1_core_t *c, int data_len, uint8_t *data)
{
	uint8_t *s;
	uint32_t i, *cigar = (uint32_t*)(data + c->l_qname);
	s = data + c->n_cigar*4 + c->l_qname + c->l_qseq + (c->l_qseq + 1)/2;
	for (i = 0; i < c->n_cigar; ++i) bam_swap_endian_4p(&cigar[i]);
	while (s < data + data_len) {
		uint8_t type;
		s += 2; // skip key
		type = toupper(*s); ++s; // skip type
		if (type == 'C' || type == 'A') ++s;
		else if (type == 'S') { bam_swap_endian_2p(s); s += 2; }
		else if (type == 'I' || type == 'F') { bam_swap_endian_4p(s); s += 4; }
		else if (type == 'D') { bam_swap_endian_8p(s); s += 8; }
		else if (type == 'Z' || type == 'H') { while (*s) ++s; ++s; }
	}
}

int bam_read1(bamFile fp, bam1_t *b)
{
	bam1_core_t *c = &b->core;
	int32_t block_len, ret, i;
	uint32_t x[8];

	if ((ret = bam_read(fp, &block_len, 4)) != 4) {
		if (ret == 0) return -1; // normal end-of-file
		else return -2; // truncated
	}
	if (bam_read(fp, x, sizeof(bam1_core_t)) != sizeof(bam1_core_t)) return -3;
	if (bam_is_be) {
		bam_swap_endian_4p(&block_len);
		for (i = 0; i < 8; ++i) bam_swap_endian_4p(x + i);
	}
	c->tid = x[0]; c->pos = x[1];
	c->bin = x[2]>>16; c->qual = x[2]>>8&0xff; c->l_qname = x[2]&0xff;
	c->flag = x[3]>>16; c->n_cigar = x[3]&0xffff;
	c->l_qseq = x[4];
	c->mtid = x[5]; c->mpos = x[6]; c->isize = x[7];
	b->data_len = block_len - sizeof(bam1_core_t);
	if (b->m_data < b->data_len) {
		b->m_data = b->data_len;
		kroundup32(b->m_data);
		b->data = (uint8_t*)realloc(b->data, b->m_data);
	}
	if (bam_read(fp, b->data, b->data_len) != b->data_len) return -4;
	b->l_aux = b->data_len - c->n_cigar * 4 - c->l_qname - c->l_qseq - (c->l_qseq+1)/2;
	if (bam_is_be) swap_endian_data(c, b->data_len, b->data);
	return 4 + block_len;
}

#define bam_get_aux(b)   ((b)->data + ((b)->core.n_cigar<<2) + (b)->core.l_qname + (((b)->core.l_qseq + 1)>>1) + (b)->core.l_qseq)

static inline int aux_type2size(uint8_t type)
{
    switch (type) {
    case 'A': case 'c': case 'C':
        return 1;
    case 's': case 'S':
        return 2;
    case 'i': case 'I': case 'f':
        return 4;
    case 'd':
        return 8;
    case 'Z': case 'H': case 'B':
        return type;
    default:
        return 0;
    }
}

static inline uint8_t* skip_aux(uint8_t* s) {
	int size = aux_type2size(*s); ++s; // skip type
	uint32_t n;
	switch (size) {
	case 'Z':
	case 'H':
		while (*s) ++s;
		return s+1;
	case 'B':
		size = aux_type2size(*s); ++s;
		memcpy(&n, s, 4); s += 4;
		return s + size * n;
	case 0:
		abort();
		break;
	default:
		return s + size;
	}
}

#define __skip_tag(s) do { \
		int type = toupper(*(s)); \
		++(s); \
		if (type == 'Z' || type == 'H') { while (*(s)) ++(s); ++(s); } \
		else if (type == 'B') (s) += 5 + bam_aux_type2size(*(s)) * (*(int32_t*)((s)+1)); \
		else (s) += bam_aux_type2size(*(s)); \
	} while(0)

uint8_t *bam_aux_get(const bam1_t *b, const char tag[2])
{
	uint8_t *s;
	int y = tag[0]<<8 | tag[1];
	s = bam1_aux(b);
	while (s < b->data + b->data_len) {
		int x = (int)s[0]<<8 | s[1];
		s += 2;
		if (x == y) return s;
		//__skip_tag(s);
		s=skip_aux(s);
	}
	return 0;
}
// s MUST BE returned by bam_aux_get()
int bam_aux_del(bam1_t *b, uint8_t *s)
{
	uint8_t *p, *aux;
	aux = bam1_aux(b);
	p = s - 2;
	//__skip_tag(s);
	s=skip_aux(s);
	memmove(p, s, b->l_aux - (s - aux));
	b->data_len -= s - p;
	b->l_aux -= s - p;
	return 0;
}

int32_t bam_aux2i(const uint8_t *s)
{
	int type;
	if (s == 0) return 0;
	type = *s++;
	if (type == 'c') return (int32_t)*(int8_t*)s;
	else if (type == 'C') return (int32_t)*(uint8_t*)s;
	else if (type == 's') return (int32_t)*(int16_t*)s;
	else if (type == 'S') return (int32_t)*(uint16_t*)s;
	else if (type == 'i' || type == 'I') return *(int32_t*)s;
	else return 0;
}

float bam_aux2f(const uint8_t *s)
{
	int type;
	type = *s++;
	if (s == 0) return 0.0;
	if (type == 'f') return *(float*)s;
	else return 0.0;
}

double bam_aux2d(const uint8_t *s)
{
	int type;
	type = *s++;
	if (s == 0) return 0.0;
	if (type == 'd') return *(double*)s;
	else return 0.0;
}

char bam_aux2A(const uint8_t *s)
{
	int type;
	type = *s++;
	if (s == 0) return 0;
	if (type == 'A') return *(char*)s;
	else return 0;
}

char *bam_aux2Z(const uint8_t *s)
{
	int type;
	type = *s++;
	if (s == 0) return 0;
	if (type == 'Z' || type == 'H') return (char*)s;
	else return 0;
}

#ifdef USE_VERBOSE_ZLIB_WRAPPERS
// Versions of gzopen, gzread and gzclose that print up error messages

gzFile bamlite_gzopen(const char *fn, const char *mode) {
	gzFile fp;
	if (strcmp(fn, "-") == 0) {
		fp = gzdopen(fileno((strstr(mode, "r"))? stdin : stdout), mode);
		if (!fp) {
			fprintf(stderr, "Couldn't open %s : %s",
					(strstr(mode, "r"))? "stdin" : "stdout",
					strerror(errno));
		}
		return fp;
	}
	if ((fp = gzopen(fn, mode)) == 0) {
		fprintf(stderr, "Couldn't open %s : %s\n", fn,
				errno ? strerror(errno) : "Out of memory");
	}
	return fp;
}

int bamlite_gzread(gzFile file, void *ptr, unsigned int len) {
	int ret = gzread(file, ptr, len);
	
	if (ret < 0) {
		int errnum = 0;
		const char *msg = gzerror(file, &errnum);
		fprintf(stderr, "gzread error: %s\n",
				Z_ERRNO == errnum ? strerror(errno) : msg);
	}
	return ret;
}

int bamlite_gzclose(gzFile file) {
	int ret = gzclose(file);
	if (Z_OK != ret) {
		fprintf(stderr, "gzclose error: %s\n",
						  Z_ERRNO == ret ? strerror(errno) : zError(ret));
	}
	
	return ret;
}
#endif /* USE_VERBOSE_ZLIB_WRAPPERS */
