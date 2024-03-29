VERSION=1.0.0
RELEASE=20201104

CC  := gcc

ARCH = $(shell uname -p)
ifeq ($(ARCH), x86_64)
ifeq (1, ${AVX2})
EXTRA_FLAGS=-mavx2 -mavx
else
EXTRA_FLAGS=-msse4.2
endif
else ifeq ($(ARCH), aarch64)
$(shell cd sse2neon; rm -rf [est]mmintrin.h ;ln -sf sse2neon.h emmintrin.h;ln -sf sse2neon.h smmintrin.h;ln -sf sse2neon.h tmmintrin.h)
EXTRA_FLAGS=-fsigned-char -mabi=lp64 -Isse2neon
endif


ifeq (1, ${DEBUG})
CFLAGS=-g3 -W -Wall -Wno-unused-but-set-variable -O0 -DDEBUG=1 -DVERSION="$(VERSION)" -DRELEASE="$(RELEASE)" -D_FILE_OFFSET_BITS=64 -D_GNU_SOURCE ${EXTRA_FLAGS}           
else
CFLAGS=-g3  -W -Wall -Wno-unused-but-set-variable -O4 -DVERSION="$(VERSION)" -DRELEASE="$(RELEASE)" -D_FILE_OFFSET_BITS=64 -D_GNU_SOURCE ${EXTRA_FLAGS}           
endif


GLIBS=-lm -lrt -lpthread -lz

PROGS=ccsx

all: $(PROGS)

src = $(wildcard *.c)
obj = $(src:%.c=%.o) 
    
$(PROGS): $(obj)
	$(CC) $(CFLAGS) $^ -o $@ $(GLIBS)

%.o: %.c
	$(CC) $(CFLAGS) -c $< -o $@ -I bsalign

.PHONY:clean
clean:
	rm -rf $(obj) $(PROGS) 

