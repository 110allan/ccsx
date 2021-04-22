VERSION=1.0.0
RELEASE=20201104

CC  := gcc

ARCH = $(shell uname -p)
ifeq ($(ARCH), x86_64)
EXTRA_FLAGS=-msse4.2
else ifeq ($(ARCH), aarch64)
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
	rm -rf $(obj) $(PROGS) bsalign/{*.o,bsalign} 

