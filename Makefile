######################################
# Makefile written by Daniel Baker   #
#          dnb@cs.jhu.edu            #
######################################


CXXSTD=c++11
CSTD=gnu99
CC=g++
GIT_VERSION := $(shell git describe --abbrev=4 --always --dirty)
FLAGS= -Wno-unused-result -Wno-unused-function -Wunreachable-code -Wall -DMASKRIPPER_VERSION=\"$(GIT_VERSION)\" -std=$(CXXSTD) -pedantic
LD= -lm -lz -lpthread
INCLUDE= -Ihtslib -I.
LIB=
INSTALL=/usr/bin/install -c

prefix = /usr/local
bindir = $(prefix)/bin
binprefix =

THREADS=12

OPT_FLAGS = -finline-functions -O3 -DNDEBUG -flto -fivopts
DB_FLAGS = -fno-inline

DLIB_SRC = dlib/cstr_util.c dlib/math_util.c dlib/vcf_util.c dlib/io_util.c dlib/bam_util.c dlib/nix_util.c \
		   dlib/bed_util.c dlib/misc_util.c

SOURCES = src/maskripper.cpp $(DLIB_SRC)

D_OBJS = $(SOURCES:.c=.dbo)
OBJS = $(SOURCES:.c=.o)
DLIB_OBJS = $(DLIB_SRC:.c=.o)

.PHONY: all clean install mostlyclean update_dlib

all: update_dlib libhts.a maskripper maskripper_db

install: all
	$(INSTALL) maskripper $(bindir)/$(binprefix)maskripper
	$(INSTALL) maskripper_db $(bindir)/$(binprefix)maskripper_db

%.o: %.cpp
	$(CC) -c $(FLAGS) $(INCLUDE) $(LIB) $(LD) $(OPT_FLAGS) $< -o $@

%.dbo: %.cpp
	$(CC) -c $(FLAGS) $(INCLUDE) $(LIB) $(LD) $(DB_FLAGS) $< -o $@


libhts.a:
	+cd htslib && echo "/* Empty config.h */" >> config.h && make -j $(THREADS) && cp libhts.a ../
maskripper_db: $(D_OBJS) libhts.a
	$(CC) $(FLAGS) $(INCLUDE) $(LIB) $(LD) $(DB_FLAGS) $(D_OBJS) libhts.a -o maskripper_db
maskripper: $(OBJS) libhts.a
	$(CC) $(FLAGS) $(INCLUDE) $(LIB) $(LD) $(OPT_FLAGS) $(OBJS) libhts.a -o maskripper

clean: mostlyclean
		cd htslib && make clean && cd ..

update_dlib:
	cd dlib && git checkout master && git pull origin master && cd ..

mostlyclean:
	rm -f *.*o && rm -f maskripper* && rm -f src/*.*o && rm -f dlib/*.*o
