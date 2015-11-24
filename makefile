#!/bin/bash

# Georgios Karagiannis
#
# Postdoctoral research associate
# Department of Mathematics, Purdue University
# 150 N. University Street
# West Lafayette, IN 47907-2067, USA
#
# Telephone: +1 765 496-1007
#
# Email: gkaragia@purdue.edu
#
# Contact email: georgios.stats@gmail.com
#
# Georgios Karagiannis © 2014  

# COMPILERS
				
CC=icc
CFLAGS=-O2
LDFLAGS=
CPPFLAGS=

FUN=cost_BNDV.c

SOURCES=pisaa.c \
	Crossover_operations.c \
	Mutation_operations.c \
	Self_adjastment_procedure.c \
	$(FUN) \
	permutrng.c \
	uniformrng.c \
	nrutil.c 
	
OBJECTS=$(SOURCES:.c=.o)

EXECUTABLE=exe

# BUILD

build: $(SOURCES) $(EXECUTABLE)
	
$(EXECUTABLE): $(OBJECTS) 
	$(CC) $(CPPFLAGS) $(CFLAGS) $(OBJECTS) $(LDFLAGS) -o $@

.c.o:
	$(CC) $(CPPFLAGS) $(CFLAGS) -c $< -o $@

# CLEAR

clean:
	rm -rf *o exe

# DETAILS

details:
	@echo CC       : $(CC)
	@echo CFLAGS   : $(CFLAGS)
	@echo LDFLAGS  : $(LDFLAGS)
	@echo FUN      : $(FUN)
	@echo CPPFLAGS : $(CPPFLAGS)

# RUN

run:
	./exe

test_run:
	./exe -ID 1 \
			-Data ./SPECT.dat \
			-Niter 1000000 \
			-Npop 5 \
			-Nsam 100 \
			-Gwarm 50000 \
			-Ghigh 1.0 \
			-Gpow 0.55 \
			-Hlow 2000.0 \
			-Hhigh 3999.0 \
			-Hsize 2000 \
			-Hzeta 0.0 \
			-Hconst 1.0 \
			-Twarm 500000 \
			-Thigh 1.09 \
			-Tlow 0.01 \
			-Tpow 0.5



