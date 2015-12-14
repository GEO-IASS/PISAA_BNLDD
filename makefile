#!/bin/bash


# --------------------------------------------------------------------------------
# 
# Copyrigtht 2014 Georgios Karagiannis
# 
# This file is part of PISAA_BNLDD.
# 
# PISAA_BNLDD is free software: you can redistribute it and/or modify
# it under the terms of the GNU General Public License as published by
# the Free Software Foundation version 2 of the License.
# 
# PISAA_BNLDD is distributed in the hope that it will be useful,
# but WITHOUT ANY WARRANTY; without even the implied warranty of
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
# GNU General Public License for more details.
# 
# You should have received a copy of the GNU General Public License
# along with PISAA_BNLDD.  If not, see <http://www.gnu.org/licenses/>.
# 
# --------------------------------------------------------------------------------

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



