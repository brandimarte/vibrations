#  *****************************************************  #
#             ** Phonon Vibration Analysis **             #
#                                                         #
#                    **  Version 1  **                    #
#                                                         #
#                         IF/USP                          #
#                                                         #
#   Advisor: Prof. Dr. Alexandre Reily Rocha              #
#                                                         #
#   Author: Pedro Brandimarte Mendonca                    #
#                                                         #
#   File: Makefile                                        #
#                                                         #
#  *****************************************************  #
#  Makefile to build 'vibranal' code.                     #
#  *****************************************************  #

CFLAGS = -g -I. -Wall -ansi -O3 -axSSE4.2
MKL    = /opt/intel/composer_xe_2013.1.119/mkl/lib
LDLIBS = $(MKL)/libmkl_intel_lp64.a $(MKL)/libmkl_sequential.a $(MKL)/libmkl_core.a
#LDLIBS = -lmkl_intel -lmkl_sequential -lmkl_core

RM = /bin/rm -f
CC = icc

#  *****************************************************  #

.c.o:
	$(CC) $(CFLAGS) -c $*.c

.c:
	make $*.o
	$(CC) $(CFLAGS) -o $* $*.o $(LDLIBS) 

#  *****************************************************  #

vibranal: Check.o Utils.o Phonon.o vibranal.o
	$(CC) $(CFLAGS) -o vibranal Check.o Utils.o Phonon.o vibranal.o $(LDLIBS) 

#  *****************************************************  #

clean:
	$(RM) *~ \#~ .\#* *.o core a.out

