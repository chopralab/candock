# put common definitions in here
CC	= g++
#~ MPICC = mpic++
PRJCFLAGS	= -static -pthread -std=c++0x -O3
#~ PRJCFLAGS	= -static -pthread -std=c++0x -O3 -DNDEBUG
#~ PRJCFLAGS	= -static -pthread -std=c++0x -O3
#~ PRJCFLAGS	= -static -std=c++0x -O3
LD	= g++
#~ LD	= mpic++
LDFLAGS	= 
AR	= ar
ARFLAGS	= 
RANLIB	= ranlib
MFLAGS = -j8
RM	= rm
ECHO	= echo

SHELL	= /bin/sh

.SILENT :
