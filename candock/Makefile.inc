# put common definitions in here
CC	= g++
#~ MPICC = mpic++
#~ PRJCFLAGS	= -g -O0 -static -pthread -std=c++0x -DNDEBUG
#~ PRJCFLAGS	= -static -pthread -std=c++0x -O3
#~ PRJCFLAGS	= -g -O0 -static -Wl,--whole-archive -lpthread -Wl,--no-whole-archive -std=c++0x -DNDEBUG
PRJCFLAGS	= -static -Wl,--whole-archive -lpthread -Wl,--no-whole-archive -std=c++0x -O3 -DNDEBUG
#~ PRJCFLAGS	= -static -Wl,--whole-archive -lpthread -Wl,--no-whole-archive -std=c++0x -O3 -DNDEBUG -m32
#~ PRJCFLAGS	= -static -Wl,--whole-archive -lpthread -Wl,--no-whole-archive -std=c++0x -O3
LD	= g++
#~ LD	= mpic++
LDFLAGS	= 
AR	= ar
ARFLAGS	= 
RANLIB	= ranlib
MFLAGS = -j4
RM	= rm
ECHO	= echo

SHELL	= /bin/sh

.SILENT :
