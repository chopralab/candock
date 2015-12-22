include ../../candock/Makefile.inc

CFLAGS	= $(PRJCFLAGS) -I..
OBJLIBS	= libhelper.a
OBJS	= help.o inout.o benchmark.o array2d.o path.o
INCL= -I ../tclap-1.2.1/include/ -I ../all3/idatmtest \
-I ../jsoncpp-src-0.6.0-rc2/include -I ../jsoncpp-src-0.6.0-rc2/src \
-I ../../lib/modeler -I ..

all : $(OBJLIBS)

help.o: help.cpp help.hpp
	$(CC) $(PRJCFLAGS) $(INCL) -c $<

inout.o: inout.cpp inout.hpp error.hpp
	$(CC) $(PRJCFLAGS) $(INCL) -c $<

benchmark.o: benchmark.cpp benchmark.hpp
	$(CC) $(PRJCFLAGS) $(INCL) -c $<

array2d.o: array2d.cpp array2d.hpp
	$(CC) $(PRJCFLAGS) $(INCL) -c $<

path.o: path.cpp path.hpp
	$(CC) $(PRJCFLAGS) $(INCL) -c $<

libhelper.a : $(OBJS)
	$(ECHO) $(AR) $(ARFLAGS) rv libhelper.a $?
	$(AR) $(ARFLAGS) rv libhelper.a $?
	$(ECHO) $(RANLIB) libhelper.a
	$(RANLIB) libhelper.a

#~ ../libsubsub.a : force_look
	#~ $(ECHO) looking into subsubdir : $(MAKE) $(MFLAGS)
	#~ cd subsubdir; $(MAKE) $(MFLAGS)
#~ 
clean :
	$(ECHO) cleaning up in subdir
	-$(RM) -f $(OBJS)
#	cd subsubdir; $(MAKE) $(MFLAGS) clean

force_look :
	true
