include ../../candock/Makefile.inc

CFLAGS	= $(PRJCFLAGS) -I..
OBJLIBS	= libkabsch.a
OBJS	= kabsch.o
INCL= -I ../tclap-1.2.1/include/ -I ../all3/idatmtest \
-I ../jsoncpp-src-0.6.0-rc2/include -I ../jsoncpp-src-0.6.0-rc2/src \
-I ../../lib/openmm -I ..

all : $(OBJLIBS)

kabsch.o: *.cpp *.hpp ../helper/debug.hpp ../helper/error.hpp
	$(CC) $(PRJCFLAGS) $(INCL) -c $<

libkabsch.a : $(OBJS)
	$(ECHO) $(AR) $(ARFLAGS) rv libkabsch.a $?
	$(AR) $(ARFLAGS) rv libkabsch.a $?
	$(ECHO) $(RANLIB) libkabsch.a
	$(RANLIB) libkabsch.a

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
