include ../../src/Makefile.inc

CFLAGS	= $(PRJCFLAGS) -I..
OBJLIBS	= libcluster.a
OBJS	= greedy.o
INCL= -I ../tclap-1.2.1/include/ -I ../all3/idatmtest \
-I ../jsoncpp-src-0.6.0-rc2/include -I ../jsoncpp-src-0.6.0-rc2/src \
-I ../../lib/modeler -I ..

all : $(OBJLIBS)

greedy.o: greedy.cpp greedy.hpp
	$(CC) $(PRJCFLAGS) $(INCL) -c $<

libcluster.a : $(OBJS)
	$(ECHO) $(AR) $(ARFLAGS) rv libcluster.a $?
	$(AR) $(ARFLAGS) rv libcluster.a $?
	$(ECHO) $(RANLIB) libcluster.a
	$(RANLIB) libcluster.a

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
