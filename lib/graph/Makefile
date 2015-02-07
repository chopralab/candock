include ../../CompBio/src/Makefile.inc

CFLAGS	= $(PRJCFLAGS) -I..
OBJLIBS	= libgraph.a
OBJS	= mcqd.o mnts.o
INCL= -I ../tclap-1.2.1/include/ -I ../all3/idatmtest \
-I ../jsoncpp-src-0.6.0-rc2/include -I ../jsoncpp-src-0.6.0-rc2/src \
-I ../../lib/openmm -I ..

all : $(OBJLIBS)

mnts.o: mnts.cpp mnts.hpp
	$(CC) $(PRJCFLAGS) $(INCL) -c $<

mcqd.o: mcqd.cpp mcqd.hpp
	$(CC) $(PRJCFLAGS) $(INCL) -c $<

libgraph.a : $(OBJS)
	$(ECHO) $(AR) $(ARFLAGS) rv libgraph.a $?
	$(AR) $(ARFLAGS) rv libgraph.a $?
	$(ECHO) $(RANLIB) libgraph.a
	$(RANLIB) libgraph.a

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
