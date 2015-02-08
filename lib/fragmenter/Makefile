include ../../candock/Makefile.inc

CFLAGS	= $(PRJCFLAGS) -I..
OBJLIBS	= libfragmenter.a
OBJS	= fragmenter.o unique.o
INCL= -I ../tclap-1.2.1/include/ -I ../all3/idatmtest \
-I ../jsoncpp-src-0.6.0-rc2/include -I ../jsoncpp-src-0.6.0-rc2/src \
-I ../../lib/openmm -I ..

all : $(OBJLIBS)

fragmenter.o: fragmenter.cpp fragmenter.hpp unique.hpp
	$(CC) $(PRJCFLAGS) $(INCL) -c $<

unique.o: unique.cpp unique.hpp
	$(CC) $(PRJCFLAGS) $(INCL) -c $<

libfragmenter.a : $(OBJS)
	$(ECHO) $(AR) $(ARFLAGS) rv libfragmenter.a $?
	$(AR) $(ARFLAGS) rv libfragmenter.a $?
	$(ECHO) $(RANLIB) libfragmenter.a
	$(RANLIB) libfragmenter.a

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
