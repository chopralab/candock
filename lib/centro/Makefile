include ../../candock/Makefile.inc

CFLAGS	= $(PRJCFLAGS) -I..
OBJLIBS	= libcentro.a
OBJS	= centroids.o
INCL= -I ../tclap-1.2.1/include/ -I ../all3/idatmtest \
-I ../jsoncpp-src-0.6.0-rc2/include -I ../jsoncpp-src-0.6.0-rc2/src \
-I ../../lib/openmm -I ..

all : $(OBJLIBS)

centroids.o: centroids.cpp centroids.hpp
	$(CC) $(PRJCFLAGS) $(INCL) -c $<

libcentro.a : $(OBJS)
	$(ECHO) $(AR) $(ARFLAGS) rv libcentro.a $?
	$(AR) $(ARFLAGS) rv libcentro.a $?
	$(ECHO) $(RANLIB) libcentro.a
	$(RANLIB) libcentro.a

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
