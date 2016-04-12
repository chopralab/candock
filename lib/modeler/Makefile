include ../../src/Makefile.inc
#~ include Makefile.inc

CFLAGS	= $(PRJCFLAGS) -I..
OBJLIBS	= libmodeler.a
OBJS	= modeler.o forcefield.o systemtopology.o topology.o
INCL= -I ../tclap-1.2.1/include/ -I ../all3/idatmtest \
-I ../jsoncpp-src-0.6.0-rc2/include -I ../jsoncpp-src-0.6.0-rc2/src \
-I ../../lib/modeler -I ..

all : $(OBJLIBS)

modeler.o: modeler.cpp modeler.hpp
	$(CC) $(PRJCFLAGS) $(INCL) -c $<

forcefield.o: forcefield.cpp forcefield.hpp
	$(CC) $(PRJCFLAGS) $(INCL) -c $<

systemtopology.o: systemtopology.cpp systemtopology.hpp
	$(CC) $(PRJCFLAGS) $(INCL) -c $<

topology.o: topology.cpp topology.hpp
	$(CC) $(PRJCFLAGS) $(INCL) -c $<

libmodeler.a : $(OBJS)
	$(ECHO) $(AR) $(ARFLAGS) rv libmodeler.a $?
	$(AR) $(ARFLAGS) rv libmodeler.a $?
	$(ECHO) $(RANLIB) libmodeler.a
	$(RANLIB) libmodeler.a

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
