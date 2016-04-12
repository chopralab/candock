include ../../src/Makefile.inc

CFLAGS	= $(PRJCFLAGS) -I..
OBJLIBS	= libdocker.a
OBJS	= dock.o gpoints.o conformations.o
INCL= -I ../tclap-1.2.1/include/ -I ../all3/idatmtest \
-I ../jsoncpp-src-0.6.0-rc2/include -I ../jsoncpp-src-0.6.0-rc2/src \
-I ../../lib/modeler -I ..

all : $(OBJLIBS)

dock.o: dock.cpp dock.hpp
	$(CC) $(PRJCFLAGS) $(INCL) -c $<

gpoints.o: gpoints.cpp gpoints.hpp
	$(CC) $(PRJCFLAGS) $(INCL) -c $<

conformations.o: conformations.cpp conformations.hpp
	$(CC) $(PRJCFLAGS) $(INCL) -c $<

libdocker.a : $(OBJS)
	$(ECHO) $(AR) $(ARFLAGS) rv libdocker.a $?
	$(AR) $(ARFLAGS) rv libdocker.a $?
	$(ECHO) $(RANLIB) libdocker.a
	$(RANLIB) libdocker.a

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
