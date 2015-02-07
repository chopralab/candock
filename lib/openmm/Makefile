include ../../CompBio/src/Makefile.inc
#~ include Makefile.inc

CFLAGS	= $(PRJCFLAGS) -I..
OBJLIBS	= libomm.a
OBJS	= omm.o forcefield.o moleculeinfo.o
INCL= -I ../tclap-1.2.1/include/ -I ../all3/idatmtest \
-I ../jsoncpp-src-0.6.0-rc2/include -I ../jsoncpp-src-0.6.0-rc2/src \
-I ../../lib/openmm -I ..

all : $(OBJLIBS)

omm.o: omm.cpp omm.hpp
	$(CC) $(PRJCFLAGS) $(INCL) -c $<

forcefield.o: forcefield.cpp forcefield.hpp
	$(CC) $(PRJCFLAGS) $(INCL) -c $<

moleculeinfo.o: moleculeinfo.cpp moleculeinfo.hpp
	$(CC) $(PRJCFLAGS) $(INCL) -c $<

libomm.a : $(OBJS)
	$(ECHO) $(AR) $(ARFLAGS) rv libomm.a $?
	$(AR) $(ARFLAGS) rv libomm.a $?
	$(ECHO) $(RANLIB) libomm.a
	$(RANLIB) libomm.a

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
