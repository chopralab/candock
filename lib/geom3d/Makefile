include ../../candock/Makefile.inc
#~ include Makefile.inc

CFLAGS	= $(PRJCFLAGS) -I..
OBJLIBS	= libgeom3d.a
OBJS	= pca.o coordinate.o geom3d.o
INCL= -I ../tclap-1.2.1/include/ -I ../all3/idatmtest \
-I ../jsoncpp-src-0.6.0-rc2/include -I ../jsoncpp-src-0.6.0-rc2/src \
-I ../../lib/openmm -I ..

all : $(OBJLIBS)

pca.o: pca.cpp pca.hpp
	$(CC) $(PRJCFLAGS) $(INCL) -c $<

coordinate.o: coordinate.cpp coordinate.hpp matrix.hpp
	$(CC) $(PRJCFLAGS) $(INCL) -c $<

geom3d.o: geom3d.cpp geom3d.hpp coordinate.hpp
	$(CC) $(PRJCFLAGS) $(INCL) -c $<

libgeom3d.a : $(OBJS)
	$(ECHO) $(AR) $(ARFLAGS) rv libgeom3d.a $?
	$(AR) $(ARFLAGS) rv libgeom3d.a $?
	$(ECHO) $(RANLIB) libgeom3d.a
	$(RANLIB) libgeom3d.a

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
