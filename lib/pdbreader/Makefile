include ../../CompBio/src/Makefile.inc

OBJLIBS	= libpdb.a
OBJS	= molecule.o bond.o pdbreader.o element.o internal.o idatm.o \
atomtype.o bondtype.o
INCL= -I ../tclap-1.2.1/include/ -I ../all3/idatmtest \
-I ../jsoncpp-src-0.6.0-rc2/include -I ../jsoncpp-src-0.6.0-rc2/src \
-I ../../lib/openmm -I ..

all : $(OBJLIBS)

pdbreader.o: pdbreader.cpp pdbreader.hpp \
	  molecule.hpp
	$(CC) $(PRJCFLAGS) $(INCL) -c $<

molecule.o: molecule.cpp molecule.hpp it.hpp \
	  element.hpp
	$(CC) $(PRJCFLAGS) $(INCL) -c $<

element.o: element.cpp element.hpp
	$(CC) $(PRJCFLAGS) $(INCL) -c $<

internal.o: internal.cpp internal.hpp \
	  molecule.hpp
	$(CC) $(PRJCFLAGS) $(INCL) -c $<

bond.o: bond.cpp bond.hpp
	$(CC) $(PRJCFLAGS) $(INCL) -c $<

atomtype.o: atomtype.cpp atomtype.hpp
	$(CC) $(PRJCFLAGS) $(INCL) -c $<

bondtype.o: bondtype.cpp bondtype.hpp
	$(CC) $(PRJCFLAGS) $(INCL) -c $<

idatm.o: idatm.cpp idatm.hpp
	$(CC) $(PRJCFLAGS) $(INCL) -c $<

libpdb.a : $(OBJS)
	$(ECHO) $(AR) $(ARFLAGS) rv libpdb.a $?
	$(AR) $(ARFLAGS) rv libpdb.a $?
	$(ECHO) $(RANLIB) libpdb.a
	$(RANLIB) libpdb.a

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
