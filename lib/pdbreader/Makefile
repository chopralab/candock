include ../../candock/Makefile.inc

OBJLIBS	= libpdb.a
OBJS	= atom.o molecule.o bond.o pdbreader.o element.o internal.o idatm.o \
atomtype.o bondtype.o hydrogens.o residue.o chain.o model.o assembly.o \
molecules.o nrset.o
INCL= -I ../tclap-1.2.1/include/ -I ../all3/idatmtest \
-I ../jsoncpp-src-0.6.0-rc2/include -I ../jsoncpp-src-0.6.0-rc2/src \
-I ../../lib/modeler -I ..

all : $(OBJLIBS)

pdbreader.o: pdbreader.cpp pdbreader.hpp \
	  nrset.hpp molecules.hpp molecule.hpp \
	  assembly.hpp model.hpp residue.hpp atom.hpp
	$(CC) $(PRJCFLAGS) $(INCL) -c $<

atom.o: atom.cpp atom.hpp it.hpp \
	  element.hpp
	$(CC) $(PRJCFLAGS) $(INCL) -c $<

residue.o: residue.cpp residue.hpp it.hpp \
	  element.hpp
	$(CC) $(PRJCFLAGS) $(INCL) -c $<

chain.o: chain.cpp chain.hpp it.hpp \
	  element.hpp
	$(CC) $(PRJCFLAGS) $(INCL) -c $<

model.o: model.cpp model.hpp it.hpp \
	  element.hpp
	$(CC) $(PRJCFLAGS) $(INCL) -c $<

assembly.o: assembly.cpp assembly.hpp it.hpp \
	  element.hpp
	$(CC) $(PRJCFLAGS) $(INCL) -c $<

molecule.o: molecule.cpp molecule.hpp it.hpp \
	  element.hpp
	$(CC) $(PRJCFLAGS) $(INCL) -c $<

molecules.o: molecules.cpp molecules.hpp it.hpp \
	  element.hpp
	$(CC) $(PRJCFLAGS) $(INCL) -c $<

nrset.o: nrset.cpp nrset.hpp it.hpp \
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

hydrogens.o: hydrogens.cpp hydrogens.hpp
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
