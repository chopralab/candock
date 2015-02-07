include Makefile.inc

DIRS	= ../../lib/pdbreader ../../lib/kabsch ../../lib/helper ../../lib/linker \
../../lib/score ../../lib/fragmenter ../../lib/graph ../../lib/kbforce ../../lib/geom3d \
../../lib/openmm ../../lib/probis ../../lib/ligands
TEST = test_minimization
EXE	= candock
OBJS_TEST	= test_minimization.o common.o
OBJS	= candock.o common.o
OBJLIBS	=  libfragmenter.a libgeom3d.a libgraph.a libhelper.a libkabsch.a \
libkbforce.a liblinker.a libomm.a libpdb.a libscore.a libprobis.a libligands.a
LIBS	= -L. -L../../lib/linker -L ../../lib/pdbreader -L ../../lib/kabsch \
-L ../../lib/helper -L ../../lib/score -L ../../lib/fragmenter -L ../../lib/graph \
-L ../../lib/kbforce -L ../../lib/geom3d -L ../../lib/openmm -L ../../lib/probis \
-L ../../lib/ligands \
-llinker -lpdb -lscore -lfragmenter -lhelper -lgraph -lkabsch \
-lgeom3d -lomm -lOpenMM_static -ldl -lc -lkbforce -lprobis -lligands -lpdb \
-lboost_regex -lboost_filesystem -lboost_system -lpthread -lboost_date_time \
-lgsl -lgslcblas -lm
INCL = -I ../../lib/tclap-1.2.1/include -I ../../lib/jsoncpp-src-0.6.0-rc2/include \
-I ../../lib/jsoncpp-src-0.6.0-rc2/src -I ../../lib -I ../../lib/openmm
test : $(TEST)
all : $(EXE)

$(TEST) : test_minimization.o common.o $(OBJLIBS)
	$(ECHO) $(LD) $(PRJCFLAGS) -o $(TEST) $(OBJS_TEST) $(LIBS)
	$(LD) $(PRJCFLAGS) -o $(TEST) $(OBJS_TEST) $(LIBS)

$(EXE) : candock.o common.o $(OBJLIBS)
	$(ECHO) $(LD) $(PRJCFLAGS) -o $(EXE) $(OBJS) $(LIBS)
	$(LD) $(PRJCFLAGS) -o $(EXE) $(OBJS) $(LIBS)

test_minimization.o: test_minimization.cpp common.cpp common.hpp \
	   opts_candock.hpp
	$(CC) $(PRJCFLAGS) $(INCL) -c $<

candock.o: candock.cpp common.cpp common.hpp \
	   opts_candock.hpp
	$(CC) $(PRJCFLAGS) $(INCL) -c $<

common.o: common.cpp common.hpp
	$(CC) $(PRJCFLAGS) $(INCL) -c $<

libligands.a : force_look
	$(ECHO) looking into subdir : $(MAKE) $(MFLAGS)
	cd ../../lib/ligands; $(MAKE) $(MFLAGS)

libprobis.a : force_look
	$(ECHO) looking into subdir : $(MAKE) $(MFLAGS)
	cd ../../lib/probis; $(MAKE) $(MFLAGS)

libpdb.a : force_look
	$(ECHO) looking into subdir : $(MAKE) $(MFLAGS)
	cd ../../lib/pdbreader; $(MAKE) $(MFLAGS)

libkabsch.a : force_look
	$(ECHO) looking into subdir : $(MAKE) $(MFLAGS)
	cd ../../lib/kabsch; $(MAKE) $(MFLAGS)

libhelper.a : force_look
	$(ECHO) looking into subdir : $(MAKE) $(MFLAGS)
	cd ../../lib/helper; $(MAKE) $(MFLAGS)

liblinker.a : force_look
	$(ECHO) looking into subdir : $(MAKE) $(MFLAGS)
	cd ../../lib/linker; $(MAKE) $(MFLAGS)

libscore.a : force_look
	$(ECHO) looking into subdir : $(MAKE) $(MFLAGS)
	cd ../../lib/score; $(MAKE) $(MFLAGS)

libfragmenter.a : force_look
	$(ECHO) looking into subdir : $(MAKE) $(MFLAGS)
	cd ../../lib/fragmenter; $(MAKE) $(MFLAGS)

libgraph.a : force_look
	$(ECHO) looking into subdir : $(MAKE) $(MFLAGS)
	cd ../../lib/graph; $(MAKE) $(MFLAGS)

libkbforce.a : force_look
	$(ECHO) looking into subdir : $(MAKE) $(MFLAGS)
	cd ../../lib/kbforce; $(MAKE) $(MFLAGS)

libgeom3d.a : force_look
	$(ECHO) looking into subdir : $(MAKE) $(MFLAGS)
	cd ../../lib/geom3d; $(MAKE) $(MFLAGS)

libomm.a: force_look
	$(ECHO) looking into subdir : $(MAKE) $(MFLAGS)
	cd ../../lib/openmm; $(MAKE) $(MFLAGS)

clean :
	$(ECHO) cleaning up in .
	-$(RM) -f $(EXE) $(TEST) $(OBJS) $(OBJS_TEST) $(OBJLIBS)
	-for d in $(DIRS); do (cd $$d; $(MAKE) clean ); done

force_look :
	true
