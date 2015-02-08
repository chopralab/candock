include Makefile.inc

LIBDIR = ../lib
DIRS	= $(LIBDIR)/pdbreader $(LIBDIR)/kabsch $(LIBDIR)/helper $(LIBDIR)/linker \
$(LIBDIR)/score $(LIBDIR)/fragmenter $(LIBDIR)/graph $(LIBDIR)/kbforce $(LIBDIR)/geom3d \
$(LIBDIR)/openmm $(LIBDIR)/probis $(LIBDIR)/ligands
EXE	= candock
TEST = test_minimization
OBJS	= candock.o common.o
OBJS_TEST	= test_minimization.o common.o
OBJLIBS	=  libfragmenter.a libgeom3d.a libgraph.a libhelper.a libkabsch.a \
libkbforce.a liblinker.a libomm.a libpdb.a libscore.a libprobis.a libligands.a
LIBS	= -L. -L$(LIBDIR)/linker -L $(LIBDIR)/pdbreader -L $(LIBDIR)/kabsch \
-L $(LIBDIR)/helper -L $(LIBDIR)/score -L $(LIBDIR)/fragmenter -L $(LIBDIR)/graph \
-L $(LIBDIR)/kbforce -L $(LIBDIR)/geom3d -L $(LIBDIR)/openmm -L $(LIBDIR)/probis \
-L $(LIBDIR)/ligands \
-llinker -lpdb -lscore -lfragmenter -lhelper -lgraph -lkabsch \
-lgeom3d -lomm -lOpenMM_static -ldl -lc -lkbforce -lprobis -lligands -lpdb \
-lboost_regex -lboost_filesystem -lboost_system -lpthread -lboost_date_time \
-lgsl -lgslcblas -lm
INCL = -I $(LIBDIR)/tclap-1.2.1/include -I $(LIBDIR)/jsoncpp-src-0.6.0-rc2/include \
-I $(LIBDIR)/jsoncpp-src-0.6.0-rc2/src -I $(LIBDIR) -I $(LIBDIR)/openmm

all : $(EXE)
test : $(TEST)

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
	cd $(LIBDIR)/ligands; $(MAKE) $(MFLAGS)

libprobis.a : force_look
	$(ECHO) looking into subdir : $(MAKE) $(MFLAGS)
	cd $(LIBDIR)/probis; $(MAKE) $(MFLAGS)

libpdb.a : force_look
	$(ECHO) looking into subdir : $(MAKE) $(MFLAGS)
	cd $(LIBDIR)/pdbreader; $(MAKE) $(MFLAGS)

libkabsch.a : force_look
	$(ECHO) looking into subdir : $(MAKE) $(MFLAGS)
	cd $(LIBDIR)/kabsch; $(MAKE) $(MFLAGS)

libhelper.a : force_look
	$(ECHO) looking into subdir : $(MAKE) $(MFLAGS)
	cd $(LIBDIR)/helper; $(MAKE) $(MFLAGS)

liblinker.a : force_look
	$(ECHO) looking into subdir : $(MAKE) $(MFLAGS)
	cd $(LIBDIR)/linker; $(MAKE) $(MFLAGS)

libscore.a : force_look
	$(ECHO) looking into subdir : $(MAKE) $(MFLAGS)
	cd $(LIBDIR)/score; $(MAKE) $(MFLAGS)

libfragmenter.a : force_look
	$(ECHO) looking into subdir : $(MAKE) $(MFLAGS)
	cd $(LIBDIR)/fragmenter; $(MAKE) $(MFLAGS)

libgraph.a : force_look
	$(ECHO) looking into subdir : $(MAKE) $(MFLAGS)
	cd $(LIBDIR)/graph; $(MAKE) $(MFLAGS)

libkbforce.a : force_look
	$(ECHO) looking into subdir : $(MAKE) $(MFLAGS)
	cd $(LIBDIR)/kbforce; $(MAKE) $(MFLAGS)

libgeom3d.a : force_look
	$(ECHO) looking into subdir : $(MAKE) $(MFLAGS)
	cd $(LIBDIR)/geom3d; $(MAKE) $(MFLAGS)

libomm.a: force_look
	$(ECHO) looking into subdir : $(MAKE) $(MFLAGS)
	cd $(LIBDIR)/openmm; $(MAKE) $(MFLAGS)

clean :
	$(ECHO) cleaning up in .
	-$(RM) -f $(EXE) $(TEST) $(OBJS) $(OBJS_TEST) $(OBJLIBS)
	-for d in $(DIRS); do (cd $$d; $(MAKE) clean ); done

force_look :
	true
