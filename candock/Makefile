include Makefile.inc

LIBDIR = ../lib
DIRS	= $(LIBDIR)/pdbreader $(LIBDIR)/kabsch $(LIBDIR)/helper $(LIBDIR)/linker \
$(LIBDIR)/score $(LIBDIR)/fragmenter $(LIBDIR)/graph $(LIBDIR)/kbforce $(LIBDIR)/geom3d \
$(LIBDIR)/probis $(LIBDIR)/ligands $(LIBDIR)/docker $(LIBDIR)/centro \
$(LIBDIR)/modeler $(LIBDIR)/cluster
EXE	= candock
POTENTIAL = genpot
TEST = test_minimization
TEST2 = test_fragmenting
TEST3 = test_bsite
TEST4 = test_clq
TEST5 = test_link
TEST6 = test_singlepoint
TEST7 = test_cofactor
OBJS	= candock.o common.o
OBJS_POTENTIAL	= genpot.o
OBJS_TEST	= test_minimization.o common.o
OBJS_TEST2	= test_fragmenting.o common.o
OBJS_TEST3	= test_bsite.o common.o
OBJS_TEST4	= test_clq.o common.o
OBJS_TEST5	= test_link.o common.o
OBJS_TEST6	= test_singlepoint.o common.o
OBJS_TEST7	= test_cofactor.o common.o
OBJLIBS	=  libfragmenter.a libgeom3d.a libgraph.a libhelper.a libkabsch.a \
libkbforce.a liblinker.a libmodeler.a libpdb.a libscore.a libprobis.a libligands.a \
libdocker.a libcentro.a libcluster.a
LIBS	= -L. -L$(LIBDIR)/linker -L $(LIBDIR)/pdbreader -L $(LIBDIR)/kabsch \
-L $(LIBDIR)/helper -L $(LIBDIR)/score -L $(LIBDIR)/fragmenter -L $(LIBDIR)/graph \
-L $(LIBDIR)/kbforce -L $(LIBDIR)/geom3d -L $(LIBDIR)/modeler -L $(LIBDIR)/probis \
-L $(LIBDIR)/ligands -L $(LIBDIR)/docker -L $(LIBDIR)/centro -L $(LIBDIR)/cluster \
-lligands -ldocker -lscore -lpdb -llinker -lfragmenter -lgraph -lkabsch -lcentro -lcluster -llinker \
-lmodeler -lhelper -lOpenMM_static -ldl -lc -lkbforce -lprobis -lgeom3d \
-lboost_regex -lboost_filesystem -lboost_system -lpthread -lboost_date_time \
-lgsl -lgslcblas -lm
INCL = -I $(LIBDIR)/tclap-1.2.1/include -I $(LIBDIR)/jsoncpp-src-0.6.0-rc2/include \
-I $(LIBDIR)/jsoncpp-src-0.6.0-rc2/src -I $(LIBDIR) -I $(LIBDIR)/modeler

all : $(EXE)
pote : $(POTENTIAL)
test : $(TEST)
test2 : $(TEST2)
test3 : $(TEST3)
test4 : $(TEST4)
test5 : $(TEST5)
test6 : $(TEST6)
test7 : $(TEST7)

$(POTENTIAL) : genpot.o $(OBJLIBS)
	$(ECHO) $(LD) $(PRJCFLAGS) -o $(POTENTIAL) $(OBJS_POTENTIAL) $(LIBS)
	$(LD) $(PRJCFLAGS) -o $(POTENTIAL) $(OBJS_POTENTIAL) $(LIBS)

$(TEST) : test_minimization.o common.o $(OBJLIBS)
	$(ECHO) $(LD) $(PRJCFLAGS) -o $(TEST) $(OBJS_TEST) $(LIBS)
	$(LD) $(PRJCFLAGS) -o $(TEST) $(OBJS_TEST) $(LIBS)

$(TEST2) : test_fragmenting.o common.o $(OBJLIBS)
	$(ECHO) $(LD) $(PRJCFLAGS) -o $(TEST2) $(OBJS_TEST2) $(LIBS)
	$(LD) $(PRJCFLAGS) -o $(TEST2) $(OBJS_TEST2) $(LIBS)

$(TEST3) : test_bsite.o common.o $(OBJLIBS)
	$(ECHO) $(LD) $(PRJCFLAGS) -o $(TEST3) $(OBJS_TEST3) $(LIBS)
	$(LD) $(PRJCFLAGS) -o $(TEST3) $(OBJS_TEST3) $(LIBS)

$(TEST4) : test_clq.o common.o $(OBJLIBS)
	$(ECHO) $(LD) $(PRJCFLAGS) -o $(TEST4) $(OBJS_TEST4) $(LIBS)
	$(LD) $(PRJCFLAGS) -o $(TEST4) $(OBJS_TEST4) $(LIBS)

$(TEST5) : test_link.o common.o $(OBJLIBS)
	$(ECHO) $(LD) $(PRJCFLAGS) -o $(TEST5) $(OBJS_TEST5) $(LIBS)
	$(LD) $(PRJCFLAGS) -o $(TEST5) $(OBJS_TEST5) $(LIBS)

$(TEST6) : test_singlepoint.o common.o $(OBJLIBS)
	$(ECHO) $(LD) $(PRJCFLAGS) -o $(TEST6) $(OBJS_TEST6) $(LIBS)
	$(LD) $(PRJCFLAGS) -o $(TEST6) $(OBJS_TEST6) $(LIBS)

$(TEST7) : test_cofactor.o common.o $(OBJLIBS)
	$(ECHO) $(LD) $(PRJCFLAGS) -o $(TEST7) $(OBJS_TEST7) $(LIBS)
	$(LD) $(PRJCFLAGS) -o $(TEST7) $(OBJS_TEST7) $(LIBS)

$(EXE) : candock.o common.o $(OBJLIBS)
	$(ECHO) $(LD) $(PRJCFLAGS) -o $(EXE) $(OBJS) $(LIBS)
	$(LD) $(PRJCFLAGS) -o $(EXE) $(OBJS) $(LIBS)

genpot.o: genpot.cpp opts_genpot.hpp
	$(CC) $(PRJCFLAGS) $(INCL) -c $<

test_cofactor.o: test_cofactor.cpp common.cpp common.hpp \
	   opts_candock.hpp
	$(CC) $(PRJCFLAGS) $(INCL) -c $<

test_singlepoint.o: test_singlepoint.cpp common.cpp common.hpp \
	   opts_candock.hpp
	$(CC) $(PRJCFLAGS) $(INCL) -c $<

test_link.o: test_link.cpp common.cpp common.hpp \
	   opts_candock.hpp
	$(CC) $(PRJCFLAGS) $(INCL) -c $<

test_clq.o: test_clq.cpp common.cpp common.hpp \
	   opts_candock.hpp
	$(CC) $(PRJCFLAGS) $(INCL) -c $<

test_bsite.o: test_bsite.cpp common.cpp common.hpp \
	   opts_candock.hpp
	$(CC) $(PRJCFLAGS) $(INCL) -c $<

test_minimization.o: test_minimization.cpp common.cpp common.hpp \
	   opts_candock.hpp
	$(CC) $(PRJCFLAGS) $(INCL) -c $<

test_fragmenting.o: test_fragmenting.cpp common.cpp common.hpp \
	   opts_candock.hpp
	$(CC) $(PRJCFLAGS) $(INCL) -c $<

candock.o: candock.cpp common.cpp common.hpp \
	   opts_candock.hpp
	$(CC) $(PRJCFLAGS) $(INCL) -c $<

common.o: common.cpp common.hpp
	$(CC) $(PRJCFLAGS) $(INCL) -c $<

libcluster.a : force_look
	$(ECHO) looking into subdir : $(MAKE) $(MFLAGS)
	cd $(LIBDIR)/cluster; $(MAKE) $(MFLAGS)

libcentro.a : force_look
	$(ECHO) looking into subdir : $(MAKE) $(MFLAGS)
	cd $(LIBDIR)/centro; $(MAKE) $(MFLAGS)

libdocker.a : force_look
	$(ECHO) looking into subdir : $(MAKE) $(MFLAGS)
	cd $(LIBDIR)/docker; $(MAKE) $(MFLAGS)

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

libmodeler.a: force_look
	$(ECHO) looking into subdir : $(MAKE) $(MFLAGS)
	cd $(LIBDIR)/modeler; $(MAKE) $(MFLAGS)

clean :
	$(ECHO) cleaning up in .
	-$(RM) -f $(EXE) $(TEST) $(TEST2) $(TEST3) $(TEST4) $(TEST5) $(OBJS) $(OBJS_TEST) $(OBJS_TEST2) $(OBJS_TEST3) $(OBJS_TEST4) $(OBJS_TEST5) $(OBJLIBS)
	-for d in $(DIRS); do (cd $$d; $(MAKE) clean ); done

force_look :
	true
