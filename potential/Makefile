include Makefile.inc

LIBDIR = ../lib
DIRS	= $(LIBDIR)/pdbreader $(LIBDIR)/helper $(LIBDIR)/score $(LIBDIR)/geom3d \
$(LIBDIR)/fragmenter $(LIBDIR)/graph
EXE	= genpot
OBJS	= genpot.o
OBJLIBS	=  libfragmenter.a libgraph.a libpdb.a libgeom3d.a libhelper.a libscore.a
LIBS	= -L. -L $(LIBDIR)/pdbreader -L $(LIBDIR)/helper -L $(LIBDIR)/score \
-L $(LIBDIR)/fragmenter -L $(LIBDIR)/graph -L $(LIBDIR)/geom3d \
-lscore -lhelper -ldl -lc -lpdb -lgeom3d -lfragmenter -lgraph \
-lboost_regex -lboost_filesystem -lboost_system -lpthread -lboost_date_time \
-lgsl -lgslcblas -lm
INCL = -I $(LIBDIR)/tclap-1.2.1/include -I $(LIBDIR)/jsoncpp-src-0.6.0-rc2/include \
-I $(LIBDIR)/jsoncpp-src-0.6.0-rc2/src -I $(LIBDIR) -I $(LIBDIR)/modeler

all : $(EXE)

$(EXE) : genpot.o $(OBJLIBS)
	$(ECHO) $(LD) $(PRJCFLAGS) -o $(EXE) $(OBJS) $(LIBS)
	$(LD) $(PRJCFLAGS) -o $(EXE) $(OBJS) $(LIBS)

genpot.o: genpot.cpp \
	   opts_genpot.hpp
	$(CC) $(PRJCFLAGS) $(INCL) -c $<

libfragmenter.a : force_look
	$(ECHO) looking into subdir : $(MAKE) $(MFLAGS)
	cd $(LIBDIR)/fragmenter; $(MAKE) $(MFLAGS)

libgraph.a : force_look
	$(ECHO) looking into subdir : $(MAKE) $(MFLAGS)
	cd $(LIBDIR)/graph; $(MAKE) $(MFLAGS)

libpdb.a : force_look
	$(ECHO) looking into subdir : $(MAKE) $(MFLAGS)
	cd $(LIBDIR)/pdbreader; $(MAKE) $(MFLAGS)

libhelper.a : force_look
	$(ECHO) looking into subdir : $(MAKE) $(MFLAGS)
	cd $(LIBDIR)/helper; $(MAKE) $(MFLAGS)

libscore.a : force_look
	$(ECHO) looking into subdir : $(MAKE) $(MFLAGS)
	cd $(LIBDIR)/score; $(MAKE) $(MFLAGS)

libgeom3d.a : force_look
	$(ECHO) looking into subdir : $(MAKE) $(MFLAGS)
	cd $(LIBDIR)/geom3d; $(MAKE) $(MFLAGS)

clean :
	$(ECHO) cleaning up in .
	-$(RM) -f $(EXE) $(OBJS) $(OBJLIBS)
	-for d in $(DIRS); do (cd $$d; $(MAKE) clean ); done

force_look :
	true
