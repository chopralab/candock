include ../../candock/Makefile.inc

CFLAGS	= $(PRJCFLAGS) -I..
OBJLIBS	= liblinker.a
OBJS	= linker.o seed.o segment.o state.o poses.o partial.o
INCL= -I ../tclap-1.2.1/include/ -I ../all3/idatmtest \
-I ../jsoncpp-src-0.6.0-rc2/include -I ../jsoncpp-src-0.6.0-rc2/src \
-I ../../lib/modeler -I ..

all : $(OBJLIBS)

linker.o: linker.cpp linker.hpp state.hpp segment.hpp seed.hpp
	$(CC) $(PRJCFLAGS) $(INCL) -c $<

seed.o: seed.cpp seed.hpp segment.hpp state.hpp
	$(CC) $(PRJCFLAGS) $(INCL) -c $<

segment.o: segment.cpp segment.hpp state.hpp
	$(CC) $(PRJCFLAGS) $(INCL) -c $<

state.o: state.cpp state.hpp segment.hpp
	$(CC) $(PRJCFLAGS) $(INCL) -c $<

poses.o: poses.cpp poses.hpp segment.hpp state.hpp
	$(CC) $(PRJCFLAGS) $(INCL) -c $<

partial.o: partial.cpp partial.hpp segment.hpp state.hpp
	$(CC) $(PRJCFLAGS) $(INCL) -c $<

liblinker.a : $(OBJS)
	$(ECHO) $(AR) $(ARFLAGS) rv liblinker.a $?
	$(AR) $(ARFLAGS) rv liblinker.a $?
	$(ECHO) $(RANLIB) liblinker.a
	$(RANLIB) liblinker.a

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
