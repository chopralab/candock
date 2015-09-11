include ../../candock/Makefile.inc

CFLAGS	= $(PRJCFLAGS) -I..
OBJLIBS	= libscore.a
OBJS	= score.o
INCL= -I ../tclap-1.2.1/include/ -I ../all3/idatmtest \
-I ../jsoncpp-src-0.6.0-rc2/include -I ../jsoncpp-src-0.6.0-rc2/src \
-I ../../lib/modeler -I ..

all : $(OBJLIBS)

score.o: score.cpp score.hpp
	$(CC) $(PRJCFLAGS) $(INCL) -c $<

libscore.a : $(OBJS)
	$(ECHO) $(AR) $(ARFLAGS) rv libscore.a $?
	$(AR) $(ARFLAGS) rv libscore.a $?
	$(ECHO) $(RANLIB) libscore.a
	$(RANLIB) libscore.a

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
