include ../../CompBio/src/Makefile.inc
#~ include Makefile.inc

#OBJLIBS	= libgenclus.a libgenlig.a libgenbio.a
OBJLIBS	= libligands.a
OBJS	= genclus.o genlig.o genbio.o jsonreader.o nosqlreader.o common.o
INCL= -I ../tclap-1.2.1/include/ -I ../all3/idatmtest \
-I ../jsoncpp-src-0.6.0-rc2/include -I ../jsoncpp-src-0.6.0-rc2/src \
-I ../../lib/openmm -I ..

all : $(OBJLIBS)

genbio.o: genbio.cpp genbio.hpp jsonreader.hpp \
	  nosqlreader.hpp opts_genbio.hpp common.hpp
	$(CC) $(PRJCFLAGS) $(INCL) -c $<

genclus.o: genclus.cpp genclus.hpp jsonreader.hpp \
	  nosqlreader.hpp opts_genclus.hpp \
	  common.hpp
	$(CC) $(PRJCFLAGS) $(INCL) -c $<

genlig.o: genlig.cpp genlig.hpp jsonreader.hpp \
	  nosqlreader.hpp opts_genlig.hpp common.hpp
	$(CC) $(PRJCFLAGS) $(INCL) -c $<

jsonreader.o: jsonreader.cpp jsonreader.hpp
	$(CC) $(PRJCFLAGS) $(INCL) -c $<

nosqlreader.o: nosqlreader.cpp nosqlreader.hpp
	$(CC) $(PRJCFLAGS) $(INCL) -c $<

common.o: common.cpp common.hpp
	$(CC) $(PRJCFLAGS) $(INCL) -c $<

libligands.a : $(OBJS)
	$(ECHO) $(AR) $(ARFLAGS) rv libligands.a $?
	$(AR) $(ARFLAGS) rv libligands.a $?
	$(ECHO) $(RANLIB) libligands.a
	$(RANLIB) libligands.a

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
