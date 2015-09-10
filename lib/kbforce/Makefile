include ../../candock/Makefile.inc

CFLAGS	= $(PRJCFLAGS) -I..
OBJLIBS	= libkbforce.a
OBJS	= KBForce.o KBForceImpl.o ReferenceKBKernelFactory.o ReferenceKBKernels.o
INCL= -I ../tclap-1.2.1/include/ -I ../all3/idatmtest -I ../jsoncpp-src-0.6.0-rc2/include \
-I ../jsoncpp-src-0.6.0-rc2/src -I .. -I ../modeler -I openmmapi/include -I platforms/reference/src \
-I platforms/reference/include -I ..

all : $(OBJLIBS)

KBForce.o: openmmapi/src/KBForce.cpp openmmapi/include/KBForce.h
	$(CC) $(PRJCFLAGS) $(INCL) -c $<

KBForceImpl.o: openmmapi/src/KBForceImpl.cpp openmmapi/include/KBKernels.h openmmapi/include/KBForce.h
	$(CC) $(PRJCFLAGS) $(INCL) -c $<

ReferenceKBKernelFactory.o: platforms/reference/src/ReferenceKBKernelFactory.cpp platforms/reference/include/ReferenceKBKernelFactory.h platforms/reference/src/ReferenceKBKernels.h
	$(CC) $(PRJCFLAGS) $(INCL) -c $<

ReferenceKBKernels.o: platforms/reference/src/ReferenceKBKernels.cpp platforms/reference/src/ReferenceKBKernels.h openmmapi/include/KBForce.h
	$(CC) $(PRJCFLAGS) $(INCL) -c $<

libkbforce.a : $(OBJS)
	$(ECHO) $(AR) $(ARFLAGS) rv libkbforce.a $?
	$(AR) $(ARFLAGS) rv libkbforce.a $?
	$(ECHO) $(RANLIB) libkbforce.a
	$(RANLIB) libkbforce.a

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
