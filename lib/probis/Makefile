include ../../candock/Makefile.inc

OBJLIBS	= libprobis.a
OBJS	= probis.o ligands.o bsite.o eelement.o probe.o \
        colors.o sphere.o atom.o geo.o grid.o \
        residue.o ligand.o molecule.o desc.o \
        schmitt.o item.o output.o kabsch.o score.o bit.o \
        clique.o subgraph.o product.o cluster.o \
        clusterlib.o args.o const.o states.o \
        motif.o denorm.o nosql.o
INCL= -I ../tclap-1.2.1/include/ -I ../all3/idatmtest \
-I ../jsoncpp-src-0.6.0-rc2/include -I ../jsoncpp-src-0.6.0-rc2/src \
-I ../../lib/modeler -I ..

all : $(OBJLIBS)

probis.o: probis.cpp probis.hpp states.h args.h \
        const.h 
	$(CC) $(PRJCFLAGS) $(INCL) -c $<

ligands.o: ligands.cc ligands.h molecule.h \
	   ligand.h residue.h clusterlib.h const.h \
           geo.h bsite.h grid.h clusterdata.h \
           kabsch.h
	$(CC) $(PRJCFLAGS) $(INCL) -c $<

bsite.o: bsite.cc bsite.h const.h geo.h ligand.h
	$(CC) $(PRJCFLAGS) $(INCL) -c $<

motif.o: motif.cc motif.h const.h
	$(CC) $(PRJCFLAGS) $(INCL) -c $<

denorm.o: denorm.cc denorm.h geo.h kabsch.h clusterdata.h \
	  residue.h const.h
	$(CC) $(PRJCFLAGS) $(INCL) -c $<

eelement.o: eelement.cc eelement.h const.h
	$(CC) $(PRJCFLAGS) $(INCL) -c $<

probe.o: probe.cc probe.h const.h sphere.h
	$(CC) $(PRJCFLAGS) $(INCL) -c $<

colors.o: colors.cc colors.h
	$(CC) $(PRJCFLAGS) $(INCL) -c $<

sphere.o: sphere.cc sphere.h
	$(CC) $(PRJCFLAGS) $(INCL) -c $<

atom.o: atom.cc atom.h const.h eelement.h colors.h sphere.h \
	geo.h
	$(CC) $(PRJCFLAGS) $(INCL) -c $<

geo.o: geo.cc geo.h const.h
	$(CC) $(PRJCFLAGS) $(INCL) -c $<

grid.o: grid.cc grid.h desc.h const.h geo.h atom.h probe.h \
	eelement.h ligand.h
	$(CC) $(PRJCFLAGS) $(INCL) -c $<

residue.o: residue.cc residue.h const.h
	$(CC) $(PRJCFLAGS) $(INCL) -c $<

ligand.o: ligand.cc ligand.h const.h molecule.h eelement.h \
          sphere.h clusterdata.h
	$(CC) $(PRJCFLAGS) $(INCL) -c $<

molecule.o: molecule.cc molecule.h desc.h const.h atom.h \
	kabsch.h eelement.h probe.h \
        biounit.h ligands.h ligand.h residue.h geo.h grid.h \
        item.h output.h clusterdata.h motif.h
	$(CC) $(PRJCFLAGS) $(INCL) -c $<

desc.o: desc.cc desc.h sphere.h product.h const.h \
        schmitt.h atom.h eelement.h grid.h probe.h
	$(CC) $(PRJCFLAGS) $(INCL) -c $<

schmitt.o: schmitt.cc schmitt.h
	$(CC) $(PRJCFLAGS) $(INCL) -c $<

item.o: item.cc item.h vert.h product.h \
        const.h
	$(CC) $(PRJCFLAGS) $(INCL) -c $<

output.o: output.cc output.h atom.h item.h molecule.h \
          product.h desc.h score.h const.h \
          eelement.h
	$(CC) $(PRJCFLAGS) $(INCL) -c $<

kabsch.o: kabsch.cc kabsch.h geo.h item.h desc.h product.h \
          const.h atom.h
	$(CC) $(PRJCFLAGS) $(INCL) -c $<

score.o: score.cc score.h kabsch.h item.h molecule.h subgraph.h \
         desc.h probe.h product.h eelement.h atom.h const.h
	$(CC) $(PRJCFLAGS) $(INCL) -c $<

bit.o: bit.cc bit.h const.h
	$(CC) $(PRJCFLAGS) $(INCL) -c $<

clique.o: clique.cc clique.h kabsch.h const.h item.h subgraph.h \
          bit.h
	$(CC) $(PRJCFLAGS) $(INCL) -c $<

subgraph.o: subgraph.cc subgraph.h item.h vert.h molecule.h \
            const.h atom.h product.h kabsch.h output.h \
            desc.h score.h
	$(CC) $(PRJCFLAGS) $(INCL) -c $<

product.o: product.cc product.h desc.h item.h atom.h eelement.h \
           score.h const.h
	$(CC) $(PRJCFLAGS) $(INCL) -c $<

cluster.o: cluster.cc cluster.h product.h subgraph.h item.h \
           const.h
	$(CC) $(PRJCFLAGS) $(INCL) -c $<

nosql.o: nosql.cc nosql.h molecule.h residue.h clusterdata.h \
         denorm.h item.h subgraph.h product.h \
         desc.h atom.h const.h
	$(CC) $(PRJCFLAGS) $(INCL) -c $<

clusterlib.o: clusterlib.cc clusterlib.h
	$(CC) $(PRJCFLAGS) $(INCL) -c $<

args.o: args.cc args.h const.h geo.h motif.h
	$(CC) $(PRJCFLAGS) $(INCL) -c $<

const.o: const.cc const.h args.h geo.h ligands.h product.h \
         grid.h probe.h desc.h clique.h subgraph.h cluster.h \
         motif.h molecule.h output.h score.h nosql.h
	$(CC) $(PRJCFLAGS) $(INCL) -c $<

states.o: states.cc states.h molecule.h grid.h desc.h \
          score.h probe.h clique.h product.h args.h const.h \
          subgraph.h cluster.h output.h ligands.h item.h \
          motif.h nosql.h
	$(CC) $(PRJCFLAGS) $(INCL) -c $<

libprobis.a : $(OBJS)
	$(ECHO) $(AR) $(ARFLAGS) rv libprobis.a $?
	$(AR) $(ARFLAGS) rv libprobis.a $?
	$(ECHO) $(RANLIB) libprobis.a
	$(RANLIB) libprobis.a

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
