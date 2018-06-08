#ifndef RENAMERULES_H
#define RENAMERULES_H

#include "candock/helper/smiles.hpp"

namespace help {

	const rename_rules special { // idatm rules for complicated groups
		//~ {{{"^O#1#1","^P|^S#2",""},{"#2","^O#3",""},{"#2","^O|^N|^S#4",""}}, {{"1:idatm=O3-"}}}, // phosphate, sulfate, N-oxide...
		{{{"^O#1#1","^Pox|^Pac|^Son|^Sac|^Sxd#2",""}}, {{"1:idatm=O3-"}}}, // resonance equivalent terminal oxygen on tetrahedral center (phosphate, sulfate, sulfone...)
		{{{"^O#1#1","N3\\+|N2\\+#2",""}}, {{"1:idatm=O3-"}}}, // amine-N-oxide, pyridine-N-oxide
		{{{"^S#1#1","Pox|Pac|P3\\+#2",""}}, {{"1:idatm=S3-"}}}, // terminal sulfur on tetrahedral center (thiophosphate)
		//~ {{{"^S#1#2","Car#2",""},{"#1","Car#3",""}}, {{"1:idatm=Sar"}}}, // aromatic sulfur
		{{{"^S#1#2#ag",".*#2",""}}, {{"1:idatm=Sar"}}}, // aromatic sulfur
		//~ {{{"^O#1#2","Car#2",""},{"#1","Car#3",""}}, {{"1:idatm=Oar+"}}}, // aromatic oxygen, formally positive (pyrylium)
		//~ {{{"^O#1#2#1ag,ag6","Car#2",""},{"#1","Car#3",""}}, {{"1:idatm=Oar+"}}}, // aromatic oxygen, formally positive (pyrylium)
		{{{"Npl#1#2#ag","Npl#2#3#ag",""},{"#1","Car#3#3",""}}, {{"1:idatm=N2"}}}, // change aromatic Npl (bonded to Npl and Car) to N2
		//~ {{{"Npl#1##1ag,ag6",".*#2",""}}, {{"1:idatm=N2"}}}, // change Npl in 6-membered aromatic ring to N2
		{{{"Npl#1#2#1ag,ag6",".*#2",""}}, {{"1:idatm=N2"}}}, // change 2-substituted Npl in 6-membered aromatic ring to N2
		//~ {{{"Npl#1#3#1ag,ag6",".*#2",""}}, {{"1:idatm=N2+"}}}, // change 3-substituted Npl in 6-membered aromatic ring to N2+
		//~ {{{"Npl#1#3,1H#1ag,ag6",".*#2",""}}, {{"1:idatm=N2"}}}, // change 3-substituted (one bondee is hydrogen) Npl in 6-membered aromatic ring to N2
		{{{"Npl#1#3,1H#1ag,ag5","Car#2##1ag,ag5",""},{"#2","Npl#3#3,0H#1ag,ag5",""}}, {{"1:idatm=N2"}}}, // change the first Npl in 5-ring Npl(-H)-Car-Npl(-C,-C,-C) to N2
		{{{"Npl#1#3,1H#1ag,ag5","^S|^Oar#2#2#1ag,ag5",""}}, {{"1:idatm=N2"}}}, // change 3-substituted Npl bound to S or O in 5-ring to N2
		{{{"Npl#1#2#1ag,ag5","^S|^Oar#2#2#1ag,ag5",""}}, {{"1:idatm=N2"}}}, // change 2-substituted Npl bound to S or O in 5-ring to N2
		{{{"N2#1#2","C2#2#3",""},{"#2","O3#3#1",""}}, {{"1:idatm=Npl"},{"3:idatm=O2"}}}, // correct peptide bond
		//~ {{{"N1$#1#2",".*#2",""}}, {{"1:idatm=N1+"}}}, // change 2-substituted N1 to N1+ (e.g. azide)
		{{{"N1$#1#2","^N#2",""},{"#1","^N#3",""}}, {{"1:idatm=N1+"}}}, // change 2-substituted N1 to N1+ (e.g. azide)
		//~ {{{"Npl#1#1",".*#2",""}}, {{"1:idatm=N1"}}}, // change 1-substituted Npl (terminal -N=N=N) to N1 (e.g. azide)
		{{{"Npl#1#1","^N#2",""},{"#2","^N#3",""}}, {{"1:idatm=N1"}}}, // change 1-substituted Npl (terminal -N=N=N) to N1 (e.g. azide)
		{{{"Npl#1#3,1H","N1|N1\\+#2",""}}, {{"1:idatm=N2+"}}}, // change 3-substituted Npl (first nitrogen in -N(-H)=N=N) to N2+ (e.g. azide)
		{{{"^S#1#4","Npl#2#2",""},{"#1","Npl#3#2",""},{"#1","^C#4",""},{"#1","^C#5",""}}, {{"2:idatm=N2"},{"3:idatm=N2"}}}, // (R-)(R-)S(=N-)(=N-)
		{{{"^P#1#3","Npl#2#2",""},{"#1","Npl#3#2",""}}, {{"2:idatm=N2"},{"3:idatm=N2"}}}, // (R-)P(=N-)(=N-)
		{{{"Npl#1#2","C2#2#3",""},{"#2","C3#3",""},{"#2","C3#4",""}}, {{"1:idatm=N2"}}}, // change 2-substituted Npl bound to isolated C2 to N2
		{{{"C2#1#3","N3#2",""},{"#1","C3#3",""},{"#1","Pox#4",""}}, {{"1:idatm=C3"}}}, // change wrongly assigned C2 (Ligand ID: POB, Atom Name: C1') to C3
		{{{"C3#1#1","C1#2",""},{"#2","^C[^1]#3",""}}, {{"1:idatm=C1"}}}, // change terminal C3 bound to C1, which in turn is bound to any C but C1 (Ligand ID: 1DJ, Atom Name: CAA) to C1
	};
	const rename_rules refine { // idatm rules for final refinement (relies on bond gaff types)
		{{{"Npl#1#3#sb,db,ag6",".*#2",""}}, {{"1:idatm=N2+"}}}, // change 3-substituted Npl in 6-membered aromatic ring to N2+
		{{{"Npl#1#3,1H#sb,db",".*#2",""}}, {{"1:idatm=N2"}}}, // change double-bonded Npl(-H) to N2 (deletes one hydrogen)
		{{{"Oar#1#2#sb,db,ag",".*#2",""}}, {{"1:idatm=Oar+"}}}, // aromatic oxygen, formally positive (pyrylium)
	};	

	const rename_rules bond_gaff_type {
		{{{"Cac#1","^O#2",""},{"#1","^O#3",""}}, {{"1,2:bond_gaff_type=DL"},{"1,3:bond_gaff_type=DL"}}},
		{{{"Ntr#1","^O#2",""},{"#1","^O#3",""}}, {{"1,2:bond_gaff_type=DL"},{"1,3:bond_gaff_type=DL"}}},
		{{{".*#1",".*#2","bo=3"}}, {{"1,2:bond_gaff_type=tb"}}},
		//~ {{{".*#1##AR1",".*#2##AR1","bo=2"}}, {{"1,2:bond_gaff_type=db"}}},
		//~ {{{".*#1##AR1",".*#2##AR2","bo=2"}}, {{"1,2:bond_gaff_type=db"}}},
		//~ {{{".*#1##AR2",".*#2##AR2","bo=2"}}, {{"1,2:bond_gaff_type=db"}}},
		{{{".*#1",".*#2","bo=2"}}, {{"1,2:bond_gaff_type=db"}}},
		//~ {{{".*#1##AR1",".*#2##AR1","bo=1"}}, {{"1,2:bond_gaff_type=sb"}}},
		//~ {{{".*#1##AR1",".*#2##AR2","bo=1"}}, {{"1,2:bond_gaff_type=sb"}}},
		//~ {{{".*#1##AR2",".*#2##AR2","bo=1"}}, {{"1,2:bond_gaff_type=sb"}}},
		{{{".*#1",".*#2","bo=1"}}, {{"1,2:bond_gaff_type=sb"}}},
	};	

	const rename_rules rotatable {
		{{{"Npl#1","C2#2",""}, {"#1","^H$#3",""}, {"#2","O2#4",""}}, {{"1,2:rota=amide,angles={180},drive_id=1"}}},
		{{{"3$#1","3$#2",""}, {"#1","^[^H]#3",""}}, {{"1,2:rota=sp3-sp3,angles={-60,60,180},drive_id=3"}}},
		{{{"3$#1","2$#2",""}, {"#1","^[^H]#3",""}, {"#2","^[^H]#4",""}}, {{"1,2:rota=sp3-sp2,angles={-150,-90,-30,30,90,150},drive_id=6"}}},
		{{{"3$#1","ar$#2",""}, {"#1","^[^H]#3",""}}, {{"1,2:rota=sp3-aromatic,angles={-150,-90,-30,30,90,150},drive_id=6"}}},
		{{{"3$#1","N3\\+#2",""}, {"#1","^[^H]#3",""}, {"#2","^[^H]#4",""}}, {{"1,2:rota=sp3-n4,angles={-60,60,180},drive_id=3"}}},
		{{{"3$#1","Npl#2",""}, {"#1","^[^H]#3",""}}, {{"1,2:rota=sp3-npl,angles={-150,-90,-30,30,90,150},drive_id=6"}}},
		{{{"3$#1","N2\\+#2",""}, {"#1","^[^H]#3",""}}, {{"1,2:rota=sp3-n2plus,angles={-150,-90,-30,30,90,150},drive_id=6"}}}, // added (issue #103)
		{{{"3$#1","1$#2",""}, {"#1","^[^H]#3",""}, {"#2","^[^H]#4",""}, {"#4","^[^H]#5",""}}, {{"1,2:rota=sp3-sp1,angles={-60,60,180},drive_id=3"}}},
		{{{"2$#1","2$#2",""}, {"#1","^[^H]#3",""}, {"#2","^[^H]#4",""}}, {{"1,2:rota=sp2-sp2,angles={0,180},drive_id=2"}}},
		{{{"2$#1","ar$#2",""}, {"#1","^[^H]#3",""}}, {{"1,2:rota=sp2-aromatic,angles={0,180},drive_id=2"}}},
		{{{"2$#1","N3\\+#2",""}, {"#1","^[^H]#3",""}, {"#2","^[^H]#4",""}}, {{"1,2:rota=sp2-n4,angles={-150,-90,-30,30,90,150},drive_id=6"}}},
		{{{"2$#1","Npl$#2",""}, {"#1","^[^H]#3",""}}, {{"1,2:rota=sp2-npl,angles={0,180},drive_id=2"}}},
		{{{"Npl$#1","ar$#2",""}, {"#1","^[^HO]#3",""}}, {{"1,2:rota=npl3-aromatic,angles={0,180},drive_id=2"}}},
		{{{"Npl$#1","Npl$#2",""}}, {{"1,2:rota=n_amide-n_amide,angles={0,180},drive_id=2"}}},
		{{{"Npl$#1","ar$#2",""}}, {{"1,2:rota=n_amide-aromatic,angles={0,180},drive_id=2"}}},
		{{{"C2$#1","2$#2",""}, {"#1","O2$#3",""}, {"#1","3$#4",""}, {"#2","^[^H]#5",""}}, {{"1,2:rota=acetyl-sp2,angles={0,180},drive_id=2"}}},
		{{{"C2#1","^N#2",""}, {"#2","^[^H]#3",""}, {"#1","^N#4",""}, {"#1","^N#5",""}}, {{"1,2:rota=ccat-n,angles={0,180},drive_id=2"}}}, // guanidinum carbon is C2 (SYBYL type is C.cat)
		{{{"C2#1","3$#2",""}, {"#2","^[^H]#3",""}, {"#1","O2#4",""}, {"#1","O3#5",""}}, {{"1,2:rota=carboxyl-sp3,angles={-90,-30,30},drive_id=62"}}},
		{{{"Cac#1","3$#2",""}, {"#2","^[^H]#3",""}}, {{"1,2:rota=carboxylate-sp3,angles={-90,-30,30},drive_id=62"}}}, // added because Chimera atom type differs from SYBYL
		{{{"O3#1","Car#2",""}, {"#1","^[^H]#3",""}}, {{"1,2:rota=arom-oxygen,angles={0,180},drive_id=2"}}},
		{{{"O3#1","2$#2",""}, {"#1","^[^H]#3",""}, {"#2","^[^H]#4",""}}, {{"1,2:rota=conj-oxygen,angles={0,180},drive_id=2"}}},
		{{{"S3#1","Car#2",""}, {"#1","^[^H]#3",""}}, {{"1,2:rota=arom-sulfur,angles={0,180},drive_id=2"}}},
		{{{"S3#1","2$#2",""}, {"#1","^[^H]#3",""}, {"#2","^[^H]#4",""}}, {{"1,2:rota=conj-sulfur,angles={0,180},drive_id=2"}}},
		{{{"C2#1","O3#2",""}, {"#1","O2#3",""}, {"#2","^[^H]#4",""}, {"#4","^[^H]#5",""}}, {{"1,2:rota=ester,angles={0,180},drive_id=2"}}},
		{{{"C2#1","S3#2",""}, {"#1","S2#3",""}, {"#2","^[^H]#4",""}, {"#4","^[^H]#5",""}}, {{"1,2:rota=sulfyl-ester,angles={0,180},drive_id=2"}}},
		{{{"O3#1","Car#2",""}, {"#1","^H$#3",""}}, {{"1,2:rota=arom-hydroxyl,angles={0,180},drive_id=2"}}},
		{{{"O3#1","2$#2",""}, {"#1","^H$#3",""}}, {{"1,2:rota=conj-hydroxyl,angles={0,180},drive_id=2"}}},
		{{{"S3#1","Car#2",""}, {"#1","^H$#3",""}}, {{"1,2:rota=arom-hydrosulfyl,angles={0,180},drive_id=2"}}},
		{{{"S3#1","2$#2",""}, {"#1","^H$#3",""}, {"#2","^[^H]#4",""}}, {{"1,2:rota=conj-hydrosulfyl,angles={0,180},drive_id=2"}}},
		{{{"Pac#1","O3#2",""}}, {{"1,2:rota=phosph-ester,angles={-90,-30,30},drive_id=62"}}},
		{{{"Sac#1","O3#2",""}}, {{"1,2:rota=sulfate-ester,angles={-90,-30,30},drive_id=62"}}},
		{{{"^S#1","Car#2",""}, {"#1","^O#3",""}, {"#1","^O#4",""}, {"#1","^O#5",""}}, {{"1,2:rota=sulfonate-arom,angles={-30,30},drive_id=63"}}},
		{{{"^S#1","2$#2",""}, {"#2","^[^H]#3",""}, {"#1","^O#4",""}, {"#1","^O#5",""}, {"#1","^O#6",""}}, {{"1,2:rota=sulfonate-sp2,angles={-30,30},drive_id=63"}}},
		{{{"^S#1","^N#2",""}, {"#1","^O#3",""}, {"#1","^O#4",""}}, {{"1,2:rota=sulfonamide,angles={-60,60,180},drive_id=3"}}},
		{{{"^S#1","^N#2",""}, {"#1","^O#3",""}, {"#1","^O#4",""}, {"#1","Car#5",""}, {"#2","^[^H]#6",""}}, {{"1,2:rota=aromatic_sulfonamide,angles={-60,60},drive_id=32"}}},
		{{{"^S#1","3$#2",""}, {"#1","^O#3",""}, {"#1","^O#4",""}, {"#2","^[^H]#5",""}}, {{"1,2:rota=sulfone-sp3,angles={-60,60,180},drive_id=3"}}},
		{{{"^S#1","2$#2",""}, {"#1","^O#3",""}, {"#1","^O#4",""}, {"#2","^[^H]#5",""}}, {{"1,2:rota=sulfone-sp2,angles={-150,-90,-30,30,90,150},drive_id=6"}}},
		{{{"^S#1","Car#2",""}, {"#1","^O#3",""}, {"#1","^O#4",""}, {"#2","^[^H]#5",""}}, {{"1,2:rota=sulfone-aromatic,angles={-120,-90,-60,60,90,120},drive_id=122"}}},
		{{{"Si#1","3$#2",""}, {"#1","^[^H]#3",""}, {"#2","^[^H]#4",""}}, {{"1,2:rota=silicon-sp3,angles={-150,-90,-30,30,90,150},drive_id=6"}}},
		{{{"Si#1","O3#2",""}, {"#1","^[^H]#3",""}}, {{"1,2:rota=silicon-oxygen,angles={-60,60,180},drive_id=3"}}},
		{{{"Car#1","Car#2",""}, {"#1","^[^H]#3",""}, {"#1","^[^H]#4",""}}, {{"1,2:rota=biphenyl,angles={-140,-40,40,140},drive_id=1000"}}}, // biphenly single bond
	};

	const rename_rules atomic_penalty_scores {
		/* 
		 * IMPORTANT: should be ordered from more to less specific within each atom & connectivity grouping
		 */
		{{{"^Ni$#1#5",".*#2",""}}, {{"1:aps={{5,1},{6,0},{7,1}}"}}}, // Ni(X5)
		{{{"^Si$#1",".*#2",""}}, {{"1:aps={{4,0}}"}}},
		{{{"^H$|^HC$|^D$|^DC$#1",".*#2",""}}, {{"1:aps={{0,64},{1,0},{2,64}}"}}},
		{{{"^F$#1",".*#2",""}}, {{"1:aps={{0,64},{1,0},{2,64}}"}}},
		{{{"^Cl$#1",".*#2",""}}, {{"1:aps={{0,64},{1,0},{2,64}}"}}},
		{{{"^Br$#1",".*#2",""}}, {{"1:aps={{0,64},{1,0},{2,64}}"}}},
		{{{"^I$#1",".*#2",""}}, {{"1:aps={{0,64},{1,0},{2,64}}"}}},

		{{{"^C#1#1","^N#2#2",""}}, {{"1:aps={{3,0},{4,1},{5,32}}"}}}, // APSC1+ : C in C#N-R
		{{{"^C#1#1",".*#2",""}}, {{"1:aps={{3,1},{4,0},{5,32}}"}}}, // C(X1)
		{{{"^C#1#3","^O|^S#2#1",""},{"#1","^O|^S#3#1",""}}, {{"1:aps={{4,32},{5,0},{6,32}}"}}}, // APSCO2 : for CO2-, CS2-, COS- etc
		{{{"^C#1",".*#2",""}}, {{"1:aps={{2,64},{3,32},{4,0},{5,32},{6,64}}"}}}, // just C ?

		{{{"^N#1#1","^N#2#2",""}}, {{"1:aps={{2,0},{3,0}}"}}}, // APSN1- : N(X1) in N=N=R
		{{{"^N#1#1",".*#2",""}}, {{"1:aps={{2,3},{3,0},{4,32}}"}}}, // N(X1)
		{{{"^N#1#2",".*#2",""}}, {{"1:aps={{2,4},{3,0},{4,2}}"}}}, // N(X2)
		{{{"^N#1#2","^N|^C#2#1",""}}, {{"1:aps={{3,1},{4,0}}"}}}, // APSN2+ : N(X2) in N=N=R
		{{{"^N#1#3","^O|^S#2#1",""},{"#1","^O|^S#3#1",""}}, {{"1:aps={{3,64},{4,32},{5,0},{6,32}}"}}}, // APSNO2 : for NO2-, NS2-, NOS- etc
		{{{"^N#1#3","^O|^S#2#1",""},{"#1","^[^OS]#3",""},{"#1","^[^OS]#4",""}}, {{"1:aps={{3,1},{4,0}}"}}}, // APSN3+
		{{{"^N#1#3",".*#2",""}}, {{"1:aps={{2,32},{3,0},{4,1},{5,2}}"}}}, // N(X3)
		{{{"^N#1#4",".*#2",""}}, {{"1:aps={{3,64},{4,0},{5,64}}"}}}, // N(X4)

		{{{"^O#1#1","^N#2#3",""},{"#2","^[^OS]#3",""},{"#2","^[^OS]#4",""}}, {{"1:aps={{1,0},{2,1}}"}}}, // APSO1- : for O in pyridine 1-oxide etc
		//~ {{{"^O#1#1","^N#2#3",""},{"#2","^O|^S#3#2,3,4,5,6",""}}, {{"1:aps={{1,0},{2,1}}"}}}, // APSO1- : for O in pyridine 1-oxide etc
		{{{"^O#1#1","^N#2#3",""},{"#2","^O|^S#3#2",""}}, {{"1:aps={{1,0},{2,1}}"}}}, // APSO1- : for O in pyridine 1-oxide etc
		{{{"^O#1#1","^N#2#3",""},{"#2","^O|^S#3#3",""}}, {{"1:aps={{1,0},{2,1}}"}}}, // APSO1- : for O in pyridine 1-oxide etc
		{{{"^O#1#1","^N#2#3",""},{"#2","^O|^S#3#4",""}}, {{"1:aps={{1,0},{2,1}}"}}}, // APSO1- : for O in pyridine 1-oxide etc
		{{{"^O#1#1","^N#2#3",""},{"#2","^O|^S#3#5",""}}, {{"1:aps={{1,0},{2,1}}"}}}, // APSO1- : for O in pyridine 1-oxide etc
		{{{"^O#1#1","^N#2#3",""},{"#2","^O|^S#3#6",""}}, {{"1:aps={{1,0},{2,1}}"}}}, // APSO1- : for O in pyridine 1-oxide etc
		{{{"^O#1#1",".*#2",""}}, {{"1:aps={{1,1},{2,0},{3,64}}"}}}, // O(X1)
		{{{"^O#1#2",".*#2",""}}, {{"1:aps={{1,32},{2,0},{3,64}}"}}}, // O(X2)

		{{{"^P#1#4","^O|^S#2#1",""},{"#1","^O|^S#3#1",""},{"#1","^[^OS]#4",""},{"#1","^[^OS]#5",""}}, {{"1:aps={{5,32},{6,0},{7,32}}"}}}, // APSPO2 : for PO2-,PS2-, POS- 
		{{{"^P#1#4","^O|^S#2#1",""},{"#1","^O|^S#3#1",""},{"#1","^O|^S#4#1",""}}, {{"1:aps={{6,32},{7,0}}"}}}, // APSPO3 : for PO3-,POS2-, PO2S-, PS3-
		{{{"^P#1#1",".*#2",""}}, {{"1:aps={{2,2},{3,0},{4,32}}"}}}, // P(X1)
		{{{"^P#1#2",".*#2",""}}, {{"1:aps={{2,4},{3,0},{4,2}}"}}}, // P(X2)
		{{{"^P#1#3","^O#2#1",""},{"#1","^O#3#1",""},{"#1","^O#4#2",""}}, {{"1:aps={{4,32},{5,0},{6,32}}"}}}, //Janez : for -O-P(=O)(=O)
		{{{"^P#1#3","^N#2#2",""},{"#1","^N#3#2",""},{"#1","^N#4#3",""}}, {{"1:aps={{4,32},{5,0},{6,32}}"}}}, //Janez : for -NH-P(=N)(=N)
		{{{"^P#1#3",".*#2",""}}, {{"1:aps={{2,32},{3,0},{4,1},{5,2}}"}}}, // P(X3)
		{{{"^P#1#4",".*#2",""}}, {{"1:aps={{3,64},{4,1},{5,0},{6,32}}"}}}, // P(X4)
		{{{"^P#1#5",".*#2",""}}, {{"1:aps={{3,64},{4,1},{5,0},{6,32}}"}}}, // Janez : P(X5) (copy of P(X4))

		{{{"^S#1#4","^O|^S#2#1",""},{"#1","^O|^S#3#1",""},{"#1","^O|^S#4#1",""},{"#1","^O|^S#5#1",""}}, {{"1:aps={{6,32},{7,0}}"}}}, // APSSO4 : for SO4-, SOS3-, SO3S- etc
		{{{"^S#1#4","^O|^S#2#1",""},{"#1","^O|^S#3#1",""},{"#1","^O|^S#4#1",""},{"#1","^[^OS]#5",""}}, {{"1:aps={{6,32},{7,0}}"}}}, // APSSO3 : for SO3-, SOS2-, SO2S- etc
		{{{"^S#1#4","^O|^S#2#1",""},{"#1","^O|^S#3#1",""},{"#1","^[^OS]#4",""},{"#1","^[^OS]#5",""}}, {{"1:aps={{6,0},{7,32}}"}}}, // APSSO2 : for SO2-, SOS-, SS2-
		{{{"^S#1#4",".*#2",""}}, {{"1:aps={{4,4},{5,2},{6,0}}"}}}, // S(X4)
		{{{"^S#1#3",".*#2",""}}, {{"1:aps={{3,1},{4,0},{5,2},{6,2}}"}}}, // S(X3)
		{{{"^S#1#2",".*#2",""}}, {{"1:aps={{1,32},{2,0},{3,32}}"}}}, // S(X2)
		{{{"^S#1#1","^N#2#3",""},{"#2","^[^OS]#3",""},{"#2","^[^OS]#4",""}}, {{"1:aps={{1,0},{2,1}}"}}}, // APSS1- : for S in pyridine 1-oxide etc
		//~ {{{"^S#1#1","^N#2#3",""},{"#2","^O|^S#3#2,3,4,5,6",""}}, {{"1:aps={{1,0},{2,1}}"}}}, // APSS1- : for S in pyridine 1-oxide etc
		{{{"^S#1#1","^N#2#3",""},{"#2","^O|^S#3#2",""}}, {{"1:aps={{1,0},{2,1}}"}}}, // APSS1- : for S in pyridine 1-oxide etc
		{{{"^S#1#1","^N#2#3",""},{"#2","^O|^S#3#3",""}}, {{"1:aps={{1,0},{2,1}}"}}}, // APSS1- : for S in pyridine 1-oxide etc
		{{{"^S#1#1","^N#2#3",""},{"#2","^O|^S#3#4",""}}, {{"1:aps={{1,0},{2,1}}"}}}, // APSS1- : for S in pyridine 1-oxide etc
		{{{"^S#1#1","^N#2#3",""},{"#2","^O|^S#3#5",""}}, {{"1:aps={{1,0},{2,1}}"}}}, // APSS1- : for S in pyridine 1-oxide etc
		{{{"^S#1#1","^N#2#3",""},{"#2","^O|^S#3#6",""}}, {{"1:aps={{1,0},{2,1}}"}}}, // APSS1- : for S in pyridine 1-oxide etc
		{{{"^S#1#1",".*#2",""}}, {{"1:aps={{1,2},{2,0},{3,64}}"}}}, // S(X1)
	};
	
	const rename_rules gaff { // GAFF atom types
		{{{"Cl#1",".*#2",""}}, {{"1:gaff=cl"}}}, // Cl is first to not mistake it with C
		{{{"Si#1",".*#2",""}}, {{"1:gaff=Si"}}},
		// 3-membered ring atom
		{{{"^C#1#4#RG3",".*#2",""}}, {{"1:gaff=cx"}}}, // 3-membered ring atom
		// 4-membered ring atom
		{{{"^C#1#4#RG4",".*#2",""}}, {{"1:gaff=cy"}}}, // 4-membered ring atom
		// sp3 C
		{{{"^C#1#4",".*#2",""}}, {{"1:gaff=c3"}}}, // sp3 C
		// C=O or C=S
		{{{"^C#1#3#2DL",XA + "#2#1",""}}, {{"1:gaff=c"}}}, // C=O or C=S
		//~ {{{"^C#1#3#1DB,0DL",XA + "#2#1",""}}, {{"1:gaff=c"}}}, // C=O or C=S
		{{{"^C#1#3#1db,0DL",XA + "#2#1",""}}, {{"1:gaff=c"}}}, // C=O or C=S
		{{{"^C#1#3#3sb",XA + "#2#1",""}}, {{"1:gaff=c"}}}, // C=O or C=S
		// sp2 C in guanidine group
		//~ {{{"^C#1#3","^N#2#3",""},{"#1","^N#3#3",""},{"#1","^N#4#3",""}}, {{"1:gaff=cz"}}}, // sp2 C in guanidine group
		{{{"^C#1#3#NG","^N#2#3",""},{"#1","^N#3#3",""},{"#1","^N#4#3",""}}, {{"1:gaff=cz"}}}, // sp2 C in guanidine group (Janez added NG)
		// pure aromatic atom that can form an aromatic single bond 
		{{{"^C#1#3#AR1,1RG6",XX + "#2##AR1",""},{"#1",XX + "#3##AR1",""},{"#1",XX + "#4##AR1",""}}, {{"1:gaff=cp"}}}, // pure aromatic atom that can form an aromatic single bond 
		// pure aromatic atom 
		{{{"^C#1#3#AR1",".*#2",""}}, {{"1:gaff=ca"}}}, // pure aromatic atom 
		// sp2 C of conjugated ring systems
		{{{"^C#1#3#sb,db,AR2","^C#2#3",""},{"#2","^C#3#3",""}}, {{"1:gaff=cc"}}}, // sp2 C of conjugated ring systems
		{{{"^C#1#3#sb,db,AR2","^C#2#3",""},{"#2","^C#3#2",""}}, {{"1:gaff=cc"}}}, // sp2 C of conjugated ring systems
		{{{"^C#1#3#sb,db,AR2","^C#2#3",""},{"#2",XB + "#3#2",""}}, {{"1:gaff=cc"}}}, // sp2 C of conjugated ring systems
		{{{"^C#1#3#sb,db,AR2",XB + "#2#2",""},{"#2",XB + "#3#2",""}}, {{"1:gaff=cc"}}}, // sp2 C of conjugated ring systems
		{{{"^C#1#3#sb,db,AR2",XB + "#2#2",""},{"#2","^C#3#2",""}}, {{"1:gaff=cc"}}}, // sp2 C of conjugated ring systems
		{{{"^C#1#3#sb,db,AR2",XB + "#2#2",""},{"#2","^C#3#3",""}}, {{"1:gaff=cc"}}}, // sp2 C of conjugated ring systems
		{{{"^C#1#3#sb,db,AR2","^C#2#3","bond_gaff_type=sb"}}, {{"1:gaff=cc"}}}, // sp2 C of conjugated ring systems
		{{{"^C#1#3#sb,db,AR2",XB + "#2#2","bond_gaff_type=sb"}}, {{"1:gaff=cc"}}}, // sp2 C of conjugated ring systems
		{{{"^C#1#3#sb,db,AR2",XD + "#2#3#db","bond_gaff_type=sb"}}, {{"1:gaff=cc"}}}, // sp2 C of conjugated ring systems
		{{{"^C#1#3#sb,db,AR2",XD + "#2#4#db","bond_gaff_type=sb"}}, {{"1:gaff=cc"}}}, // sp2 C of conjugated ring systems
		{{{"^C#1#3#sb,db,AR3","^C#2#3",""},{"#2","^C#3#3",""}}, {{"1:gaff=cc"}}}, // sp2 C of conjugated ring systems
		{{{"^C#1#3#sb,db,AR3","^C#2#3",""},{"#2","^C#3#2",""}}, {{"1:gaff=cc"}}}, // sp2 C of conjugated ring systems
		{{{"^C#1#3#sb,db,AR3","^C#2#3",""},{"#2",XB + "#3#2",""}}, {{"1:gaff=cc"}}}, // sp2 C of conjugated ring systems
		{{{"^C#1#3#sb,db,AR3",XB + "#2#2",""},{"#2",XB + "#3#2",""}}, {{"1:gaff=cc"}}}, // sp2 C of conjugated ring systems
		{{{"^C#1#3#sb,db,AR3",XB + "#2#2",""},{"#2","^C#3#2",""}}, {{"1:gaff=cc"}}}, // sp2 C of conjugated ring systems
		{{{"^C#1#3#sb,db,AR3",XB + "#2#2",""},{"#2","^C#3#3",""}}, {{"1:gaff=cc"}}}, // sp2 C of conjugated ring systems
		{{{"^C#1#3#sb,db,AR3","^C#2#3","bond_gaff_type=sb"}}, {{"1:gaff=cc"}}}, // sp2 C of conjugated ring systems
		{{{"^C#1#3#sb,db,AR3",XB + "#2#2","bond_gaff_type=sb"}}, {{"1:gaff=cc"}}}, // sp2 C of conjugated ring systems
		{{{"^C#1#3#sb,db,AR3",XD + "#2#3#db","bond_gaff_type=sb"}}, {{"1:gaff=cc"}}}, // sp2 C of conjugated ring systems
		{{{"^C#1#3#sb,db,AR3",XD + "#2#4#db","bond_gaff_type=sb"}}, {{"1:gaff=cc"}}}, // sp2 C of conjugated ring systems

		{{{"^C#1#3#sb,db,AR4","^C#2#3",""},{"#2","^C#3#3",""}}, {{"1:gaff=cc"}}}, // sp2 C of conjugated ring systems
		{{{"^C#1#3#sb,db,AR4","^C#2#3",""},{"#2","^C#3#2",""}}, {{"1:gaff=cc"}}}, // sp2 C of conjugated ring systems
		{{{"^C#1#3#sb,db,AR4","^C#2#3",""},{"#2",XB + "#3#2",""}}, {{"1:gaff=cc"}}}, // sp2 C of conjugated ring systems
		{{{"^C#1#3#sb,db,AR4",XB + "#2#2",""},{"#2",XB + "#3#2",""}}, {{"1:gaff=cc"}}}, // sp2 C of conjugated ring systems
		{{{"^C#1#3#sb,db,AR4",XB + "#2#2",""},{"#2","^C#3#2",""}}, {{"1:gaff=cc"}}}, // sp2 C of conjugated ring systems
		{{{"^C#1#3#sb,db,AR4",XB + "#2#2",""},{"#2","^C#3#3",""}}, {{"1:gaff=cc"}}}, // sp2 C of conjugated ring systems
		{{{"^C#1#3#sb,db,AR4","^C#2#3","bond_gaff_type=sb"}}, {{"1:gaff=cc"}}}, // sp2 C of conjugated ring systems
		{{{"^C#1#3#sb,db,AR4",XB + "#2#2","bond_gaff_type=sb"}}, {{"1:gaff=cc"}}}, // sp2 C of conjugated ring systems
		{{{"^C#1#3#sb,db,AR4",XD + "#2#3#db","bond_gaff_type=sb"}}, {{"1:gaff=cc"}}}, // sp2 C of conjugated ring systems
		{{{"^C#1#3#sb,db,AR4",XD + "#2#4#db","bond_gaff_type=sb"}}, {{"1:gaff=cc"}}}, // sp2 C of conjugated ring systems

		{{{"^C#1#3#sb,db,AR2",".*#2",""}}, {{"1:gaff=cc"}}}, // sp2 C of conjugated ring systems
		{{{"^C#1#3#sb,db,AR3",".*#2",""}}, {{"1:gaff=cc"}}}, // sp2 C of conjugated ring systems

		{{{"^C#1#3#sb,db,AR4",".*#2",""}}, {{"1:gaff=cc"}}}, // sp2 C of conjugated ring systems

		// sp2 C of conjugated chain systems
		//~ {{{"^C#1#3#sb,db","^C#2#3","bond_gaff_type=SB"}}, {{"1:gaff=ce"}}}, // sp2 C of conjugated chain systems
		//~ {{{"^C#1#3#sb,db","^C#2#2","bond_gaff_type=SB"}}, {{"1:gaff=ce"}}}, // sp2 C of conjugated chain systems
		//~ {{{"^C#1#3#sb,db",XB + "#2#2","bond_gaff_type=SB"}}, {{"1:gaff=ce"}}}, // sp2 C of conjugated chain systems
		//~ {{{"^C#1#3#sb,db",XD + "#2#3#db","bond_gaff_type=SB"}}, {{"1:gaff=ce"}}}, // sp2 C of conjugated chain systems
		//~ {{{"^C#1#3#sb,db",XD + "#2#4#db","bond_gaff_type=SB"}}, {{"1:gaff=ce"}}}, // sp2 C of conjugated chain systems
		{{{"^C#1#3#sb,db","^C#2#3","bond_gaff_type=sb"}}, {{"1:gaff=ce"}}}, // sp2 C of conjugated chain systems
		{{{"^C#1#3#sb,db","^C#2#2","bond_gaff_type=sb"}}, {{"1:gaff=ce"}}}, // sp2 C of conjugated chain systems
		{{{"^C#1#3#sb,db",XB + "#2#2","bond_gaff_type=sb"}}, {{"1:gaff=ce"}}}, // sp2 C of conjugated chain systems
		{{{"^C#1#3#sb,db",XD + "#2#3#db","bond_gaff_type=sb"}}, {{"1:gaff=ce"}}}, // sp2 C of conjugated chain systems
		{{{"^C#1#3#sb,db",XD + "#2#4#db","bond_gaff_type=sb"}}, {{"1:gaff=ce"}}}, // sp2 C of conjugated chain systems
		//sp2 carbon in a 3-membered ring
		{{{"^C#1#3#RG3",".*#2",""}}, {{"1:gaff=cu"}}}, //sp2 carbon in a 3-membered ring
		//sp2 carbon in a 4-membered ring
		{{{"^C#1#3#RG4",".*#2",""}}, {{"1:gaff=cv"}}}, //sp2 carbon in a 4-membered ring
		// other sp2 C
		{{{"^C#1#3",".*#2",""}}, {{"1:gaff=c2"}}}, // other sp2 C
		// sp C of conjugated systems
		//~ {{{"^C#1#2#sb,tb","^C#2#2","bond_gaff_type=SB"}}, {{"1:gaff=cg"}}}, // sp C of conjugated systems
		//~ {{{"^C#1#2#sb,tb","^C#2#3","bond_gaff_type=SB"}}, {{"1:gaff=cg"}}}, // sp C of conjugated systems
		//~ {{{"^C#1#2#sb,tb","^N#2#1","bond_gaff_type=SB"}}, {{"1:gaff=cg"}}}, // sp C of conjugated systems
		//~ {{{"^C#1#2#sb,tb",XB + "#2#2","bond_gaff_type=SB"}}, {{"1:gaff=cg"}}}, // sp C of conjugated systems
		{{{"^C#1#2#sb,tb","^C#2#2","bond_gaff_type=sb"}}, {{"1:gaff=cg"}}}, // sp C of conjugated systems
		{{{"^C#1#2#sb,tb","^C#2#3","bond_gaff_type=sb"}}, {{"1:gaff=cg"}}}, // sp C of conjugated systems
		{{{"^C#1#2#sb,tb","^N#2#1","bond_gaff_type=sb"}}, {{"1:gaff=cg"}}}, // sp C of conjugated systems
		{{{"^C#1#2#sb,tb",XB + "#2#2","bond_gaff_type=sb"}}, {{"1:gaff=cg"}}}, // sp C of conjugated systems
		// other sp C 
		{{{"^C#1#2",".*#2",""}}, {{"1:gaff=c1"}}}, // other sp C
		{{{"^C#1#1",".*#2",""}}, {{"1:gaff=c1"}}}, // other sp C
		{{{"^H#1#1","^N#2",""}}, {{"1:gaff=hn"}}}, // H on N
		{{{"^H#1#1","^O#2",""}}, {{"1:gaff=ho"}}}, // H on O (hydroxyl group)
		{{{"^H#1#1","^S#2",""}}, {{"1:gaff=hs"}}}, // H on S
		{{{"^H#1#1","^P#2",""}}, {{"1:gaff=hp"}}}, // H on P
		{{{"^H#1#1","^C#2",""},{"#2","^N#3#4",""}}, {{"1:gaff=hx"}}},
		{{{"^H#1#1","^O#2",""},{"#2","^H#3#1",""}}, {{"1:gaff=hw"}}},
		{{{"^H#1#1","^C#2#4",""},{"#2",EW + "#3",""},{"#2",EW + "#4",""},{"#2",EW + "#5",""}}, {{"1:gaff=h3"}}}, // H bonded to aliphatic carbon with 3 electrwd. group
		{{{"^H#1#1","^C#2#4",""},{"#2",EW + "#3",""},{"#2",EW + "#4",""}}, {{"1:gaff=h2"}}}, // H bonded to aliphatic carbon with 2 electrwd. group
		{{{"^H#1#1","^C#2#4",""},{"#2",EW + "#3",""}}, {{"1:gaff=h1"}}}, // H bonded to aliphatic carbon with 1 electrwd. group
		{{{"^H#1#1","^C#2#4",""}}, {{"1:gaff=hc"}}}, // H bonded to aliphatic carbon without electrwd. group
		{{{"^H#1#1","^C#2#3",""},{"#2",EW + "#3",""},{"#2",EW + "#4",""}}, {{"1:gaff=h5"}}}, // H bonded to non-sp3 carbon with 2 electrwd. group
		{{{"^H#1#1","^C#2#3",""},{"#2",EW + "#3",""}}, {{"1:gaff=h4"}}}, // H bonded to non-sp3 carbon with 1 electrwd. group
		{{{"^H#1#1",".*#2",""}}, {{"1:gaff=ha"}}}, // H bonded to aromatic carbon
		{{{"F#1",".*#2",""}}, {{"1:gaff=f"}}},
		//~ {{{"Cl#1",".*#2",""}}, {{"1:gaff=cl"}}},
		{{{"Br#1",".*#2",""}}, {{"1:gaff=br"}}},
		{{{"I#1",".*#2",""}}, {{"1:gaff=i"}}},
		// sp2 P of conjugated ring systems
		{{{"^P#1#2#sb,db,AR2","^C#2#3",""},{"#2","^C#3#3",""}}, {{"1:gaff=pc"}}}, // sp2 P of conjugated ring systems
		{{{"^P#1#2#sb,db,AR2","^C#2#3",""},{"#2","^C#3#2",""}}, {{"1:gaff=pc"}}}, // sp2 P of conjugated ring systems
		{{{"^P#1#2#sb,db,AR2","^C#2#3",""},{"#2",XB + "#3#2",""}}, {{"1:gaff=pc"}}}, // sp2 P of conjugated ring systems
		{{{"^P#1#2#sb,db,AR2",XB + "#2#2",""},{"#2","^C#3#3",""}}, {{"1:gaff=pc"}}}, // sp2 P of conjugated ring systems
		{{{"^P#1#2#sb,db,AR2",XB + "#2#2",""},{"#2","^C#3#2",""}}, {{"1:gaff=pc"}}}, // sp2 P of conjugated ring systems
		{{{"^P#1#2#sb,db,AR2",XB + "#2#2",""},{"#2",XB + "#3#2",""}}, {{"1:gaff=pc"}}}, // sp2 P of conjugated ring systems
		{{{"^P#1#2#sb,db,AR2","^C#2#3","bond_gaff_type=sb"}}, {{"1:gaff=pc"}}}, // sp2 P of conjugated ring systems
		{{{"^P#1#2#sb,db,AR2","^C#2#2","bond_gaff_type=sb"}}, {{"1:gaff=pc"}}}, // sp2 P of conjugated ring systems
		{{{"^P#1#2#sb,db,AR2",XB + "#2#2","bond_gaff_type=sb"}}, {{"1:gaff=pc"}}}, // sp2 P of conjugated ring systems
		{{{"^P#1#2#sb,db,AR2",XD + "#2#3#db","bond_gaff_type=sb"}}, {{"1:gaff=pc"}}}, // sp2 P of conjugated ring systems
		{{{"^P#1#2#sb,db,AR2",XD + "#2#4#db","bond_gaff_type=sb"}}, {{"1:gaff=pc"}}}, // sp2 P of conjugated ring systems
		{{{"^P#1#2#sb,db,AR3","^C#2#3",""},{"#2","^C#3#3",""}}, {{"1:gaff=pc"}}}, // sp2 P of conjugated ring systems
		{{{"^P#1#2#sb,db,AR3","^C#2#3",""},{"#2","^C#3#2",""}}, {{"1:gaff=pc"}}}, // sp2 P of conjugated ring systems
		{{{"^P#1#2#sb,db,AR3","^C#2#3",""},{"#2",XB + "#3#2",""}}, {{"1:gaff=pc"}}}, // sp2 P of conjugated ring systems
		{{{"^P#1#2#sb,db,AR3",XB + "#2#2",""},{"#2","^C#3#3",""}}, {{"1:gaff=pc"}}}, // sp2 P of conjugated ring systems
		{{{"^P#1#2#sb,db,AR3",XB + "#2#2",""},{"#2","^C#3#2",""}}, {{"1:gaff=pc"}}}, // sp2 P of conjugated ring systems
		{{{"^P#1#2#sb,db,AR3",XB + "#2#2",""},{"#2",XB + "#3#2",""}}, {{"1:gaff=pc"}}}, // sp2 P of conjugated ring systems
		{{{"^P#1#2#sb,db,AR3","^C#2#3","bond_gaff_type=sb"}}, {{"1:gaff=pc"}}}, // sp2 P of conjugated ring systems
		{{{"^P#1#2#sb,db,AR3","^C#2#2","bond_gaff_type=sb"}}, {{"1:gaff=pc"}}}, // sp2 P of conjugated ring systems
		{{{"^P#1#2#sb,db,AR3",XB + "#2#2","bond_gaff_type=sb"}}, {{"1:gaff=pc"}}}, // sp2 P of conjugated ring systems
		{{{"^P#1#2#sb,db,AR3",XD + "#2#3#db","bond_gaff_type=sb"}}, {{"1:gaff=pc"}}}, // sp2 P of conjugated ring systems
		{{{"^P#1#2#sb,db,AR3",XD + "#2#4#db","bond_gaff_type=sb"}}, {{"1:gaff=pc"}}}, // sp2 P of conjugated ring systems
		// Sp2 P in pure aromatic systems
		{{{"^P#1#2#AR1",".*#2",""}}, {{"1:gaff=pb"}}}, // Sp2 P in pure aromatic systems
		// sp2 P of conjugated chain systems
		{{{"^P#1#2#sb,db","^C#2#3","bond_gaff_type=sb"}}, {{"1:gaff=pe"}}}, // sp2 P of conjugated chain systems
		//~ {{{"^P#1#2#sb,db","^C#2#2","bond_gaff_type=SB"}}, {{"1:gaff=pe"}}}, // sp2 P of conjugated chain systems
		//~ {{{"^P#1#2#sb,db","^C#2#2","bond_gaff_type=SB"}}, {{"1:gaff=pe"}}}, // sp2 P of conjugated chain systems
		//~ {{{"^P#1#2#sb,db",XA + "#2#1","bond_gaff_type=SB"}}, {{"1:gaff=pe"}}}, // sp2 P of conjugated chain systems
		//~ {{{"^P#1#2#sb,db",XB + "#2#2","bond_gaff_type=SB"}}, {{"1:gaff=pe"}}}, // sp2 P of conjugated chain systems
		//~ {{{"^P#1#2#sb,db",XD + "#2#3#DB","bond_gaff_type=SB"}}, {{"1:gaff=pe"}}}, // sp2 P of conjugated chain systems
		//~ {{{"^P#1#2#sb,db",XD + "#2#4#DB","bond_gaff_type=SB"}}, {{"1:gaff=pe"}}}, // sp2 P of conjugated chain systems
		{{{"^P#1#2#sb,db","^C#2#2","bond_gaff_type=sb"}}, {{"1:gaff=pe"}}}, // sp2 P of conjugated chain systems
		{{{"^P#1#2#sb,db","^C#2#2","bond_gaff_type=sb"}}, {{"1:gaff=pe"}}}, // sp2 P of conjugated chain systems
		{{{"^P#1#2#sb,db",XA + "#2#1","bond_gaff_type=sb"}}, {{"1:gaff=pe"}}}, // sp2 P of conjugated chain systems
		{{{"^P#1#2#sb,db",XB + "#2#2","bond_gaff_type=sb"}}, {{"1:gaff=pe"}}}, // sp2 P of conjugated chain systems
		{{{"^P#1#2#sb,db",XD + "#2#3#db","bond_gaff_type=sb"}}, {{"1:gaff=pe"}}}, // sp2 P of conjugated chain systems
		{{{"^P#1#2#sb,db",XD + "#2#4#db","bond_gaff_type=sb"}}, {{"1:gaff=pe"}}}, // sp2 P of conjugated chain systems
		{{{"^P#1#2",".*#2",""}}, {{"1:gaff=p2"}}}, // Phosphorous with two connected atoms
		{{{"^P#1#1",".*#2",""}}, {{"1:gaff=p2"}}}, // Phosphorous with two connected atoms
		{{{"^P#1#3#db",XB + "#2","bond_gaff_type=sb"}}, {{"1:gaff=px"}}}, // Special p4 in conjugated systems (conj. P, 3 subst)
		{{{"^P#1#3#db","^C#2#3","bond_gaff_type=sb"}}, {{"1:gaff=px"}}}, // Special p4 in conjugated systems (conj. P, 3 subst)
		{{{"^P#1#3#db",XD + "#2#3#db","bond_gaff_type=sb"}}, {{"1:gaff=px"}}}, // Special p4 in conjugated systems (conj. P, 3 subst)
		{{{"^P#1#3#db",XD + "#2#4#db","bond_gaff_type=sb"}}, {{"1:gaff=px"}}}, // Special p4 in conjugated systems (conj. P, 3 subst)
		{{{"^P#1#3#db",XA + "#2#1",""}}, {{"1:gaff=p4"}}}, // Phosphorous with three connected atoms, such as O=P(CH3)2
		{{{"^P#1#3",".*#2",""}}, {{"1:gaff=p3"}}}, // Phosphorous with three connected atoms, such as PH3
		{{{"^P#1#4#db",XB + "#2#2","bond_gaff_type=sb"}}, {{"1:gaff=py"}}}, // Special p5 in conjugated systems (conj. P, 4 subst.)
		{{{"^P#1#4#db","^C#2#3","bond_gaff_type=sb"}}, {{"1:gaff=py"}}}, // Special p5 in conjugated systems (conj. P, 4 subst.)
		{{{"^P#1#4#db",XD + "#2#3#db","bond_gaff_type=sb"}}, {{"1:gaff=py"}}}, // Special p5 in conjugated systems (conj. P, 4 subst.)
		{{{"^P#1#4#db",XD + "#2#4#db","bond_gaff_type=sb"}}, {{"1:gaff=py"}}}, // Special p5 in conjugated systems (conj. P, 4 subst.)
		{{{"^P#1#4",".*#2",""}}, {{"1:gaff=p5"}}}, // Phosphate with four connected atoms, such as O=P(OH)3
		{{{"^P#1#5",".*#2",""}}, {{"1:gaff=p5"}}}, // Phosphate with four connected atoms, such as O=P(OH)3
		{{{"^P#1#6",".*#2",""}}, {{"1:gaff=p5"}}}, // Phosphate with four connected atoms, such as O=P(OH)3

		//~ {{{"^N#1#3","^C#2#3",""},{"#2",XA + "#2#1",""}}, {{"1:gaff=n"}}}, // Sp2 nitrogen in amide groups N-CO, N-SO2, N-PO
		{{{"^N#1#3","^C#2#3",""},{"#2",XA + "#3#1",""}}, {{"1:gaff=n"}}}, // Sp2 nitrogen in amide groups N-CO, N-SO2, N-PO
		{{{"^N#1#3","^S#2#4",""},{"#2",XA + "#3#1",""}}, {{"1:gaff=n"}}}, // Sp2 nitrogen in amide groups N-CO, N-SO2, N-PO
		{{{"^N#1#3","^P#2#4",""},{"#2",XA + "#3#1",""}}, {{"1:gaff=n"}}}, // Sp2 nitrogen in amide groups N-CO, N-SO2, N-PO

		{{{"^N#1#4",".*#2",""}}, {{"1:gaff=n4"}}}, // Sp3 N with four connected atoms
		{{{"^N#1#3","^O#2#1",""},{"#1","^O#3#1",""}}, {{"1:gaff=no"}}}, // Nitro N
		
		//~ {{{"^N#1#3#AR1|AR2|AR3",".*#2",""}}, {{"1:gaff=na"}}}, // Sp2 N with three connected atoms
		//~ {{{"^N#1#3",XX + "#2##AR1|AR2|AR3",""}}, {{"1:gaff=nh"}}}, // Amine N connected to one or more aromatic rings
		{{{"^N#1#3#AR1",".*#2",""}}, {{"1:gaff=na"}}}, // Sp2 N with three connected atoms
		{{{"^N#1#3#AR2",".*#2",""}}, {{"1:gaff=na"}}}, // Sp2 N with three connected atoms
		{{{"^N#1#3#AR3",".*#2",""}}, {{"1:gaff=na"}}}, // Sp2 N with three connected atoms

		{{{"^N#1##NG","^C#2##NG",""},{"#2","^N#3",""},{"#2","^N#4",""}}, {{"1:gaff=nh"}}}, // Janez : Amine N in guanidino group
		{{{"^N#1#3",XX + "#2##AR1",""}}, {{"1:gaff=nh"}}}, // Amine N connected to one or more aromatic rings
		{{{"^N#1#3",XX + "#2##AR2",""}}, {{"1:gaff=nh"}}}, // Amine N connected to one or more aromatic rings
		{{{"^N#1#3",XX + "#2##AR3",""}}, {{"1:gaff=nh"}}}, // Amine N connected to one or more aromatic rings

		//~ {{{"^N#1#3","^C#2#3#DB",""}}, {{"1:gaff=nh"}}}, // Amine N connected to one or more aromatic rings
		//~ {{{"^N#1#3","^N#2#2#DB",""}}, {{"1:gaff=nh"}}}, // Amine N connected to one or more aromatic rings
		//~ {{{"^N#1#3","^P#2#2#DB",""}}, {{"1:gaff=nh"}}}, // Amine N connected to one or more aromatic rings
		{{{"^N#1#3","^C#2#3#AR1",""}}, {{"1:gaff=nh"}}}, // Amine N connected to one or more aromatic rings
		{{{"^N#1#3","^C#2#3#AR2",""}}, {{"1:gaff=nh"}}}, // Amine N connected to one or more aromatic rings
		{{{"^N#1#3","^C#2#3#AR3",""}}, {{"1:gaff=nh"}}}, // Amine N connected to one or more aromatic rings
		{{{"^N#1#3","^N#2#2#AR1",""}}, {{"1:gaff=nh"}}}, // Amine N connected to one or more aromatic rings
		{{{"^N#1#3","^N#2#2#AR2",""}}, {{"1:gaff=nh"}}}, // Amine N connected to one or more aromatic rings
		{{{"^N#1#3","^N#2#2#AR3",""}}, {{"1:gaff=nh"}}}, // Amine N connected to one or more aromatic rings
		{{{"^N#1#3","^P#2#2#AR1",""}}, {{"1:gaff=nh"}}}, // Amine N connected to one or more aromatic rings
		{{{"^N#1#3","^P#2#2#AR2",""}}, {{"1:gaff=nh"}}}, // Amine N connected to one or more aromatic rings
		{{{"^N#1#3","^P#2#2#AR3",""}}, {{"1:gaff=nh"}}}, // Amine N connected to one or more aromatic rings
		{{{"^N#1#3",".*#2",""}}, {{"1:gaff=n3"}}}, // Sp3 N with three connected atoms
		{{{"^N#1#2#AR1",".*#2",""}}, {{"1:gaff=nb"}}}, // Sp2 N in pure aromatic systems (aromatic nitrogen)
		// sp2 N of conjugated ring systems
		{{{"^N#1#2#sb,db,AR2","^C#2#3",""},{"#2","^C#3#3",""}}, {{"1:gaff=nc"}}}, // sp2 N of conjugated ring systems
		{{{"^N#1#2#sb,db,AR2","^C#2#3",""},{"#2","^C#3#2",""}}, {{"1:gaff=nc"}}}, // sp2 N of conjugated ring systems
		{{{"^N#1#2#sb,db,AR2","^C#2#3",""},{"#2",XB + "#3#2",""}}, {{"1:gaff=nc"}}}, // sp2 N of conjugated ring systems
		{{{"^N#1#2#sb,db,AR2",XB + "#2#2",""},{"#2","^C#3#3",""}}, {{"1:gaff=nc"}}}, // sp2 N of conjugated ring systems
		{{{"^N#1#2#sb,db,AR2",XB + "#2#2",""},{"#2","^C#3#2",""}}, {{"1:gaff=nc"}}}, // sp2 N of conjugated ring systems
		{{{"^N#1#2#sb,db,AR2",XB + "#2#2",""},{"#2",XB + "#3#2",""}}, {{"1:gaff=nc"}}}, // sp2 N of conjugated ring systems
		{{{"^N#1#2#sb,db,AR2","^C#2#3","bond_gaff_type=sb"}}, {{"1:gaff=nc"}}}, // sp2 N of conjugated ring systems
		{{{"^N#1#2#sb,db,AR2",XB + "#2#2","bond_gaff_type=sb"}}, {{"1:gaff=nc"}}}, // sp2 N of conjugated ring systems
		{{{"^N#1#2#sb,db,AR2",XD + "#2#3#db","bond_gaff_type=sb"}}, {{"1:gaff=nc"}}}, // sp2 N of conjugated ring systems
		{{{"^N#1#2#sb,db,AR2",XD + "#2#4#db","bond_gaff_type=sb"}}, {{"1:gaff=nc"}}}, // sp2 N of conjugated ring systems

		{{{"^N#1#2#sb,db,AR3","^C#2#3",""},{"#2","^C#3#3",""}}, {{"1:gaff=nc"}}}, // sp2 N of conjugated ring systems
		{{{"^N#1#2#sb,db,AR3","^C#2#3",""},{"#2","^C#3#2",""}}, {{"1:gaff=nc"}}}, // sp2 N of conjugated ring systems
		{{{"^N#1#2#sb,db,AR3","^C#2#3",""},{"#2",XB + "#3#2",""}}, {{"1:gaff=nc"}}}, // sp2 N of conjugated ring systems
		{{{"^N#1#2#sb,db,AR3",XB + "#2#2",""},{"#2","^C#3#3",""}}, {{"1:gaff=nc"}}}, // sp2 N of conjugated ring systems
		{{{"^N#1#2#sb,db,AR3",XB + "#2#2",""},{"#2","^C#3#2",""}}, {{"1:gaff=nc"}}}, // sp2 N of conjugated ring systems
		{{{"^N#1#2#sb,db,AR3",XB + "#2#2",""},{"#2",XB + "#3#2",""}}, {{"1:gaff=nc"}}}, // sp2 N of conjugated ring systems
		{{{"^N#1#2#sb,db,AR3","^C#2#3","bond_gaff_type=sb"}}, {{"1:gaff=nc"}}}, // sp2 N of conjugated ring systems
		{{{"^N#1#2#sb,db,AR3",XB + "#2#2","bond_gaff_type=sb"}}, {{"1:gaff=nc"}}}, // sp2 N of conjugated ring systems
		{{{"^N#1#2#sb,db,AR3",XD + "#2#3#db","bond_gaff_type=sb"}}, {{"1:gaff=nc"}}}, // sp2 N of conjugated ring systems
		{{{"^N#1#2#sb,db,AR3",XD + "#2#4#db","bond_gaff_type=sb"}}, {{"1:gaff=nc"}}}, // sp2 N of conjugated ring systems

		{{{"^N#1#2#sb,db,AR4","^C#2#3",""},{"#2","^C#3#3",""}}, {{"1:gaff=nc"}}}, // sp2 N of conjugated ring systems
		{{{"^N#1#2#sb,db,AR4","^C#2#3",""},{"#2","^C#3#2",""}}, {{"1:gaff=nc"}}}, // sp2 N of conjugated ring systems
		{{{"^N#1#2#sb,db,AR4","^C#2#3",""},{"#2",XB + "#3#2",""}}, {{"1:gaff=nc"}}}, // sp2 N of conjugated ring systems
		{{{"^N#1#2#sb,db,AR4",XB + "#2#2",""},{"#2","^C#3#3",""}}, {{"1:gaff=nc"}}}, // sp2 N of conjugated ring systems
		{{{"^N#1#2#sb,db,AR4",XB + "#2#2",""},{"#2","^C#3#2",""}}, {{"1:gaff=nc"}}}, // sp2 N of conjugated ring systems
		{{{"^N#1#2#sb,db,AR4",XB + "#2#2",""},{"#2",XB + "#3#2",""}}, {{"1:gaff=nc"}}}, // sp2 N of conjugated ring systems
		{{{"^N#1#2#sb,db,AR4","^C#2#3","bond_gaff_type=sb"}}, {{"1:gaff=nc"}}}, // sp2 N of conjugated ring systems
		{{{"^N#1#2#sb,db,AR4",XB + "#2#2","bond_gaff_type=sb"}}, {{"1:gaff=nc"}}}, // sp2 N of conjugated ring systems
		{{{"^N#1#2#sb,db,AR4",XD + "#2#3#db","bond_gaff_type=sb"}}, {{"1:gaff=nc"}}}, // sp2 N of conjugated ring systems
		{{{"^N#1#2#sb,db,AR4",XD + "#2#4#db","bond_gaff_type=sb"}}, {{"1:gaff=nc"}}}, // sp2 N of conjugated ring systems

		{{{"^N#1#2#RG5","^N#2#2#RG5",""},{"#2","^N#3#2#RG5",""},{"#3","^N#4#2#RG5",""}}, {{"1:gaff=nc"}}}, // sp2 N of conjugated ring systems
		{{{"^N#1#2#RG5","^N#2#2#RG5",""},{"#1","^N#3#2#RG5",""},{"#3","^N#4#2#RG5",""}}, {{"1:gaff=nc"}}}, // sp2 N of conjugated ring systems
		// sp2 N of conjugated chain systems
		//~ {{{"^N#1#2#sb,db","^C#2#2","bond_gaff_type=SB"}}, {{"1:gaff=ne"}}}, // sp2 N of conjugated chain systems
		//~ {{{"^N#1#2#sb,db","^C#2#3","bond_gaff_type=SB"}}, {{"1:gaff=ne"}}}, // sp2 N of conjugated chain systems
		//~ {{{"^N#1#2#sb,db",XA + "#2#1","bond_gaff_type=SB"}}, {{"1:gaff=ne"}}}, // sp2 N of conjugated chain systems
		//~ {{{"^N#1#2#sb,db",XB + "#2#2","bond_gaff_type=SB"}}, {{"1:gaff=ne"}}}, // sp2 N of conjugated chain systems
		//~ {{{"^N#1#2#sb,db",XD + "#2#3#db","bond_gaff_type=SB"}}, {{"1:gaff=ne"}}}, // sp2 N of conjugated chain systems
		//~ {{{"^N#1#2#sb,db",XD + "#2#4#db","bond_gaff_type=SB"}}, {{"1:gaff=ne"}}}, // sp2 N of conjugated chain systems
		{{{"^N#1#2#sb,db","^C#2#2","bond_gaff_type=sb"}}, {{"1:gaff=ne"}}}, // sp2 N of conjugated chain systems
		{{{"^N#1#2#sb,db","^C#2#3","bond_gaff_type=sb"}}, {{"1:gaff=ne"}}}, // sp2 N of conjugated chain systems
		{{{"^N#1#2#sb,db",XA + "#2#1","bond_gaff_type=sb"}}, {{"1:gaff=ne"}}}, // sp2 N of conjugated chain systems
		{{{"^N#1#2#sb,db",XB + "#2#2","bond_gaff_type=sb"}}, {{"1:gaff=ne"}}}, // sp2 N of conjugated chain systems
		{{{"^N#1#2#sb,db",XD + "#2#3#db","bond_gaff_type=sb"}}, {{"1:gaff=ne"}}}, // sp2 N of conjugated chain systems
		{{{"^N#1#2#sb,db",XD + "#2#4#db","bond_gaff_type=sb"}}, {{"1:gaff=ne"}}}, // sp2 N of conjugated chain systems
		{{{"^N#1#2#2db",".*#2",""}}, {{"1:gaff=n1"}}}, // Sp1 N
		{{{"^N#1#2#tb,sb",".*#2",""}}, {{"1:gaff=n1"}}}, // Sp1 N
		{{{"^N#1#2",".*#2",""}}, {{"1:gaff=n2"}}}, // aliphatic Sp2 N with two connected atoms
		{{{"^N#1#1",".*#2",""}}, {{"1:gaff=n1"}}}, // Sp1 N
		{{{"^O#1#1",".*#2",""}}, {{"1:gaff=o"}}}, // Oxygen with one connected atom (Sp2 oxygen in C=O, COO-)
		{{{"^O#1#2,1H",".*#2",""}}, {{"1:gaff=oh"}}}, // Oxygen in hydroxyl group
		{{{"^O#1#2,2H",".*#2",""}}, {{"1:gaff=oh"}}}, // Oxygen in hydroxyl group
		{{{"^O#1#3,1H",".*#2",""}}, {{"1:gaff=oh"}}}, // Oxygen in hydroxyl group
		{{{"^O#1#3,2H",".*#2",""}}, {{"1:gaff=oh"}}}, // Oxygen in hydroxyl group
		{{{"^O#1#3,3H",".*#2",""}}, {{"1:gaff=oh"}}}, // Oxygen in hydroxyl group
		{{{"^O#1#2",".*#2",""}}, {{"1:gaff=os"}}}, // Sp3 oxygen in ethers and esters
		{{{"^O#1#3",".*#2",""}}, {{"1:gaff=os"}}}, // Sp3 oxygen in ethers and esters
		{{{"^O#1",".*#2",""}}, {{"1:gaff=os"}}}, // Sp3 oxygen in ethers and esters
		{{{"^S#1#1",".*#2",""}}, {{"1:gaff=s"}}}, // S with one connected atom
		//~ {{{"^S#1#2#DB",".*#2",""}}, {{"1:gaff=s2"}}}, // S with two connected atoms, involved at least one double bond (sp2 S (P=S, C=S, ...) - old???) 
		{{{"^S#1#2#db",".*#2",""}}, {{"1:gaff=s2"}}}, // S with two connected atoms, involved at least one double bond (sp2 S (P=S, C=S, ...) - old???) 

		//~ {{{"^S#1#2#TB",".*#2",""}}, {{"1:gaff=s2"}}}, // S with two connected atoms, involved at least one double bond (sp2 S (P=S, C=S, ...) - old???) 
		{{{"^S#1#2#tb",".*#2",""}}, {{"1:gaff=s2"}}}, // S with two connected atoms, involved at least one double bond (sp2 S (P=S, C=S, ...) - old???) 

		{{{"^S#1#2,1H",".*#2",""}}, {{"1:gaff=sh"}}}, // Sp3 S in thiol groups
		{{{"^S#1#2,2H",".*#2",""}}, {{"1:gaff=sh"}}}, // Sp3 S in thiol groups
		{{{"^S#1#2",".*#2",""}}, {{"1:gaff=ss"}}}, // Sp3 S in thio-ester and thio-ether (sp3 sulfur in —SR and S—S)
		{{{"^S#1#3#db",XB + "#2","bond_gaff_type=sb"}}, {{"1:gaff=sx"}}}, // Special s4 in conjugated systems
		{{{"^S#1#3#db","^C#2#3","bond_gaff_type=sb"}}, {{"1:gaff=sx"}}}, // Special s4 in conjugated systems
		{{{"^S#1#3#db",XD + "#2#3#db","bond_gaff_type=sb"}}, {{"1:gaff=sx"}}}, // Special s4 in conjugated systems
		{{{"^S#1#3#db",XD + "#2#4#db","bond_gaff_type=sb"}}, {{"1:gaff=sx"}}}, // Special s4 in conjugated systems
		{{{"^S#1#3",".*#2",""}}, {{"1:gaff=s4"}}}, // S with three connected atoms
		{{{"^S#1#4#db",XB + "#2","bond_gaff_type=sb"}}, {{"1:gaff=sy"}}}, // Special s6 in conjugated systems
		{{{"^S#1#4#db","^C#2#3","bond_gaff_type=sb"}}, {{"1:gaff=sy"}}}, // Special s6 in conjugated systems
		{{{"^S#1#4#db",XD + "#2#3#db","bond_gaff_type=sb"}}, {{"1:gaff=sy"}}}, // Special s6 in conjugated systems
		{{{"^S#1#4#db",XD + "#2#4#db","bond_gaff_type=sb"}}, {{"1:gaff=sy"}}}, // Special s6 in conjugated systems
		{{{"^S#1#4",".*#2",""}}, {{"1:gaff=s6"}}}, // S with four (or five, or six) connected atoms
		{{{"^S#1#5",".*#2",""}}, {{"1:gaff=s6"}}}, // S with four (or five, or six) connected atoms
		{{{"^S#1#6",".*#2",""}}, {{"1:gaff=s6"}}}, // S with four (or five, or six) connected atoms
		{{{"He#1",".*#2",""}}, {{"1:gaff=He"}}},
		{{{"Li#1",".*#2",""}}, {{"1:gaff=Li"}}},
		{{{"Be#1",".*#2",""}}, {{"1:gaff=Be"}}},
		{{{"B#1",".*#2",""}}, {{"1:gaff=B"}}},
		{{{"Ne#1",".*#2",""}}, {{"1:gaff=Ne"}}},
		{{{"Na#1",".*#2",""}}, {{"1:gaff=Na"}}},
		{{{"Mg#1",".*#2",""}}, {{"1:gaff=Mg"}}},
		{{{"Al#1",".*#2",""}}, {{"1:gaff=Al"}}},
		//~ {{{"Si#1",".*#2",""}}, {{"1:gaff=Si"}}},
		{{{"Ar#1",".*#2",""}}, {{"1:gaff=Ar"}}},
		{{{"K#1",".*#2",""}}, {{"1:gaff=K"}}},
		{{{"Ca#1",".*#2",""}}, {{"1:gaff=Ca"}}},
		{{{"Sc#1",".*#2",""}}, {{"1:gaff=Sc"}}},
		{{{"Ti#1",".*#2",""}}, {{"1:gaff=Ti"}}},
		{{{"V#1",".*#2",""}}, {{"1:gaff=V"}}},
		{{{"Cr#1",".*#2",""}}, {{"1:gaff=Cr"}}},
		{{{"Mn#1",".*#2",""}}, {{"1:gaff=Mn"}}},
		{{{"Fe#1",".*#2",""}}, {{"1:gaff=Fe"}}},
		{{{"Co#1",".*#2",""}}, {{"1:gaff=Co"}}},
		{{{"Ni#1",".*#2",""}}, {{"1:gaff=Ni"}}},
		{{{"Cu#1",".*#2",""}}, {{"1:gaff=Cu"}}},
		{{{"Zn#1",".*#2",""}}, {{"1:gaff=Zn"}}},
		{{{"Ga#1",".*#2",""}}, {{"1:gaff=Ga"}}},
		{{{"Ge#1",".*#2",""}}, {{"1:gaff=Ge"}}},
		{{{"As#1",".*#2",""}}, {{"1:gaff=As"}}},
		{{{"Se#1",".*#2",""}}, {{"1:gaff=Se"}}},
		{{{"Kr#1",".*#2",""}}, {{"1:gaff=Kr"}}},
		{{{"Rb#1",".*#2",""}}, {{"1:gaff=Rb"}}},
		{{{"Sr#1",".*#2",""}}, {{"1:gaff=Sr"}}},
		{{{"Y#1",".*#2",""}}, {{"1:gaff=Y"}}},
		{{{"Zr#1",".*#2",""}}, {{"1:gaff=Zr"}}},
		{{{"Nb#1",".*#2",""}}, {{"1:gaff=Nb"}}},
		{{{"Mo#1",".*#2",""}}, {{"1:gaff=Mo"}}},
		{{{"Tc#1",".*#2",""}}, {{"1:gaff=Tc"}}},
		{{{"Ru#1",".*#2",""}}, {{"1:gaff=Ru"}}},
		{{{"Rh#1",".*#2",""}}, {{"1:gaff=Rh"}}},
		{{{"Pd#1",".*#2",""}}, {{"1:gaff=Pd"}}},
		{{{"Ag#1",".*#2",""}}, {{"1:gaff=Ag"}}},
		{{{"Cd#1",".*#2",""}}, {{"1:gaff=Cd"}}},
		{{{"In#1",".*#2",""}}, {{"1:gaff=In"}}},
		{{{"Sn#1",".*#2",""}}, {{"1:gaff=Sn"}}},
		{{{"Sb#1",".*#2",""}}, {{"1:gaff=Sb"}}},
		{{{"Te#1",".*#2",""}}, {{"1:gaff=Te"}}},
		{{{"Xe#1",".*#2",""}}, {{"1:gaff=Xe"}}},
		{{{"Cs#1",".*#2",""}}, {{"1:gaff=Cs"}}},
		{{{"Ba#1",".*#2",""}}, {{"1:gaff=Ba"}}},
		{{{"La#1",".*#2",""}}, {{"1:gaff=La"}}},
		{{{"Ce#1",".*#2",""}}, {{"1:gaff=Ce"}}},
		{{{"Pr#1",".*#2",""}}, {{"1:gaff=Pr"}}},
		{{{"Nd#1",".*#2",""}}, {{"1:gaff=Nd"}}},
		{{{"Pm#1",".*#2",""}}, {{"1:gaff=Pm"}}},
		{{{"Sm#1",".*#2",""}}, {{"1:gaff=Sm"}}},
		{{{"Eu#1",".*#2",""}}, {{"1:gaff=Eu"}}},
		{{{"Gd#1",".*#2",""}}, {{"1:gaff=Gd"}}},
		{{{"Tb#1",".*#2",""}}, {{"1:gaff=Tb"}}},
		{{{"Dy#1",".*#2",""}}, {{"1:gaff=Dy"}}},
		{{{"Ho#1",".*#2",""}}, {{"1:gaff=Ho"}}},
		{{{"Er#1",".*#2",""}}, {{"1:gaff=Er"}}},
		{{{"Tm#1",".*#2",""}}, {{"1:gaff=Tm"}}},
		{{{"Yb#1",".*#2",""}}, {{"1:gaff=Yb"}}},
		{{{"Lu#1",".*#2",""}}, {{"1:gaff=Lu"}}},
		{{{"Hf#1",".*#2",""}}, {{"1:gaff=Hf"}}},
		{{{"Ta#1",".*#2",""}}, {{"1:gaff=Ta"}}},
		{{{"W#1",".*#2",""}}, {{"1:gaff=W"}}},
		{{{"Re#1",".*#2",""}}, {{"1:gaff=Re"}}},
		{{{"Os#1",".*#2",""}}, {{"1:gaff=Os"}}},
		{{{"Ir#1",".*#2",""}}, {{"1:gaff=Ir"}}},
		{{{"Pt#1",".*#2",""}}, {{"1:gaff=Pt"}}},
		{{{"Au#1",".*#2",""}}, {{"1:gaff=Au"}}},
		{{{"Hg#1",".*#2",""}}, {{"1:gaff=Hg"}}},
		{{{"Tl#1",".*#2",""}}, {{"1:gaff=Tl"}}},
		{{{"Pb#1",".*#2",""}}, {{"1:gaff=Pb"}}},
		{{{"Bi#1",".*#2",""}}, {{"1:gaff=Bi"}}},
		{{{"Po#1",".*#2",""}}, {{"1:gaff=Po"}}},
		{{{"At#1",".*#2",""}}, {{"1:gaff=At"}}},
		{{{"Rn#1",".*#2",""}}, {{"1:gaff=Rn"}}},
		{{{"Fr#1",".*#2",""}}, {{"1:gaff=Fr"}}},
		{{{"Ra#1",".*#2",""}}, {{"1:gaff=Ra"}}},
		{{{"Ac#1",".*#2",""}}, {{"1:gaff=Ac"}}},
		{{{"Th#1",".*#2",""}}, {{"1:gaff=Th"}}},
		{{{"Pa#1",".*#2",""}}, {{"1:gaff=Pa"}}},
		{{{"U#1",".*#2",""}}, {{"1:gaff=U"}}},
		{{{"Np#1",".*#2",""}}, {{"1:gaff=Np"}}},
		{{{"Pu#1",".*#2",""}}, {{"1:gaff=Pu"}}},
		{{{"Am#1",".*#2",""}}, {{"1:gaff=Am"}}},
		{{{"Cm#1",".*#2",""}}, {{"1:gaff=Cm"}}},
		{{{"Bk#1",".*#2",""}}, {{"1:gaff=Bk"}}},
		{{{"Cf#1",".*#2",""}}, {{"1:gaff=Cf"}}},
		{{{"Es#1",".*#2",""}}, {{"1:gaff=Es"}}},
		{{{"Fm#1",".*#2",""}}, {{"1:gaff=Fm"}}},
		{{{"Md#1",".*#2",""}}, {{"1:gaff=Md"}}},
		{{{"No#1",".*#2",""}}, {{"1:gaff=No"}}},
		{{{"Lr#1",".*#2",""}}, {{"1:gaff=Lr"}}},
		{{{"Rf#1",".*#2",""}}, {{"1:gaff=Rf"}}},
		{{{"Db#1",".*#2",""}}, {{"1:gaff=Db"}}},
		{{{"Sg#1",".*#2",""}}, {{"1:gaff=Sg"}}},
		{{{"Bh#1",".*#2",""}}, {{"1:gaff=Bh"}}},
		{{{"Hs#1",".*#2",""}}, {{"1:gaff=Hs"}}},
		{{{"Mt#1",".*#2",""}}, {{"1:gaff=Mt"}}},
		{{{"Ds#1",".*#2",""}}, {{"1:gaff=Ds"}}},
		{{{"LP#1",".*#2",""}}, {{"1:gaff=LP"}}},
		{{{"lp#1",".*#2",""}}, {{"1:gaff=lp"}}},
		{{{"DU#1",".*#2",""}}, {{"1:gaff=DU"}}},
	};
	
	const map<const string, const string> gaff_flip {
		{"cc", "cd"},
		{"ce", "cf"},
		{"nc", "nd"},
		{"ne", "nf"},
		{"pe", "pf"},
		{"cp", "cq"},
	};
	const set<string> gaff_group_1 {
		{"cc"},
		{"ce"},
		{"nc"},
		{"ne"},
		{"pe"},
		{"cp"},
	};
	const set<string> gaff_group_2 {
		{"cd"},
		{"cf"},
		{"nd"},
		{"nf"},
		{"pf"},
		{"cq"},
	};

	const map<string, vector<string>> gaff_replacement {
		{"cl",{"cl"}},
		{"Si",{"Si"}},
		{"cx",{"cx","c3"}},
		{"cy",{"cy","c3"}},
		{"c3",{"c3"}},
		{"c",{"c","c2"}},
		{"cz",{"cz"}},
		{"cp",{"cp","ca","cc","cd","c2"}},
		{"cq",{"cq","ca","cc","cd","c2"}}, // fix issue #118
		{"ca",{"ca","cp","cc","cd","c2"}},
		{"cc",{"cc","ca","c2"}},
		{"cd",{"cd","ca","c2"}},
		{"ce",{"ce","c2"}},
		{"cf",{"cf","c2"}},
		{"cu",{"cu","c2"}},
		{"cv",{"cv","c2"}},
		{"c2",{"c2","ce","cf"}},
		{"cg",{"cg","c1"}},
		{"c1",{"c1","cg"}},
		{"hn",{"hn"}},
		{"ho",{"ho"}},
		{"hs",{"hs"}},
		{"hp",{"hp"}},
		{"hx",{"hx"}},
		{"hw",{"hw"}},
		{"h3",{"h3"}},
		{"h2",{"h2"}},
		{"h1",{"h1"}},
		{"hc",{"hc"}},
		{"h5",{"h5"}},
		{"h4",{"h4"}},
		{"ha",{"ha"}},
		{"f",{"f"}},
		{"br",{"br"}},
		{"i",{"i"}},
		{"pc",{"pc","pb","p2"}},
		{"pd",{"pd","pb","p2"}},
		{"pe",{"pe","p2"}},
		{"pf",{"pf","p2"}},
		{"px",{"px","p4","p3"}},
		{"p4",{"p4","px","p3"}},
		{"p3",{"p3","px","p4"}},
		{"py",{"py","p5"}},
		{"p5",{"p5","py"}},
		{"n",{"n","n3"}},
		{"n4",{"n4","n3"}},
		{"no",{"no","n3"}},
		{"na",{"na","n3"}},
		{"nh",{"nh","n3"}},
		{"n",{"n","n3"}},
		{"n3",{"n3","nh"}},
		{"nb",{"nb","nc","nd","n2"}},
		{"nc",{"nc","nb","n2"}},
		{"nd",{"nd","nb","n2"}},
		{"ne",{"ne","n2"}},
		{"nf",{"nf","n2"}},
		{"n1",{"n1"}},
		{"n2",{"n2"}},
		{"o",{"o"}},
		{"oh",{"oh"}},
		{"os",{"os"}},
		{"s",{"s"}},
		{"s2",{"s2","sh"}},
		{"sh",{"sh","s2"}},
		{"ss",{"ss","s2"}},
		{"sx",{"sx","s4"}},
		{"s4",{"s4","sx"}},
		{"sy",{"sy","s6"}},
		{"s6",{"s6","sy"}},
	};

}
#endif
