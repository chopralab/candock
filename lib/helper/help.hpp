#ifndef HELP_H
#define HELP_H
#include <string>
#include <set>
#include <sstream>
#include <vector>
#include <memory>
#include <map>
#include <tuple>
#include "error.hpp"
#include "debug.hpp"
using namespace std;

namespace help {
	string memusage(const string&);
	std::tuple<double, double, double> gnuplot(const double &x1, const double &x2, const string &datapoints);
	
	struct edge { 
		string atom_property1;
		string atom_property2;
		string bond_property;
	};
	typedef vector<edge> smiles;
	ostream& operator<<(ostream& os, const smiles& edges);
	struct rename_rule {
		smiles pattern;
		vector<string> rule;
	};
	typedef vector<rename_rule> rename_rules;

	ostream& operator<<(ostream& os, const rename_rule& rule);
	
	enum IdatmGeometry {
		Ion=0, Single=1, Linear=2, Planar=3, Tetrahedral=4, TrigonalBipyramidal=5
	};
	struct IdatmEntry {
		IdatmGeometry geometry;
		int substituents;
		string description;
	};
	typedef map<const string, const IdatmEntry> IdatmInfoMap;

	template<class Set1, class Set2> 
	bool is_disjoint(const Set1 &set1, const Set2 &set2) {
		if(set1.empty() || set2.empty()) return true;
		typename Set1::const_iterator it1 = set1.begin(), 
			it1End = set1.end();
		typename Set2::const_iterator it2 = set2.begin(), 
			it2End = set2.end();
		if(*it1 > *set2.rbegin() || *it2 > *set1.rbegin()) return true;
		while(it1 != it1End && it2 != it2End) {
			if(*it1 == *it2) return false;
			if(*it1 < *it2) { it1++; }
			else { it2++; }
		}
		return true;
	}

	template <typename T>
	string to_string ( T num ) {
		ostringstream ss;
		ss << num;
		return ss.str();
	}
	static string replace_str_char(string str, const string& replace, char ch) {
		size_t found = str.find_first_of(replace);
		while (found != string::npos) { // While our position in the sting is in range.
			str[found] = ch; // Change the character at position.
			found = str.find_first_of(replace, found+1); // Relocate again.
		}
		return str; // return our new string.
	}
	static vector<string> ssplit(const string &source, const char*delimiter = " ", bool keepEmpty = false) {
		vector<string> results;
		size_t prev = 0;
		size_t next = 0;
		while ((next = source.find_first_of(delimiter, prev)) != string::npos){
			if (keepEmpty || (next - prev != 0)) {
				results.push_back(source.substr(prev, next - prev));
			}
			prev = next + 1;
		}
		if (prev < source.size()) {
			results.push_back(source.substr(prev));
		}
		return results;
	}
	
	const map<const string, string> sybyl {
		{"C.3", "C"},
		{"C.2", "C"},
		{"C.ar", "C"},
		{"C.1", "C"},
		{"N.3", "N"},
		{"N.2", "N"},
		{"N.1", "N"},
		{"O.3", "O"},
		{"O.2", "O"},
		{"S.3", "S"},
		{"N.ar", "N"},
		{"P.3", "P"},
		{"H", "H"},
		{"Br", "Br"},
		{"Cl", "Cl"},
		{"F", "F"},
		{"I", "I"},
		{"S.2", "S"},
		{"N.pl3", "N"},
		{"LP", "LP"},
		{"Na", "Na"},
		{"K", "K"},
		{"Ca", "Ca"},
		{"Li", "Li"},
		{"Al", "Al"},
		{"Du", "Du"},
		{"Du.C", "Du"},
		{"Si", "Si"},
		{"N.am", "N"},
		{"S.o", "S"},
		{"S.o2", "S"},
		{"N.4", "N"},
		{"O.co2", "O"},
		{"C.cat", "C"},
		{"H.spc", "H"},
		{"O.spc", "O"},
		{"H.t3p", "H"},
		{"O.t3p", "O"},
		{"ANY", "ANY"},
		{"HEV", "HEV"},
		{"HET", "HET"},
		{"HAL", "HAL"},
		{"Mg", "Mg"},
		{"Cr.oh", "Cr"},
		{"Cr.th", "Cr"},
		{"Se", "Se"},
		{"Fe", "Fe"},
		{"Cu", "Cu"},
		{"Zn", "Zn"},
		{"Sn", "Sn"},
		{"Mo", "Mo"},
		{"Mn", "Mn"},
		{"Co.oh", "Co"}
	};
		
	const map<string, const char> one_letter {
		{"ALA", 'A'},
		{"CYS", 'C'},
		{"ASP", 'D'},
		{"GLU", 'E'},
		{"PHE", 'F'},
		{"GLY", 'G'},
		{"HIS", 'H'},
		{"ILE", 'I'},
		{"LYS", 'K'},
		{"LEU", 'L'},
		{"MET", 'M'},
		{"ASN", 'N'},
		{"PRO", 'P'},
		{"GLN", 'Q'},
		{"ARG", 'R'},
		{"SER", 'S'},
		{"THR", 'T'},
		{"VAL", 'V'},
		{"TRP", 'W'},
		{"TYR", 'Y'}
	}; 
	static char get_one_letter(string s) { 
		auto it = one_letter.find(s); 
		if (it == one_letter.end()) { 
			dbgmsg("warn : residues name " << s << 
				" not found in to_one_letter..."); 
			return ' '; 
		} 
		return it->second; 
	}
	const set<string> amino_acids {
		{"PPP"}, // artificial for geo
		{"ALA"},
		{"ARG"},
		{"ASN"},
		{"ASP"},
		{"CYS"},
		{"GLN"},
		{"GLU"},
		{"GLY"},
		{"HIS"},
		{"ILE"},
		{"LEU"},
		{"LYS"},
		{"MET"},
		{"PHE"},
		{"PRO"},
		{"SER"},
		{"THR"},
		{"TRP"},
		{"TYR"},
		{"VAL"}
	};
	const set<string> nucleic_acids {
		{"NNN"}, // artificial for geo
		{"A"},
		{"T"},
		{"G"},
		{"C"},
		{"U"},
		{"DA"},
		{"DT"},
		{"DG"},
		{"DC"},
		{"DU"}
	};
	const set<string> ions {
		{"LI"},
		{"BE"},
		{"NA"},
		{"MG"},
		{"AL"},
		{"K"},
		{"CA"},
		{"CR"},
		{"MN"},
		{"FE"},
		{"FE2"},
		{"CO"},
		{"NI"},
		{"CU"},
		{"CU1"},
		{"ZN"},
		{"PD"},
		{"PB"},
		{"RU"},
		{"AG"},
		{"CD"},
		{"SN"},
		{"SB"},
		{"SR"},
		{"CS"},
		{"BA"},
		{"PT"},
		{"AU"},
		{"HG"},
		{"UNK"},
		{"BR"},
		{"CL"},
		{"F"},
		{"IOD"},
		{"D8U"}
	};
	const map<const string, const char> heavy_atoms {
		{"N", 'N'},
		{"CA", 'C'},
		{"C", 'C'},
		{"O", 'O'},
		{"CB", 'C'},
		{"CG", 'C'},
		{"CG1", 'C'},
		{"CG2", 'C'},
		{"CD", 'C'},
		{"CD1", 'C'},
		{"CD2", 'C'},
		{"CE", 'C'},
		{"CE1", 'C'},
		{"CE2", 'C'},
		{"CE3", 'C'},
		{"CZ", 'C'},
		{"CZ2", 'C'},
		{"CZ3", 'C'},
		{"CH2", 'C'},
		{"ND1", 'N'},
		{"ND2", 'N'},
		{"NE", 'N'},
		{"NE1", 'N'},
		{"NE2", 'N'},
		{"NZ", 'N'},
		{"NH1", 'N'},
		{"NH2", 'N'},
		{"OG", 'O'},
		{"OG1", 'O'},
		{"OD1", 'O'},
		{"OD2", 'O'},
		{"OE1", 'O'},
		{"OE11", 'O'},
		{"OE12", 'O'},
		{"OE21", 'O'},
		{"OE22", 'O'},
		{"OE2", 'O'},
		{"OH", 'O'},
		{"OXT", 'O'},
		{"SD", 'S'},
		{"SG", 'S'},
		{"C2", 'C'},
		{"C4", 'C'},
		{"C5", 'C'},
		{"C5M", 'C'},
		{"C6", 'C'},
		{"C8", 'C'},
		{"N1", 'N'},
		{"N2", 'N'},
		{"N3", 'N'},
		{"N4", 'N'},
		{"N6", 'N'},
		{"N7", 'N'},
		{"N9", 'N'},
		{"O1P", 'O'},
		{"O2", 'O'},
		{"O2P", 'O'},
		{"O4", 'O'},
		{"O6", 'O'},
		{"P", 'P'},
		{"C1*", 'C'},
		{"C2*", 'C'},
		{"C3*", 'C'},
		{"C4*", 'C'},
		{"C5*", 'C'},
		{"O3*", 'O'},
		{"O4*", 'O'},
		{"O5*", 'O'},
		{"C1\'", 'C'},
		{"C2\'", 'C'},
		{"C3\'", 'C'},
		{"C4\'", 'C'},
		{"C5\'", 'C'},
		{"C7", 'C'},
		{"O3\'", 'O'},
		{"O4\'", 'O'},
		{"O5\'", 'O'},
		{"OP1", 'O'},
		{"OP2", 'O'}
	};
	const map<const string, const char> protein_hydrogen_atoms {
		{"H", 'H'},
		{"HA", 'H'},
		{"HA1", 'H'},
		{"HA2", 'H'},
		{"HB", 'H'},
		{"HB1", 'H'},
		{"HB2", 'H'},
		{"HB3", 'H'},
		{"HD1", 'H'},
		{"HD11", 'H'},
		{"HD12", 'H'},
		{"HD13", 'H'},
		{"HD2", 'H'},
		{"HD21", 'H'},
		{"HD22", 'H'},
		{"HD23", 'H'},
		{"HE", 'H'},
		{"HE1", 'H'},
		{"HE2", 'H'},
		{"HE21", 'H'},
		{"HE22", 'H'},
		{"HE3", 'H'},
		{"HG", 'H'},
		{"HG1", 'H'},
		{"HG11", 'H'},
		{"HG12", 'H'},
		{"HG13", 'H'},
		{"HG2", 'H'},
		{"HG21", 'H'},
		{"HG22", 'H'},
		{"HG23", 'H'},
		{"HH", 'H'},
		{"HH11", 'H'},
		{"HH12", 'H'},
		{"HH2", 'H'},
		{"HH21", 'H'},
		{"HH22", 'H'},
		{"HOCA", 'H'},
		{"HZ", 'H'},
		{"HZ1", 'H'},
		{"HZ2", 'H'},
		{"HZ3", 'H'},
		{"1HZ", 'H'},
		{"2HZ", 'H'},
		{"3HZ", 'H'},
		{"1HD2", 'H'},
		{"2HD2", 'H'},
		{"1HE2", 'H'},
		{"2HE2", 'H'},
		{"1HH1", 'H'},
		{"2HH1", 'H'},
		{"1HH2", 'H'},
		{"2HH2", 'H'}
	};
	const map<const string, const char> dna_hydrogen_atoms {
		{"D1", 'D'},
		{"D21", 'D'},
		{"D22", 'D'},
		{"D3", 'D'},
		{"D41", 'D'},
		{"D42", 'D'},
		{"D61", 'D'},
		{"D62", 'D'},
		{"D8", 'D'},
		{"DO3\'", 'D'},
		{"DO5\'", 'D'},
		{"H1\'", 'H'},
		{"H2", 'H'},
		{"H2\'", 'H'},
		{"H2\'\'", 'H'},
		{"H3\'", 'H'},
		{"H4\'", 'H'},
		{"H5", 'H'},
		{"H5\'", 'H'},
		{"H5\'\'", 'H'},
		{"H6", 'H'},
		{"H71", 'H'},
		{"H72", 'H'},
		{"H73", 'H'},
		{"H8", 'H'}
	};
	const set<string> metals {
		{"Li"},
		{"Be"},
		{"B"}, 
		{"Na"},
		{"Mg"},
		{"Al"},
		{"Si"},
		{"K"}, 
		{"Ca"},
		{"Sc"},
		{"Ti"},
		{"V"}, 
		{"Cr"},
		{"Mn"},
		{"Fe"},
		{"Co"},
		{"Ni"},
		{"Cu"},
		{"Zn"},
		{"Ga"},
		{"Ge"},
		{"As"},
		{"Rb"},
		{"Sr"},
		{"Y"}, 
		{"Zr"},
		{"Nb"},
		{"Mo"},
		{"Tc"},
		{"Ru"},
		{"Rh"},
		{"Pd"},
		{"Ag"},
		{"Cd"},
		{"In"},
		{"Sn"},
		{"Sb"},
		{"Te"},
		{"Cs"},
		{"Ba"},
		{"La"},
		{"Ce"},
		{"Pr"},
		{"Nd"},
		{"Pm"},
		{"Sm"},
		{"Eu"},
		{"Gd"},
		{"Tb"},
		{"Dy"},
		{"Ho"},
		{"Er"},
		{"Tm"},
		{"Yb"},
		{"Lu"},
		{"Hf"},
		{"Ta"},
		{"W"}, 
		{"Re"},
		{"Os"},
		{"Ir"},
		{"Pt"},
		{"Au"},
		{"Hg"},
		{"Tl"},
		{"Pb"},
		{"Bi"},
		{"Po"},
		{"Fr"},
		{"Ra"},
		{"Ac"},
		{"Th"},
		{"Pa"},
		{"U"}, 
		{"Np"},
		{"Pu"},
		{"Am"},
		{"Cm"},
		{"Bk"},
		{"Cf"},
		{"Es"},
		{"Fm"},
		{"Md"},
		{"No"},
		{"Lr"}
	};
	const map<const string, const pair<string,string>> non_specific_binders {
		{"BE7", {"(4-CARBOXYPHENYL)(CHLORO)MERCURY", "C7 H5 O2 CL1 HG1"}},
		{"MRD", {"(4R)-2-METHYLPENTANE-2,4-DIOL", "C6 H14 O2"}},
		{"MHA", {"(CARBAMOYLMETHYL-CARBOXYMETHYL-AMINO)-ACETIC ACID", "C6 H10 N2 O5"}},
		{"BU3", {"(R,R)-2,3-BUTANEDIOL", "C4 H10 O2"}},
		{"EDO", {"1,2-ETHANEDIOL", "C2 H6 O2"}},
		{"PGO", {"1,2-PROPANEDIOL", "C3 H8 O2"}},
		{"BU2", {"1,3-BUTANEDIOL", "C4 H10 O2"}},
		{"PDO", {"1,3-PROPANEDIOL", "C3 H8 O2"}},
		{"BU1", {"1,4-BUTANEDIOL", "C4 H10 O2"}},
		{"PG6", {"1-(2-METHOXY-ETHOXY)-2-{2-[2-(2-METHOXY-ETHOXY]-ETHOXY}-ETHANE", "C12 H26 O6"}},
		{"1BO", {"1-BUTANOL", "C4 H10 O1"}},
		{"PE7", {"1-DEOXY-1-THIO-HEPTAETHYLENE GLYCOL", "C14 H30 O7 S1"}},
		{"PG5", {"1-METHOXY-2-[2-(2-METHOXY-ETHOXY]-ETHANE", "C8 H18 O4"}},
		{"TFP", {"10-[3-(4-METHYL-PIPERAZIN-1-YL)-PROPYL]-2-TRIFLUOROMETHYL-10H-PHENOTHIAZINE", "C21 H24 N3 F3 S1"}},
		{"DHD", {"2,4-DIHYDROXY-5-OXO-HEXA-2,3-DIENOIC ACID", "C5 H4 O6"}},
		{"PEU", {"2,5,8,11,14,17,20,23,26,29,32,35,38,41,44,47,50,53,56,59,62,65,68,71,74,77,80-HEPTACOSAOXADOOCTACONTAN-82-OL", "C55 H112 O28"}},
		{"TRS", {"2-AMINO-2-HYDROXYMETHYL-PROPANE-1,3-DIOL", "C4 H12 N1 O3"}},
		{"TAU", {"2-AMINOETHANESULFONIC ACID", "C2 H7 N1 O3 S1"}},
		{"SBT", {"2-BUTANOL", "C4 H10 O1"}},
		{"SAL", {"2-HYDROXYBENZOIC ACID", "C7 H6 O3"}},
		{"MPD", {"2-METHYL-2,4-PENTANEDIOL", "C6 H14 O2"}},
		{"IOH", {"2-PROPANOL, ISOPROPANOL", "C3 H8 O1"}},
		{"IPA", {"2-PROPANOL, ISOPROPANOL", "C3 H8 O1"}},
		{"PGE", {"2-[2-(2-HYDROXY-ETHOXY)-ETHOXY]-ETHANOL", "C6 H14 O4"}},
		{"PIG", {"2-[2-(2-HYDROXY-ETHOXY)-ETHOXY]-ETHANOL", "C6 H14 O4"}},
		{"B3P", {"2-[3-(2-HYDROXY-1,1-DIHYDROXYMETHYL-ETHYLAMINO)-PROPYLAMINO]-2-HYDROXYMETHYL-PROPANE-1,3-DIOL", "C11 H26 N2 O6"}},
		{"BTB", {"2-[BIS-(2-HYDROXY-ETHYL)-AMINO]-2-HYDROXYMETHYL- PROPANE-1,3-DIOL", "C8 H19 N1 O5"}},
		{"NHE", {"2-[N-CYCLOHEXYLAMINO]ETHANE SULFONIC ACID", "C8 H17 N1 O3 S1"}},
		{"C8E", {"2-{2-[2-(2-OCTYLOXYETHOXY)-ETHOXYL]-ETHOXY}ETHANOL", "C16 H34 O5"}},
		{"OTE", {"2-{2-[2-(2-OCTYLOXYETHOXY)-ETHOXYL]-ETHOXY}ETHANOL", "C16 H34 O5"}},
		{"PE4", {"2-{2-[2-(2-{2-[2-(2-ETHOXY-ETHOXY)-ETHOXY]-ETHOXY}-ETHOXY)-ETHOXY]-ETHOXY}-ETHANOL", "C16 H34 O8"}},
		{"XPE", {"3,6,9,12,15,18,21,24,27-NONAOXANONACOSANE-1,29-DIOL", "C20 H42 O11"}},
		{"PE8", {"3,6,9,12,15,18,21-HEPTAOXATRICOSANE-1,23-DIOL", "C16 H34 O9"}},
		{"P33", {"3,6,9,12,15,18-HEXAOXAICOSANE-1,20-DIOL", "C14 H30 O8"}},
		{"N8E", {"3,6,9,12,15-PENTAOXATRICOSAN-1-OL", "C18 H38 O6"}},
		{"2OS", {"3-N-OCTANOYLSUCROSE", "C20 H36 O12"}},
		{"1PS", {"3-PYRIDINIUM-1-YLPROPANE-1-SULFONATE", "C8 H11 N1 O3 S1"}},
		{"CPS", {"3-[(3-CHOLAMIDOPROPYL)DIMETHYLAMMONIO]-1-PROPANESULFONATE", "C32 H58 N2 O7 S1"}},
		{"DMX", {"3-[BENZYL(DIMETHYL)AMMONIO]PROPANE-1-SULFONATE", "C12 H19 N1 O3 S1"}},
		{"MPO", {"3[N-MORPHOLINE]PROPANE SULFONIC ACID", "C7 H15 N1 O4 S1"}},
		{"GCD", {"4,5-DEHYDRO-D-GLUCURONIC ACID", "C6 H8 O7"}},
		{"IDT", {"4,5-DEHYDRO-L-IDURONIC ACID", "C6 H8 O7"}},
		{"DXG", {"4-DEOXYGLUCARATE", "C6 H8 O7"}},
		{"CM5", {"5-CYCLOHEXYL-1-PENTYL-BETA-D-MALTOSIDE", "C23 H42 O11"}},
		{"ACA", {"6-AMINOHEXANOIC ACID", "C6 H13 N1 O2"}},
		{"ACT", {"ACETATE ION", "C2 H3 O2"}},
		{"ACN", {"ACETONE", "C3 H6 O1"}},
		{"CCN", {"ACETONITRILE", "C2 H3 N1"}},
		{"AGC", {"ALPHA-D-GLUCOSE", "C6 H12 O6"}},
		{"GLC", {"ALPHA-D-GLUCOSE", "C6 H12 O6"}},
		{"MAN", {"ALPHA-D-MANNOSE", "C6 H12 O6"}},
		{"DR6", {"ALPHA-[4-(1,1,3,3 - TETRAMETHYLBUTYL)PHENYL]- OMEGA-HYDROXY-POLY(OXY-1,2-ETHANEDIYL", "C74 H142 O31"}},
		{"AL", {"ALUMINUM ION", "AL1"}},
		{"NH4", {"AMMONIUM ION", "H4 N1"}},
		{"AZI", {"AZIDE ION", "N3"}},
		{"BNG", {"B-NONYLGLUCOSIDE", "C15H30O6"}},
		{"BOG", {"B-OCTYLGLUCOSIDE", "C14 H28 O6"}},
		{"BA", {"BARIUM ION", "BA1"}},
		{"LAK", {"BETA-D-GALACTOPYRANOSYL-1-6-BETA-D-GLUCOPYRANOSE", "C12 H22 O11"}},
		{"BGC", {"BETA-D-GLUCOSE", "C6 H12 O6"}},
		{"BMA", {"BETA-D-MANNOSE", "C6 H12 O6"}},
		{"BCN", {"BICINE", "C6 H13 N1 O4"}},
		{"BR", {"BROMIDE ION", "BR1"}},
		{"BR", {"BROMO GROUP", "BR1"}},
		{"BRO", {"BROMO GROUP", "BR1"}},
		{"CAC", {"CACODYLATE ION", "C2 H6 O2 AS1 1-"}},
		{"CD", {"CADMIUM ION", "CD1"}},
		{"CA", {"CALCIUM ION", "CA1"}},
		{"CBX", {"CARBOXY GROUP", "C1 H1 O2"}},
		{"FMT", {"CARBOXY GROUP", "C1 H1 O2"}},
		{"ACY", {"CARBOXYMETHYL GROUP", "C2 H3 O2"}},
		{"CBM", {"CARBOXYMETHYL GROUP", "C2 H3 O2"}},
		{"CM", {"CARBOXYMETHYL GROUP", "C2 H3 O2"}},
		{"CS", {"CESIUM ION", "CS1"}},
		{"CL", {"CHLORIDE ION", "CL1"}},
		{"CL", {"CHLORO GROUP", "CL1"}},
		{"CLO", {"CHLORO GROUP", "CL1"}},
		{"FCL", {"CITRATE ANION", "C6 H5 O7 3-"}},
		{"CIT", {"CITRIC ACID", "C6 H8 O7"}},
		{"CO", {"COBALT (II) ION", "CO1"}},
		{"3CO", {"COBALT (III) ION", "CO1"}},
		{"NCO", {"COBALT HEXAMMINE ION", "H18 N6 CO1"}},
		{"CU1", {"COPPER (I) ION", "CU1"}},
		{"CU", {"COPPER (II) ION", "CU1"}},
		{"CN", {"CYANIDE GROUP", "C1 N1"}},
		{"CYN", {"CYANIDE GROUP", "C1 N1"}},
		{"CYN", {"CYANIDE ION", "C1 N1"}},
		{"MA4", {"CYCLOHEXYL-HEXYL-BETA-D-MALTOSIDE", "C24 H44 O11"}},
		{"BTC", {"CYSTEINE", "C3 H7 N1 O2 S1"}},
		{"CYS", {"CYSTEINE", "C3 H7 N1 O2 S1"}},
		{"TAR", {"D(-)-TARTARIC ACID", "C4 H6 O6"}},
		{"GLO", {"D-GLUCOSE IN LINEAR FORM", "C6 H12 O6"}},
		{"MTL", {"D-MANNITOL", "C6 H14 O6"}},
		{"DPR", {"D-PROLINE", "C5 H9 N1 O2"}},
		{"SOR", {"D-SORBITOL", "C6 H14 O6"}},
		{"SYL", {"D-XYLITOL", "C7 H16 O5"}},
		{"DMU", {"DECYL-BETA-D-MALTOPYRANOSIDE", "C22 H42 O11"}},
		{"DDQ", {"DECYLAMINE-N,N-DIMETHYL-N-OXIDE", "C12 H27 N1 O1"}},
		{"DMS", {"DIMETHYL SULFOXIDE", "C2 H6 O1 S1"}},
		{"DMF", {"DIMETHYLFORMAMIDE", "C3 H7 N1 O1"}},
		{"DIO", {"DIOXANE", "C4 H8 O2"}},
		{"DOX", {"DIOXANE", "C4 H8 O2"}},
		{"12P", {"DODECAETHYLENE GLYCOL", "C24 H50 O13"}},
		{"SDS", {"DODECYL SULFATE", "C12 H26 O4 S1"}},
		{"LMT", {"DODECYL-BETA-D-MALTOSIDE", "C24 H46 O11"}},
		{"EOH", {"ETHANOL", "C2 H6 O1"}},
		{"EEE", {"ETHYL ACETATE", "C4 H8 O2"}},
		{"EDO", {"ETHYLENE GLYCOL", "C2 H6 O2"}},
		{"EGL", {"ETHYLENE GLYCOL", "C2 H6 O2"}},
		{"FE2", {"FE (II) ION", "FE1"}},
		{"FE", {"FE (III) ION", "FE1"}},
		{"F", {"FLUORIDE ION", "F1"}},
		{"F", {"FLUORO GROUP", "F1"}},
		{"FLO", {"FLUORO GROUP", "F1"}},
		{"TRT", {"FRAGMENT OF TRITON X-100", "C21 H36 O4"}},
		{"CYS", {"FREE CYSTEINE", "C3 H7 N1 O2 S1"}},
		{"FCY", {"FREE CYSTEINE", "C3 H7 N1 O2 S1"}},
		{"FRU", {"FRUCTOSE", " C6 H12 O6"}},
		{"GBL", {"GAMMA-BUTYROLACTONE", "C4 H6 O2"}},
		{"GLC", {"GLUCOSE", "C6 H12 O6"}},
		{"GOL", {"GLYCEROL", "C3 H8 O3"}},
		{"GLY", {"GLYCINE", "C2 H5 N1 O2"}},
		{"GPX", {"GUANOSINE 5'-DIPHOSPHATE 2':3'-CYCLIC MONOPHOSPHATE", "C10 H14 N5 O13 P3"}},
		{"HTO", {"HEPTANE-1,2,3-TRIOL", "C7 H16 O3"}},
		{"HTG", {"HEPTYL 1-THIOHEXOPYRANOSIDE", "C13 H26 O5 S1"}},
		{"B7G", {"HEPTYL-BETA-D-GLUCOPYRANOSIDE", "C13 H26 O6"}},
		{"C10", {"HEXAETHYLENE GLYCOL MONODECYL ETHER", "C22 H46 O7"}},
		{"16D", {"HEXANE-1,6-DIAMINE", "C6 H16 N2"}},
		{"HEZ", {"HEXANE-1,6-DIOL", "C6 H14 O2"}},
		{"IOD", {"IODIDE ION", "I1"}},
		{"IDO", {"IODO GROUP", "I1"}},
		{"IOD", {"IODO GROUP", "I1"}},
		{"ICI", {"ISOCITRIC ACID", "C6 H8 O7"}},
		{"ICT", {"ISOCITRIC ACID", "C6 H8 O7"}},
		{"IPA", {"ISOPROPYL ALCOHOL", "C3 H8 01"}},
		{"TLA", {"L(+)-TARTARIC ACID", "C4 H6 O6"}},
		{"LAT", {"LACTOSE", "C12 H22 O11"}},
		{"LBT", {"LACTOSE", "C12 H22 O11"}},
		{"LDA", {"LAURYL DIMETHYLAMINE-N-OXIDE", "C14 H31 N1 O1"}},
		{"PB", {"LEAD (II) ION", "PB1"}},
		{"LI", {"LITHIUM ION", "LI1"}},
		{"MG", {"MAGNESIUM ION", "MG1"}},
		{"MN", {"MANGANESE (II) ION", "MN1"}},
		{"MN3", {"MANGANESE (III) ION", "MN1"}},
		{"HG", {"MERCURY (II) ION", "HG1"}},
		{"MRY", {"MESO-ERYTRHITOL", "C4 H10 O4"}},
		{"MOH", {"METHANOL", "C1 H4 O1"}},
		{"BEQ", {"N-(CARBOXYMETHYL)-N,N-DIMETHYL-3-[(1-OXODODECYL)AMINO]-1-PROPANAMINIUM INNER SALT", "C19 H38 N2 O3"}},
		{"C15", {"N-DODECYL-N,N-DIMETHYL-3-AMMONIO-1-PROPANESULFONATE", "C17 H38 N1 O3 S1"}}, 
		{"MG8", {"N-OCTANOYL-N-METHYLGLUCAMINE", "C15 H31 N1 O6"}},
		{"POL", {"N-PROPANOL", "C3 H8 O1"}},
		{"NI", {"NICKEL (II) ION", "NI1"}},
		{"3NI", {"NICKEL (III) ION", "NI1"}},
		{"NO3", {"NITRATE ION", "N1 O3"}},
		{"JEF", {"O-(O-(2-AMINOPROPYL)-O'-(2-METHOXYETHYL)POLYPROPYLENEGLYCOL 500)", "C31 H65 N1 O10"}},
		{"P4C", {"O-ACETALDEHYDYL-HEXAETHYLENE GLYCOL", "C14 H28 O8"}},
		{"CE1", {"O-DODECANYL OCTAETHYLENE GLYCOL", "C28 H58 O9"}},
		{"DIA", {"OCTANE 1,8-DIAMINE", "C8 H20 N2"}},
		{"CXE", {"PENTAETHYLENE GLYCOL MONODECYL ETHER", "C20 H42 O6"}},
		{"IPH", {"PHENOL", "C6 H6 O1"}},
		{"PIN", {"PIPERAZINE-N,N'-BIS(2-ETHANESULFONIC ACID)", "C8 H18 N2 O6 S2"}},
		{"15P", {"POLYETHYLENE GLYCOL (N=34)", "C69 H140 O35"}},
		{"K", {"POTASSIUM ION", "K1"}},
		{"CRY", {"PROPANE-1,2,3-TRIOL", "C3 H8 O3"}},
		{"GOL", {"PROPANE-1,2,3-TRIOL", "C3 H8 O3"}},
		{"PGR", {"R-1,2-PROPANEDIOL", "C3 H8 O2"}},
		{"RB", {"RUBIDIUM ION", "RB1"}},
		{"PGO", {"S-1,2-PROPANEDIOL", "C3 H8 O2"}},
		{"PGQ", {"S-1,2-PROPANEDIOL", "C3 H8 O2"}},
		{"AG", {"SILVER ION", "AG1 1+"}},
		{"NA", {"SODIUM ION", "NA1"}},
		{"SPD", {"SPERMIDINE", "C7 H19 N3"}},
		{"SPK", {"SPERMINE (FULLY PROTONATED FORM)", "C10 H30 N4 4+"}},
		{"SPM", {"SPERMINE", "C10 H26 N4"}},
		{"SR", {"STRONTIUM ION", "SR1"}},
		{"SUC", {"SUCROSE", "C12 H22 O11"}},
		{"SO4", {"SULFATE ANION", "O4 S1"}},
		{"SUL", {"SULFATE ANION", "O4 S1"}},
		{"SO4", {"SULFATE ION", "O4 S1"}},
		{"TBU", {"TERTIARY-BUTYL ALCOHOL", "C4 H10 O1"}},
		{"TMA", {"TETRAMETHYLAMMONIUM ION", "C4 H12 N1 1+"}},
		{"TEP", {"THEOPHYLLINE", "C7 H8 N4 O2"}},
		{"SCN", {"THIOCYANATE ION", "C1 N1 S1"}},
		{"TRE", {"TREHALOSE", "C12 H22 O11"}},
		{"PGE", {"TRIETHYLENE GLYCOL", "C6 H14 O4"}},
		{"ETF", {"TRIFLUOROETHANOL", "C2 H3 O1 F3"}},
		{"144", {"TRIS-HYDROXYMETHYL-METHYL-AMMONIUM", "C4 H12 N1 O3"}},
		{"UMQ", {"UNDECYL-MALTOSIDE", "C23 H44 O11"}},
		{"URE", {"UREA", "C1 H4 N2 O1"}},
		{"YT3", {"YTTRIUM (III) ION", "Y1"}},
		{"Y1", {"YTTRIUM ION", "Y1"}},
		{"ZN2", {"ZINC ION ON 3-FOLD CRYSTAL AXIS", "ZN1"}},
		{"ZN", {"ZINC ION ON 3-FOLD CRYSTAL AXIS", "ZN1"}},
		{"ZN", {"ZINC ION", "ZN1"}},
		{"PEG", {"DI(HYDROXYETHYL)ETHER", ""}},
		{"BEN", {"BENZAMIDINE", ""}},
		{"2PE", {"NONAETHYLENE GLYCOL", ""}},
		{"P6G", {"HEXAETHYLENE GLYCOL", ""}},
		{"1PE", {"PENTAETHYLENE GLYCOL", ""}},
		{"DOD", {"DEUTERATED WATER", ""}},
		{"ACE", {"ACETYL GROUP", ""}},
		{"BEZ", {"BENZOIC ACID", ""}},
		{"XYP", {"BETA-D-XYLOPYRANOSE", ""}},
		{"NAG", {"N-ACETYL-D-GLUCOSAMINE", ""}},
		{"NDG", {"2-(ACETYLAMINO)-2-DEOXY-A-D-GLUCOPYRANOSE", ""}},
		{"MSE", {"SELENOMETHIONINE", ""}},
		{"PO4", {"PHOSPHATE ION", ""}},
		{"VO4", {"VANADATE ION", ""}},
		{"EPE", {"4-(2-HYDROXYETHYL)-1-PIPERAZINE ETHANESULFONIC ACID", ""}},
		{"MES", {"2-(N-MORPHOLINO)-ETHANESULFONIC ACID", ""}},
		{"MXE", {"2-METHOXYETHANOL", ""}},
		{"XE", {"XENON", ""}},
	};
	const map<const string, const map<string,string>> standard_residues {
		{"ALA",{{"OXT","O2-"},{"CA","C3"},{"CB","C3"},{"C","C2"},{"N","Npl"},{"O","O2"},{"H","H"},{"HA","HC"},{"HB1","HC"},{"HB2","HC"},{"HB3","HC"}}},
		{"ARG",{{"OXT","O2-"},{"CA","C3"},{"CB","C3"},{"C","C2"},{"CD","C3"},{"CG","C3"},{"CZ","C2"},{"NE","Ng+"},{"NH1","Ng+"},{"NH2","Ng+"},{"N","Npl"},{"O","O2"},{"H","H"},{"HA","HC"},{"HB2","HC"},{"HB3","HC"},{"HG2","HC"},{"HG3","HC"},{"HD2","HC"},{"HD3","HC"},{"HE","H"},{"HH11","H"},{"HH12","H"},{"HH21","H"},{"HH22","H"}}},
		{"ASN",{{"OXT","O2-"},{"CA","C3"},{"CB","C3"},{"C","C2"},{"CG","C2"},{"ND2","Npl"},{"N","Npl"},{"OD1","O2"},{"O","O2"},{"H","H"},{"HA","HC"},{"HB2","HC"},{"HB3","HC"},{"HD21","H"},{"HD22","H"}}},

		{"ASP",{{"OXT","O2-"},{"CA","C3"},{"CB","C3"},{"C","C2"},{"CG","Cac"},{"N","Npl"},{"OD1","O2-"},{"OD2","O2-"},{"O","O2"},{"H","H"},{"HA","HC"},{"HB2","HC"},{"HB3","HC"}}},

		{"CYS",{{"OXT","O2-"},{"CA","C3"},{"CB","C3"},{"C","C2"},{"N","Npl"},{"O","O2"},{"SG","S3"},{"H","H"},{"HA","HC"},{"HB2","HC"},{"HB3","HC"},{"HG","H"}}},
		{"CYX",{{"OXT","O2-"},{"CA","C3"},{"CB","C3"},{"C","C2"},{"N","Npl"},{"O","O2"},{"SG","S3"},{"H","H"},{"HA","HC"},{"HB2","HC"},{"HB3","HC"}}},
		{"GLN",{{"OXT","O2-"},{"CA","C3"},{"CB","C3"},{"C","C2"},{"CD","C2"},{"CG","C3"},{"NE2","Npl"},{"N","Npl"},{"OE1","O2"},{"O","O2"},{"H","H"},{"HA","HC"},{"HB2","HC"},{"HB3","HC"},{"HG2","HC"},{"HG3","HC"},{"HE21","H"},{"HE22","H"}}},
		{"GLU",{{"OXT","O2-"},{"CA","C3"},{"CB","C3"},{"C","C2"},{"CD","Cac"},{"CG","C3"},{"N","Npl"},{"OE1","O2-"},{"OE2","O2-"},{"O","O2"},{"H","H"},{"HA","HC"},{"HB2","HC"},{"HB3","HC"},{"HG2","HC"},{"HG3","HC"}}},
		{"GLY",{{"OXT","O2-"},{"CA","C3"},{"C","C2"}, {"N","Npl"},{"O","O2"},{"H","H"},{"HA2","HC"},{"HA3","HC"}}},
		{"HID",{{"OXT","O2-"},{"CA","C3"},{"CB","C3"},{"C","C2"},{"CD2","Car"},{"CE1","Car"},{"CG","Car"},{"ND1","Npl"},{"NE2","N2"},{"N","Npl"},{"O","O2"},{"H","H"},{"HA","HC"},{"HB2","HC"},{"HB3","HC"},{"HD1","H"},{"HE1","HC"},{"HD2","HC"}}},
		{"HIE",{{"OXT","O2-"},{"CA","C3"},{"CB","C3"},{"C","C2"},{"CD2","Car"},{"CE1","Car"},{"CG","Car"},{"ND1","N2"},{"NE2","Npl"},{"N","Npl"},{"O","O2"},{"H","H"},{"HA","HC"},{"HB2","HC"},{"HB3","HC"},{"HE1","HC"},{"HE2","H"},{"HD2","HC"}}},
		{"HIP",{{"OXT","O2-"},{"CA","C3"},{"CB","C3"},{"C","C2"},{"CD2","C2"},{"CE1","C2"},{"CG","C2"},{"ND1","Npl"},{"NE2","Npl"},{"N","Npl"},{"O","O2"},{"H","H"},{"HA","HC"},{"HB2","HC"},{"HB3","HC"},{"HD1","H"},{"HE1","HC"},{"HE2","H"},{"HD2","HC"}}},
		{"HIS",{{"OXT","O2-"},{"CA","C3"},{"CB","C3"},{"C","C2"},{"CD2","Car"},{"CE1","Car"},{"CG","Car"},{"ND1","Npl"},{"NE2","N2"},{"N","Npl"},{"O","O2"},{"H","H"},{"HN","H"},{"HA","HC"},{"HB1","HC"},{"HB2","HC"},{"HB3","HC"},{"HD1","H"},{"HD2","HC"},{"HE1","HC"},{"HE2","H"}}},
		{"ILE",{{"OXT","O2-"},{"CA","C3"},{"CB","C3"},{"C","C2"},{"CD1","C3"},{"CG1","C3"},{"CG2","C3"},{"N","Npl"},{"O","O2"},{"H","H"},{"HA","HC"},{"HB","HC"},{"HG21","HC"},{"HG22","HC"},{"HG23","HC"},{"HG12","HC"},{"HG13","HC"},{"HD11","HC"},{"HD12","HC"},{"HD13","HC"}}},
		{"LEU",{{"OXT","O2-"},{"CA","C3"},{"CB","C3"},{"C","C2"},{"CD1","C3"},{"CD2","C3"},{"CG","C3"},{"N","Npl"},{"O","O2"},{"H","H"},{"HA","HC"},{"HB2","HC"},{"HB3","HC"},{"HG","HC"},{"HD11","HC"},{"HD12","HC"},{"HD13","HC"},{"HD21","HC"},{"HD22","HC"},{"HD23","HC"}}},
		{"LYS",{{"OXT","O2-"},{"CA","C3"},{"CB","C3"},{"C","C2"},{"CD","C3"},{"CE","C3"},{"CG","C3"},{"N","Npl"},{"NZ","N3+"},{"O","O2"},{"H","H"},{"HA","HC"},{"HB2","HC"},{"HB3","HC"},{"HG2","HC"},{"HG3","HC"},{"HD2","HC"},{"HD3","HC"},{"HE2","HC"},{"HE3","HC"},{"HZ1","H"},{"HZ2","H"},{"HZ3","H"}}},
		{"MET",{{"OXT","O2-"},{"CA","C3"},{"CB","C3"},{"C","C2"},{"CE","C3"},{"CG","C3"},{"N","Npl"},{"O","O2"},{"SD","S3"},{"H","H"},{"HA","HC"},{"HB2","HC"},{"HB3","HC"},{"HG2","HC"},{"HG3","HC"},{"HE1","HC"},{"HE2","HC"},{"HE3","HC"}}},
		{"MSE",{{"OXT","O2-"},{"CA","C3"},{"CB","C3"},{"C","C2"},{"CE","C3"},{"CG","C3"},{"N","Npl"},{"O","O2"},{"SE","Se"},{"H","H"},{"HA","HC"},{"HB2","HC"},{"HB3","HC"},{"HG2","HC"},{"HG3","HC"},{"HE1","HC"},{"HE2","HC"},{"HE3","HC"}}},
		{"PHE",{{"OXT","O2-"},{"CA","C3"},{"CB","C3"},{"C","C2"},{"CD1","Car"},{"CD2","Car"},{"CE1","Car"},{"CE2","Car"},{"CG","Car"},{"CZ","Car"},{"N","Npl"},{"O","O2"},{"H","H"},{"HA","HC"},{"HB2","HC"},{"HB3","HC"},{"HD1","HC"},{"HE1","HC"},{"HZ","HC"},{"HE2","HC"},{"HD2","HC"}}},
		{"PRO",{{"OXT","O2-"},{"CA","C3"},{"CB","C3"},{"C","C2"},{"CD","C3"},{"CG","C3"},{"N","Npl"},{"O","O2"},{"HD2","HC"},{"HD3","HC"},{"HG2","HC"},{"HG3","HC"},{"HB2","HC"},{"HB3","HC"},{"HA","HC"}}},
		{"SER",{{"OXT","O2-"},{"CA","C3"},{"CB","C3"},{"C","C2"},{"N","Npl"},{"OG","O3"},{"O","O2"},{"H","H"},{"HA","HC"},{"HB2","HC"},{"HB3","HC"},{"HG","H"}}},
		{"THR",{{"OXT","O2-"},{"CA","C3"},{"CB","C3"},{"C","C2"},{"CG2","C3"},{"N","Npl"},{"OG1","O3"},{"O","O2"},{"H","H"},{"HA","HC"},{"HB","HC"},{"HG21","HC"},{"HG22","HC"},{"HG23","HC"},{"HG1","H"}}},
		{"TRP",{{"OXT","O2-"},{"CA","C3"},{"CB","C3"},{"C","C2"},{"CD1","Car"},{"CD2","Car"},{"CE2","Car"},{"CE3","Car"},{"CG","Car"},{"CH2","Car"},{"CZ2","Car"},{"CZ3","Car"},{"NE1","Npl"},{"N","Npl"},{"O","O2"},{"H","H"},{"HA","HC"},{"HB2","HC"},{"HB3","HC"},{"HD1","HC"},{"HE1","H"},{"HZ2","HC"},{"HH2","HC"},{"HZ3","HC"},{"HE3","HC"}}},
		{"TYR",{{"OXT","O2-"},{"CA","C3"},{"CB","C3"},{"C","C2"},{"CD1","Car"},{"CD2","Car"},{"CE1","Car"},{"CE2","Car"},{"CG","Car"},{"CZ","Car"},{"N","Npl"},{"OH","O3"},{"O","O2"},{"H","H"},{"HA","HC"},{"HB2","HC"},{"HB3","HC"},{"HD1","HC"},{"HE1","HC"},{"HH","H"},{"HE2","HC"},{"HD2","HC"}}},
		{"VAL",{{"OXT","O2-"},{"CA","C3"},{"CB","C3"},{"C","C2"},{"CG1","C3"},{"CG2","C3"},{"N","Npl"},{"O","O2"},{"H","H"},{"HA","HC"},{"HB","HC"},{"HG11","HC"},{"HG12","HC"},{"HG13","HC"},{"HG21","HC"},{"HG22","HC"},{"HG23","HC"}}},

		{"HOH",{{"O","O3"},{"H1","H"},{"H2","H"}}},
	};
	const IdatmInfoMap infoMap {
	        { "Car", { Planar, 3, "aromatic carbon" } },
	        { "C3", { Tetrahedral, 4, "sp3-hybridized carbon" } },
	        { "C2", { Planar, 3, "sp2-hybridized carbon" } },
	        { "C1", { Linear, 2, "sp-hybridized carbon bonded to 2 other atoms" } },
	        { "C1-", { Linear, 1, "sp-hybridized carbon bonded to 1 other atom" } },
	        { "Cac", { Planar, 3, "carboxylate carbon" } },
	        { "N3+", { Tetrahedral, 4, "sp3-hybridized nitrogen, formal positive charge" } },
	        { "N3", { Tetrahedral, 3, "sp3-hybridized nitrogen, neutral" } },
	        { "Npl", { Planar, 3, "sp2-hybridized nitrogen, not double bonded" } },
	        { "N2+", { Planar, 3, "sp2-hybridized nitrogen, double bonded, formal positive charge" } },
	        { "N2", { Planar, 2, "sp2-hybridized nitrogen, double bonded" } },
	        { "N1+", { Linear, 2, "sp-hybridized nitrogen bonded to 2 other atoms" } },
	        { "N1", { Linear, 1, "sp-hybridized nitrogen bonded to 1 other atom" } },
	        { "Ntr", { Planar, 3, "nitro nitrogen" } },
	        { "Ng+", { Planar, 3, "guanidinium/amidinium nitrogen, partial positive charge" } },
	        { "O3", { Tetrahedral, 2, "sp3-hybridized oxygen" } },
	        { "O3-", { Tetrahedral, 1, "phosphate or sulfate oxygen sharing formal negative charge" } },
	        { "Oar", { Planar, 2, "aromatic oxygen" } },
	        { "Oar+", { Planar, 2, "aromatic oxygen, formal positive charge" } },
	        { "O2", { Planar, 1, "sp2-hybridized oxygen" } },
	        { "O2-", { Planar, 1, "carboxylate oxygen sharing formal negative charge; nitro group oxygen" } },
	        { "O1", { Linear, 1, "sp-hybridized oxygen" } },
	        { "O1+", { Linear, 1, "sp-hybridized oxygen, formal positive charge" } },
	        { "S3+", { Tetrahedral, 3, "sp3-hybridized sulfur, formal positive charge" } },
	        { "S3", { Tetrahedral, 2, "sp3-hybridized sulfur, neutral" } },
	        { "S3-", { Tetrahedral, 1, "thiophosphate sulfur, sharing formal negative charge" } },
	        { "S2", { Planar, 1, "sp2-hybridized sulfur" } },
	        { "Sar", { Planar, 2, "aromatic sulfur" } },
	        { "Sac", { Tetrahedral, 4, "sulfate, sulfonate, or sulfamate sulfur" } },
	        { "Son", { Tetrahedral, 4, "sulfone sulfur" } },
	        { "Sxd", { Tetrahedral, 3, "sulfoxide sulfur" } },
	        { "S", { Tetrahedral, 4, "other sulfur" } },
	        { "Pac", { Tetrahedral, 4, "phosphate phosphorus" } },
	        { "Pox", { Tetrahedral, 4, "P-oxide phosphorus" } },
	        { "P3+", { Tetrahedral, 4, "sp3-hybridized phosphorus, formal positive charge" } },
	        { "P", { TrigonalBipyramidal, 5, "other phosphorus" } }, // Janez : see exception in compute_hydrogen (molecule.cpp)
	        { "HC", { Single, 1, "hydrogen bonded to carbon" } },
	        { "H", { Single, 1, "other hydrogen" } },
	        { "DC", { Single, 1, "deuterium bonded to carbon" } },
	        { "D", { Single, 1, "other deuterium" } },
	        { "F", { Single, 1, "fluoride" } },
	        { "Cl", { Single, 1, "chloride" } },
	        { "Br", { Single, 1, "bromium" } },
	        { "I", { Single, 1, "iodide" } },
	        { "Si", { Tetrahedral, 4, "silicon" } }, // Janez : added for CHEMBL95846 and alike
	        { "Mg", { Ion, 0, "magnesium" } },
	        { "Mn", { Ion, 0, "manganese" } },
	        { "Zn", { Ion, 0, "zinc" } },
	};
	const char* const idatm_unmask[] {	
		"Ac",
		"Ag",
		"Al",
		"Am",
		"Ar",
		"As",
		"At",
		"Au",
		"B",
		"Ba",
		"Be",
		"Bh",
		"Bi",
		"Bk",
		"Br",
		"C",
		"C1",
		"C1-",
		"C2",
		"C3",
		"Ca",
		"Cac",
		"Car",
		"Cd",
		"Ce",
		"Cf",
		"Cl",
		"Cm",
		"Co",
		"Cr",
		"Cs",
		"Cu",
		"D",
		"Db",
		"DC",
		"Ds",
		"Dy",
		"Er",
		"Es",
		"Eu",
		"F",
		"Fe",
		"Fm",
		"Fr",
		"Ga",
		"Gd",
		"Ge",
		"H",
		"HC",
		"He",
		"Hf",
		"Hg",
		"Ho",
		"Hs",
		"I",
		"In",
		"Ir",
		"K",
		"Kr",
		"La",
		"Li",
		"Lr",
		"Lu",
		"Lw",
		"Md",
		"Mg",
		"Mn",
		"Mo",
		"Mt",
		"N",
		"N1",
		"N1+",
		"N2",
		"N2+",
		"N3",
		"N3+",
		"Na",
		"Nb",
		"Nd",
		"Ne",
		"Ng+",
		"Ni",
		"No",
		"Nox",
		"Np",
		"Npl",
		"Ntr",
		"O",
		"O1",
		"O1+",
		"O2",
		"O2-",
		"O3",
		"O3-",
		"Oar",
		"Oar+",
		"Os",
		"P",
		"P3+",
		"Pa",
		"Pac",
		"Pb",
		"Pd",
		"Pm",
		"Po",
		"Pox",
		"Pr",
		"Pt",
		"Pu",
		"Ra",
		"Rb",
		"Re",
		"Rf",
		"Rh",
		"Rn",
		"Ru",
		"S",
		"S2",
		"S3",
		"S3-",
		"S3+",
		"Sac",
		"Sar",
		"Sb",
		"Sc",
		"Se",
		"Sg",
		"Si",
		"Sm",
		"Sn",
		"Son",
		"Sr",
		"Sxd",
		"Ta",
		"Tb",
		"Tc",
		"Te",
		"Th",
		"Ti",
		"Tl",
		"Tm",
		"U",
		"V",
		"W",
		"Xe",
		"Y",
		"Yb",
		"Zn",
		"Zr",
		"???"
	};

	const map<pair<string, string>, double> repulsion_idx {
		{{"F", "Sar"}, 3.6}, 
		{{"Npl", "Npl"}, 2.5}, 
		{{"Npl", "S3"}, 2.8},
		{{"O2", "O2-"}, 2.2},
		{{"O2-", "S3"}, 2.9},
		{{"Cac", "S3"}, 3.4},
		{{"Cl", "O2-"}, 2.8},
		{{"N2", "O2-"}, 2.8},
		{{"N2", "S3"}, 3.0},
		{{"N3+", "N3+"}, 2.5},
		{{"N3+", "Ng+"}, 2.9},
		{{"N3+", "O2-"}, 2.5},
		{{"N3+", "O3"}, 2.5},
		{{"N3+", "S3"}, 3.0},
		{{"Ng+", "Ng+"}, 2.7}
	};

	const map<const string, const int> idatm_mask {	
		{"Ac",  0},
		{"Ag",  1},
		{"Al",  2},
		{"Am",  3},
		{"Ar",  4},
		{"As",  5},
		{"At",  6},
		{"Au",  7},
		{"B",   8},
		{"Ba",  9},
		{"Be",  10},
		{"Bh",  11},
		{"Bi",  12},
		{"Bk",  13},
		{"Br",  14},
		{"C",   15},
		{"C1",  16},
		{"C1-", 17},
		{"C2",  18},
		{"C3",  19},
		{"Ca",  20},
		{"Cac", 21},
		{"Car", 22},
		{"Cd",  23},
		{"Ce",  24},
		{"Cf",  25},
		{"Cl",  26},
		{"Cm",  27},
		{"Co",  28},
		{"Cr",  29},
		{"Cs",  30},
		{"Cu",  31},
		{"D",   32},
		{"Db",  33},
		{"DC",  34},
		{"Ds",  35},
		{"Dy",  36},
		{"Er",  37},
		{"Es",  38},
		{"Eu",  39},
		{"F",   40},
		{"Fe",  41},
		{"Fm",  42},
		{"Fr",  43},
		{"Ga",  44},
		{"Gd",  45},
		{"Ge",  46},
		{"H",   47},
		{"HC",  48},
		{"He",  49},
		{"Hf",  50},
		{"Hg",  51},
		{"Ho",  52},
		{"Hs",  53},
		{"I",   54},
		{"In",  55},
		{"Ir",  56},
		{"K",   57},
		{"Kr",  58},
		{"La",  59},
		{"Li",  60},
		{"Lr",  61},
		{"Lu",  62},
		{"Lw",  63},
		{"Md",  64},
		{"Mg",  65},
		{"Mn",  66},
		{"Mo",  67},
		{"Mt",  68},
		{"N",   69},
		{"N1",  70},
		{"N1+", 71},
		{"N2",  72},
		{"N2+", 73},
		{"N3",  74},
		{"N3+", 75},
		{"Na",  76},
		{"Nb",  77},
		{"Nd",  78},
		{"Ne",  79},
		{"Ng+", 80},
		{"Ni",  81},
		{"No",  82},
		{"Nox", 83},
		{"Np",  84},
		{"Npl", 85},
		{"Ntr", 86},
		{"O",   87},
		{"O1",  88},
		{"O1+", 89},
		{"O2",  90},
		{"O2-", 91},
		{"O3",  92},
		{"O3-", 93},
		{"Oar", 94},
		{"Oar+",95},
		{"Os",  96},
		{"P",   97},
		{"P3+", 98},
		{"Pa",  99},
		{"Pac", 100},
		{"Pb",  101},
		{"Pd",  102},
		{"Pm",  103},
		{"Po",  104},
		{"Pox", 105},
		{"Pr",  106},
		{"Pt",  107},
		{"Pu",  108},
		{"Ra",  109},
		{"Rb",  110},
		{"Re",  111},
		{"Rf",  112},
		{"Rh",  113},
		{"Rn",  114},
		{"Ru",  115},
		{"S",   116},
		{"S2",  117},
		{"S3",  118},
		{"S3-", 119},
		{"S3+", 120},
		{"Sac", 121},
		{"Sar", 122},
		{"Sb",  123},
		{"Sc",  124},
		{"Se",  125},
		{"Sg",  126},
		{"Si",  127},
		{"Sm",  128},
		{"Sn",  129},
		{"Son", 130},
		{"Sr",  131},
		{"Sxd", 132},
		{"Ta",  133},
		{"Tb",  134},
		{"Tc",  135},
		{"Te",  136},
		{"Th",  137},
		{"Ti",  138},
		{"Tl",  139},
		{"Tm",  140},
		{"U",   141},
		{"V",   142},
		{"W",   143},
		{"Xe",  144},
		{"Y",   145},
		{"Yb",  146},
		{"Zn",  147},
		{"Zr",  148},
		{"???",  149},
	};
	
	
	const double vdw_radius[] = {	
		2, // Ac
		1.72, // Ag
		2, // Al
		2, // Am
		1.88, // Ar
		1.85, // As
		2, // At
		1.66, // Au
		2, // B
		2, // Ba
		2, // Be
		2, // Bh
		2, // Bi
		2, // Bk
		1.85, // Br
		1.7, // C
		1.7, // C1
		1.7, // C1-
		1.7, // C2
		1.7, // C3
		2, // Ca
		1.7, // Cac
		1.7, // Car
		1.58, // Cd
		2, // Ce
		2, // Cf
		1.75, // Cl
		2, // Cm
		2, // Co
		2, // Cr
		2, // Cs
		1.4, // Cu
		1.2, // D
		2, // Db
		1.2, // DC
		2, // Ds
		2, // Dy
		2, // Er
		2, // Es
		2, // Eu
		1.47, // F
		2, // Fe
		2, // Fm
		2, // Fr
		1.87, // Ga
		2, // Gd
		2, // Ge
		1.2, // H
		1.2, // HC
		1.4, // He
		2, // Hf
		1.55, // Hg
		2, // Ho
		2, // Hs
		1.98, // I
		1.93, // In
		2, // Ir
		2.75, // K
		2.02, // Kr
		2, // La
		1.82, // Li
		2, // Lr
		2, // Lu
		2, // Lw
		2, // Md
		1.73, // Mg
		2, // Mn
		2, // Mo
		2, // Mt
		1.55, // N
		1.55, // N1
		1.55, // N1+
		1.55, // N2
		1.55, // N2+
		1.55, // N3
		1.55, // N3+
		2.27, // Na
		2, // Nb
		2, // Nd
		1.54, // Ne
		1.55, // Ng+
		1.63, // Ni
		2, // No
		1.55, // Nox
		2, // Np
		1.55, // Npl
		1.55, // Ntr
		1.52, // O
		1.52, // O1
		1.52, // O1+
		1.52, // O2
		1.52, // O2-
		1.52, // O3
		1.52, // O3-
		1.52, // Oar
		1.52, // Oar+
		2, // Os
		1.8, // P
		1.8, // P3+
		2, // Pa
		1.8, // Pac
		2.02, // Pb
		1.63, // Pd
		2, // Pm
		2, // Po
		1.8, // Pox
		2, // Pr
		1.72, // Pt
		2, // Pu
		2, // Ra
		2, // Rb
		2, // Re
		2, // Rf
		2, // Rh
		2, // Rn
		2, // Ru
		1.8, // S
		1.8, // S2
		1.8, // S3
		1.8, // S3-
		1.8, // S3+
		1.8, // Sac
		1.8, // Sar
		2, // Sb
		2, // Sc
		1.9, // Se
		2, // Sg
		2.1, // Si
		2, // Sm
		2.17, // Sn
		1.8, // Son
		2, // Sr
		1.8, // Sxd
		2, // Ta
		2, // Tb
		2, // Tc
		2.06, // Te
		2, // Th
		2, // Ti
		1.96, // Tl
		2, // Tm
		1.86, // U
		2, // V
		2, // W
		2.16, // Xe
		2, // Y
		2, // Yb
		1.39, // Zn
		2 // Zr
	};
	
	const string EW = "^N|^O|F|Cl|Br"; // electron withdrawing atoms
	const string XX = "^C|^N|^O|^S|^P";
	const string XA = "^O|^S";
	const string XB = "^N|^P";
	const string XC = "F|Cl|Br|I";
	const string XD = "^S|^P";

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

	const map<string, vector<string>> gaff_replacement {
		{"cl",{"cl"}},
		{"Si",{"Si"}},
		{"cx",{"cx","c3"}},
		{"cy",{"cy","c3"}},
		{"c3",{"c3"}},
		{"c",{"c","c2"}},
		{"cz",{"cz"}},
		{"cp",{"cp","ca","cc","cd","c2"}},
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

	vector<vector<string>> get_replacement(const vector<string> &initial);

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

};
#endif
