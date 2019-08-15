/* Copyright (c) 2016-2019 Chopra Lab at Purdue University, 2013-2016 Janez Konc at National Institute of Chemistry and Samudrala Group at University of Washington
 *
 * This program is free for educational and academic use
 *
 * This program is free software; you can redistribute it and/or modify
 * it under the terms of the GNU General Public License as published by
 * the Free Software Foundation version 3 of the License.
 *
 * This program is distributed in the hope that it will be useful,
 * but WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the
 * GNU General Public License for more details.
 */

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

	string replace_str_char(string str, const string& replace, char ch);
	vector<string> ssplit(const string &source, const char*delimiter = " ", bool keepEmpty = false);
	
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
		{"Du", "H"}, // dummy atom changed to element H (fix issue #112)
		{"Du.C", "C"}, // dummyC atom changed to element C (fix issue #112)
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
	char get_one_letter(string s);
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
		{"VAL"},
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

        extern const map<const string, const pair<string,string>> non_specific_binders;

        extern const map<const string, const map<string,string>> standard_residues;

        extern const map<const string, const map<string,string>> cofactor_residues;

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
                { "Nox", { Tetrahedral, 4, "N-oxide amine" } },
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
	        { "Ca", { Ion, 0, "calcium" } },
	        { "Na", { Ion, 0, "sodium" } },
                { "K",  { Ion, 0, "potassium" } },
                { "Fe", { Ion, 0, "iron" } },
                { "Co", { Ion, 0, "cobolt"} },
                { "Cu", { Ion, 0, "copper"} },
                { "Ni", { Ion, 0, "nickel"} }
	};

	const IdatmEntry &get_info_map(const string &name);
	
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

	vector<vector<string>> get_replacement(const vector<string> &initial);

};
#endif
