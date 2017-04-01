#include "element.hpp"
#include "helper/debug.hpp"
#include <algorithm>
using namespace std;

namespace Molib {
        const vector<string> Element::symbols {
		"LP",  "H", "He", "Li", "Be",  "B",  "C",  "N",  "O",
		 "F", "Ne", "Na", "Mg", "Al", "Si",  "P",  "S", "Cl",
		"Ar",  "K", "Ca", "Sc", "Ti",  "V", "Cr", "Mn", "Fe",
		"Co", "Ni", "Cu", "Zn", "Ga", "Ge", "As", "Se", "Br",
		"Kr", "Rb", "Sr",  "Y", "Zr", "Nb", "Mo", "Tc", "Ru",
		"Rh", "Pd", "Ag", "Cd", "In", "Sn", "Sb", "Te",  "I",
		"Xe", "Cs", "Ba", "La", "Ce", "Pr", "Nd", "Pm", "Sm",
		"Eu", "Gd", "Tb", "Dy", "Ho", "Er", "Tm", "Yb", "Lu",
		"Hf", "Ta",  "W", "Re", "Os", "Ir", "Pt", "Au", "Hg",
		"Tl", "Pb", "Bi", "Po", "At", "Rn", "Fr", "Ra", "Ac",
		"Th", "Pa",  "U", "Np", "Pu", "Am", "Cm", "Bk", "Cf",
		"Es", "Fm", "Md", "No", "Lr", "Rf", "Db", "Sg", "Bh",
		"Hs", "Mt", "Uun", "Uuu", "Uub"
	};
	const vector<double> covalent {
		0.0 /* */, 0.23 /* H */, 0.0 /* He */, 0.68 /* Li */,
		0.35 /* Be */, 0.83 /* B */, 0.68 /* C */, 0.68 /* N */,
		0.68 /* O */, 0.64 /* F */, 0.0 /* Ne */, 0.97 /* Na */,
		1.10 /* Mg */, 1.35 /* Al */, 1.20 /* Si */, 1.05 /* P */,
		1.02 /* S */, 0.99 /* Cl */, 0.0 /* Ar */, 1.33 /* K */,
		0.99 /* Ca */, 1.44 /* Sc */, 1.47 /* Ti */, 1.33 /* V */,
		1.35 /* Cr */, 1.35 /* Mn */, 1.34 /* Fe */, 1.33 /* Co */,
		1.50 /* Ni */, 1.52 /* Cu */, 1.45 /* Zn */, 1.22 /* Ga */,
		1.17 /* Ge */, 1.21 /* As */, 1.22 /* Se */, 1.21 /* Br */,
		0.0 /* Kr */, 1.47 /* Rb */, 1.12 /* Sr */, 1.78 /* Y */,
		1.56 /* Zr */, 1.48 /* Nb */, 1.47 /* Mo */, 1.35 /* Tc */,
		1.40 /* Ru */, 1.45 /* Rh */, 1.50 /* Pd */, 1.59 /* Ag */,
		1.69 /* Cd */, 1.63 /* In */, 1.46 /* Sn */, 1.46 /* Sb */,
		1.47 /* Te */, 1.40 /* I */, 0.0 /* Xe */, 1.67 /* Cs */,
		1.34 /* Ba */, 1.87 /* La */, 1.83 /* Ce */, 1.82 /* Pr */,
		1.81 /* Nd */, 1.80 /* Pm */, 1.80 /* Sm */, 1.99 /* Eu */,
		1.79 /* Gd */, 1.76 /* Tb */, 1.75 /* Dy */, 1.74 /* Ho */,
		1.73 /* Er */, 1.72 /* Tm */, 1.94 /* Yb */, 1.72 /* Lu */,
		1.57 /* Hf */, 1.43 /* Ta */, 1.37 /* W */, 1.35 /* Re */,
		1.37 /* Os */, 1.32 /* Ir */, 1.50 /* Pt */, 1.50 /* Au */,
		1.70 /* Hg */, 1.55 /* Tl */, 1.54 /* Pb */, 1.54 /* Bi */,
		1.68 /* Po */, 0.0 /* At */, 0.0 /* Rn */, 0.0 /* Fr */,
		1.90 /* Ra */, 1.88 /* Ac */, 1.79 /* Th */, 1.61 /* Pa */,
		1.58 /* U */, 1.55 /* Np */, 1.53 /* Pu */, 1.51 /* Am */,
	};
        const vector<double> masses {
            0.00000 /*    */, 1.00800 /* H  */, 4.00200 /* He */, 6.94000 /* Li */, 
            9.01200 /* Be */, 10.8100 /* B  */, 12.0110 /* C  */, 14.0070 /* N  */, 
            15.9990 /* O  */, 18.9980 /* F  */, 20.1797 /* Ne */, 22.9890 /* Na */, 
            24.3050 /* Mg */, 26.9810 /* Al */, 28.0850 /* Si */, 30.9730 /* P  */, 
            32.0600 /* S  */, 35.4500 /* Cl */, 39.9480 /* Ar */, 39.0980 /* K  */, 
            40.0780 /* Ca */, 44.9550 /* Sc */, 47.8670 /* Ti */, 50.9415 /* V  */, 
            51.9961 /* Cr */, 54.9380 /* Mn */, 55.8450 /* Fe */, 58.9330 /* Co */, 
            58.6934 /* Ni */, 63.5460 /* Cu */, 65.3800 /* Zn */, 69.7230 /* Ga */, 
            72.6300 /* Ge */, 74.9210 /* As */, 78.9710 /* Se */, 79.9040 /* Br */, 
            83.7980 /* Kr */, 85.4678 /* Rb */, 87.6200 /* Sr */, 88.9050 /* Y  */, 
            91.2240 /* Zr */, 92.9060 /* Nb */, 95.9500 /* Mo */, 97.0000 /* Tc */, 
            101.070 /* Ru */, 102.905 /* Rh */, 106.420 /* Pd */, 107.868 /* Ag */, 
            112.414 /* Cd */, 114.818 /* In */, 118.710 /* Sn */, 121.760 /* Sb */, 
            127.600 /* Te */, 126.904 /* I  */, 131.293 /* Xe */, 132.905 /* Cs */, 
            137.327 /* Ba */, 138.905 /* La */, 140.116 /* Ce */, 140.907 /* Pr */, 
            144.242 /* Nd */, 145.000 /* Pm */, 150.360 /* Sm */, 151.964 /* Eu */, 
            157.250 /* Gd */, 158.925 /* Tb */, 162.500 /* Dy */, 164.930 /* Ho */, 
            167.259 /* Er */, 168.934 /* Tm */, 173.045 /* Yb */, 174.967 /* Lu */, 
            178.490 /* Hf */, 180.947 /* Ta */, 183.840 /* W  */, 186.207 /* Re */, 
            190.230 /* Os */, 192.217 /* Ir */, 195.084 /* Pt */, 196.966 /* Au */, 
            200.592 /* Hg */, 204.380 /* Tl */, 207.200 /* Pb */, 208.980 /* Bi */, 
            209.000 /* Po */, 210.000 /* At */, 222.000 /* Rn */, 223.000 /* Fr */, 
            226.000 /* Ra */, 227.000 /* Ac */, 232.038 /* Th */, 231.035 /* Pa */, 
            238.028 /* U  */, 237.000 /* Np */, 244.000 /* Pu */, 243.000 /* Am */, 
        };
	ostream &operator<<(ostream &os, const Element &a) {
		os << a.name();
		return os;
	}

	double Element::bondRadius(Element a) {
		if (a.number() < 0 || a.number() >= Element::NumCovalent)
			return 0.0;
		else
			return covalent[a.number()];
	}
	double Element::bondLength(Element a0, Element a1) {
		return bondRadius(a0) + bondRadius(a1);
	}
	Element::AS Element::atomicNumber(const string &name) {
		string symbol;
		dbgmsg(name);
		if (name == "")
			return LonePair;
		//~ string n;
		//~ std::remove_copy_if(name.begin(), name.end(), std::back_inserter( n ), [](const char &c){ return isdigit(c)||c=='\''; });
		//~ dbgmsg(n << " " << n.size());
		//~ symbol.push_back(n[0]);
		//~ if (n.size() > 1) {
			 //~ if(isupper(n[1])) {
				//~ symbol.push_back(tolower(n[1]));
			//~ } else {
				//~ symbol.push_back(n[1]);
			//~ }
		//~ }
		//~ dbgmsg(symbol);
		
		//~ string n;
		//~ std::remove_copy_if(name.begin(), name.end(), std::back_inserter( n ), [](const char &c){ return isdigit(c)||c=='\''; });
		//~ dbgmsg(n << " " << n.size());
		symbol.push_back(name[0]);
		if (name.size() > 1 && isalpha(name[1])) {
			symbol.push_back(tolower(name[1]));
		}
		dbgmsg(symbol);
	
		if (symbol.size() == 1)
			switch (symbol[0]) {
			case 'H': return H;
			case 'D': return D;	// deuterium
			case 'T': return T;	// tritium
			case 'B': return B;
			case 'C': return C;
			case 'N': return N;
			case 'O': return O;
			case 'F': return F;
			case 'P': return P;
			case 'S': return S;
			case 'K': return K;
			case 'V': return V;
			case 'Y': return Y;
			case 'I': return I;
			case 'W': return W;
			case 'U': return U;
			default: return LonePair;
			}
		for (vector<string>::const_iterator e = symbols.begin() + 1; e != symbols.end(); e++)
			if (*e == symbol)
				return AS(e - symbols.begin());
		for (vector<string>::const_iterator e = symbols.begin() + 1; e != symbols.end(); e++)
			if (symbol.find(*e) != string::npos)
				return AS(e - symbols.begin());
		return LonePair;
	}
	double Element::mass() const {
            return masses[number()];
        }

}
