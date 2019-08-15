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

#ifndef ELEMENT_H
#define ELEMENT_H
#include <iostream>
#include <vector>

namespace Molib {
	class Element {
	public:
		static const std::vector<std::string> symbols;
		// Atomic Symbols:
		enum AS {
			LonePair, H, D = 1, T = 1, He,
			Li, Be, B, C, N, O, F, Ne,
			Na, Mg, Al, Si, P, S, Cl, Ar,
			K, Ca, Sc, Ti, V, Cr, Mn, Fe, Co, Ni, Cu, Zn,
			Ga, Ge, As, Se, Br, Kr,
			Rb, Sr, Y, Zr, Nb, Mo, Tc, Ru, Rh, Pd, Ag, Cd,
			In, Sn, Sb, Te, I, Xe,
			Cs, Ba, La,
			Ce, Pr, Nd, Pm, Sm, Eu, Gd, Tb, Dy, Ho, Er, Tm, Yb, Lu,
			Hf, Ta, W, Re, Os, Ir, Pt, Au, Hg, Tl, Pb, Bi, Po, At, Rn,
			Fr, Ra, Ac,
			Th, Pa, U, Np, Pu, Am, Cm, Bk, Cf, Es, Fm, Md, No, Lr,
			Rf, Db, Sg, Bh, Hs, Mt, Uun, Uuu, Uub
		};
		static const int NumSymbols = 115;
		static const int NumCovalent = 96;
		static double bondRadius(Element);
		static double bondLength(Element, Element);
	private:
		AS atomicNumber(const std::string &name);
		AS as; // atomic number
	public:
		Element(const std::string &name) : as(atomicNumber(name)) {}
		Element(int i) : as(AS(i)) {}
		Element(AS a) : as(a) {}
		Element(const Element &other) { as = other.as; }
		const std::string name() const {
			if (as >= NumSymbols)
				return "??";
			return symbols[as];
		}
		int     number() const { return int(as); }
		double  mass() const;
		long    hash() const { return number(); }
		const Element& operator=(const Element &a) { this->as = a.as; return *this; }
		bool    operator==(const Element &a) const { return as == a.as; }
		bool    operator!=(const Element &a) const { return as != a.as; }
		bool    operator<(const Element &a) const { return as < a.as; }
		bool    operator<=(const Element &a) const { return as <= a.as; }
		bool    operator>(const Element &a) const { return as > a.as; }
		bool    operator>=(const Element &a) const { return as >= a.as; }
		friend std::ostream& operator<< (std::ostream&, const Element&);
	};
}  // namespace Molib
#endif
