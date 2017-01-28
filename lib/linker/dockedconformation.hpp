#ifndef DOCKEDCONFORMATION_H
#define DOCKEDCONFORMATION_H
#include "pdbreader/molecule.hpp"
#include <memory>

namespace Linker {

	class DockedConformation {
	public:
		typedef vector<DockedConformation> Vec;
	private:
		unique_ptr<Molib::Molecule> __ligand;
		unique_ptr<Molib::Molecule> __receptor;
		double __energy;
	public:
		DockedConformation() : __ligand(nullptr), __receptor(nullptr), __energy(0) {}
		DockedConformation(Molib::Molecule ligand, Molib::Molecule receptor, double energy) 
			: __ligand(unique_ptr<Molib::Molecule>(new Molib::Molecule(ligand))), 
			__receptor(unique_ptr<Molib::Molecule>(new Molib::Molecule(receptor))), __energy(energy) {}
		Molib::Molecule& get_ligand() { return *__ligand; }
		const Molib::Molecule& get_ligand() const { return *__ligand; }
		Molib::Molecule& get_receptor() { return *__receptor; }
		const Molib::Molecule& get_receptor() const { return *__receptor; }
		double get_energy() { return __energy; }
		double get_energy() const { return __energy; }

		static void sort(DockedConformation::Vec &v);

        friend bool operator<( const DockedConformation& first, const DockedConformation& second ) { return first.get_energy() < second.get_energy(); }
		friend ostream& operator<<(ostream& os, const DockedConformation &conf);
	};
	
};

#endif
