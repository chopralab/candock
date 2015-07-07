#include "segment.hpp"
//~ #include "linker.hpp"
#include "pdbreader/molecule.hpp"
#include "helper/benchmark.hpp"
#include "helper/help.hpp"

namespace Molib {
	ostream& operator<< (ostream& stream, const Segment& s) {
		stream << "Segment(" << s.__name << ") : atom numbers = ";
		for (auto &pa : s.__atoms) stream << pa->atom_number() << " ";
		return stream;
	}
	const Atom& Segment::adjacent_in_segment(const Atom &atom, 
		const Atom &forbidden) const { 
		for (auto &adj : atom) {
			if (&adj != &forbidden && has_atom(adj)) 
				return adj; 
		}
		throw Error("die : couldn't find adjacent in segment");
	}
	
	Segment::Graph Segment::create_graph(const Molecule &molecule) {
		dbgmsg("Create segment graph ...");
		const Model &model = molecule.first().first();
		vector<unique_ptr<Segment>> vertices;
		for (auto &kv : model.get_rigid()) { // make vertices (segments) of a graph
			const string &nm = kv.first;
			const set<AtomSet> &atom_sets = kv.second;
			dbgmsg(nm);
			for (auto &atoms : atom_sets) {
				vertices.push_back(unique_ptr<Segment>(new Segment(atoms)));
				// If it has <=3 atoms (1 core and 2 join atoms) it's a linker, else it's a seed, 
				// meaning it was rigidly docked. Seeds are identified by non-empty names.
				if (atoms.size() > 3)
					vertices.back()->set_name(nm);
			}
		}
		// connect segments
		for (int i = 0; i < vertices.size(); ++i) {
			Segment &s1 = *vertices[i];
			for (int j = i + 1; j < vertices.size(); ++j) {
				Segment &s2 = *vertices[j];
				auto inter = Glib::intersection(s1.get_atoms(), s2.get_atoms());
				if (inter.size() == 2) {
					s1.add(&s2);
					s2.add(&s1);
					auto &atom1 = **inter.begin();
					auto &atom2 = **inter.rbegin();
					int num_bonds = 0; 
					for (auto &adj : atom1) {
						if (s1.has_atom(adj))
							num_bonds++;
					}
					if (num_bonds == 1) {
						dbgmsg("atom " << atom1 << " belongs to segment " << s2);
						dbgmsg("atom " << atom2 << " belongs to segment " << s1);
						s1.set_bond(s2, Bond(&atom2, &atom1));
						s2.set_bond(s1, Bond(&atom1, &atom2));
					} else {
						dbgmsg("atom " << atom1 << " belongs to segment " << s1);
						dbgmsg("atom " << atom2 << " belongs to segment " << s2);
						s1.set_bond(s2, Bond(&atom1, &atom2));
						s2.set_bond(s1, Bond(&atom2, &atom1));
					}
				}
			}
		}
		return Segment::Graph(std::move(vertices), true, false);
	}

};
