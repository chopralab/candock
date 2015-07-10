#include "segment.hpp"
//~ #include "linker.hpp"
#include "pdbreader/molecule.hpp"
#include "helper/benchmark.hpp"
#include "helper/help.hpp"

namespace Molib {
	ostream& operator<< (ostream& stream, const Segment& s) {
		stream << "Segment(" << s.__seed_id << ") : atom numbers = ";
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
		dbgmsg(model.get_rigid());
		for (auto &fragment : model.get_rigid()) { // make vertices (segments) of a graph
			dbgmsg(fragment.get_all());
			vertices.push_back(unique_ptr<Segment>(new Segment(fragment.get_all(), fragment.get_seed_id())));
		}
		// connect segments
		for (int i = 0; i < vertices.size(); ++i) {
			Segment &s1 = *vertices[i];
			for (int j = i + 1; j < vertices.size(); ++j) {
				Segment &s2 = *vertices[j];
				auto inter = Glib::intersection(s1.get_atoms(), s2.get_atoms());
				dbgmsg(s1.get_atoms().size() << " " << s2.get_atoms().size() << " " << inter.size());
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
