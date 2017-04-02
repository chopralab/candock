#include <memory>
#include <iostream>
#include <set>
#include <map>
#include <string>
#include <vector>
#include <cstdlib>
#include "geom3d/geom3d.hpp"
#include "atom.hpp"
#include "residue.hpp"
#include "chain.hpp"
#include "model.hpp"
#include "assembly.hpp"
using namespace std;

namespace Molib {
	ostream& operator<< (ostream& stream, const Assembly& a) {
		stream << setw(6) << left << "REMARK   6 " << a.name() << " " << a.number() << endl;
		for (auto &model : a) { 
			stream << model;
		}
		stream << "REMARK   6 END" << endl;
		return stream;
	}

	void Assembly::init_bio(const Assembly &asym, map<int, Geom3D::Matrix> &matrices, const set<char> &chains) {
		for (auto &i : matrices) {
			int matrix_number = i.first;
			Geom3D::Matrix &matrix = i.second;
			Model &last = add(new Model(matrix_number));
			last.init_bio(asym.first(), matrix, chains); // asymmetric unit has just one model - the 1-st
		}
	}

	void Assembly::rotate(const Geom3D::Matrix &rota, const bool inverse) {
		for (auto &model : *this) {
			model.rotate(rota, inverse);
		}
	}

	Atom::Vec Assembly::get_atoms(const string &chain_ids, const Residue::res_type &rest,
		const int model_number) const {
		Atom::Vec atoms;
		for (auto &model : *this) {
			if (model_number == -1 || model.number() == model_number) {
				auto ret = model.get_atoms(chain_ids, rest);
				atoms.insert(atoms.end(), ret.begin(), ret.end());
			}
		}
		return atoms;
	}


};
