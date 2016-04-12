#include "common.hpp"
#include "helper/inout.hpp"
#include "pdbreader/grid.hpp"
#include "pdbreader/molecule.hpp"
#include "helper/benchmark.hpp"
#include "helper/path.hpp"
#include "geom3d/matrix.hpp"
#include "geom3d/pca.hpp"
#include "geom3d/geom3d.hpp"
#include "kabsch/kabsch.hpp"
#include "score/score.hpp"
#include "pdbreader/pdbreader.hpp"
#include "pdbreader/atom.hpp"
#include <gsl/gsl_matrix_double.h>
#include <gsl/gsl_statistics_double.h>
#include <gsl/gsl_eigen.h>
#include <gsl/gsl_blas.h>

namespace common {

	void change_residue_name(Molib::Molecule &ligand, const string &resn) {
		for (auto &presidue : ligand.get_residues()) {
			presidue->set_resn(resn);
		}
	}

	void change_residue_name(Molib::Molecule &ligand, std::mutex &mtx, int &ligand_cnt) {
		lock_guard<std::mutex> guard(mtx);
		++ligand_cnt;
		for (auto &presidue : ligand.get_residues()) {
			presidue->set_resn("ligand_" + help::to_string(ligand_cnt));
		}
	}

	/* Part7a stuff
	 * 
	 */
	Molib::NRset read_top_seeds_files(const Molib::Molecule &ligand, const string &top_seeds_dir, const string &top_seeds_file, const double top_percent) {
		Molib::NRset top_seeds;
		const Molib::Model &model = ligand.first().first();
		for (auto &fragment : model.get_rigid()) { // iterate over seeds
			if (fragment.is_seed()) {

				dbgmsg("reading top_seeds_file for seed id = " << fragment.get_seed_id());
				Molib::PDBreader pdb(Path::join(Path::join(top_seeds_dir, help::to_string(fragment.get_seed_id())), 
					top_seeds_file), Molib::PDBreader::all_models);
					
				Molib::Molecules all = pdb.parse_molecule();

				const int sz = (int) (all.size() * top_percent);
				dbgmsg("taking " << sz << " top seeds for seed " << fragment.get_seed_id());
				
				// deleting excess top seeds
				for (int i = all.size() - 1; i > sz; --i) all.erase(i);
				
				dbgmsg("number of top seeds left = " << all.size());
				
				Molib::Molecules &last = top_seeds.add(new Molib::Molecules(all));

				if (last.empty()) {
					throw Error("die : there are no docked conformations for seed " + help::to_string(fragment.get_seed_id()));
				}
			}
		}
		return top_seeds;
	}

	void create_mols_from_seeds(set<int> &added, Molib::Molecules &seeds, const Molib::Molecules &mols) {
		for (auto &molecule : mols)
		for (auto &assembly : molecule)
		for (auto &model : assembly) {
			for (auto &fragment : model.get_rigid()) { // iterate over seeds
				if (fragment.is_seed()) {
					dbgmsg("considering to add " << fragment.get_seed_id());
					if (!added.count(fragment.get_seed_id())) { // take seeds that haven't been docked already
						dbgmsg("added " << fragment.get_seed_id());
						added.insert(fragment.get_seed_id());
						// add to new molecules
						Molib::Molecule &seed = seeds.add(new Molib::Molecule(help::to_string(fragment.get_seed_id())));
						Molib::Assembly &a = seed.add(new Molib::Assembly(0));
						Molib::Model &mod = a.add(new Molib::Model(1));
						Molib::Chain &c = mod.add(new Molib::Chain('X'));
						Molib::Residue &r = c.add(new Molib::Residue("XXX", 1, ' ', Molib::Residue::hetero));
						for (const Molib::Atom *atom : fragment.get_all()) {
							Molib::Atom &at = r.add(new Molib::Atom(*atom));
							dbgmsg("added atom = " << at);
						}
						seed.regenerate_bonds(molecule);
					}
				}
			}
		}
	}
//~ #ifndef NDEBUG
	void create_mols_from_fragments(set<int> &added, Molib::Molecules &seeds, const Molib::Molecules &mols) {
		for (auto &molecule : mols)
		for (auto &assembly : molecule)
		for (auto &model : assembly) {
			for (auto &fragment : model.get_rigid()) { // iterate over seeds
				if (fragment.is_seed()) {
					dbgmsg("considering to add " << fragment.get_seed_id());
					if (!added.count(fragment.get_seed_id())) { // take seeds that haven't been docked already
						dbgmsg("added " << fragment.get_seed_id());
						added.insert(fragment.get_seed_id());
						// add to new molecules
						Molib::Molecule &seed = seeds.add(new Molib::Molecule(help::to_string(fragment.get_seed_id())));
						Molib::Assembly &a = seed.add(new Molib::Assembly(0));
						Molib::Model &mod = a.add(new Molib::Model(1));
						Molib::Chain &c = mod.add(new Molib::Chain('X'));
						Molib::Residue &r = c.add(new Molib::Residue("XXX", 1, ' ', Molib::Residue::hetero));
						for (const Molib::Atom *atom : fragment.get_all()) {
							Molib::Atom &at = r.add(new Molib::Atom(*atom));
							dbgmsg("added atom = " << at);
						}
						seed.regenerate_bonds(molecule);
					}
				}
			}
		}
	}
//~ #endif

};
