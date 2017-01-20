#include "common.hpp"
#include "helper/inout.hpp"
#include "pdbreader/grid.hpp"
#include "pdbreader/molecule.hpp"
#include "helper/benchmark.hpp"
#include "helper/path.hpp"
#include "helper/grep.hpp"
#include "score/score.hpp"
#include "pdbreader/pdbreader.hpp"
#include "pdbreader/atom.hpp"

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
			presidue->set_resn("ligand_" + std::to_string(ligand_cnt));
		}
	}

	/* Part7a stuff
	 * 
	 */
	Molib::NRset read_top_seeds_files(const std::set<std::string> &seeds, const string &top_seeds_dir, const string &top_seeds_file, const double top_percent) {
		Molib::NRset top_seeds;

		boost::regex regex;
		regex.assign("REMARK   5 MOLECULE ", boost::regex_constants::basic);

		for ( auto &fragment : seeds ) {

			boost::filesystem::path file_to_read(top_seeds_dir);
			file_to_read /= fragment;
			file_to_read /= top_seeds_file;
			std::ifstream file( file_to_read.c_str() );

			const size_t number_of_seeds = Grep::count_matches(file, regex);
			const int sz = static_cast<int> (number_of_seeds * top_percent);
			dbgmsg("taking " << sz << " top seeds for seed " << fragment);

			// Add one in case the user is silly enough to select a top_percent of 0.000
			Molib::PDBreader pdb( file_to_read.string(), Molib::PDBreader::all_models, sz + 1 );

			dbgmsg("reading top_seeds_file for seed id = " << fragment);
			Molib::Molecules all = pdb.parse_molecule();

			dbgmsg("number of top seeds left = " << all.size());

			Molib::Molecules &last = top_seeds.add(new Molib::Molecules(all));

			if (last.empty()) {
				throw Error("die : there are no docked conformations for seed " + fragment);
			}
		}

		return top_seeds;
	}

	Molib::NRset read_top_seeds_files(const Molib::Molecule &ligand, const string &top_seeds_dir, const string &top_seeds_file, const double top_percent) {
		std::set<string> seeds_to_read;
		const Molib::Model &model = ligand.first().first();
		for (auto &fragment : model.get_rigid()) { // iterate over seeds
			if (fragment.is_seed()) {
				seeds_to_read.insert(std::to_string(fragment.get_seed_id()));
			}
		}
		return read_top_seeds_files(seeds_to_read, top_seeds_dir, top_seeds_file, top_percent);
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
						Molib::Molecule &seed = seeds.add(new Molib::Molecule(std::to_string(fragment.get_seed_id())));
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
};
