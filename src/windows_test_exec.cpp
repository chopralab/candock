#include "program/cmdlnopts.hpp"
#include "program/fragmentligands.hpp"
#include "program/target.hpp"

#include "version.hpp"
using namespace std;
using namespace Program;

////////////////// TEST SINGLEPOINT ENERGY CALCULATION OF COMPLEX ///////////////////////////

int main(int argc, char* argv[]) {
        try {

                help::Options::set_options(new Program::CmdLnOpts(
                    argc, argv, Program::CmdLnOpts::STARTING |
                                Program::CmdLnOpts::LIG_FRAMGENT| 
                                Program::CmdLnOpts::SCORING));

                Benchmark main_timer;
                main_timer.display_time("Starting");

                cout << Version::get_banner()   <<
                        Version::get_version()  <<
                        Version::get_run_info() <<
                        help::Options::get_options()->configuration_file() << endl;

                Molib::PDBreader drpdb(cmdl.get_string_option("receptor"), Molib::PDBreader::first_model | Molib::PDBreader::skip_hetatm, 1);
                Molib::PDBreader dlpdb(cmdl.get_string_option("prep"), Molib::PDBreader::first_model | Molib::PDBreader::skip_atom, 1);

                Molib::Molecules receptor, ligand;
                if (drpdb.parse_molecule(receptor) && dlpdb.parse_molecule(ligand)) {

                        receptor.compute_idatm_type()
                                .compute_hydrogen()
                                .compute_bond_order()
                                .compute_bond_gaff_type()
                                .refine_idatm_type()
                                .erase_hydrogen()  // needed because refine changes connectivities
                                .compute_hydrogen()   // needed because refine changes connectivities
                                .compute_ring_type()
                                .compute_gaff_type()
                                .compute_rotatable_bonds() // relies on hydrogens being assigned
                                .erase_hydrogen();

                        Molib::Atom::Grid gridrec (receptor.get_atoms());

                        Molib::Score score (cmdl.get_string_option("ref"), cmdl.get_string_option("comp"),
                                            cmdl.get_string_option("func"),cmdl.get_int_option("cutoff"),
                                            cmdl.get_double_option("step"));

                        score.define_composition( receptor.get_idatm_types(),
                                                  ligand.get_idatm_types())
                             .process_distributions_file(cmdl.get_string_option("dist"))
                             .compile_scoring_function()
                             .parse_objective_function(cmdl.get_string_option("obj_dir"), cmdl.get_double_option("scale"));

                        OMMIface::ForceField ffield;

                        ffield.parse_gaff_dat_file(cmdl.get_string_option("gaff_dat"))
                                .add_kb_forcefield(score, cmdl.get_double_option("step"))
                                .parse_forcefield_file(cmdl.get_string_option("amber_xml"))
                                .parse_forcefield_file(cmdl.get_string_option("water_xml"));

                        receptor[0].prepare_for_mm(ffield, gridrec);

                        const double energy = score.non_bonded_energy(gridrec, ligand[0]);

                        receptor.clear();
                        ligand.clear();

                        cout << "Energy is " << energy << endl;
                }

        } catch (exception& e) {
                cerr << e.what() << endl;
        }
        return 0;
}
