#include "program/cmdlnopts.hpp"
#include "program/target.hpp"

#include "score/score.hpp"
#include "modeler/forcefield.hpp"
#include "modeler/modeler.hpp"

#include "version.hpp"
#include "drm/drm.hpp"

#include "fileout/fileout.hpp"

using namespace std;
using namespace Program;

////////////////// TEST SINGLEPOINT ENERGY CALCULATION OF COMPLEX ///////////////////////////

int main(int argc, char* argv[]) {
        try {
                if(!drm::check_drm(Version::get_install_path() + "/.candock")) {
                    throw logic_error("CANDOCK has expired. Please contact your CANDOCK distributor to get a new version.");
                }

                Inout::Logger::set_all_stderr(true);

                help::Options::set_options(new Program::CmdLnOpts(
                    argc, argv, Program::CmdLnOpts::STARTING |
                                Program::CmdLnOpts::FORCE_FIELD| 
                                Program::CmdLnOpts::SCORING));

                cerr << Version::get_banner()   <<
                        Version::get_version()  <<
                        Version::get_run_info();
                cerr << help::Options::get_options()->configuration_file() << endl;

                Molib::Molecules receptor_mols;
                Molib::Molecules ligand_mols;

                Parser::FileParser drpdb (cmdl.get_string_option("receptor"),
                                          Parser::pdb_read_options::protein_poses_only |
                                          Parser::pdb_read_options::all_models);

                drpdb.parse_molecule(receptor_mols);

                // Check to see if the user set the ligand option, use the receptor as a last resort.
                const std::string &ligand_file = Inout::file_size(cmdl.get_string_option("ligand")) == 0 ?
                                                    cmdl.get_string_option("receptor") :
                                                    cmdl.get_string_option("ligand");

                Parser::FileParser dlpdb (ligand_file,
                                          Parser::pdb_read_options::docked_poses_only |
                                          Parser::pdb_read_options::skip_atom |
                                          Parser::pdb_read_options::all_models);

                dlpdb.parse_molecule(ligand_mols);

                Score::Score score(cmdl.get_string_option("ref"), cmdl.get_string_option("comp"),
                                   cmdl.get_string_option("func"),cmdl.get_int_option("cutoff"),
                                   cmdl.get_double_option("step"));

                score.define_composition(receptor_mols.get_idatm_types(),
                                         ligand_mols.get_idatm_types())
                     .process_distributions_file(cmdl.get_string_option("dist"))
                     .compile_scoring_function()
                     .parse_objective_function(cmdl.get_string_option("obj_dir"), cmdl.get_double_option("scale"), 1501);

                OMMIface::ForceField ffield;

                ffield.parse_gaff_dat_file(cmdl.get_string_option("gaff_dat"))
                      .add_kb_forcefield(score, cmdl.get_double_option("step"), 15)
                      .parse_forcefield_file(cmdl.get_string_option("amber_xml"))
                      .parse_forcefield_file(cmdl.get_string_option("water_xml"));

                if ( ! cmdl.get_string_option("gaff_heme").empty() ) {
                        dbgmsg( "Adding " << cmdl.get_string_option("gaff_heme") << endl);
                        ffield.parse_gaff_dat_file(cmdl.get_string_option("gaff_heme"));
                }

                OMMIface::SystemTopology::loadPlugins();

                for (size_t i = 0; i < receptor_mols.size(); ++i) {
                        
                        Molib::Molecule& protein = receptor_mols[i];
                        Molib::Molecule& ligand  = ligand_mols[i];

                        Molib::Atom::Grid gridrec(protein.get_atoms());                        
                        protein.prepare_for_mm(ffield, gridrec);

                        ffield.insert_topology(protein);
                        ffield.insert_topology(ligand);

                        OMMIface::Modeler modeler (ffield, cmdl.get_string_option("fftype"),   cmdl.get_int_option("cutoff"),
                                                           cmdl.get_double_option("mini_tol"), cmdl.get_int_option("max_iter"),
                                                           cmdl.get_int_option("update_freq"), cmdl.get_double_option("pos_tol"), false,
                                                           cmdl.get_double_option("dynamic_step_size"),
                                                           cmdl.get_double_option("temperature"), cmdl.get_double_option("friction") );

                        modeler.add_topology (protein.get_atoms());
                        modeler.add_topology (ligand.get_atoms());

                        modeler.init_openmm();

                        modeler.add_crds (protein.get_atoms(), protein.get_crds());
                        modeler.add_crds (ligand.get_atoms(), ligand.get_crds());

                        modeler.init_openmm_positions();

                        modeler.minimize_state ();

                        // init with minimized coordinates
                        Molib::Molecule minimized_receptor (protein,modeler.get_state (protein.get_atoms()));
                        Molib::Molecule minimized_ligand   (ligand, modeler.get_state (ligand.get_atoms()));

                        minimized_receptor.undo_mm_specific();

                        Molib::Atom::Grid gridrec_min (minimized_receptor.get_atoms());
                        const double energy = score.non_bonded_energy (gridrec_min, minimized_ligand);

                        Fileout::print_complex_pdb (std::cout, minimized_ligand, minimized_receptor, energy, std::nan(""), i + 1, 0xFFFFFF, std::nan(""));

                        ffield.erase_topology(ligand);

                }

        } catch (exception& e) {
                cerr << e.what() << endl;
        }
        return 0;
}

