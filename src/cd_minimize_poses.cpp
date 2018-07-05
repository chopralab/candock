#include "candock/program/cmdlnopts.hpp"
#include "candock/program/target.hpp"

#include "candock/score/kbff.hpp"
#include "candock/modeler/forcefield.hpp"
#include "candock/modeler/modeler.hpp"

#include "version.hpp"
#include "candock/drm/drm.hpp"

#include "candock/fileout/fileout.hpp"

#include "candock/modeler/systemtopology.hpp"

using namespace std;
using namespace candock;
using namespace candock::Program;

////////////////// TEST SINGLEPOINT ENERGY CALCULATION OF COMPLEX ///////////////////////////

int main(int argc, char *argv[])
{
        try
        {
                if (!drm::check_drm(Version::get_install_path() + "/.candock"))
                {
                        throw logic_error("CANDOCK has expired. Please contact your CANDOCK distributor to get a new version.");
                }

                Inout::Logger::set_all_stderr(true);

                help::Options::set_options(new Program::CmdLnOpts(
                    argc, argv, Program::CmdLnOpts::STARTING | Program::CmdLnOpts::FORCE_FIELD | Program::CmdLnOpts::SCORING));

                cerr << Version::get_banner() << Version::get_version() << Version::get_run_info();
                cerr << help::Options::get_options()->configuration_file() << endl;

                molib::Molecules receptor_mols;
                molib::Molecules ligand_mols;

                Parser::FileParser drpdb(cmdl.get_string_option("receptor"),
                                         Parser::pdb_read_options::protein_poses_only |
                                             Parser::pdb_read_options::all_models);

                drpdb.parse_molecule(receptor_mols);

                // Check to see if the user set the ligand option, use the receptor as a last resort.
                const std::string &ligand_file = Inout::file_size(cmdl.get_string_option("ligand")) == 0 ? cmdl.get_string_option("receptor") : cmdl.get_string_option("ligand");

                Parser::FileParser dlpdb(ligand_file,
                                         Parser::pdb_read_options::docked_poses_only |
                                             Parser::pdb_read_options::skip_atom |
                                             Parser::pdb_read_options::all_models);

                dlpdb.parse_molecule(ligand_mols);

                score::KBFF score(cmdl.get_string_option("ff_ref"),  cmdl.get_string_option("ff_comp"),
                                   cmdl.get_string_option("ff_func"), cmdl.get_int_option("ff_cutoff"),
                                   cmdl.get_double_option("step"));

                score.define_composition(receptor_mols.get_idatm_types(),
                                         ligand_mols.get_idatm_types())
                    .process_distributions_file(cmdl.get_string_option("dist"))
                    .compile_scoring_function();

                if (cmdl.get_string_option("obj_dir").empty())
                {
                        score.compile_objective_function(cmdl.get_double_option("scale"));
                }
                else
                {
                        score.parse_objective_function(cmdl.get_string_option("obj_dir"),
                                                       cmdl.get_double_option("scale"),
                                                       15);
                }

                OMMIface::ForceField ffield;

                ffield.parse_gaff_dat_file(cmdl.get_string_option("gaff_dat"))
                .add_kb_forcefield(score, cmdl.get_double_option("dist_cutoff"))
                .parse_forcefield_file(cmdl.get_string_option("amber_xml"))
                .parse_forcefield_file(cmdl.get_string_option("water_xml"));

                if (!cmdl.get_string_option("gaff_heme").empty())
                {
                        dbgmsg("Adding " << cmdl.get_string_option("gaff_heme") << endl);
                        ffield.parse_gaff_dat_file(cmdl.get_string_option("gaff_heme"));
                }

                OMMIface::SystemTopology::loadPlugins();

                for (size_t i = 0; i < receptor_mols.size(); ++i)
                {

                        molib::Molecule &protein = receptor_mols[i];
                        molib::Molecule &ligand = ligand_mols[i];

                        molib::Atom::Grid gridrec(protein.get_atoms());
                        protein.prepare_for_mm(ffield, gridrec);

                        ffield.insert_topology(protein);
                        ffield.insert_topology(ligand);

                        OMMIface::Modeler modeler(
                                ffield,cmdl.get_string_option("fftype"),
                                cmdl.get_double_option("mini_tol"),
                                cmdl.get_int_option("max_iter"),
                                cmdl.get_double_option("pos_tol")
                        );

                        modeler.add_topology(protein.get_atoms());
                        modeler.add_topology(ligand.get_atoms());

                        modeler.init_openmm(cmdl.get_string_option("platform"),
                                            cmdl.get_string_option("precision"),
                                            cmdl.get_string_option("accelerators"),
                                            OMMIface::SystemTopology::integrator_type::none);

                        modeler.add_crds(protein.get_atoms(), protein.get_crds());
                        modeler.add_crds(ligand.get_atoms(), ligand.get_crds());

                        modeler.unmask(ligand.get_atoms());
                        modeler.unmask(protein.get_atoms());

                        modeler.init_openmm_positions();

                        modeler.minimize_state();

                        // init with minimized coordinates
                        molib::Molecule minimized_receptor(protein, modeler.get_state(protein.get_atoms()));
                        molib::Molecule minimized_ligand(ligand, modeler.get_state(ligand.get_atoms()));

                        minimized_receptor.undo_mm_specific();

                        molib::Atom::Grid gridrec_min(minimized_receptor.get_atoms());
                        const double energy = score.non_bonded_energy(gridrec_min, minimized_ligand);
                        modeler.mask(protein.get_atoms());
                        const double potential = modeler.potential_energy();

                        fileout::print_complex_pdb(
                                std::cout,
                                minimized_ligand,
                                minimized_receptor,
                                energy,
                                potential,
                                i + 1,
                                0xFFFFFF,
                                std::nan("")
                        );

                        ffield.erase_topology(ligand);

                        // delete *system;
                }
        }
        catch (exception &e)
        {
                cerr << e.what() << endl;
        }
        return 0;
}
