#include "candock/program/cmdlnopts.hpp"
#include "candock/program/target.hpp"

#include "candock/score/score.hpp"
#include "candock/modeler/forcefield.hpp"
#include "candock/modeler/modeler.hpp"

#include "version.hpp"
#include "candock/drm/drm.hpp"

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

                molib::Molecules starting_mol;

                parser::FileParser drpdb(cmdl.get_string_option("receptor"),
                                         parser::pdb_read_options::first_model);

                drpdb.parse_molecule(starting_mol);

                score::KBFF score(cmdl.get_string_option("ff_ref"), cmdl.get_string_option("ff_comp"),
                                   cmdl.get_string_option("ff_func"), cmdl.get_int_option("ff_cutoff"),
                                   cmdl.get_double_option("step"));

                score.define_composition(starting_mol.get_idatm_types(),
                                         starting_mol.get_idatm_types())
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

                molib::Molecule &starting = starting_mol[0];

                molib::Atom::Grid gridrec(starting.get_atoms());
                starting_mol[0].prepare_for_mm(ffield, gridrec);

                ffield.insert_topology(starting);

                OMMIface::Modeler modeler(ffield,
                        cmdl.get_string_option("fftype"),
                        cmdl.get_double_option("mini_tol"),
                        cmdl.get_int_option("max_iter"),
                        false,
                        cmdl.get_double_option("dynamic_step_size"),
                        cmdl.get_double_option("temperature"),
                        cmdl.get_double_option("friction")
                );

                modeler.add_topology(starting.get_atoms());

                if (cmdl.get_string_option("integrator") == "verlet")
                {
                        modeler.init_openmm(cmdl.get_string_option("platform"),
                                            cmdl.get_string_option("precision"),
                                            cmdl.get_string_option("accelerators"),
                                            OMMIface::SystemTopology::integrator_type::verlet);
                }
                else if (cmdl.get_string_option("integrator") == "langevin")
                {
                        modeler.init_openmm(cmdl.get_string_option("platform"),
                                            cmdl.get_string_option("precision"),
                                            cmdl.get_string_option("accelerators"),
                                            OMMIface::SystemTopology::integrator_type::langevin);
                }
                else if (cmdl.get_string_option("integrator") == "brownian")
                {
                        modeler.init_openmm(cmdl.get_string_option("platform"),
                                            cmdl.get_string_option("precision"),
                                            cmdl.get_string_option("accelerators"),
                                            OMMIface::SystemTopology::integrator_type::brownian);
                }
                else
                {
                        throw Error("You need to select a valid integrator for dynamics");
                }

                modeler.add_crds(starting.get_atoms(), starting.get_crds());

                modeler.init_openmm_positions();

                modeler.set_num_steps_to_run(cmdl.get_int_option("dynamic_steps"));
                modeler.dynamics();

                // init with minimized coordinates
                molib::Molecule final_state(starting, modeler.get_state(starting.get_atoms()));

                final_state.undo_mm_specific();

                std::cout << final_state;
        }
        catch (exception &e)
        {
                cerr << e.what() << endl;
        }
        return 0;
}
