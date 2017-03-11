#include "cmdlnopts.hpp"

#include <boost/program_options/errors.hpp>

#include <iostream>
#include <fstream>
#include <string>

namespace po = boost::program_options;

namespace Program {

	void CmdLnOpts::__init (int argc, char *argv[], int opts_to_parse) {
 
		__program_name = argv[0];

		try {

			std::string config_file;

			po::options_description generic ("Generic options");
			generic.add_options()
			("help,h", "Show this help")
			("quiet,q", "Quiet mode (default is verbose)")
			("conifg,c", po::value<std::string> (&config_file)->default_value (""), "Configuration File")
			("ncpu",     po::value<int> (&__ncpu)             ->default_value(-1),
			 "Number of CPUs to use concurrently (use -1 to use all CPUs)")
			;

			po::options_description starting_inputs ("Starting Input Files");
			starting_inputs.add_options()
			("receptor", po::value<std::string> ()->default_value("receptor.pdb"),
			 "Receptor filename")
			("ligand",   po::value<std::string> ()  ->default_value("ligands.mol2"),
			 "Ligand filename")
			;

			po::options_description probis_options ("Probis (binding site indentification) Options");
			probis_options.add_options()
			("bslib", po::value<std::string> ()->default_value ("bslibdb/bslib.txt"),
			 "Read binding sites library from this file")
			("names", po::value<std::string> () ->default_value ("bslibdb/data/names"),
			 "Directory with ligand names")
			("bio",   po::value<std::string> ()   ->default_value ("bslibdb/data/bio"),
			 "Directory with ProBiS-ligands bio database")
			("nosql", po::value<std::string> ()->default_value ("probis.nosql"),
			 "NoSql-formatted ProBiS alignments output file")
			("json",  po::value<std::string> () ->default_value ("probis.json"),
			 "Json-formatted ProBiS alignments output file")
			("jsonwl", po::value<std::string> ()->default_value ("probis_with_ligands.json"),
			 "Json-formatted ProBiS alignments with transposed ligands output file")
			("lig_clus_file", po::value<std::string> ()->default_value ("ligand_clusters.pdb"),
			 "Ligand clusters found by ProBiS are outputted to this file")
			("z_scores_file", po::value<std::string> ()->default_value ("z_scores.pdb"),
			 "Binding site z-scores are outputted to this file")
			("probis_min_z_score", po::value<double> ()->default_value (2.5, "2.5"),
			 "Minimium z-score of ligands to be considered in clustering")
			("probis_min_pts", po::value<int> ()->default_value (10),
			 "The minimum number of points (for predicted ligands) required to form a cluster")
			("probis_clus_rad", po::value<double> ()->default_value (3.0, "3.0"),
			 "Cluster radius for predicted ligands by probis")
			("centro_clus_rad", po::value<double> ()->default_value (3.0, "3.0"),
			 "Cluster radius for centroid centers")
			("centroid", po::value<std::string> ()->default_value ("site.cen"),
			 "Filename for reading and writing centroids")
			("neighb", po::value<bool> ()    ->default_value (false, "false"),
			 "Allow only ligands that are in the similar regions according to REMARKs")
			("num_bsites", po::value<int> ()->default_value (3),
			 "Maximum number of predicted (or given) binding sites to consider for docking")
                        ("srf_file", po::value<std::string> ()->default_value("probis.srf"),
                         "File for storing the protein surface calculated by probis.")
			;

			po::options_description ligand_fragmention_options ("Ligand Fragmention Options");
			ligand_fragmention_options.add_options()
			("seeds", po::value<std::string> ()->default_value ("seeds.txt"),
			 "Read unique seeds from this file, if it exists, and append new unique seeds if found")
			("prep",  po::value<std::string> () ->default_value ("prepared_ligands.pdb"),
			 "Prepared small molecule(s) are outputted to this filename")
			("seeds_pdb", po::value<std::string> ()->default_value("seeds.pdb"),
			 "File to save full seeds into.")
			("max_num_ligands", po::value<int> ()->default_value (10),
			 "Maximum number of ligands to read in one chunk")
			;

			po::options_description frag_dock_options ("Fragment Docking Options");
			frag_dock_options.add_options()
			("top_seeds_dir",  po::value<std::string> ()->default_value (""),
			 "Directory for saving top docked seeds")
			("top_seeds_file",   po::value<std::string> ()->default_value ("top_seeds.pdb"),
			 "Top seeds output file")
			("gridpdb_hcp",    po::value<std::string> ()->default_value ("gridpdb_hcp.pdb"),
			 "Grid pdb hcp file for output")
			("max_frag_radius", po::value<double> ()->default_value (16.0, "16.0"),
			 "Maximum fragment radius for creating the initial rotamers")
			("grid",           po::value<double> ()->default_value (0.375, "0.375"),
			 "Grid spacing")
			("num_univec",     po::value<int> ()->default_value (256),
			 "Number of unit vectors evenly distributed on a sphere for conformation generation")
			("conf_spin",      po::value<double> ()->default_value (10, "10"),
			 "Spin degrees for conformation generation")
			("excluded",       po::value<double> ()->default_value (0.8, "0.8"),
			 "Excluded radius")
			("interatomic",    po::value<double> ()->default_value (8.0, "8.0"),
			 "Maximum interatomic distance")
			("clus_rad",       po::value<double> ()->default_value (2.0, "2.0"),
			 "Cluster radius for docked seeds")
			("clusterfile",    po::value<std::string> ()->default_value ("clustered_seeds.txt"),
			 "Clustered representative docked-seed conformations output file")
			;

			po::options_description scoring_options ("Scoring Function Arguments");
			scoring_options.add_options()
			("dist",    po::value<std::string> ()->default_value ("data/csd_complete_distance_distributions.txt"),
			 "Select one of the interatomic distance distribution file(s) provided with this program")
			("ref",     po::value<std::string> ()->default_value ("mean"),
			 "Normalization method for the reference state ('mean' is averaged over all atom type pairs, whereas 'cumulative' is a summation for atom type pairs)")
			("comp",    po::value<std::string> ()->default_value ("reduced"),
			 "Atom types used in calculating reference state 'reduced' or 'complete' ('reduced' includes only those atom types present in the specified receptor and small molecule, whereas 'complete' includes all atom types)")
			("func",    po::value<std::string> ()->default_value ("radial"),
			 "Function for calculating scores 'radial' or 'normalized_frequency'")
			("cutoff",  po::value<int> ()->default_value (6),
			 "Cutoff length [4-15].")
			("step",    po::value<double> ()->default_value (0.01, "0.01"),
			 "Step for spline generation of non-bonded knowledge-based potential [0.0-1.0]")
			("scale",   po::value<double> ()->default_value (10.0, "10.0"),
			 "Scale non-bonded forces and energy for knowledge-based potential [0.0-1000.0]")
			("obj_dir", po::value<std::string> ()->default_value ("obj"),
			 "Output directory for objective function and derivatives")
			("potential_file", po::value<std::string> ()->default_value ("potentials.txt"),
			 "Output file for potentials and derivatives")
			;

			po::options_description force_field_min ("Forcefield and Minimization Options");
			force_field_min.add_options()
			("amber_xml", po::value<std::string> ()->default_value ("data/amber10.xml"),
			 "Receptor XML parameters (and topology) input file")
			("water_xml", po::value<std::string> ()->default_value ("data/tip3p.xml"),
			 "Water XML parameters (and topology) input file")
			("gaff_dat",  po::value<std::string> () ->default_value ("data/gaff.dat"),
			 "Gaff DAT forcefield input file")
			("gaff_xml",  po::value<std::string> () ->default_value ("data/gaff.xml"),
			 "Gaff XML forcefield and ligand topology output file")
			("fftype"   , po::value<std::string> ()->default_value ("kb"),
			 "Forcefield to use 'kb' (knowledge-based) or 'phy' (physics-based)")
			("pos_tol",   po::value<double> ()->default_value (0.00000000001, "0.00000000001"),
			 "Minimization position tolerance in Angstroms - only for KB")
			("mini_tol",  po::value<double> ()->default_value (0.0001),
			 "Minimization tolerance")
			("max_iter",  po::value<int> ()->default_value (10),
			 "Maximum iterations for minimization during linking")
			("max_iter_final", po::value<int> ()->default_value (100),
			 "Maximum iterations for final minimization")
			("update_freq", po::value<int> ()->default_value (10),
			 "Update non-bond frequency")
			;

			po::options_description linking_step ("Fragment Linking Options");
			linking_step.add_options()
			("docked_dir",      po::value<std::string> ()->default_value("docked"),
			 "Docked ligands output directory")
			("iterative",       po::value<bool> ()->default_value (false, "false")->implicit_value(true),
			 "Enable iterative minimization during linking")
			("cuda",       po::value<bool> ()->default_value (false, "false")->implicit_value(true),
			 "Enable cuda iterative linker during linking")
			("top_percent",     po::value<double> ()->default_value (0.05, "0.05"),
			 "Top percent of each docked seed to extend to full molecule")
			("max_clique_size", po::value<int> ()->default_value (3),
			 "Maximum clique size for initial partial conformations generation")
			("spin",            po::value<int> ()->default_value (60),
			 "Spin degrees to rotate ligand. Allowed values are 5, 10, 15, 20, 30, 60, 90")
			("clash_coeff",     po::value<double> ()->default_value (0.75, "0.75"),
			 "Clash coefficient for determining whether two atoms clash by eq. dist12 s< C * (vdw1 + vdw2)")
			("tol_seed_dist",   po::value<double> ()->default_value (2.0, "2.0"),
			 "Tolerance on seed distance in-between linking")
			("lower_tol_seed_dist", po::value<double> ()->default_value (2.0, "2.0"),
			 "Lower tolerance on seed distance for getting initial conformations of docked fragments")
			("upper_tol_seed_dist", po::value<double> ()->default_value (2.0, "2.0"),
			 "Upper tolerance on seed distance for getting initial conformations of docked fragments")
			("max_possible_conf",   po::value<int> ()->default_value (-1),
			 "Maximum number of possible conformations to link (-1 means unlimited)")
			("link_iter",           po::value<int> ()->default_value (1000),
			 "Maximum iterations for linking procedure")
			("docked_clus_rad",     po::value<double> ()->default_value (2.0, "2.0"),
			 "Cluster radius between docked ligand conformations")
			("max_allow_energy",    po::value<double> ()->default_value (0.0, "0.0"),
			 "Maximum allowed energy for seed conformations")
			("max_num_possibles",   po::value<int> ()->default_value (200000),
			 "Maximum number of possibles conformations considered for clustering")
			;

			po::options_description design_step ("Automated Design Options");
			design_step.add_options()
			("target_dir",          po::value<std::string>()->default_value("")->implicit_value("targets"),
			 "Directory containing PDB files. These are docked against and labeled as targets. ")
			("antitarget_dir",      po::value<std::string>()->default_value("")->implicit_value("atargets"),
			 "Directory containing PDB files. These are docked against and labeled as antitargets")
			("target_linking",      po::value<bool>()->default_value(true,"true"),
			 "Should the ligands be linked for target")
			("antitarget_linking",  po::value<bool>()->default_value(true,"true"),
			 "Shoutd the ligands be linked for antitargets")
			("fragment_bag",        po::value<std::string>()->default_value("")->implicit_value("fragment_bag.mol2"),
			 "Additional fragments to be added to seeds.pdb")
			("fragment_mol",        po::value<std::string>()->default_value("")->implicit_value("fragment_mol.mol2"),
			 "Additional fragments to be added to seeds.pdb without rotatable bonds being cut.")
			("seeds_to_add",        po::value<int>()->default_value(50),
			 "Number of seeds from seeds.pdb to be considered for addition to the ligands in prepared_ligands.pdb")
			("seeds_to_avoid",      po::value<int>()->default_value(50),
			 "Number of seeds from seeds.pdb to be considered for removal from determined from seeds_to_add")
			("seeds_till_good",     po::value<int>()->default_value(-1),
			 "Number of times a seed must be present in the top_seeds for targets until it is considered for addition")
			("seeds_till_bad",      po::value<int>()->default_value(-1),
			 "Number of times a seed must be present in the top_seeds for antitargets until it is removed from the good list")
			("force_seed",          po::value< std::vector<std::string>>()->multitoken()->default_value(std::vector<std::string>(),""),
			 "Force addition of a certain seed from seeds.pdb. Multiple seeds can be given")
			("add_single_atoms",    po::value< std::vector<std::string>>()->multitoken()->default_value(std::vector<std::string>(),""),
			 "Change hydrogens to given atoms. Multiple atoms can be given.")
			("change_terminal_atom",po::value< std::vector<std::string>>()->multitoken()->default_value(std::vector<std::string>(),""),
			 "Change non-hydrogen atoms that terminate chains to given atoms. Multiple atoms can be given.")
			;

			po::options_description config_options;
			config_options.add(starting_inputs)
			              .add(probis_options)
			              .add(ligand_fragmention_options)
			              .add(frag_dock_options)
			              .add(scoring_options)
			              .add(force_field_min)
			              .add(linking_step)
			              .add(design_step);

			po::options_description cmdln_options;
			cmdln_options.add (generic);
			cmdln_options.add (config_options);

			po::store (po::parse_command_line (argc, argv, cmdln_options), __vm);
			po::notify (__vm);

			if (__vm.count ("help")) {

				po::options_description print_options;
				print_options.add(generic);

				if (opts_to_parse & STARTING ) {
					print_options.add(starting_inputs);
				}

				if (opts_to_parse & PROBIS ) {
					print_options.add(probis_options);
				} 

				if (opts_to_parse & LIG_FRAMGENT) {
					print_options.add(ligand_fragmention_options);
				} 

				if (opts_to_parse & FRAG_DOCKING) {
					print_options.add(frag_dock_options);
				}

				if (opts_to_parse & SCORING) {
					print_options.add (scoring_options);
				} 

				if (opts_to_parse & FORCE_FIELD ) {
					print_options.add (force_field_min);
				} 

				if (opts_to_parse & LINKING ) {
					print_options.add (linking_step);
				}

				if (opts_to_parse & DESIGN ) {
					print_options.add (design_step);
				}

				std::cout << print_options << std::endl;
				exit (0);
			}

			if (!config_file.empty()) {
				std::ifstream config_stream (config_file.c_str());

				if (!config_stream) {
					throw Error ("Unable to open conifguration file!\n");
				}

				store (parse_config_file (config_stream, config_options), __vm);
				notify (__vm);
			}

			po::store(po::parse_environment(config_options,"CANDOCK_"), __vm);
			po::notify(__vm);

			if (__ncpu == -1) {
				__ncpu = thread::hardware_concurrency();
			}

			const string &fftype = get_string_option("fftype");
			if (fftype != "kb" && fftype != "phy") {
				throw po::validation_error (po::validation_error::invalid_option_value,
				                            "fftype", fftype);
			}

			const string &ref_state = get_string_option("ref");
			if (ref_state != "mean" && ref_state != "cumulative") {
				throw po::validation_error (po::validation_error::invalid_option_value,
				                            "ref", ref_state);
			}

			const string &comp = get_string_option("comp");
			if (comp != "reduced" && comp != "complete") {
				throw po::validation_error (po::validation_error::invalid_option_value,
				                            "comp", comp);
			}

			const string &func = get_string_option("func");
			if (func != "radial" && func != "normalized_frequency") {
				throw po::validation_error (po::validation_error::invalid_option_value,
				                            "func", func);
			}

			const int cutoff = get_int_option("cutoff");
			if (cutoff < 4 || cutoff > 15) {
				throw po::validation_error (po::validation_error::invalid_option_value,
				                            "cutoff", std::to_string (cutoff));
			}

			const double step = get_double_option("step");
			if (step < 0.0 || step > 1.0) {
				throw po::validation_error (po::validation_error::invalid_option_value,
				                            "step", std::to_string (step));
			}

			const double scale = get_double_option("scale");
			if (scale < 0.0 || scale > 1000.0) {
				throw po::validation_error (po::validation_error::invalid_option_value,
				                            "scale", std::to_string (scale));
			}

			const int spin_degrees = get_int_option("spin");
			switch (spin_degrees) {
			case  5:
			case 10:
			case 15:
			case 20:
			case 30:
			case 60:
			case 90:
				break;

			default:
				throw po::validation_error (po::validation_error::invalid_option_value,
				                            "spin", std::to_string (spin_degrees));
			}

		} catch (po::error &e) {
			std::cerr << "error: " << e.what() << std::endl;
			std::cerr << "Please see the help (-h) for more information" << std::endl;
			throw Error ("die: arguments error");
		}
	}

        const std::string& CmdLnOpts::get_string_option(const std::string& option) const {
                return __vm[option].as<std::string>();
        }

        int CmdLnOpts::get_int_option           (const std::string& option) const {
                return __vm[option].as<int>();
        }

        double CmdLnOpts::get_double_option     (const std::string& option) const {
                return __vm[option].as<double>();
        }

        bool CmdLnOpts::get_bool_option         (const std::string& option) const {
                return __vm[option].as<bool>();
        }

        const std::vector<std::string>& CmdLnOpts::get_string_vector (const std::string& option) const {
                return __vm[option].as<std::vector<std::string>>();
        }

        std::ostream & operator<<(std::ostream &stream, const CmdLnOpts &cmdl_) {
                for ( const auto& a : cmdl_.__vm ) {
                        stream << std::setw(22)<< a.first << " = ";
                        if        ( auto v = boost::any_cast<std::string>(&a.second.value()) ) {
                                stream << std::setw(47) << *v;
                        } else if ( auto v = boost::any_cast<int>(&a.second.value()) ) {
                                stream << std::setw(47) << *v;
                        } else if ( auto v = boost::any_cast<double>(&a.second.value()) ) {
                                stream << std::setw(47) << *v;
                        } else if ( auto v = boost::any_cast<bool>(&a.second.value()) ) {
                                stream << std::setw(47) << *v;
                        } else if ( auto v = boost::any_cast<std::vector<std::string>>(&a.second.value()) ) {
                                std::string combination("");
                                for ( const auto& s : *v ) {
                                        combination += s;
                                        combination += ", ";
                                }
                                stream << std::setw(47) << combination;
                        }

                        stream << std::endl;
                }

                return stream;
        }
}
