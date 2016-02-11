#ifndef OPTS_CANDOCK_H
#define OPTS_CANDOCK_H
#include <tclap/CmdLine.h>
#include <boost/regex.hpp>
#include <boost/algorithm/string.hpp>
#include <boost/filesystem.hpp>
#include <boost/date_time/posix_time/posix_time.hpp>
#include <boost/algorithm/string/split.hpp>
#include <boost/algorithm/string.hpp>
#include <boost/asio/ip/host_name.hpp>
#include "helper/error.hpp"
#include <thread>
using namespace std;

class CmdLnOpts {
	string __receptor_file;
	string __ligand_file;
	string __seeds_file;
	string __prep_file;
	
	string __bslib_file;
	string __nosql_file;
	string __json_file;
	string __json_with_ligs_file;

	string __top_seeds_dir;

	string __names_dir;
	bool __neighb;
	double __probis_clus_rad;
	int __probis_min_pts;
	double __probis_min_z_score;
	
	
	string __bio_dir;
	string __lig_clus_file;
	string __z_scores_file;
	double __centro_clus_rad;
	
	string __centroid_in_file;
	string __centroid_out_file;

	string __gridpdb_hcp_file;
	double __max_frag_radius;

	string __ref_state;
	string __comp;
	string __rad_or_raw;
	int __dist_cutoff;
	string __distributions_file;
	double __step_non_bond;
	double __scale_non_bond;

	string __potential_file;
	string __obj_dir;

	string __cluster_file;
	string __top_seeds_file;

	int __max_top_seeds;
	
	int __max_num_ligands;
	
	string __gaff_dat_file;
	string __gaff_xml_file;

	int __spin_degrees;
	double __clash_coeff;
	double __tol_seed_dist;
	double __lower_tol_seed_dist;
	double __upper_tol_seed_dist;
	int __max_possible_conf;
	int __link_iter;

	string __docked_file;
	
	double __docked_clus_rad;
	double __max_allow_energy;
	bool __iterative;
	int __max_num_possibles;
	
	string __amber_xml_file;
	string __water_xml_file;
	string __fftype;
	double __tolerance;
	int __max_iterations;
	int __max_iterations_final;
	int __update_freq;
	double __position_tolerance;

	double __top_percent;
	int __max_clique_size;
	
	int __num_bsites;
	double __grid_spacing;
	int __num_univec;
	double __conf_spin;
	double __excluded_radius;
	double __max_interatomic_distance;

	double __clus_rad;

	int __ncpu;
	
	bool __quiet;
	string __program_name;
	string __version;
public:
	CmdLnOpts() : __quiet(false), __version("1.00") {}
	void init(int argc, char* argv[]) {
		__program_name = argv[0];
		try {			
			TCLAP::CmdLine cmd("Command description message", ' ', __version);
			TCLAP::SwitchArg quietSwitch("q","quiet","Quiet mode (default is verbose)", cmd, false);

			TCLAP::ValueArg<string> receptor_fileArg("","receptor","Receptor filename",true,"default","string", cmd);
			TCLAP::ValueArg<string> ligand_fileArg("","ligand","Ligand filename",true,"default","string", cmd);
			TCLAP::ValueArg<string> seeds_fileArg("","seeds","Read unique seeds from this file, if it exists, and append new \
				unique seeds if found",false,"seeds.txt","string", cmd);
			TCLAP::ValueArg<string> prep_fileArg("","prep","Prepared small molecule(s) are outputted to this filename",
				false,"prepared_ligands.pdb","string", cmd);

			TCLAP::ValueArg<string> bslib_fileArg("","bslib","Read binding sites library from this file (default is bslibdb/bslib.txt)",
				false,"bslibdb/bslib.txt","string", cmd);
			TCLAP::ValueArg<string> nosql_fileArg("","nosql","NoSql-formatted ProBiS alignments output file (default is probis.nosql)",false,"probis.nosql","string", cmd);
			TCLAP::ValueArg<string> json_fileArg("","json","Json-formatted ProBiS alignments output file (default is probis.json",false,"probis.json","string", cmd);
			TCLAP::ValueArg<string> json_with_ligs_fileArg("","jsonwl","Json-formatted ProBiS alignments with transposed ligands output file (default is probis_with_ligands.json",false,"probis_with_ligands.json","string", cmd);

			TCLAP::ValueArg<string> top_seeds_dirArg("","top_seeds_dir","Directory for saving top docked seeds (default is top_seeds)",false,"top_seeds","string", cmd);

			TCLAP::ValueArg<string> names_dirArg("","names","Directory with ligand names (default is bslibdb/names)",false,"bslibdb/names","string", cmd);
			TCLAP::SwitchArg neighbSwitch("","neighb","Allow only ligands that are in the similar regions according to REMARKs (not enabled by default)", cmd, false);

			TCLAP::ValueArg<double> probis_clus_radArg("","probis_clus_rad","Cluster radius for predicted ligands by probis (default is 3.0)",false,3.0,"double", cmd);
			TCLAP::ValueArg<int> probis_min_ptsArg("","probis_min_pts","The minimum number of points (for predicted ligands) required to form a cluster (default is 10)",false,10,"int", cmd);
			TCLAP::ValueArg<double> probis_min_z_scoreArg("","probis_min_z_score","Minimium z-score of ligands to be considered in clustering (default is 2.5)",false,2.5,"double", cmd);

			TCLAP::ValueArg<string> bio_dirArg("","bio","Directory with ProBiS-ligands bio database (default is bslibdb/bio)",false,"bslibdb/bio","string", cmd);
			TCLAP::ValueArg<string> lig_clus_fileArg("","lig_clus_file","Ligand clusters found by ProBiS are outputted to this file",false,"ligand_clusters.pdb","string", cmd);
			TCLAP::ValueArg<string> z_scores_fileArg("","z_scores_file","Binding site z-scores are outputted to this file",false,"z_scores.pdb","string", cmd);

			TCLAP::ValueArg<double> centro_clus_radArg("","centro_clus_rad","Cluster radius for centroid centers (default is 3.0)",false,3.0,"double", cmd);

			TCLAP::ValueArg<string> centroid_in_fileArg("","centro_in","Filename for reading centroids",false,"","string", cmd);
			TCLAP::ValueArg<string> centroid_out_fileArg("","centro_out","Filename for outputting calculated centroids",false,"","string", cmd);

			TCLAP::ValueArg<string> gridpdb_hcp_fileArg("","gridpdb_hcp","Grid pdb hcp file for output",false,"gridpdb_hcp.pdb","string", cmd);
			
			TCLAP::ValueArg<double> max_frag_radiusArg("","max_frag_radius","Maximum fragment radius for creating the initial rotamers (default is 16.0)",false,16.0,"double", cmd);

			vector<string> allowedRef{"mean","cumulative"};
			TCLAP::ValuesConstraint<string> allowedValsRef( allowedRef );
			TCLAP::ValueArg<string> ref_stateArg("","ref", "Normalization method for the reference state ('mean' (default) is averaged over all atom type pairs, whereas 'cumulative' is a summation for atom type pairs)",false,"mean",&allowedValsRef, cmd);
			vector<string> allowedComp{"reduced","complete"};
			TCLAP::ValuesConstraint<string> allowedValsComp( allowedComp );
			TCLAP::ValueArg<string> compArg("","comp", "Atom types used in calculating reference state 'reduced' (default) or 'complete' \
				('reduced' includes only those atom types present in the specified receptor and small molecule, whereas	'complete' \
				includes all atom types)",false,"reduced",&allowedValsComp, cmd);
			vector<string> allowedRad{"radial","normalized_frequency"};
			TCLAP::ValuesConstraint<string> allowedValsRad( allowedRad );
			TCLAP::ValueArg<string> rad_or_rawArg("","func","Function for calculating scores 'radial' (default) or \
				'normalized_frequency'",false,"radial",&allowedValsRad, cmd);
			vector<int> allowedDistCutoff{4,5,6,7,8,9,10,11,12,13,14,15};
			TCLAP::ValuesConstraint<int> allowedValsDistCutoff( allowedDistCutoff );
			
			TCLAP::ValueArg<int> dist_cutoffArg("","cutoff","Cutoff length (default is 6)",false,6,&allowedValsDistCutoff, cmd);
			TCLAP::ValueArg<string> distributions_fileArg("","dist","Select one of the interatomic distance distribution \
				file(s) provided with this script",false,"data/csd_complete_distance_distributions.txt","string", cmd);
			TCLAP::ValueArg<double> step_non_bondArg("","step","Step for spline generation of non-bonded knowledge-based \
				potential [0.0-1.0] (default is 0.01)",false,0.01,"double", cmd);
			TCLAP::ValueArg<double> scale_non_bondArg("","scale","Scale non-bonded forces and energy for knowledge-based \
				potential [0.0-1000.0] (default is 10.0)",false,10.0,"double", cmd);

			TCLAP::ValueArg<string> potential_fileArg("","potential_file","Output file for potentials and derivatives",false,"potentials.txt","string", cmd);
			TCLAP::ValueArg<string> obj_dirArg("","obj_dir","Output directory for objective function and derivatives",false,"obj","string", cmd);

			TCLAP::ValueArg<string> cluster_fileArg("","clusterfile","Clustered representative docked-seed conformations \
				output file",false,"clustered_seeds.txt","string", cmd);
			TCLAP::ValueArg<string> top_seeds_fileArg("","topseedsfile","Top seeds \
				output file",false,"top_seeds.pdb","string", cmd);

			TCLAP::ValueArg<int> max_top_seedsArg("","max_top_seeds","Maximum number of top seeds to read for linking (default is 20)",false,20,"int", cmd);
			
			TCLAP::ValueArg<int> max_num_ligandsArg("","max_num_ligands","Maximum number of ligands to read in one chunk (default is 10)",false,10,"int", cmd);

			TCLAP::ValueArg<string> gaff_dat_fileArg("","gaff_dat","Gaff DAT forcefield input file",false,"data/gaff.dat","string", cmd);
			TCLAP::ValueArg<string> gaff_xml_fileArg("","gaff_xml","Gaff XML forcefield and ligand topology output file",false,"gaff.xml","string", cmd);

			vector<int> allowedSpinDegrees{5,10,15,20,30,60,90};
			TCLAP::ValuesConstraint<int> allowedValsSpinDegrees( allowedSpinDegrees );
			TCLAP::ValueArg<int> spin_degreesArg("","spin","Spin degrees to rotate ligand (default 60)",false,60,&allowedValsSpinDegrees, cmd);

			TCLAP::ValueArg<double> clash_coeffArg("","clash_coeff","Clash coefficient for determining whether two atoms clash by eq. dist12 s< C * (vdw1 + vdw2) (default is 0.75)",false,0.75,"double", cmd);
			TCLAP::ValueArg<double> tol_seed_distArg("","tol_seed_dist","Tolerance on seed distance in-between linking (default is 2.0)",false,2.0,"double", cmd);
			TCLAP::ValueArg<double> lower_tol_seed_distArg("","lower_tol_seed_dist","Lower tolerance on seed distance for getting initial conformations of docked fragments (default is 2.0)",false,2.0,"double", cmd);
			TCLAP::ValueArg<double> upper_tol_seed_distArg("","upper_tol_seed_dist","Upper tolerance on seed distance for getting initial conformations of docked fragments (default is 2.0)",false,2.0,"double", cmd);

			TCLAP::ValueArg<int> max_possible_confArg("","max_possible_conf","Maximum number of possible conformations to link (default is 20, -1 means unlimited)",false,20,"int", cmd);
			TCLAP::ValueArg<int> link_iterArg("","link_iter","Maximum iterations for linking procedure (default is 1000)",false,1000,"int", cmd);

			TCLAP::ValueArg<string> docked_fileArg("","docked_file","Docked ligands output filename",false,"docked.pdb","string", cmd);
			
			TCLAP::ValueArg<double> docked_clus_radArg("","docked_clus_rad","Cluster radius between docked ligand conformations (default is 2.0)",false,2.0,"double", cmd);
			TCLAP::ValueArg<double> max_allow_energyArg("","max_allow_energy","Maximum allowed energy for seed conformations (default is 0.0)",false,0.0,"double", cmd);
			TCLAP::SwitchArg iterativeSwitch("","iterative","Enable iterative minimization during linking (not enabled by default)", cmd, false);
			TCLAP::ValueArg<int> max_num_possiblesArg("","max_num_possibles","Maximum number of possibles conformations considered for clustering (default is 200000)",false,200000,"int", cmd);

			TCLAP::ValueArg<string> amber_xml_fileArg("","amber_xml","Receptor XML parameters (and topology) input file",false,"data/amber10.xml","string", cmd);
			TCLAP::ValueArg<string> water_xml_fileArg("","water_xml","Water XML parameters (and topology) input file",false,"data/tip3p.xml","string", cmd);
			vector<string> allowedFf{"kb","phy"};
			TCLAP::ValuesConstraint<string> allowedValsFf( allowedFf );
			TCLAP::ValueArg<string> fftypeArg("","fftype","Forcefield to use 'kb' (knowledge-based, default) or 'phy' (physics-based)",false,"kb",&allowedValsFf, cmd);
			TCLAP::ValueArg<double> toleranceArg("","mini_tol","Minimization tolerance (default is 0.0001)",false,0.0001,"double", cmd);
			TCLAP::ValueArg<int> max_iterationsArg("","max_iter","Maximum iterations for minimization during linking (default is 100)",false,100,"int", cmd);
			TCLAP::ValueArg<int> max_iterations_finalArg("","max_iter_final","Maximum iterations for final minimization (default is 100)",false,100,"int", cmd);
			TCLAP::ValueArg<int> update_freqArg("","update_freq","Update non-bond frequency (default is 10)",false,10,"int", cmd);
			TCLAP::ValueArg<double> position_toleranceArg("","pos_tol","Minimization position tolerance in Angstroms - only for KB (default is 0.00000000001)",false, 0.00000000001,"double", cmd);
			
			TCLAP::ValueArg<double> top_percentArg("","top_percent","Top percent of each docked seed to extend to full molecule - only applied to iterative minimization (default is 0.05)",false, 0.05,"double", cmd);
			TCLAP::ValueArg<int> max_clique_sizeArg("","max_clique_size","Maximum clique size for initial partial conformations generation (default is 3)",false,3,"int", cmd);

			TCLAP::ValueArg<int> num_bsitesArg("","num_bsites","Maximum number of predicted (or given) binding sites to consider for docking (default is 3)",false,3,"int", cmd);
			
			TCLAP::ValueArg<double> grid_spacingArg("","grid","Grid spacing (default is 0.5)",false,0.5,"double", cmd);
			TCLAP::ValueArg<int> num_univecArg("","num_univec","Number of unit vectors evenly distributed on a sphere for conformation generation (default is 256)",false,256,"int", cmd);
			TCLAP::ValueArg<double> conf_spinArg("","conf_spin","Spin degrees for conformation generation (default is 10)",false,10,"double", cmd);

			TCLAP::ValueArg<double> excluded_radiusArg("","excluded","Excluded radius (default is 0.8)",false,0.8,"double", cmd);
			TCLAP::ValueArg<double> max_interatomic_distanceArg("","interatomic","Maximum interatomic distance (default is 8.0)",false,8.0,"double", cmd);

			TCLAP::ValueArg<double> clus_radArg("","clus_rad","Cluster radius for docked seeds (default is 2.0)",false,2.0,"double", cmd);

			TCLAP::ValueArg<int> ncpuArg("","ncpu","Number of CPUs to use concurrently (default is 1)",false,1,"int", cmd);

			cmd.parse( argc, argv );
			__quiet = quietSwitch.getValue();

			__receptor_file = receptor_fileArg.getValue();
			__ligand_file = ligand_fileArg.getValue();
			__seeds_file = seeds_fileArg.getValue();
			__prep_file = prep_fileArg.getValue();
			__bslib_file = bslib_fileArg.getValue();
			__nosql_file = nosql_fileArg.getValue();
			__json_file = json_fileArg.getValue();
			__json_with_ligs_file = json_with_ligs_fileArg.getValue();
			__top_seeds_dir = top_seeds_dirArg.getValue();
			__names_dir = names_dirArg.getValue();
			__neighb = neighbSwitch.getValue();
			__probis_clus_rad = probis_clus_radArg.getValue();
			__probis_min_pts = probis_min_ptsArg.getValue();
			__probis_min_z_score = probis_min_z_scoreArg.getValue();
			__bio_dir = bio_dirArg.getValue();
			__lig_clus_file = lig_clus_fileArg.getValue();
			__z_scores_file = z_scores_fileArg.getValue();

			__centro_clus_rad = centro_clus_radArg.getValue();

			__centroid_in_file = centroid_in_fileArg.getValue();
			__centroid_out_file = centroid_out_fileArg.getValue();
			__gridpdb_hcp_file = gridpdb_hcp_fileArg.getValue();
			__max_frag_radius = max_frag_radiusArg.getValue();
			
			__ref_state = ref_stateArg.getValue();
			__comp = compArg.getValue();
			__rad_or_raw = rad_or_rawArg.getValue();
			__dist_cutoff = dist_cutoffArg.getValue();
			__distributions_file = distributions_fileArg.getValue();
			__step_non_bond = step_non_bondArg.getValue();
			__scale_non_bond = scale_non_bondArg.getValue();

			__potential_file = potential_fileArg.getValue();
			__obj_dir = obj_dirArg.getValue();

			__cluster_file = cluster_fileArg.getValue();
			__top_seeds_file = top_seeds_fileArg.getValue();
			
			__max_top_seeds = max_top_seedsArg.getValue();
			
			__max_num_ligands = max_num_ligandsArg.getValue();
			__gaff_dat_file = gaff_dat_fileArg.getValue();
			__gaff_xml_file = gaff_xml_fileArg.getValue();
			__spin_degrees = spin_degreesArg.getValue();
			
			__clash_coeff = clash_coeffArg.getValue();
			__tol_seed_dist = tol_seed_distArg.getValue();
			__lower_tol_seed_dist = lower_tol_seed_distArg.getValue();
			__upper_tol_seed_dist = upper_tol_seed_distArg.getValue();
			__max_possible_conf = max_possible_confArg.getValue();
			__link_iter = link_iterArg.getValue();
			__docked_file = docked_fileArg.getValue();
			
			__docked_clus_rad = docked_clus_radArg.getValue();
			__max_allow_energy = max_allow_energyArg.getValue();
			__iterative = iterativeSwitch.getValue();
			__max_num_possibles = max_num_possiblesArg.getValue();
			
			__amber_xml_file = amber_xml_fileArg.getValue();
			__water_xml_file = water_xml_fileArg.getValue();
			__fftype = fftypeArg.getValue();
			__tolerance = toleranceArg.getValue();
			__max_iterations = max_iterationsArg.getValue();
			__max_iterations_final = max_iterations_finalArg.getValue();
			__update_freq = update_freqArg.getValue();
			__position_tolerance = position_toleranceArg.getValue();

			__top_percent = top_percentArg.getValue();
			__max_clique_size = max_clique_sizeArg.getValue();
			
			__num_bsites = num_bsitesArg.getValue();
			__grid_spacing = grid_spacingArg.getValue();
			__num_univec = num_univecArg.getValue();
			__conf_spin = conf_spinArg.getValue();
			__excluded_radius = excluded_radiusArg.getValue();
			__max_interatomic_distance = max_interatomic_distanceArg.getValue();

			__clus_rad = clus_radArg.getValue();

			__ncpu = ncpuArg.getValue();
		} 
		catch (TCLAP::ArgException &e) { 
			cerr << "error: " << e.error() << " for arg " << e.argId() << endl; 
			throw Error("die: arguments error\n");
		}		
	}
	void display_time(string what) {
		cout << what << " on " << boost::posix_time::to_simple_string(boost::posix_time::second_clock::local_time()) << "\n";
		cout << "running " << __program_name << " version " << __version << " on hostname " << boost::asio::ip::host_name() << "\n";
	}
	// interface
	bool quiet() const { return __quiet; }

	string receptor_file() const { return __receptor_file; }
	string ligand_file() const { return __ligand_file; }
	string seeds_file() const { return __seeds_file; }
	string prep_file() const { return __prep_file; }
	string bslib_file() const { return __bslib_file; }
	string nosql_file() const { return __nosql_file; }
	string json_file() const { return __json_file; }
	string json_with_ligs_file() const { return __json_with_ligs_file; }
	string top_seeds_dir() const { return __top_seeds_dir; }
	string names_dir() const { return __names_dir; }
	bool neighb() const { return __neighb; }
	double probis_clus_rad() const { return __probis_clus_rad; }
	int probis_min_pts() const { return __probis_min_pts; }
	double probis_min_z_score() const { return __probis_min_z_score; }
	string bio_dir() const { return __bio_dir; }
	string lig_clus_file() const { return __lig_clus_file; }
	string z_scores_file() const { return __z_scores_file; }
	double centro_clus_rad() const { return __centro_clus_rad; }
	string centroid_in_file() const { return __centroid_in_file; }
	string centroid_out_file() const { return __centroid_out_file; }
	string gridpdb_hcp_file() const { return __gridpdb_hcp_file; }
	double max_frag_radius() const { return __max_frag_radius; }
	
	string ref_state() const { return __ref_state; }
	string comp() const { return __comp; }
	string rad_or_raw() const { return __rad_or_raw; }
	int dist_cutoff() const { return __dist_cutoff; }
	string distributions_file() const { return __distributions_file; }
	double step_non_bond() const { return __step_non_bond; }
	double scale_non_bond() const { return __scale_non_bond; }

	string potential_file() const { return __potential_file; }
	string obj_dir() const { return __obj_dir; }

	string cluster_file() const { return __cluster_file; }
	string top_seeds_file() const { return __top_seeds_file; }

	int max_top_seeds() const { return __max_top_seeds; }

	int max_num_ligands() const { return __max_num_ligands; }
	string gaff_dat_file() const { return __gaff_dat_file; }
	string gaff_xml_file() const { return __gaff_xml_file; }
	int spin_degrees() const { return __spin_degrees; }
	
	double clash_coeff() const { return __clash_coeff; }
	double tol_seed_dist() const { return __tol_seed_dist; }
	double lower_tol_seed_dist() const { return __lower_tol_seed_dist; }
	double upper_tol_seed_dist() const { return __upper_tol_seed_dist; }
	int max_possible_conf() const { return __max_possible_conf; }
	int link_iter() const { return __link_iter; }
	string docked_file() const { return __docked_file; }
	
	double docked_clus_rad() const { return __docked_clus_rad; }
	double max_allow_energy() const { return __max_allow_energy; }
	bool iterative() const { return __iterative; }
	int max_num_possibles() const { return __max_num_possibles; }
	
	string amber_xml_file() const { return __amber_xml_file; }
	string water_xml_file() const { return __water_xml_file; }
	string fftype() const { return __fftype; }
	double tolerance() const { return __tolerance; }
	int max_iterations() const { return __max_iterations; }
	int max_iterations_final() const { return __max_iterations_final; }
	int update_freq() const { return __update_freq; }
	double position_tolerance() const { return __position_tolerance; }
	
	double top_percent() const { return __top_percent; }
	int max_clique_size() const { return __max_clique_size; }
	
	int num_bsites() const { return __num_bsites; }
	double grid_spacing() const { return __grid_spacing; }
	int num_univec() const { return __num_univec; }
	double conf_spin() const { return __conf_spin; }
	double excluded_radius() const { return __excluded_radius; }
	double max_interatomic_distance() const { return __max_interatomic_distance; }

	double clus_rad() const { return __clus_rad; }

	int ncpu() const { return __ncpu; }
	
	friend ostream& operator<< (ostream& stream, const CmdLnOpts &cmdl) {
		unsigned int n = thread::hardware_concurrency();
		stream << endl << "Detected support for " << n << " concurrent threads." << endl;
		return stream;
	}

};
extern CmdLnOpts cmdl;
#endif
