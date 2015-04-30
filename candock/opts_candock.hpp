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
	string __receptor_chain_id;
	string __seeds_file;
	string __prep_file;
	
	string __bslib_file;
	string __nosql_file;
	string __json_file;
	string __json_with_ligs_file;
	string __geo_dir;
	string __names_dir;
	bool __neighb;
	double __probis_clus_rad;
	int __probis_min_pts;
	double __probis_min_z_score;
	
	
	string __bio_dir;
	string __lig_clus_file;
	
	string __centroid_file;

	string __gridpdb_hcp_file;
	//~ string __gridxyz_hcp_file;

	string __ref_state;
	string __comp;
	string __rad_or_raw;
	int __dist_cutoff;
	string __distributions_file;
	double __step_non_bond;
	double __scale_non_bond;

	string __egrid_file;
	string __docked_seeds_file;
	//~ string __score_file;

	//~ string __rmsd_file;
	string __cluster_file;
	string __top_seeds_file;

	int __max_num_ligands;
	
	string __gaff_dat_file;
	string __gaff_xml_file;

	int __spin_degrees;
	double __tol_dist;
	double __tol_max_coeff;
	double __tol_min_coeff;
	int __max_possible_conf;
	int __link_iter;

	string __docked_ligands_file;
	//~ string __docked_rmsd_file;
	double __docked_clus_rad;
	//~ string __docked_score_file;
	int __docked_max_num_clus;
	int __docked_min_pts;

	string __amber_xml_file;
	string __fftype;
	double __tolerance;
	int __max_iterations;
	
	string __mini_ligands_file;
	string __energy_file;
	
	double __def_radial_check;
	//~ double __max_radial_check;
	int __num_bsites;
	double __grid_spacing;
	double __excluded_radius;
	double __max_interatomic_distance;

	double __top_percent;

	double __clus_rad;
	int __min_pts;
	int __max_num_clus;
	int __max_seeds_to_cluster;

	int __num_iter;

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
			TCLAP::ValueArg<string> receptor_chain_idArg("","receptor_chain_id","Chain id(s) of the receptor (default is A)",
				false,"A","string", cmd);
			TCLAP::ValueArg<string> seeds_fileArg("","seeds","Read unique seeds from this file, if it exists, and append new \
				unique seeds if found",false,"seeds.txt","string", cmd);
			TCLAP::ValueArg<string> prep_fileArg("","prep","Prepared small molecule(s) are outputted to this filename",
				false,"prepared_ligands.pdb","string", cmd);

			TCLAP::ValueArg<string> bslib_fileArg("","bslib","Read binding sites library from this file (default is data/probis_ligands/bslib.txt)",
				false,"data/probis_ligands/bslib.txt","string", cmd);
			TCLAP::ValueArg<string> nosql_fileArg("","nosql","NoSql-formatted ProBiS alignments output file (default is probis.nosql)",false,"probis.nosql","string", cmd);
			TCLAP::ValueArg<string> json_fileArg("","json","Json-formatted ProBiS alignments output file (default is probis.json",false,"probis.json","string", cmd);
			TCLAP::ValueArg<string> json_with_ligs_fileArg("","jsonwl","Json-formatted ProBiS alignments with transposed ligands output file (default is probis_with_ligands.json",false,"probis_with_ligands.json","string", cmd);
			TCLAP::ValueArg<string> geo_dirArg("","geo","Directory with ProBiS-ligands geo database (default is data/probis_ligands/geo)",false,"data/probis_ligands/geo","string", cmd);
			TCLAP::ValueArg<string> names_dirArg("","names","Directory with ligand names (default is data/probis_ligands/names)",false,"data/probis_ligands/names","string", cmd);
			TCLAP::SwitchArg neighbSwitch("","neighb","Allow only ligands that are in the similar regions according to REMARKs (not enabled by default)", cmd, false);

			TCLAP::ValueArg<double> probis_clus_radArg("","probis_clus_rad","Cluster radius for predicted ligands by probis (default is 2.0)",false,2.0,"double", cmd);
			TCLAP::ValueArg<int> probis_min_ptsArg("","probis_min_pts","The minimum number of points (for predicted ligands) required to form a cluster (default is 10)",false,10,"int", cmd);
			TCLAP::ValueArg<double> probis_min_z_scoreArg("","probis_min_z_score","Minimium z-score of ligands to be considered in clustering (default is 2.5)",false,2.5,"double", cmd);

			TCLAP::ValueArg<string> bio_dirArg("","bio","Directory with ProBiS-ligands bio database (default is data/probis_ligands/bio)",false,"data/probis_ligands/bio","string", cmd);
			TCLAP::ValueArg<string> lig_clus_fileArg("","lig_clus_file","Ligand clusters found by ProBiS are outputted to this file",false,"ligand_clusters.pdb","string", cmd);

			TCLAP::ValueArg<string> centroid_fileArg("","centroid","Centroid filename",false,"","string", cmd);

			TCLAP::ValueArg<string> gridpdb_hcp_fileArg("","gridpdb_hcp","Grid pdb hcp file for output",false,"gridpdb_hcp.pdb","string", cmd);
			//~ TCLAP::ValueArg<string> gridxyz_hcp_fileArg("","gridxyz_hcp","Grid xyz hcp file for output",false,"gridxyz_hcp.xyz","string", cmd);

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
			//~ TCLAP::ValueArg<int> dist_cutoffArg("","cutoff","Cutoff length (default is 6)",false,6,&allowedValsDistCutoff, cmd);
			TCLAP::ValueArg<int> dist_cutoffArg("","cutoff","Cutoff length (default is 8)",false,8,&allowedValsDistCutoff, cmd);
			TCLAP::ValueArg<string> distributions_fileArg("","dist","Select one of the interatomic distance distribution \
				file(s) provided with this script",false,"data/csd_complete_distance_distributions.txt","string", cmd);
			TCLAP::ValueArg<double> step_non_bondArg("","step","Step for spline generation of non-bonded knowledge-based \
				potential [0.0-1.0] (default is 0.1)",false,0.1,"double", cmd);
			TCLAP::ValueArg<double> scale_non_bondArg("","scale","Scale non-bonded forces and energy for knowledge-based \
				potential [0.0-1.0] (default is 0.01)",false,0.01,"double", cmd);

			TCLAP::ValueArg<string> egrid_fileArg("","egrid", "Energy grid output file",false,"energy_grid.txt","string", cmd);
			TCLAP::ValueArg<string> docked_seeds_fileArg("","docked_seeds_file", "Docked seeds output file",false,
				"docked_seeds.pdb","string", cmd);
			//~ TCLAP::ValueArg<string> score_fileArg("","scorefile", "Scores for docked seeds output file",false,
				//~ "scored_docked_seeds.txt","string", cmd);

			//~ TCLAP::ValueArg<string> rmsd_fileArg("","rmsdfile","RMSDs between all docked-seed conformations output \
				//~ file",false,"all_all_rmsd.txt","string", cmd);
			TCLAP::ValueArg<string> cluster_fileArg("","clusterfile","Clustered representative docked-seed conformations \
				output file",false,"clustered_seeds.txt","string", cmd);
			TCLAP::ValueArg<string> top_seeds_fileArg("","topseedsfile","Top seeds \
				output file",false,"top_seeds.pdb","string", cmd);

			TCLAP::ValueArg<int> max_num_ligandsArg("","max_num_ligands","Maximum number of ligands to read in one chunk (default is 1000)",false,1000,"int", cmd);

			TCLAP::ValueArg<string> gaff_dat_fileArg("","gaff_dat","Gaff DAT forcefield input file",false,"data/gaff.dat","string", cmd);
			TCLAP::ValueArg<string> gaff_xml_fileArg("","gaff_xml","Gaff XML forcefield and ligand topology output file",false,"gaff.xml","string", cmd);

			vector<int> allowedSpinDegrees{5,10,15,20,30};
			TCLAP::ValuesConstraint<int> allowedValsSpinDegrees( allowedSpinDegrees );
			TCLAP::ValueArg<int> spin_degreesArg("","spin","Spin degrees to rotate ligand (default 15)",false,15,&allowedValsSpinDegrees, cmd);
			TCLAP::ValueArg<double> tol_distArg("","tol_dist","Allowed distance between last path- and goal- state in linker (default is 3.0)",false,3.0,"double", cmd);

			TCLAP::ValueArg<double> tol_max_coeffArg("","tol_max_coeff","Multiply maximum linker length (upper bound) (default is 1.3)",false,1.3,"double", cmd);
			TCLAP::ValueArg<double> tol_min_coeffArg("","tol_min_coeff","Multiply maximum linker length (lower bound) (default is 0.7)",false,0.7,"double", cmd);

			TCLAP::ValueArg<int> max_possible_confArg("","max_possible_conf","Maximum number of possible conformations to link (default is -1 [no limit])",false,-1,"int", cmd);
			TCLAP::ValueArg<int> link_iterArg("","link_iter","Maximum iterations for linking procedure (default is 1000)",false,1000,"int", cmd);

			TCLAP::ValueArg<string> docked_ligands_fileArg("","docked_ligands","Docked ligands output filename",false,"docked_ligands.pdb","string", cmd);
			//~ TCLAP::ValueArg<string> docked_rmsd_fileArg("","docked_rmsd_file","RMSDs between all docked ligand conformations output filename",false,"docked_ligands_rmsd.txt","string", cmd);
			TCLAP::ValueArg<double> docked_clus_radArg("","docked_clus_rad","Cluster radius between docked ligand conformations (default is 2.0)",false,2.0,"double", cmd);
			//~ TCLAP::ValueArg<string> docked_score_fileArg("","docked_scorefile", "Docked ligands scores (before minimization) output file",false,
				//~ "docked_ligands_scores.txt","string", cmd);
			TCLAP::ValueArg<int> docked_max_num_clusArg("","docked_max_num_clus","Maximum number of clustered docked ligands to consider further for minimization (default is 10)",false,10,"int", cmd);
			TCLAP::ValueArg<int> docked_min_ptsArg("","docked_min_pts","The minimum number of points (docked ligands) required to form a cluster (default is 3)",false,3,"int", cmd);

			TCLAP::ValueArg<string> amber_xml_fileArg("","amber_xml","Receptor XML parameters (and topology) input file",false,"data/amber10.xml","string", cmd);
			vector<string> allowedFf{"kb","phy"};
			TCLAP::ValuesConstraint<string> allowedValsFf( allowedFf );
			TCLAP::ValueArg<string> fftypeArg("","fftype","Forcefield to use 'kb' (knowledge-based, default) or 'phy' (physics-based)",false,"kb",&allowedValsFf, cmd);
			TCLAP::ValueArg<double> toleranceArg("","mini_tol","Minimization tolerance (default is 0.1)",false,0.1,"double", cmd);
			TCLAP::ValueArg<int> max_iterationsArg("","max_iter","Maximum iterations for minimization (default is 0 - meaning minimize until convergence is reached)",false,0,"int", cmd);

			TCLAP::ValueArg<string> mini_ligands_fileArg("","mini_ligands","Docked & minimized ligands output filename (default minimized.pdb)",false,"minimized.pdb","string", cmd);
			TCLAP::ValueArg<string> energy_fileArg("","energy","Energies of minimized ligands output filename (default energies.txt)",false,"energies.txt","string", cmd);

			TCLAP::ValueArg<double> def_radial_checkArg("","radial","Radial check - the 2*r of the sphere around the binding site central point (default is 10.0)",false,10.0,"double", cmd);
			//~ TCLAP::ValueArg<double> max_radial_checkArg("","max_radial","Maximum radial check - the 2*r of the sphere around the binding site central point (default is 20.0)",false,20.0,"double", cmd);
			TCLAP::ValueArg<int> num_bsitesArg("","num_bsites","Maximum number of predicted (or given) binding sites to consider for docking (default is 3)",false,3,"int", cmd);

			//~ TCLAP::ValueArg<double> grid_spacingArg("","grid","Grid spacing (default is 0.375)",false,0.375,"double", cmd);
			TCLAP::ValueArg<double> grid_spacingArg("","grid","Grid spacing (default is 2.0)",false,2.0,"double", cmd);
			TCLAP::ValueArg<double> excluded_radiusArg("","excluded","Excluded radius (default is 0.8)",false,0.8,"double", cmd);
			TCLAP::ValueArg<double> max_interatomic_distanceArg("","interatomic","Maximum interatomic distance (default is 8.0)",false,8.0,"double", cmd);

			TCLAP::ValueArg<double> top_percentArg("","top_percent","Percent of top scores to keep (default is 0.80)",false,0.80,"double", cmd);

			TCLAP::ValueArg<double> clus_radArg("","clus_rad","Cluster radius for docked seeds (default is 1.0)",false,1.0,"double", cmd);
			TCLAP::ValueArg<int> min_ptsArg("","min_pts","The minimum number of points (docked seeds) required to form a cluster (default is 5)",false,5,"int", cmd);
			TCLAP::ValueArg<int> max_num_clusArg("","max_num_clus","Maximum number of clustered seeds to consider further (default is 500)",false,500,"int", cmd);
			TCLAP::ValueArg<int> max_seeds_to_clusterArg("","max_seeds_to_cluster","Maximum number of seeds to cluster (default is 2000)",false,2000,"int", cmd);

			TCLAP::ValueArg<int> num_iterArg("","num_iter","Number of iterations for maximum weight clique algorithm (default is 1000)",false,1000,"int", cmd);

			TCLAP::ValueArg<int> ncpuArg("","ncpu","Number of CPUs to use concurrently (default is 1)",false,1,"int", cmd);

			cmd.parse( argc, argv );
			__quiet = quietSwitch.getValue();

			__receptor_file = receptor_fileArg.getValue();
			__ligand_file = ligand_fileArg.getValue();
			__receptor_chain_id = receptor_chain_idArg.getValue();
			__seeds_file = seeds_fileArg.getValue();
			__prep_file = prep_fileArg.getValue();
			__bslib_file = bslib_fileArg.getValue();
			__nosql_file = nosql_fileArg.getValue();
			__json_file = json_fileArg.getValue();
			__json_with_ligs_file = json_with_ligs_fileArg.getValue();
			__geo_dir = geo_dirArg.getValue();
			__names_dir = names_dirArg.getValue();
			__neighb = neighbSwitch.getValue();
			__probis_clus_rad = probis_clus_radArg.getValue();
			__probis_min_pts = probis_min_ptsArg.getValue();
			__probis_min_z_score = probis_min_z_scoreArg.getValue();
			__bio_dir = bio_dirArg.getValue();
			__lig_clus_file = lig_clus_fileArg.getValue();
			__centroid_file = centroid_fileArg.getValue();
			__gridpdb_hcp_file = gridpdb_hcp_fileArg.getValue();
			//~ __gridxyz_hcp_file = gridxyz_hcp_fileArg.getValue();
			__ref_state = ref_stateArg.getValue();
			__comp = compArg.getValue();
			__rad_or_raw = rad_or_rawArg.getValue();
			__dist_cutoff = dist_cutoffArg.getValue();
			__distributions_file = distributions_fileArg.getValue();
			__step_non_bond = step_non_bondArg.getValue();
			__scale_non_bond = scale_non_bondArg.getValue();
			__egrid_file = egrid_fileArg.getValue();
			__docked_seeds_file = docked_seeds_fileArg.getValue();
			//~ __score_file = score_fileArg.getValue();
			//~ __rmsd_file = rmsd_fileArg.getValue();
			__cluster_file = cluster_fileArg.getValue();
			__top_seeds_file = top_seeds_fileArg.getValue();
			__max_num_ligands = max_num_ligandsArg.getValue();
			__gaff_dat_file = gaff_dat_fileArg.getValue();
			__gaff_xml_file = gaff_xml_fileArg.getValue();
			__spin_degrees = spin_degreesArg.getValue();
			__tol_dist = tol_distArg.getValue();
			__tol_max_coeff = tol_max_coeffArg.getValue();
			__tol_min_coeff = tol_min_coeffArg.getValue();
			__max_possible_conf = max_possible_confArg.getValue();
			__link_iter = link_iterArg.getValue();
			__docked_ligands_file = docked_ligands_fileArg.getValue();
			//~ __docked_rmsd_file = docked_rmsd_fileArg.getValue();
			__docked_clus_rad = docked_clus_radArg.getValue();
			//~ __docked_score_file = docked_score_fileArg.getValue();
			__docked_max_num_clus = docked_max_num_clusArg.getValue();
			__docked_min_pts = docked_min_ptsArg.getValue();
			__amber_xml_file = amber_xml_fileArg.getValue();
			__fftype = fftypeArg.getValue();
			__tolerance = toleranceArg.getValue();
			__max_iterations = max_iterationsArg.getValue();
			__mini_ligands_file = mini_ligands_fileArg.getValue();
			__energy_file = energy_fileArg.getValue();
			__def_radial_check = def_radial_checkArg.getValue();
			//~ __max_radial_check = max_radial_checkArg.getValue();
			__num_bsites = num_bsitesArg.getValue();
			__grid_spacing = grid_spacingArg.getValue();
			__excluded_radius = excluded_radiusArg.getValue();
			__max_interatomic_distance = max_interatomic_distanceArg.getValue();
			__top_percent = top_percentArg.getValue();
			__clus_rad = clus_radArg.getValue();
			__min_pts = min_ptsArg.getValue();
			__max_num_clus = max_num_clusArg.getValue();
			__max_seeds_to_cluster = max_seeds_to_clusterArg.getValue();
			__num_iter = num_iterArg.getValue();
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
	string receptor_chain_id() const { return __receptor_chain_id; }
	string seeds_file() const { return __seeds_file; }
	string prep_file() const { return __prep_file; }
	string bslib_file() const { return __bslib_file; }
	string nosql_file() const { return __nosql_file; }
	string json_file() const { return __json_file; }
	string json_with_ligs_file() const { return __json_with_ligs_file; }
	string geo_dir() const { return __geo_dir; }
	string names_dir() const { return __names_dir; }
	bool neighb() const { return __neighb; }
	double probis_clus_rad() const { return __probis_clus_rad; }
	int probis_min_pts() const { return __probis_min_pts; }
	double probis_min_z_score() const { return __probis_min_z_score; }
	string bio_dir() const { return __bio_dir; }
	string lig_clus_file() const { return __lig_clus_file; }
	string centroid_file() const { return __centroid_file; }
	string gridpdb_hcp_file() const { return __gridpdb_hcp_file; }
	//~ string gridxyz_hcp_file() const { return __gridxyz_hcp_file; }
	string ref_state() const { return __ref_state; }
	string comp() const { return __comp; }
	string rad_or_raw() const { return __rad_or_raw; }
	int dist_cutoff() const { return __dist_cutoff; }
	string distributions_file() const { return __distributions_file; }
	double step_non_bond() const { return __step_non_bond; }
	double scale_non_bond() const { return __scale_non_bond; }
	string egrid_file() const { return __egrid_file; }
	string docked_seeds_file() const { return __docked_seeds_file; }
	//~ string score_file() const { return __score_file; }
	//~ string rmsd_file() const { return __rmsd_file; }
	string cluster_file() const { return __cluster_file; }
	string top_seeds_file() const { return __top_seeds_file; }
	int max_num_ligands() const { return __max_num_ligands; }
	string gaff_dat_file() const { return __gaff_dat_file; }
	string gaff_xml_file() const { return __gaff_xml_file; }
	int spin_degrees() const { return __spin_degrees; }
	double tol_dist() const { return __tol_dist; }
	double tol_max_coeff() const { return __tol_max_coeff; }
	double tol_min_coeff() const { return __tol_min_coeff; }
	int max_possible_conf() const { return __max_possible_conf; }
	int link_iter() const { return __link_iter; }
	string docked_ligands_file() const { return __docked_ligands_file; }
	//~ string docked_rmsd_file() const { return __docked_rmsd_file; }
	double docked_clus_rad() const { return __docked_clus_rad; }
	//~ string docked_score_file() const { return __docked_score_file; }
	int docked_max_num_clus() const { return __docked_max_num_clus; }
	int docked_min_pts() const { return __docked_min_pts; }
	string amber_xml_file() const { return __amber_xml_file; }
	string fftype() const { return __fftype; }
	double tolerance() const { return __tolerance; }
	int max_iterations() const { return __max_iterations; }
	string mini_ligands_file() const { return __mini_ligands_file; }
	string energy_file() const { return __energy_file; }
	double def_radial_check() const { return __def_radial_check; }
	//~ double max_radial_check() const { return __max_radial_check; }
	int num_bsites() const { return __num_bsites; }
	double grid_spacing() const { return __grid_spacing; }
	double excluded_radius() const { return __excluded_radius; }
	double max_interatomic_distance() const { return __max_interatomic_distance; }
	double top_percent() const { return __top_percent; }
	double clus_rad() const { return __clus_rad; }
	int min_pts() const { return __min_pts; }
	int max_num_clus() const { return __max_num_clus; }
	int max_seeds_to_cluster() const { return __max_seeds_to_cluster; }
	int num_iter() const { return __num_iter; }
	int ncpu() const { return __ncpu; }
	
	friend ostream& operator<< (ostream& stream, const CmdLnOpts &cmdl) {
		unsigned int n = thread::hardware_concurrency();
		stream << endl << "Detected support for " << n << " concurrent threads." << endl;
		return stream;
	}

};
extern CmdLnOpts cmdl;
#endif
