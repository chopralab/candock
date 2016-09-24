#ifndef OPTS_CANDOCK_H
#define OPTS_CANDOCK_H

#include <boost/program_options.hpp>
#include <boost/date_time/posix_time/posix_time.hpp>
#include <boost/asio/ip/host_name.hpp>
#include "helper/error.hpp"
#include "version.h"
#include <thread>

namespace Program {

	class CmdLnOpts {
		boost::program_options::variables_map __vm;

		std::string __receptor_file;
		std::string __ligand_file;
		std::string __seeds_file;
		std::string __prep_file;

		std::string __bslib_file;
		std::string __nosql_file;
		std::string __json_file;
		std::string __json_with_ligs_file;

		std::string __top_seeds_dir;

		std::string __names_dir;
		bool __neighb;
		double __probis_clus_rad;
		int __probis_min_pts;
		double __probis_min_z_score;


		std::string __bio_dir;
		std::string __lig_clus_file;
		std::string __z_scores_file;
		double __centro_clus_rad;

		std::string __centroid_file;

		std::string __gridpdb_hcp_file;
		double __max_frag_radius;

		std::string __ref_state;
		std::string __comp;
		std::string __rad_or_raw;
		int __dist_cutoff;
		std::string __distributions_file;
		double __step_non_bond;
		double __scale_non_bond;

		std::string __potential_file;
		std::string __obj_dir;

		std::string __cluster_file;
		std::string __top_seeds_file;

		int __max_num_ligands;

		std::string __gaff_dat_file;
		std::string __gaff_xml_file;

		int __spin_degrees;
		double __clash_coeff;
		double __tol_seed_dist;
		double __lower_tol_seed_dist;
		double __upper_tol_seed_dist;
		int __max_possible_conf;
		int __link_iter;

		std::string __docked_dir;

		double __docked_clus_rad;
		double __max_allow_energy;
		bool __iterative;
		int __max_num_possibles;

		std::string __amber_xml_file;
		std::string __water_xml_file;
		std::string __fftype;
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

		std::string __program_name;
		std::string __version;
		std::string __git_version;
	public:
		CmdLnOpts() : __quiet (false), __version (""), __git_version ("") {
			std::stringstream ss;

			ss << CANDOCK_MAJOR_VERSION << "." << CANDOCK_MINOR_VERSION << "."
			   << CANDOCK_TWEAK_VERSION;

			__version = ss.str();

			std::stringstream ss2;

			ss2 << CANDOCK_GIT_REFERENCE << " on branch " << CANDOCK_GIT_MYCBRANCH;

			__git_version = ss2.str();
		}

		enum CMDLN_OPTS_GROUPS {
			STARTING     = 1 << 0,
			PROBIS       = 1 << 1,
			LIG_FRAMGENT = 1 << 2,
			FRAG_DOCKING = 1 << 3,
			SCORING      = 1 << 4,
			FORCE_FIELD  = 1 << 5,
			LINKING      = 1 << 6,
			ALL_OPTIONS  = 0xFF
		};

		void init (int argc, char *argv[], int opts_to_parse = ALL_OPTIONS);
		void display_time (std::string what) {
			cout << what << " on " << boost::posix_time::to_simple_string (boost::posix_time::second_clock::local_time()) << "\n";
			cout << "running " << __program_name << " on hostname " << boost::asio::ip::host_name() << "\n";
			cout << "version " << __version << "\n";
			cout << "buildid " << __git_version << "\n";
		}
		// interface

		std::string get_string_option(std::string option) const;
		bool        get_bool_option(std::string option) const;
		int         get_int_option(std::string option) const;
		double      get_double_option (std::string option) const;

		bool quiet() const {
			return __quiet;
		}

		std::string program_name() const {
			return __program_name;
		}

		std::string receptor_file() const {
			return __receptor_file;
		}
		std::string ligand_file() const {
			return __ligand_file;
		}
		std::string seeds_file() const {
			return __seeds_file;
		}
		std::string prep_file() const {
			return __prep_file;
		}
		std::string bslib_file() const {
			return __bslib_file;
		}
		std::string nosql_file() const {
			return __nosql_file;
		}
		std::string json_file() const {
			return __json_file;
		}
		std::string json_with_ligs_file() const {
			return __json_with_ligs_file;
		}
		std::string top_seeds_dir() const {
			return __top_seeds_dir;
		}
		std::string names_dir() const {
			return __names_dir;
		}
		bool neighb() const {
			return __neighb;
		}
		double probis_clus_rad() const {
			return __probis_clus_rad;
		}
		int probis_min_pts() const {
			return __probis_min_pts;
		}
		double probis_min_z_score() const {
			return __probis_min_z_score;
		}
		std::string bio_dir() const {
			return __bio_dir;
		}
		std::string lig_clus_file() const {
			return __lig_clus_file;
		}
		std::string z_scores_file() const {
			return __z_scores_file;
		}
		double centro_clus_rad() const {
			return __centro_clus_rad;
		}
		std::string centroid_file() const {
			return __centroid_file;
		}
		std::string gridpdb_hcp_file() const {
			return __gridpdb_hcp_file;
		}
		double max_frag_radius() const {
			return __max_frag_radius;
		}

		std::string ref_state() const {
			return __ref_state;
		}
		std::string comp() const {
			return __comp;
		}
		std::string rad_or_raw() const {
			return __rad_or_raw;
		}
		int dist_cutoff() const {
			return __dist_cutoff;
		}
		std::string distributions_file() const {
			return __distributions_file;
		}
		double step_non_bond() const {
			return __step_non_bond;
		}
		double scale_non_bond() const {
			return __scale_non_bond;
		}

		std::string potential_file() const {
			return __potential_file;
		}
		std::string obj_dir() const {
			return __obj_dir;
		}

		std::string cluster_file() const {
			return __cluster_file;
		}
		std::string top_seeds_file() const {
			return __top_seeds_file;
		}

		int max_num_ligands() const {
			return __max_num_ligands;
		}
		std::string gaff_dat_file() const {
			return __gaff_dat_file;
		}
		std::string gaff_xml_file() const {
			return __gaff_xml_file;
		}
		int spin_degrees() const {
			return __spin_degrees;
		}

		double clash_coeff() const {
			return __clash_coeff;
		}
		double tol_seed_dist() const {
			return __tol_seed_dist;
		}
		double lower_tol_seed_dist() const {
			return __lower_tol_seed_dist;
		}
		double upper_tol_seed_dist() const {
			return __upper_tol_seed_dist;
		}
		int max_possible_conf() const {
			return __max_possible_conf;
		}
		int link_iter() const {
			return __link_iter;
		}
		std::string docked_dir() const {
			return __docked_dir;
		}

		double docked_clus_rad() const {
			return __docked_clus_rad;
		}
		double max_allow_energy() const {
			return __max_allow_energy;
		}
		bool iterative() const {
			return __iterative;
		}
		int max_num_possibles() const {
			return __max_num_possibles;
		}

		std::string amber_xml_file() const {
			return __amber_xml_file;
		}
		std::string water_xml_file() const {
			return __water_xml_file;
		}
		std::string fftype() const {
			return __fftype;
		}
		double tolerance() const {
			return __tolerance;
		}
		int max_iterations() const {
			return __max_iterations;
		}
		int max_iterations_final() const {
			return __max_iterations_final;
		}
		int update_freq() const {
			return __update_freq;
		}
		double position_tolerance() const {
			return __position_tolerance;
		}

		double top_percent() const {
			return __top_percent;
		}
		int max_clique_size() const {
			return __max_clique_size;
		}

		int num_bsites() const {
			return __num_bsites;
		}
		double grid_spacing() const {
			return __grid_spacing;
		}
		int num_univec() const {
			return __num_univec;
		}
		double conf_spin() const {
			return __conf_spin;
		}
		double excluded_radius() const {
			return __excluded_radius;
		}
		double max_interatomic_distance() const {
			return __max_interatomic_distance;
		}

		double clus_rad() const {
			return __clus_rad;
		}

		int ncpu() const {
			return __ncpu;
		}

		friend std::ostream &operator<< (std::ostream &stream, const CmdLnOpts &cmdl) {
			unsigned int n = std::thread::hardware_concurrency();
			stream << endl << "Detected support for " << n << " concurrent threads."
			       << " Using " << cmdl.ncpu() << " threads." << endl;

			for ( const auto& a : cmdl.__vm ) {
				stream << std::setw(24)<< a.first;

				if ( auto v = boost::any_cast<std::string>(&a.second.value()) )
					stream << std::setw(50) << *v;
				else if ( auto v = boost::any_cast<int>(&a.second.value()) )
					stream << std::setw(50) << *v;
				else if ( auto v = boost::any_cast<double>(&a.second.value()) )
					stream << std::setw(50) << *v;
				else if ( auto v = boost::any_cast<bool>(&a.second.value()) )
					stream << std::setw(50) << *v;

				stream << std::endl;
			}

			return stream;
		}

	};
}

#endif
