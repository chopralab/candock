#include "findcentroids.hpp"

#include <boost/filesystem.hpp>

#include "probis/probis.hpp"
#include "helper/inout.hpp"
#include "ligands/genclus.hpp"
#include "molib/molecules.hpp"

#include "helper/path.hpp"
#include "options.hpp"

#ifdef _WINDOWS
#include <direct.h>
#endif

namespace Program {

        FindCentroids::FindCentroids(const std::string &filename, const std::string &chain_ids, const std::string& out_dir) :
                     __filename(filename), __chain_ids(chain_ids), __out_dir(out_dir) {
                if (cmdl.get_string_option("centroid").empty()) {
                        __centroid_file = Path::join(out_dir, "site.cen");
                } else {
                        __centroid_file = cmdl.get_string_option("centroid");
                }
        }

        bool FindCentroids::__can_read_from_files( ) {
                return Inout::file_size( __centroid_file ) > 0;
        }

        void FindCentroids::__read_from_files( ) {
                log_note << "Reading " << cmdl.get_int_option("num_bsites") << 
                    " binding sites from " << __centroid_file << endl;
                __result = Centro::set_centroids( __centroid_file, cmdl.get_int_option("num_bsites"));
        }

        void FindCentroids::__continue_from_prev( ) {

                log_step << "Running PROBIS for receptor in file: " << __filename + ".pdb" << endl;

                // Creates an empty nosql file for probis local structural alignments
                Inout::output_file("", Path::join(__out_dir, cmdl.get_string_option("nosql")));

                // PROBIS is a bit needy and requires the directory 'bslibdb' to be in the current path
                // To make this work properly, we change directories to the directory with this directory
                // and change back. This is effectively a dirty hack around a probis problem...
                boost::filesystem::path original_file(__filename);

                boost::filesystem::path cwd = boost::filesystem::current_path();
                original_file = absolute(original_file);

                boost::filesystem::path bslibdb (cmdl.get_string_option("bslib"));
                bslibdb = bslibdb.parent_path();

#ifdef _WINDOWS
                int chdir_error = _wchdir( bslibdb.c_str() );
#else
                int chdir_error = chdir( bslibdb.c_str() );
#endif
                if ( chdir_error != 0) {
                        throw Error ("Unable to change into bslibdb dir: " + bslibdb.string() + " because " + strerror(chdir_error)) ;
                }

                probis::compare_against_bslib(original_file.string(),
                        (cwd / __out_dir / cmdl.get_string_option("srf_file")).string(),
                        __chain_ids,
                        "bslibdb/bslib.txt", cmdl.ncpu(),
                        (cwd / __out_dir / cmdl.get_string_option("nosql")).string(),
                        (cwd / __out_dir / cmdl.get_string_option("json")).string() );

#ifdef _WINDOWS
                chdir_error = _wchdir( cwd.c_str() );
#else
                chdir_error = chdir( cwd.c_str() );
#endif
                if ( chdir_error != 0) {
                        throw Error ("Unable to change into original dir: " + bslibdb.string() + " because " + strerror(chdir_error)) ;
                }

                genclus::generate_clusters_of_ligands(
                        Path::join(__out_dir, cmdl.get_string_option("json")),
                        Path::join(__out_dir, cmdl.get_string_option("jsonwl")),
                        cmdl.get_string_option("bio"),
                        cmdl.get_string_option("names"),
                        cmdl.get_bool_option("neighb"),
                        cmdl.get_double_option("probis_clus_rad"),
                        cmdl.get_int_option("probis_min_pts"),
                        cmdl.get_double_option("probis_min_z_score"));

                auto binding_sites = 
                        genlig::generate_binding_site_prediction(
                                Path::join(__out_dir, cmdl.get_string_option("jsonwl")), 
                                cmdl.get_string_option("bio"),
                                cmdl.get_int_option("num_bsites"));

                Inout::output_file(binding_sites.first,  Path::join(__out_dir, cmdl.get_string_option("lig_clus_file")));
                Inout::output_file(binding_sites.second, Path::join(__out_dir, cmdl.get_string_option("z_scores_file")));

                __result = Centro::set_centroids(binding_sites.first, cmdl.get_double_option("centro_clus_rad"));
                Inout::output_file(__result, __centroid_file); // probis local structural alignments
        }
}
