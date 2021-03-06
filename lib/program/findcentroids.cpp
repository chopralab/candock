/* This is findcentroids.cpp and is part of CANDOCK
 * Copyright (c) 2016-2019 Chopra Lab at Purdue University, 2013-2016 Janez Konc at National Institute of Chemistry and Samudrala Group at University of Washington
 *
 * This program is free for educational and academic use
 *
 * This program is free software; you can redistribute it and/or modify
 * it under the terms of the GNU General Public License as published by
 * the Free Software Foundation version 3 of the License.
 *
 * This program is distributed in the hope that it will be useful,
 * but WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the
 * GNU General Public License for more details.
 */

#include "candock/program/findcentroids.hpp"
#include "candock/program/options.hpp"

#include <boost/filesystem.hpp>

std::ostream& operator<<(std::ostream& os, const std::map<int, double>& bscores) {
    for (auto& kv : bscores) {
        const int& cluster_number = kv.first;
        const double& z_score = kv.second;
        os << cluster_number << " " << z_score << std::endl;
    }
    return os;
}

#include "candock/ligands/genclus.hpp"
#include "candock/ligands/genlig.hpp"
#include "candock/probis/probis.hpp"
#include "statchem/fileio/inout.hpp"
#include "statchem/molib/molecules.hpp"

#include "statchem/helper/logger.hpp"
#include "statchem/helper/path.hpp"

#ifdef _WINDOWS
#include <direct.h>
#endif

using namespace std;

namespace statchem {
namespace molib {

ostream& operator<<(ostream& os, const map<int, statchem::molib::Molecules>& bsites) {
    for (auto& kv : bsites) {
        const int& cluster_number = kv.first;
        const statchem::molib::Molecules& ligands = kv.second;
        os << "REMARK  99 ______________________ BEGINNING CLUSTER #"
           << cluster_number << " ______________________" << endl;
        os << ligands;
    }
    return os;
}
}
}

namespace candock {

namespace Program {

using namespace statchem;

FindCentroids::FindCentroids(const std::string& filename,
                             const std::string& chain_ids,
                             const std::string& out_dir)
    : __filename(filename), __chain_ids(chain_ids), __out_dir(out_dir) {
    if (cmdl.get_string_option("centroid").empty()) {
        __centroid_file = Path::join(out_dir, "site.cen");
    } else {
        __centroid_file = cmdl.get_string_option("centroid");
    }
}

bool FindCentroids::__can_read_from_files() {
    return fileio::file_size(__centroid_file) > 0;
}

void FindCentroids::__read_from_files() {
    log_note << "Reading " << cmdl.get_int_option("num_bsites")
             << " binding sites from " << __centroid_file << std::endl;
    __result = centro::set_centroids(__centroid_file,
                                     cmdl.get_int_option("num_bsites"));
}

void FindCentroids::__continue_from_prev() {
    log_step << "Running PROBIS for receptor in file: " << __filename
             << std::endl;

    // Creates an empty nosql file for probis local structural alignments
    fileio::output_file("",
                        Path::join(__out_dir, cmdl.get_string_option("nosql")));

    // PROBIS is a bit needy and requires the directory 'bslibdb' to be in the
    // current path
    // To make this work properly, we change directories to the directory with
    // this directory
    // and change back. This is effectively a dirty hack around a probis
    // problem...
    boost::filesystem::path original_file(__filename);

    boost::filesystem::path cwd = boost::filesystem::current_path();
    original_file = absolute(original_file);

    boost::filesystem::path bslibdb(cmdl.get_string_option("bslib"));
    bslibdb = bslibdb.parent_path();

    boost::filesystem::path protein_dir(__out_dir);

#ifdef _WINDOWS
    int chdir_error = _wchdir(bslibdb.c_str());
#else
    int chdir_error = chdir(bslibdb.c_str());
#endif
    if (chdir_error != 0) {
        throw Error("Unable to change into bslibdb dir: " + bslibdb.string() +
                    " because " + strerror(chdir_error));
    }

    boost::filesystem::path p = cwd / protein_dir;

    probis::compare_against_bslib(
        original_file.string(),
        (p / cmdl.get_string_option("srf_file")).string(), __chain_ids,
        "bslibdb/bslib.txt", cmdl.ncpu(),
        (p / cmdl.get_string_option("nosql")).string(),
        (p / cmdl.get_string_option("json")).string());

#ifdef _WINDOWS
    chdir_error = _wchdir(cwd.c_str());
#else
    chdir_error = chdir(cwd.c_str());
#endif
    if (chdir_error != 0) {
        throw Error("Unable to change into original dir: " + bslibdb.string() +
                    " because " + strerror(chdir_error));
    }

    genclus::generate_clusters_of_ligands(
        Path::join(protein_dir.string(), cmdl.get_string_option("json")),
        Path::join(protein_dir.string(), cmdl.get_string_option("jsonwl")),
        cmdl.get_string_option("bio"), cmdl.get_string_option("names"),
        cmdl.get_bool_option("neighb"),
        cmdl.get_double_option("probis_clus_rad"),
        cmdl.get_int_option("probis_min_pts"),
        cmdl.get_double_option("probis_min_z_score"));

    auto binding_sites = genlig::generate_binding_site_prediction(
        Path::join(protein_dir.string(), cmdl.get_string_option("jsonwl")),
        cmdl.get_string_option("bio"), cmdl.get_int_option("num_bsites"));

    fileio::output_file(binding_sites.first,
                        Path::join(protein_dir.string(),
                                   cmdl.get_string_option("lig_clus_file")));
    fileio::output_file(binding_sites.second,
                        Path::join(protein_dir.string(),
                                   cmdl.get_string_option("z_scores_file")));

    __result = centro::set_centroids(binding_sites.first,
                                     cmdl.get_double_option("centro_clus_rad"));
    fileio::output_file(__result,
                        __centroid_file);  // probis local structural alignments
}
}  // namespace Program
}  // namespace candock
