/* This is dockfragments.cpp and is part of CANDOCK
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

#include "candock/program/dockfragments.hpp"

#include <boost/filesystem/path.hpp>

#include "statchem/helper/grep.hpp"
#include "statchem/fileio/inout.hpp"
#include "statchem/helper/path.hpp"
#include "statchem/helper/logger.hpp"
#include "statchem/molib/atomtype.hpp"
#include "candock/program/options.hpp"

namespace candock {
namespace Program {

DockFragments::DockFragments(const FindCentroids& found_centroids,
                             const FragmentLigands& fragmented_ligands,
                             const score::Score& score,
                             const molib::Atom::Grid& gridrec,
                             const std::string& name)
    : __found_centroids(found_centroids),
      __fragmented_ligands(fragmented_ligands),
      __score(score),
      __gridrec(gridrec),
      __name(name) {
    if (!cmdl.get_string_option("top_seeds_dir").empty()) {
        __top_seeds_location = cmdl.get_string_option("top_seeds_dir");
        return;
    }

    boost::filesystem::path p(__name);
    p = p / "top_seeds";
    __top_seeds_location = p.string();
}

bool DockFragments::__can_read_from_files() {
    const molib::Molecules& all_seeds = __fragmented_ligands.seeds();

    // No early return so that we have the ability to redock missing seeds
    for (const auto& seed : all_seeds) {
        boost::filesystem::path p(__top_seeds_location);
        p = p / seed.name();
        p += ".pdb";
        if (fileio::file_size(p.string()) <= 0) return false;
    }

    return true;
}

void DockFragments::__read_from_files() {
    log_step << "All seeds are present in "
             << cmdl.get_string_option("top_seeds_dir") << " for " << __name
             << ". Docking of fragments skipped." << std::endl;
}

void DockFragments::__dock_fragment(int start, const docker::Gpoints& gpoints,
                                    const docker::Gpoints& gpoints0) {
    // iterate over docked seeds and dock unique seeds
    for (size_t j = start; j < __fragmented_ligands.seeds().size();
         j += cmdl.ncpu()) {
        boost::filesystem::path p(__top_seeds_location);
        p = p / __fragmented_ligands.seeds()[j].name();
        p += (".pdb");

        try {
            if (fileio::file_size(p.string()) > 0) {
                log_note << "Skipping docking of seed: "
                         << __fragmented_ligands.seeds()[j].name()
                         << " because it is already docked!" << std::endl;
                continue;
            } else {
                log_note << "Docking seed: "
                         << __fragmented_ligands.seeds()[j].name() << std::endl;
            }
            dbgmsg(__fragmented_ligands.seeds()[j]);
            /* Compute all conformations of this seed with the center
             * atom fixed on coordinate origin using maximum clique algorithm
             *
             */
            docker::Conformations conf(
                __fragmented_ligands.seeds()[j], gpoints0,
                cmdl.get_double_option("conf_spin"),  // degrees
                cmdl.get_int_option("num_univec")     // number of unit vectors
                );

#ifndef NDEBUG
            fileio::output_file(
                conf,
                "conf_" + __fragmented_ligands.seeds()[j].name() + ".pdb");
#endif
/* Dock this seed's conformations to the entire grid by moving them
 * over all gridpoints and probe where they clash with the receptor:
 * cluster docked conformations based on rmsd and energy and output
 * only best-scored cluster representatives
 *
 */
#ifndef NDEBUG
            docker::Dock dock(gpoints, conf, __fragmented_ligands.seeds()[j],
                              __score, __gridrec,
                              cmdl.get_double_option("clus_rad"));
#else
            docker::Dock dock(gpoints, conf, __fragmented_ligands.seeds()[j],
                              cmdl.get_double_option("clus_rad"));
#endif

            dock.run();

            fileio::output_file(dock.get_docked(),
                               p.string());  // output docked & clustered seeds
        } catch (std::exception& e) {
            log_warning << "skipping seed due to : " << e.what() << std::endl;
            fileio::output_file(e.what(), p.string());
        }
    }
}

void DockFragments::__continue_from_prev() {
    log_step << "Docking fragments into: " << __top_seeds_location << std::endl;

    /* Create gridpoints for ALL centroids representing one or more binding
     * sites
     *
     */
    docker::Gpoints gpoints = get_gridhcp();

    /* Create a zero centered centroid with 10 A radius (max fragment
     * radius) for getting all conformations of each seed
     *
     */
    docker::Gpoints gpoints0(cmdl.get_double_option("grid"),
                             cmdl.get_double_option("max_frag_radius"));

    /* Create template grids using ProBiS-ligands algorithm
     * WORK IN PROGESS WORK IN PROGESS WORK IN PROGESS WORK IN PROGESS
     * WORK IN PROGESS WORK IN PROGESS WORK IN PROGESS WORK IN PROGESS
     * WORK IN PROGESS WORK IN PROGESS WORK IN PROGESS WORK IN PROGESS
     */
    std::vector<std::thread> threads;

    for (int i = 0; i < cmdl.ncpu(); ++i) {
        threads.push_back(std::thread(
            [&, this, i] { __dock_fragment(i, gpoints, gpoints0); }));
    }
    for (auto& thread : threads) {
        thread.join();
    }

    log_step << "Done with fragment docking" << std::endl;
}

std::vector<std::pair<double, std::string>> DockFragments::get_best_seeds()
    const {
    const molib::Molecules& all_seeds = __fragmented_ligands.seeds();

    std::vector<std::pair<double, std::string>> seed_score_map;

    for (const auto& seed : all_seeds) {
        // FIXME: Do not read from disk here
        dbgmsg("Reading: " << seed.name() << endl);
        try {
            boost::filesystem::path p(__top_seeds_location);
            p = p / seed.name();
            p += ".pdb";
            parser::FileParser spdb(p.string(), parser::first_model, 1);

            molib::Molecules seed_molec = spdb.parse_molecule();

            seed_molec.compute_hydrogen();
            std::tuple<double, size_t, size_t, size_t> frag_tuple =
                molib::AtomType::determine_lipinski(seed_molec.get_atoms());

            if (std::get<2>(frag_tuple) > 1 || std::get<3>(frag_tuple) > 1) {
                continue;
            }

            seed_score_map.push_back(
                {std::stod(seed_molec.first().name()), seed.name()});
        } catch (Error e) {
            log_warning << "Skipping seed " << seed.name() << " in " << __name
                        << " because " << e.what() << std::endl;
        }
    }

    std::sort(seed_score_map.begin(), seed_score_map.end());
    return seed_score_map;
}

molib::NRset DockFragments::get_negative_seeds(
    const std::set<std::string>& seeds, const double max_value) const {
    molib::NRset top_seeds;

    std::regex regex("REMARK   5 MOLECULE (\\S*)");

    for (auto& fragment : seeds) {
        boost::filesystem::path file_to_read(__top_seeds_location);
        file_to_read /= fragment;
        file_to_read += ".pdb";

        std::ifstream file(file_to_read.c_str());

        const size_t number_of_seeds =
            grep::find_first_case_greater_than(file, regex, max_value);

        parser::FileParser pdb(file_to_read.string(), parser::all_models,
                               number_of_seeds);

        dbgmsg("reading top_seeds_file for seed id = " << fragment);
        molib::Molecules all = pdb.parse_molecule();

        dbgmsg("number of top seeds left = " << all.size());

        molib::Molecules& last = top_seeds.add(new molib::Molecules(all));

        if (last.empty()) {
            throw Error("die : there are no docked conformations for seed " +
                        fragment);
        }
    }

    return top_seeds;
}

molib::NRset DockFragments::get_top_seeds(const std::set<std::string>& seeds,
                                          const double top_percent) const {
    molib::NRset top_seeds;

    std::regex regex;
    regex.assign("REMARK   5 MOLECULE ", std::regex_constants::basic);

    for (auto& fragment : seeds) {
        boost::filesystem::path file_to_read(__top_seeds_location);
        file_to_read /= fragment;
        file_to_read += ".pdb";

        std::ifstream file(file_to_read.c_str());

        const size_t number_of_seeds = grep::count_matches(file, regex);
        const int sz = static_cast<int>(number_of_seeds * top_percent);
        dbgmsg("taking " << sz << " top seeds for seed " << fragment);

        // Add one in case the user is silly enough to select a top_percent of
        // 0.000
        parser::FileParser pdb(file_to_read.string(), parser::all_models,
                               sz + 1);

        dbgmsg("reading top_seeds_file for seed id = " << fragment);
        molib::Molecules all = pdb.parse_molecule();

        dbgmsg("number of top seeds left = " << all.size());

        molib::Molecules& last = top_seeds.add(new molib::Molecules(all));

        if (last.empty()) {
            throw Error("die : there are no docked conformations for seed " +
                        fragment);
        }
    }

    return top_seeds;
}

molib::NRset DockFragments::get_seeds(const molib::Molecule& ligand,
                                      const double top_percent) const {
    std::set<std::string> seeds_to_read;
    const molib::Model& model = ligand.first().first();
    for (auto& fragment : model.get_rigid()) {  // iterate over seeds
        if (fragment.is_seed()) {
            seeds_to_read.insert(std::to_string(fragment.get_seed_id()));
        }
    }
    return top_percent <= 0 ? get_negative_seeds(seeds_to_read, top_percent)
                            : get_top_seeds(seeds_to_read, top_percent);
}

docker::Gpoints DockFragments::get_gridhcp() {
    docker::Gpoints gpoints(__score, __fragmented_ligands.ligand_idatm_types(),
                            __found_centroids.centroids(), __gridrec,
                            cmdl.get_double_option("grid"),
                            cmdl.get_int_option("cutoff"),
                            cmdl.get_double_option("excluded"),
                            cmdl.get_double_option("interatomic"));
    return gpoints;
}
}
}
