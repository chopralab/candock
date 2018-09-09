#include "candock/ligands/genclus.hpp"
#include <stdlib.h> /* srand, rand */
#include <time.h>   /* time */
#include <algorithm>
#include <exception>
#include <iostream>
#include <typeinfo>
#include "candock/cluster/optics.hpp"
#include "statchem/geometry/matrix.hpp"
#include "statchem/helper/debug.hpp"
#include "statchem/helper/error.hpp"
#include "statchem/helper/help.hpp"
#include "statchem/fileio/inout.hpp"
#include "statchem/helper/path.hpp"
#include "candock/ligands/common.hpp"
#include "candock/ligands/jsonreader.hpp"
#include "candock/ligands/nosqlreader.hpp"
#include "statchem/molib/grid.hpp"
#include "statchem/molib/nrset.hpp"
#include "statchem/parser/fileparser.hpp"
using namespace std;

namespace candock {
namespace genclus {

using namespace statchem;

static geometry::Matrix make_matrix(const Json::Value& rota,
                                    const Json::Value& trans) {
    geometry::Matrix result;

    for (int row = 0; row < 3; ++row) {
        auto row0 = std::make_tuple(rota[row][0].asDouble(),
                                    rota[row][1].asDouble(),
                                    rota[row][2].asDouble(),
                                    trans[row].asDouble());

        result.set_row(row, row0);
    }

    return result;
}

string remove_special_characters(string str) {
    str.erase(std::remove_if(str.begin(), str.end(),
                             [](const char& c) {
                                 return c == ':' || c == '"' || c == '\\';
                             }),
              str.end());
    return str;
}
void add_to_json(JsonReader& jr,
                 cluster::Clusters<molib::Molecule> ligand_clusters,
                 string cluster_name, const string& names_dir,
                 const bool for_gclus) {
    map<string, unique_ptr<Json::Value>> ligands;
    for (auto& kv : ligand_clusters) {
        int cluster_id =
            (kv.first == -1 ? 0 : kv.first);  // unclustered is cluster #0
        const molib::Molecule& molecule = *kv.second;
        const molib::Molecules& mols =
            molecule.br();  // mols.name() is the nr-pdb name
        const molib::Assembly& assembly = molecule.first();
        const molib::Model& model = assembly.first();
        const molib::Chain& chain = model.first();
        const molib::Residue& residue = chain.first();
        dbgmsg(mols.name());
        if (ligands.find(mols.name()) == ligands.end()) {
            ligands.insert(make_pair(
                mols.name(),
                unique_ptr<Json::Value>(new Json::Value(Json::arrayValue))));
        }
        vector<string> lig_name;
        try {
            fileio::read_file(
                Path::join(names_dir,
                           ((residue.rest() == molib::Residue::protein ||
                             residue.rest() == molib::Residue::nucleic)
                                ? (molecule.name().substr(0, 4) +
                                   molecule.name().substr(4, 1))
                                : residue.resn())),
                lig_name);
        } catch (exception& e) {
            dbgmsg(e.what() << " ... skipping ... ");
        }

        ligands[mols.name()]->append(
            to_string(cluster_id) + ":" + residue.resn() + ":" +
            to_string(residue.resi()) + ":" + chain.chain_id() + ":" +
            to_string(model.number()) + ":" + to_string(assembly.number()) +
            ":" + molecule.name() + ":" + mols.name() + ":" +
            (lig_name.empty() ? "-" : remove_special_characters(lig_name[0])) +
            ":" + (help::non_specific_binders.find(residue.resn()) ==
                           help::non_specific_binders.end()
                       ? "S"
                       : "N"));
        dbgmsg(residue.resn() << ":" << to_string(cluster_id)
                              << " molecule_name = " << molecule.name()
                              << " mols_name = " << mols.name());
    }
    for (auto& i : ligands) {
        const string& mols_name = i.first;
        dbgmsg("TEST mols_name = " << mols_name);
        Json::Value& vec = *(i.second);
        JsonReader::iterator it = jr.end();
        if (for_gclus) {
            if (mols_name.size() != 5)
                throw Error("die : format should be PdbId[no-space]ChainID");
            const string pdb_id = mols_name.substr(0, 4);
            const string chain_id = mols_name.substr(4, 1);
            it = jr.find(
                {make_pair("pdb_id", pdb_id), make_pair("chain_id", chain_id)});
            dbgmsg("pdb_id = " << pdb_id << " chain_id = " << chain_id);
        } else {
            it = jr.find({make_pair("pdb_id", mols_name)});
        }
        if (it != jr.end()) {
            dbgmsg("found json value " << mols_name);
            (*it)[cluster_name] = vec;
            dbgmsg("after found json value " << mols_name);
        }
    }
}
void squeeze_proteins_nucleic(molib::Molecules& mols) {
    for (auto& molecule : mols) {
        for (auto& assembly : molecule) {
            for (auto& model : assembly) {
                for (auto& chain : model) {
                    geometry::Coordinate::Vec crdp, crdn;
                    // remove protein and nucleic ligands
                    for (auto& residue : chain)
                        residue.set_crd();  // they really need to be set!!!
                    chain.remove_if([&chain, &model, &crdp,
                                     &crdn](const molib::Residue& r) {
                        if (r.rest() == molib::Residue::protein) {
                            crdp.push_back(r.crd());
                            return true;
                        }
                        if (r.rest() == molib::Residue::nucleic) {
                            crdn.push_back(r.crd());
                        }
                        return false;
                    });
                    // add PPP or NNN instead
                    if (!crdp.empty()) {
                        molib::Residue& residue = chain.add(new molib::Residue(
                            "PPP", 1, ' ', molib::Residue::protein));
                        residue.add(new molib::Atom(
                            1, "CA", geometry::compute_geometric_center(crdp),
                            help::idatm_mask.at("???")));
                    }
                    if (!crdn.empty()) {
                        molib::Residue& residue = chain.add(new molib::Residue(
                            "NNN", 1, ' ', molib::Residue::nucleic));
                        residue.add(new molib::Atom(
                            1, "P", geometry::compute_geometric_center(crdn),
                            help::idatm_mask.at("???")));
                    }
                }
            }
        }
    }
}

void remove_ligands_not_in_aligned_region(molib::Molecules& mols,
                                          Json::Value aligned_residues) {
    set<molib::Residue::res_tuple2> a =
        common_ligands::json_to_set(aligned_residues).second;
    for (auto& molecule : mols) {
        for (auto& assembly : molecule) {
            for (auto& model : assembly) {
                for (auto& chain : model) {
                    // remove hetero, ion and water residues that are not close
                    // to the query
                    chain.remove_if([&a, &chain,
                                     &model](const molib::Residue& r) {
                        map<molib::Residue::res_tuple2,
                            molib::Residue::res_tuple2>
                            respair;
                        model.remarks(7, molib::Residue::res_tuple2(
                                             chain.chain_id(), r.resn(),
                                             r.resi(), r.ins_code()),
                                      respair);
                        set<molib::Residue::res_tuple2> b;
                        for (auto& i : respair) {
                            b.insert(i.first);
                        }
                        dbgmsg("LIGAND " << chain.chain_id() << ":" << r.resn()
                                         << ":" << r.resi() << ":"
                                         << r.ins_code());
#ifndef NDEBUG
                        for (auto& i : b) {
                            dbgmsg("B " << get<0>(i) << ":" << get<1>(i) << ":"
                                        << get<2>(i) << ":" << get<3>(i));
                        }
#endif
                        vector<molib::Residue::res_tuple2> v(b.size());
                        vector<molib::Residue::res_tuple2>::iterator it =
                            set_intersection(a.begin(), a.end(), b.begin(),
                                             b.end(), v.begin());
                        v.resize(it - v.begin());
#ifndef NDEBUG
                        for (auto& iv : v) {
                            dbgmsg("INTERSECTION "
                                   << get<0>(iv) << ":" << get<1>(iv) << ":"
                                   << get<2>(iv) << ":" << get<3>(iv));
                        }
#endif
                        // delete if less than 2 residues are in intersection
                        const bool is_prot_res =
                            r.rest() == molib::Residue::protein && v.size() < 5;
                        const bool is_nuca_res =
                            r.rest() == molib::Residue::nucleic && v.size() < 4;
                        const bool is_hetr_res =
                            r.rest() == molib::Residue::hetero && v.size() < 3;
                        const bool is_ionc_res =
                            r.rest() == molib::Residue::ion && v.size() < 3;
                        if (is_prot_res || is_nuca_res || is_hetr_res ||
                            is_ionc_res) {
                            return true;
                        }
                        return false;  // don't delete
                    });
                }
            }
        }
    }
}

/**
 * Calculate shortest distance between any two molecules.
 *
 */
cluster::PairwiseDistances<molib::Molecule> create_pairwise_distances(
    const molib::Molecule::Vec& ligands, const double eps) {
    cluster::PairwiseDistances<molib::Molecule> pairwise_distances;
    molib::Atom::Vec atoms;

    for (auto& pmolecule : ligands)
        for (auto& patom : pmolecule->get_atoms()) {
            atoms.push_back(patom);
            dbgmsg("&molecule = " << &const_cast<molib::Molecule&>(
                                         patom->br().br().br().br().br())
                                  << " molecule = "
                                  << const_cast<molib::Molecule&>(
                                         patom->br().br().br().br().br())
                                         .name()
                                  << " (pmolecule = " << pmolecule
                                  << " should be " << pmolecule->name()
                                  << ") atom = " << *patom);
        }

    molib::Atom::Grid grid(atoms);

    for (auto& patom : atoms) {
        molib::Molecule& molecule1 =
            const_cast<molib::Molecule&>(patom->br().br().br().br().br());
        for (auto& pneighbor :
             grid.get_neighbors_including_self(patom->crd(), eps)) {
            molib::Molecule& molecule2 = const_cast<molib::Molecule&>(
                pneighbor->br().br().br().br().br());
            const double distance = patom->crd().distance(pneighbor->crd());
            if (pairwise_distances[&molecule1][&molecule2] < 0.0000000001 ||
                distance < pairwise_distances[&molecule1][&molecule2]) {
                pairwise_distances[&molecule1][&molecule2] = distance;
            }
        }
    }
    return pairwise_distances;
}

/**
 * Split read molecules into smaller pieces (i.e., protein chain A and its
 * ligands
 * also chain A are outputted as separate molecules)
 */
molib::NRset split_into_molecules(molib::NRset& nrset) {
    molib::NRset cnrset;
    for (auto& mols : nrset) {
        molib::Molecules& cmols = cnrset.add(new molib::Molecules(mols.name()));
        for (auto& molecule : mols) {
            for (auto& assembly : molecule) {
                for (auto& model : assembly) {
                    for (auto& chain : model) {
                        molib::Residue::res_type prev_rest(
                            molib::Residue::notassigned);
                        for (auto& residue : chain) {
                            if (prev_rest != residue.rest() ||
                                residue.rest() == molib::Residue::water ||
                                residue.rest() == molib::Residue::ion ||
                                residue.rest() == molib::Residue::hetero) {
                                molib::Molecule& cmolecule = cmols.add(
                                    new molib::Molecule(molecule.name()));
                                molib::Assembly& cassembly = cmolecule.add(
                                    new molib::Assembly(assembly.number()));
                                molib::Model& cmodel = cassembly.add(
                                    new molib::Model(model.number()));
                                cmodel.add(new molib::Chain(chain.chain_id()));
                            }
                            cmols.last().first().first().first().add(
                                new molib::Residue(residue));
                            prev_rest = residue.rest();
                        }
                        chain.clear();
                    }
                }
            }
        }
    }
    return cnrset;
}

void generate_clusters_of_ligands(
    const string& json_file, const string& json_with_ligs_file,
    const string& bio_dir, const string& names_dir, const bool neighb,
    const double hetero_clus_rad, const int hetero_min_pts,
    const double min_z_score, const bool for_gclus) {
    JsonReader jr;
    molib::NRset nrset;
    map<string, double> scr;
    // the aligned files from json
    jr.parse_JSON(json_file);
    for (auto& d : jr.root()) {
        /* if something goes wrong, e.g., a pdb file is not found,
         * don't exit, just skip that file..
         */
        try {
            const string pdb_id = d["pdb_id"].asString();
            const string chain_ids = d["chain_id"].asString();
            const double z_score =
                d["alignment"][0]["scores"]["z_score"].asDouble();
            if (z_score < min_z_score) continue;

            const string pdb_file = Path::join(
                bio_dir, (for_gclus ? (pdb_id + chain_ids) : pdb_id) + ".pdb");

            dbgmsg(pdb_id + " " + chain_ids);

            parser::FileParser pr(pdb_file,
                                  parser::all_models | parser::sparse_macromol);
            molib::Molecules& mols =
                nrset.add(new molib::Molecules(pr.parse_molecule()));

            squeeze_proteins_nucleic(mols);

            for (auto& molecule : mols) {
                scr[molecule.name()] = z_score;
            }

            dbgmsg("from genclus : " << mols.name());
            if (neighb)
                remove_ligands_not_in_aligned_region(
                    mols,
                    d["alignment"][0]["aligned_residues"]);  // remove ligands
                                                             // that are not
                                                             // near aligned
                                                             // residues
            mols.rotate(
                make_matrix(d["alignment"][0]["rotation_matrix"],
                            d["alignment"][0]["translation_vector"]),
                true);  // inverse rotation
        } catch (exception& e) {
            nrset.erase(nrset.size() -
                        1);  // delete the last mols which is empty
            log_warning << e.what() << " ... skipping ... " << endl;
        }
    }
#ifndef NDEBUG
    for (auto& mols : nrset) {
        dbgmsg("molecules name = " << mols.name());
    }
#endif
    molib::NRset splitted = split_into_molecules(nrset);

    cluster::MapD<molib::Molecule> scores;
    for (auto& mols : splitted)
        for (auto& molecule : mols) scores[&molecule] = scr.at(molecule.name());

    molib::Molecule::Vec protein_ligands =
        splitted.get_molecules(molib::Residue::protein);
    molib::Molecule::Vec nucleic_ligands =
        splitted.get_molecules(molib::Residue::nucleic);
    molib::Molecule::Vec hetero_ligands =
        splitted.get_molecules(molib::Residue::hetero);
    molib::Molecule::Vec ion_ligands =
        splitted.get_molecules(molib::Residue::ion);

#ifndef NDEBUG
    for (auto& pmol : protein_ligands) {
        dbgmsg("protein_ligand = " << pmol->name());
        for (auto& patom : pmol->get_atoms()) {
            dbgmsg(*patom);
        }
    }
#endif

    cluster::PairwiseDistances<molib::Molecule> pd_protein =
        create_pairwise_distances(protein_ligands, 4.1);
    cluster::Optics<molib::Molecule> op(pd_protein, scores, 4.1, 3);

    cluster::PairwiseDistances<molib::Molecule> pd_nucleic =
        create_pairwise_distances(nucleic_ligands, 4.1);
    cluster::Optics<molib::Molecule> on(pd_nucleic, scores, 4.1, 3);

    //~ cluster::PairwiseDistances<molib::Molecule> pd_hetero =
    //create_pairwise_distances(hetero_ligands, 3.1);
    cluster::PairwiseDistances<molib::Molecule> pd_hetero =
        create_pairwise_distances(hetero_ligands, hetero_clus_rad + 0.1);
    cluster::Optics<molib::Molecule> oh(pd_hetero, scores,
                                        hetero_clus_rad + 0.1, hetero_min_pts);

    cluster::PairwiseDistances<molib::Molecule> pd_ion =
        create_pairwise_distances(ion_ligands, 3.1);
    cluster::Optics<molib::Molecule> oi(pd_ion, scores, 3.1, 4);

    auto protein_clusters = op.extract_dbscan(4.0);
    auto nucleic_clusters = on.extract_dbscan(4.0);
    auto hetero_clusters = oh.extract_dbscan(hetero_clus_rad, 100, true);
    auto ion_clusters = oi.extract_dbscan(3.0);

    add_to_json(jr, protein_clusters.first, "protein", names_dir, for_gclus);
    add_to_json(jr, nucleic_clusters.first, "nucleic", names_dir, for_gclus);
    add_to_json(jr, hetero_clusters.first, "hetero", names_dir, for_gclus);
    add_to_json(jr, ion_clusters.first, "ion", names_dir, for_gclus);
    fileio::file_open_put_stream(json_with_ligs_file,
                                stringstream(jr.output_json()));
}
}
}
