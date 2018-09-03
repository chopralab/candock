#include "candock/ligands/genbio.hpp"
#include <exception>
#include <iostream>
#include <typeinfo>
#include "statchem/geometry/matrix.hpp"
#include "statchem/helper/debug.hpp"
#include "statchem/helper/error.hpp"
#include "statchem/helper/help.hpp"
#include "statchem/fileio/inout.hpp"
#include "statchem/helper/path.hpp"
#include "statchem/helper/logger.hpp"
#include "candock/ligands/common.hpp"
#include "candock/ligands/jsonreader.hpp"
#include "candock/ligands/nosqlreader.hpp"
#include "statchem/molib/grid.hpp"
#include "statchem/molib/molecules.hpp"
#include "statchem/parser/fileparser.hpp"
using namespace std;

namespace candock {
namespace genbio {

using namespace statchem;

void remove_chains(molib::Molecules& mols,
                   const vector<string>& aligned_chains) {
    int i = 0;
    for (auto& molecule : mols) {
        for (auto& assembly : molecule) {
            for (char chain_id : aligned_chains[i]) {
                if (assembly.first().has_chain(chain_id)) {
                    assembly.first().chain(chain_id).remove_if(
                        [](const molib::Residue& r) {
                            return r.rest() == molib::Residue::protein;
                        });  // remove all protein residues
                    if (assembly.first().chain(chain_id).empty())
                        assembly.first().erase(chain_id);
                }
            }
        }
        i++;
    }
}
void remove_assemblies(
    molib::Molecule& molecule,
    const string
        chain_ids) {  // remove assemblies that don't have aligned chain_ids
    for (size_t i = 0; i < molecule.size(); i++) {
        molib::Assembly& assembly = molecule[i];
        size_t sz = 0;
        for (char chain_id : chain_ids) {
            if (!assembly.first().has_chain(chain_id)) {
                sz++;
            }
        }
        if (sz == chain_ids.size()) {  // if the first model in assembly doesn't
                                       // have any chain_ids, remove it!
            molecule.erase(i);
        }
    }
}

void remove_asymmetric_unit(molib::Molecule& molecule) {  // remove asymmetric
                                                          // unit if there are
                                                          // biological
                                                          // assembiles
    if (molecule.size() > 1) {
        if (molecule.has_element(0)) molecule.erase(0);
    }
}

molib::Atom::Grid set_grid_query(molib::Molecule& query_molecule,
                                 const string query_chain_ids) {
    // set grid with query chain atoms
    return molib::Atom::Grid(
        query_molecule.get_atoms(query_chain_ids, molib::Residue::protein, 1));
}
void find_neighbor_residues(
    molib::Molecule& molecule, molib::Atom::Grid& grid,
    const map<molib::Residue::res_tuple2, molib::Residue::res_tuple2>*
        aligned_residues = nullptr) {
    for (auto& assembly : molecule) {
        for (auto& model : assembly) {
            for (auto& chain : model) {
                for (auto& residue : chain) {
                    const string lresn =
                        (residue.rest() == molib::Residue::protein
                             ? "PPP"
                             : (residue.rest() == molib::Residue::nucleic
                                    ? "NNN"
                                    : residue.resn()));
                    const int lresi =
                        ((residue.rest() == molib::Residue::protein ||
                          residue.rest() == molib::Residue::nucleic)
                             ? 1
                             : residue.resi());
                    const string lig_code =
                        lresn + ":" + std::to_string(lresi) + ":" +
                        chain.chain_id() + ":" +
                        std::to_string(model.number()) + ":" +
                        std::to_string(assembly.number()) + ":" +
                        molecule.name() + ":" + molecule.br().name();
                    for (auto& atom : residue) {
                        vector<molib::Atom*> neighbors =
                            grid.get_neighbors(atom.crd(), 4.0);
                        for (molib::Atom* query_atom : neighbors) {
                            const molib::Residue& qr = query_atom->br();
                            const molib::Chain& qc = query_atom->br().br();
                            const molib::Residue::res_tuple2& r1 =
                                molib::Residue::res_tuple2(
                                    qc.chain_id(),
                                    std::to_string(
                                        help::get_one_letter(qr.resn())),
                                    qr.resi(), qr.ins_code());
                            if (aligned_residues &&
                                aligned_residues->find(r1) !=
                                    aligned_residues->end()) {
                                const molib::Residue::res_tuple2& r2 =
                                    aligned_residues->find(r1)->second;
                                model.set_remark(
                                    7,
                                    molib::Residue::res_tuple2(
                                        chain.chain_id(), lresn, lresi,
                                        (residue.rest() ==
                                                     molib::Residue::protein ||
                                                 residue.rest() ==
                                                     molib::Residue::nucleic
                                             ? ' '
                                             : residue.ins_code())),
                                    make_pair(r1, r2));
                            } else if (!aligned_residues) {
                                model.set_remark(
                                    7,
                                    molib::Residue::res_tuple2(
                                        chain.chain_id(), lresn, lresi,
                                        (residue.rest() ==
                                                     molib::Residue::protein ||
                                                 residue.rest() ==
                                                     molib::Residue::nucleic
                                             ? ' '
                                             : residue.ins_code())),
                                    make_pair(r1, r1));
                            }
                        }
                    }
                }
            }
        }
    }
}
void remove_not_ligands(molib::Molecule& molecule, molib::Atom::Grid& grid) {
    // find atoms in aligned molecules that are neighbors of query chains
    for (auto& assembly : molecule) {
        for (auto& model : assembly) {
            for (auto& chain : model) {
                // remove hetero, ion and water residues that are not close to
                // the query
                chain.remove_if([&grid](const molib::Residue& r) {
                    if (r.rest() == molib::Residue::hetero ||
                        r.rest() == molib::Residue::ion ||
                        r.rest() == molib::Residue::water) {
                        for (auto& atom : r) {
                            if (!grid.get_neighbors(atom.crd(), 4.0).empty()) {
                                return false;  // close enough, don't delete
                            }
                        }
                        return true;
                    }
                    return false;  // don't delete protein or nucleic
                });  // if atom is less than 4.0 A from some query atom, then
                     // this residue is a ligand of query, and is to keep
                bool keep_chain = false;
                for (auto& residue : chain) {
                    if (residue.rest() == molib::Residue::protein ||
                        residue.rest() == molib::Residue::nucleic) {
                        for (auto& atom : residue) {
                            // if atom is less than 4.0 A from some query atom,
                            // then this chain is a ligand of query, and is to
                            // keep
                            if (!grid.get_neighbors(atom.crd(), 4.0).empty()) {
                                keep_chain = true;
                                break;
                            }
                        }
                        if (keep_chain) break;
                    }
                }
                // if no protein or nucleic atom of this chain is a ligand of
                // the query, delete this chain
                if (!keep_chain) model.erase(chain.chain_id());
            }
        }
    }
}
void generate_biological_assemblies(const string& models, const bool hydrogens,
                                    const string& pdb_dirname,
                                    const string& qpdb_file, const string& qcid,
                                    const bool neighb, const bool rnolig,
                                    const string& bio, const bool ralch,
                                    const string& json_file,
                                    const string& bio_file, const bool noalch,
                                    const string& mols_name, const bool rasym) {
    JsonReader jr;
    molib::Atom::Grid query_grid;
    // the query file
    molib::Molecules mols;
    if (!mols_name.empty()) mols.set_name(mols_name);
    vector<string> aligned_chains;
    if (!qpdb_file.empty()) {
        parser::FileParser pr(
            Path::join(pdb_dirname, qpdb_file),
            (models == "all" ? parser::all_models : parser::first_model) |
                (hydrogens ? parser::hydrogens : 0));

        pr.parse_molecule(mols);

        mols.set_name(boost::filesystem::path(qpdb_file)
                          .stem()
                          .string());  // for probis-server, qpdb_file needs to
                                       // have pdb_id and chain_id
        mols.last().set_name(
            boost::filesystem::path(qpdb_file).stem().string());
        query_grid = genbio::set_grid_query(mols.last(), qcid);

        if (bio != "none")
            mols.last().init_bio(bio == "all" ? molib::Molecule::all_bio
                                              : molib::Molecule::first_bio);
        if (ralch)
            genbio::remove_assemblies(
                mols.last(),
                qcid);  // remove assemblies that don't contain query chains
        if (rnolig)
            genbio::remove_not_ligands(
                mols.last(),
                query_grid);  // remove chains that are not near query chains
        if (neighb) genbio::find_neighbor_residues(mols.last(), query_grid);
        aligned_chains.push_back(qcid);
    }
    // the aligned files from json
    if (!json_file.empty()) {
        jr.parse_JSON(json_file);
        for (auto& d : jr.root()) {
            try {  // if something goes wrong, e.g., pdb file is not found,
                   // don't exit..
                const string pdb_id = d["pdb_id"].asString();
                const string chain_ids = d["chain_id"].asString();
                const string pdb_file =
                    Path::join(pdb_dirname, pdb_id + ".pdb");

                dbgmsg(pdb_file << " " << chain_ids);

                parser::FileParser pr(pdb_file,
                                      (models == "all" ? parser::all_models
                                                       : parser::first_model) |
                                          (hydrogens ? parser::hydrogens : 0));
                pr.parse_molecule(mols);
                // distinguish between two use-cases
                if (pdb_id.size() == 4)
                    mols.last().set_name(pdb_id +
                                         chain_ids);  // e.g., for probis-server
                else
                    mols.last().set_name(pdb_id);  // e.g., for bslib where the
                                                   // pdb_id represents ligand
                                                   // binding site

                if (bio != "none")
                    mols.last().init_bio(
                        bio == "all"
                            ? molib::Molecule::all_bio
                            : molib::Molecule::first_bio);  // make biounits
                                                            // coordinates
                                                            // (optional)
                //mols.last().rotate(
                    //geometry::Matrix(d["alignment"][0]["rotation_matrix"],
                                     //d["alignment"][0]["translation_vector"]),
                    //true);  // inverse rotation
std::cout << "AASDF" << std::endl;
                if (ralch)
                    genbio::remove_assemblies(mols.last(),
                                              chain_ids);  // remove assemblies
                                                           // that don't contain
                                                           // aligned chains
                if (!qpdb_file.empty() && rnolig)
                    genbio::remove_not_ligands(
                        mols.last(), query_grid);  // remove chains that are not
                                                   // near query chains
                if (!qpdb_file.empty() && neighb) {
                    map<molib::Residue::res_tuple2, molib::Residue::res_tuple2>
                        res = common_ligands::json_to_map(
                            d["alignment"][0]["aligned_residues"]);
                    genbio::find_neighbor_residues(mols.last(), query_grid,
                                                   &res);
                }
                aligned_chains.push_back(chain_ids);
            } catch (exception& e) {
                log_warning << e.what() << " ... skipping ... " << endl;
            }
        }
    }
    // now it's safe to remove query chains (they are used in query_grid)
    if (noalch)
        genbio::remove_chains(mols, aligned_chains);  // remove query chains
    // remove asymmetric units when there are bio assemblies available
    if (rasym) {
        for (auto& molecule : mols) {
            genbio::remove_asymmetric_unit(molecule);
        }
    }
    // output bio file
    if (!bio_file.empty())
        fileio::output_file(mols, bio_file);  // output rotated bio assemblies
                                             // (or asymmetric units if NO
                                             // bioassembly)
}
};
}
