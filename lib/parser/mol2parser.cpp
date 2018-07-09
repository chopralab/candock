#include "candock/parser/fileparser.hpp"

using namespace candock::molib;

namespace candock {
namespace parser {
        void FileParser::Mol2Parser::parse_molecule (Molecules &mols) {
                vector<string> mol2_raw;
                {
                        std::lock_guard<std::mutex> gaurd(__concurrent_read_mtx);
                        Inout::read_stream (__stream, mol2_raw, __num_occur, "@<TRIPOS>MOLECULE");
                }
                bool found_molecule = false, found_assembly = false, found_model = false;
                map<const Model *, map<const int, Atom *>> atom_number_to_atom;

                for (size_t i = 0; i < mol2_raw.size(); ++i) {
                        const string &line = mol2_raw[i];

                        if (line.find ("@<TRIPOS>ATOM") != string::npos) {
                                __generate_molecule (mols, found_molecule, "");
                                __generate_assembly (mols, found_assembly, 0, "ASYMMETRIC UNIT");
                                __generate_model (mols, found_model, 1);
                                Chain *chain = nullptr;
                                Residue *residue = nullptr;
                                std::set<std::string> used_names;

                                for (i = i + 1; i < mol2_raw.size(); ++i) {
                                        const string &atom_line = mol2_raw[i];
                                        dbgmsg (atom_line);

                                        if (atom_line.empty() || atom_line[0] == '@') {
                                                break;
                                        }

                                        stringstream ss (atom_line);
                                        int atom_id, subst_id;
                                        string atom_name, sybyl_type, subst_name;
                                        double x, y, z, charge;
                                        ss >> atom_id >> atom_name >> x >> y >> z >> sybyl_type
                                           >> subst_id >> subst_name >> charge;

                                        subst_id = 1; // fix value to prevent multiresidue small molecules (issue #110)
                                        subst_name = "LIG"; // prevent idatm typing for standard residues (issue #111)

                                        const string element = help::sybyl.count (sybyl_type) ? help::sybyl.at (sybyl_type) : "";

                                        subst_name = subst_name.substr (0, 3); // truncate longer than 3-lett (pde5a bug)

                                        const bool hydrogen = (element == "H" || (atom_name.size() == 1 && atom_name.at (0) == 'H'));

                                        dbgmsg ("atom_id = " << atom_id << " atom_name = " << atom_name
                                                << " x = " << x << " y = " << y << " z = " << z
                                                << " sybyl_type = " << sybyl_type << " subst_id = " << subst_id
                                                << " subst_name = " << subst_name << " charge = " << charge
                                                << " element = " << element << " hydrogen = " << hydrogen);
                                        geometry::Coordinate crd (x, y, z);

                                        if ( (__hm & hydrogens) || !hydrogen) {
                                                Residue::res_type rest (Residue::hetero);

                                                if (!__giant_molecule || rest != Residue::protein || atom_name == "CA") {
                                                        Model &model = mols.last().last().last();

                                                        if (!model.has_chain ('A')) {
                                                                chain = &model.add (new Chain ('A'));
                                                        }

                                                        if (!chain->has_residue (Residue::res_pair (subst_id, ' '))) {
                                                                residue = &chain->add (new Residue (subst_name, subst_id, ' ', rest));
                                                        }

                                                        if (used_names.find(atom_name) != used_names.end()) {
                                                                atom_name += std::to_string(used_names.size());
                                                        }

                                                        used_names.insert(atom_name);

                                                        Atom &a = residue->add (new Atom (atom_id,
                                                                                          atom_name,
                                                                                          crd,
                                                                                          help::idatm_mask.at ("???"),
                                                                                          element,
                                                                                          sybyl_type
                                                                                         ));

                                                        atom_number_to_atom[&model][atom_id] = &a;
                                                }
                                        }
                                }

                                --i;
                        } else if (line.find ("@<TRIPOS>MOLECULE") != string::npos) {
                                if (i + 1 >= mol2_raw.size()) {
                                        throw Error ("die : wrong format of mol2 file");
                                }

                                const string &name = mol2_raw[i + 1];
                                dbgmsg (name);
                                found_molecule = false;
                                found_assembly = false;
                                found_model = false;
                                __generate_molecule (mols, found_molecule, name);
                                __generate_assembly (mols, found_assembly, 0, "ASYMMETRIC UNIT");
                                __generate_model (mols, found_model, 1);
                        } else if (line.find ("@<TRIPOS>BOND") != string::npos) {
                                Molecule &molecule = mols.last(); // the CONECT words are valid for the whole molecule

                                for (i = i + 1; i < mol2_raw.size(); ++i) {
                                        const string &bond_line = mol2_raw[i];

                                        if (bond_line.empty() || bond_line[0] == '@') {
                                                break;
                                        }

                                        int bond_id, origin_atom_id, target_atom_id;
                                        string bond_type;
                                        stringstream ss (bond_line);
                                        ss >> bond_id >> origin_atom_id >> target_atom_id >> bond_type;

                                        for (auto &assembly : molecule) {
                                                for (auto &model : assembly) {
                                                        auto it1 = atom_number_to_atom[&model].find (origin_atom_id);
                                                        auto it2 = atom_number_to_atom[&model].find (target_atom_id);

                                                        if (it1 != atom_number_to_atom[&model].end()
                                                                        && it2 != atom_number_to_atom[&model].end()) {
                                                                Atom &a1 = *it1->second;
                                                                Atom &a2 = *it2->second;
                                                                a1.connect (a2);
                                                        }
                                                }
                                        }
                                }

                                connect_bonds (get_bonds_in (molecule.get_atoms()));
                                --i;
                        }
                }
        }
}
}
