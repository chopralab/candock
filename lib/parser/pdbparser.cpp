/* Copyright (c) 2016-2019 Chopra Lab at Purdue University, 2013-2016 Janez Konc at National Institute of Chemistry and Samudrala Group at University of Washington
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

#include "fileparser.hpp"

#include <boost/regex.hpp>
#include <boost/algorithm/string.hpp>
#include <boost/filesystem.hpp>
#include <boost/date_time/posix_time/posix_time.hpp>
#include <boost/algorithm/string/split.hpp>
#include <boost/algorithm/string.hpp>
#include <boost/asio/ip/host_name.hpp>

using namespace std;
using namespace Molib;

namespace Parser {

        void FileParser::PdbParser::parse_molecule (Molecules &mols) {
                vector<string> pdb_raw;
                dbgmsg ("num_occur = " << __num_occur);
                Inout::read_file (__molecule_file, pdb_raw, __pos,
                                  Inout::panic, __num_occur, "REMARK   5 MOLECULE");
                boost::smatch m;
                set<char> bio_chain;
                int biomolecule_number = -1;
                bool found_molecule = false, found_assembly = false, found_model = false;
                map<const Model *, map<const int, Atom *>> atom_number_to_atom;
                Chain *chain = nullptr;
                Residue *residue = nullptr;
                bool ter_found=false;

                char chain_id_replacement = 'A';

                for (string &line : pdb_raw) {

                        if ( (__hm & skip_hetatm) && line.compare (0, 6, "HETATM") == 0) {
                                continue;
                        }

                        if ( (__hm & skip_atom) && line.compare (0, 4, "ATOM")   == 0) {
                                continue;
                        }

                        if (line.compare (0, 4, "ATOM") == 0 || line.compare (0, 6, "HETATM") == 0) {

                                if ( (__hm & docked_poses_only) && (! ter_found)) {
                                        continue;
                                }

                                if ( (__hm & protein_poses_only)&& (  ter_found)) {
                                        continue;
                                }

                                __generate_molecule (mols, found_molecule, "");
                                __generate_assembly (mols, found_assembly, 0, "ASYMMETRIC UNIT");
                                __generate_model (mols, found_model, 1);
                                string record_name = boost::algorithm::trim_copy (line.substr (0,6));
                                int atom_number = stoi (line.substr (6, 5));
                                dbgmsg ("atom_number = " << atom_number);
                                string atom_name = boost::algorithm::trim_copy (line.substr (12,4));
                                char alt_loc = line.at (16);
                                string resn = boost::algorithm::trim_copy (line.substr (17,3));
                                char chain_id = line.at (21);
                                                                
                                int resi = stoi (line.substr (22, 4));
                                char ins_code = line.at (26);
                                double x_coord = stof (line.substr (30, 8));
                                double y_coord = stof (line.substr (38, 8));
                                double z_coord = stof (line.substr (46, 8));
                                Geom3D::Coordinate crd (x_coord, y_coord, z_coord);
                                string element = line.size() > 77 ? boost::algorithm::trim_copy (line.substr (76,2)) : "";

                                if (element.empty()) {
                                        if (std::isdigit (line.at (12))) {
                                                element = line.at (13);
                                        } else {
                                                element = boost::algorithm::trim_copy (line.substr (12,2));
                                        }
                                }

                                string idatm_type = line.size() > 80 ? boost::algorithm::trim_copy (line.substr (80, 5)) : "???";
                                idatm_type = help::idatm_mask.count (idatm_type) ? idatm_type : "???";
                                string gaff_type = line.size() > 85 ? boost::algorithm::trim_copy (line.substr (85, 5)) : "???";
                                string rest_of_line = line.size() > 90 ? boost::algorithm::trim_copy (line.substr (90)) : "";
                                const bool hydrogen = (element == "H" || (atom_name.size() == 1 && atom_name.at (0) == 'H'));
                                dbgmsg ("hydrogen = " << boolalpha << hydrogen);

                                if ( (__hm & hydrogens) || !hydrogen) {
                                        if (alt_loc == ' ' || alt_loc == 'A') {
                                                if (!mols.last().is_modified (Residue::res_tuple2 (chain_id, resn, resi, ins_code))) {
                                                        Residue::res_type rest;

                                                        if (help::amino_acids.find (resn) != help::amino_acids.end()) {
                                                                rest = Residue::protein;
                                                        } else if (help::nucleic_acids.find (resn) != help::nucleic_acids.end()) {
                                                                rest = Residue::nucleic;
                                                        } else if (help::ions.find (resn) != help::ions.end()) {
                                                                rest = Residue::ion;
                                                        } else if (resn == "HOH" || resn == "WAT" || resn == "WAM" || resn == "WAU") {
                                                                resn = "HOH";
                                                                rest = Residue::water;
                                                        } else {
                                                                // if nothing known then hetero...
                                                                rest = Residue::hetero;
                                                        }

                                                        if (rest == Residue::water && __hm & skip_atom) {
                                                                continue;
                                                        }

                                                        if ( (__hm & sparse_macromol) &&
                                                                        ( (rest == Residue::protein && atom_name != "CA")
                                                                          || (rest == Residue::nucleic && atom_name != "P"))) {
                                                                continue;
                                                        }

                                                        if (!__giant_molecule || rest != Residue::protein || atom_name == "CA") {
                                                                Model &model = mols.last().last().last();

                                                                if ( chain_id == ' ' ) {
                                                                        chain_id = chain_id_replacement;
                                                                        if ( model.has_chain(chain_id) 
                                                                                && model.chain(chain_id).size() > 0
                                                                                && model.chain(chain_id).last().resi() > resi ) {
                                                                                chain_id = ++chain_id_replacement;
                                                                        }
                                                                }
                                                                
                                                                // chains may alternate : A B A ...
                                                                if (!model.has_chain (chain_id)) {
                                                                        chain = &model.add (new Chain (chain_id));
                                                                } else {
                                                                        chain = &model.chain (chain_id);
                                                                }

                                                                if (!chain->has_residue (Residue::res_pair (resi, ins_code))) {
                                                                        residue = &chain->add (new Residue (resn, resi, ins_code, rest));
                                                                }

                                                                Atom &a = residue->add (new Atom (atom_number,
                                                                                                  atom_name,
                                                                                                  crd,
                                                                                                  help::idatm_mask.at (idatm_type),
                                                                                                  element
                                                                                                 ));

                                                                if (gaff_type != "???") {
                                                                        a.set_gaff_type (gaff_type);
                                                                }

                                                                if (rest_of_line != "") {
                                                                        a.set_members (rest_of_line);
                                                                }

                                                                atom_number_to_atom[&model][atom_number] = &a;
                                                        }
                                                }
                                        }
                                }
                        } else if (line.compare (0, 14, "REMARK   4 NRP") == 0) {
                                string name = "";

                                if (boost::regex_search (line, m, boost::regex ("REMARK   4 NRPDB\\s+(.*)"))) {
                                        if (m[1].matched) {
                                                name = m[1].str();
                                                dbgmsg (name);
                                        }
                                }

                                mols.set_name (name);
                        } else if (line.compare (0, 14, "REMARK   5 MOL") == 0) {
                                string name = "";

                                if (boost::regex_search (line, m, boost::regex ("REMARK   5 MOLECULE\\s+(.*)"))) {
                                        if (m[1].matched) {
                                                name = m[1].str();
                                                dbgmsg (name);
                                        }
                                }

                                found_molecule = false;
                                __generate_molecule (mols, found_molecule, name);
                        } else if (line.compare (0, 14, "REMARK   6 BIO") == 0 || line.compare (0, 14, "REMARK   6 ASY") == 0) {
                                int number = 0;
                                string name = "";

                                if (boost::regex_search (line, m, boost::regex ("REMARK   6 (\\S+ \\S+)\\s+(\\d+)"))) {
                                        if (m[1].matched && m[2].matched) {
                                                name = m[1].str();
                                                number = stoi (m[2].str());
                                        }
                                }

                                found_assembly = false;
                                __generate_assembly (mols, found_assembly, number, name);
                        } else if (line.compare (0, 5, "MODEL") == 0) {
                                if (boost::regex_search (line, m, boost::regex ("MODEL\\s+(\\d+)"))) {
                                        if (m[1].matched) {
                                                
                                                if ( __hm & skip_atom || __hm & protein_poses_only) {
                                                        found_molecule = false;
                                                        found_assembly = false;
                                                }

                                                ter_found = false;

                                                __generate_molecule (mols, found_molecule, "");
                                                __generate_assembly (mols, found_assembly, 0, "ASYMMETRIC UNIT");

                                                if (__hm & all_models) {
                                                        found_model = false;
                                                }
                                                __generate_model (mols, found_model, stoi (m[1].str()));
                                        }
                                } else {
                                        throw Error ("die : model number not found in pdb");
                                }
                        } else if (line.compare (0, 6, "ENDMDL") == 0) {
                                if (__hm & first_model) {
                                        return;
                                }
                        } else if (line.compare (0, 3, "TER") == 0) {
                                ter_found = true;
                        } else if (line.compare (0, 10, "REMARK 350") == 0) {
                                if (boost::regex_search (line, m, boost::regex ("REMARK 350 BIOMOLECULE:\\s*(\\d+)"))) {
                                        if (m[1].matched) {
                                                biomolecule_number = stoi (m[1].str());

                                                 // if there were no REMARKs...
                                                if (biomolecule_number == 1) {
                                                        __generate_molecule (mols, found_molecule, "");
                                                }
                                        }

                                        bio_chain.clear();
                                } else if (boost::regex_search (line, m, boost::regex ("REMARK 350 APPLY THE FOLLOWING TO CHAINS:\\s?(\\S{1})?,?\\s?(\\S{1})?,?\\s?(\\S{1})?,?\\s?(\\S{1})?,?\\s?(\\S{1})?,?\\s?(\\S{1})?,?\\s?(\\S{1})?,?\\s?(\\S{1})?,?\\s?(\\S{1})?,?"))
                                                || boost::regex_search (line, m, boost::regex ("REMARK 350                    AND CHAINS:\\s?(\\S{1})?,?\\s?(\\S{1})?,?\\s?(\\S{1})?,?\\s?(\\S{1})?,?\\s?(\\S{1})?,?\\s?(\\S{1})?,?\\s?(\\S{1})?,?\\s?(\\S{1})?,?\\s?(\\S{1})?,?"))) {
                                        if (m.empty()) {
                                                throw Error ("die : problem in REMARK 350 APPLY");
                                        }

                                        for (auto &el : m) {
                                                if (el != *m.begin() && el.matched) {
                                                        bio_chain.insert (el.str().at (0));
                                                }
                                        }
                                } else if (boost::regex_search (line, m, boost::regex ("REMARK 350   BIOMT(\\d+)\\s+(\\d+)\\s+(\\S+)\\s+(\\S+)\\s+(\\S+)\\s+(\\S+)"))) {
                                        if (m.empty()) {
                                                throw Error ("die : problem in REMARK 350   BIOMT");
                                        }

                                        if (m[1].matched && m[2].matched && m[3].matched && m[4].matched && m[5].matched && m[6].matched) {
                                                if (biomolecule_number == -1) {
                                                        throw Error ("die : missing 'REMARK 350 BIOMOLECULE'");
                                                }

                                                int row = stoi (m[1].str());
                                                int rn = stoi (m[2].str());
                                                double a = stod (m[3].str());
                                                double b = stod (m[4].str());
                                                double c = stod (m[5].str());
                                                double w = stod (m[6].str());

                                                if (rn < 15) {
                                                        if (row == 1) {
                                                                mols.last().add_bio_chain (biomolecule_number, bio_chain);
                                                        }

                                                        mols.last().add_bio_rota (biomolecule_number, rn, row - 1, make_tuple (a, b, c, w));
                                                } else {
                                                         // if giant molecule, read only CA
                                                        __giant_molecule = true;
                                                }
                                        }
                                }
                        } else if (line.compare (0, 6, "MODRES") == 0) {
                                string resn = line.substr (24, 3);

                                if (help::amino_acids.find (resn) != help::amino_acids.end()) {
                                        string resn_mod = line.substr (12, 3);
                                        char chain_id = line.at (16);
                                        int resi = stoi (line.substr (18, 4));
                                        char ins_code = line.at (22);
                                        __generate_molecule (mols, found_molecule, "");

                                        // modified 'standard' residue gets saved,e.g., MET --> MSE, glycosylated...
                                        mols.last().add_modified (Residue::res_tuple2 (chain_id, resn_mod, resi, ins_code));
                                }
                        } else if (line.compare (0, 4, "SITE") == 0 && boost::regex_search (line, m, boost::regex ("SITE\\s+\\d+\\s+(\\S+)\\s+\\S+(\\s+\\S+\\s?\\S{1}\\s{0,3}\\d{1,4}\\S?)?(\\s+\\S+\\s?\\S{1}\\s{0,3}\\d{1,4}\\S?)?(\\s+\\S+\\s?\\S{1}\\s{0,3}\\d{1,4}\\S?)?(\\s+\\S+\\s?\\S{1}\\s{0,3}\\d{1,4}\\S?)?"))) {
                                if (m.empty()) {
                                        throw Error ("die : problem in SITE");
                                }

                                string site_name = (m[1].matched ? m[1].str() : "");

                                for (auto &el : m) {
                                        if (el != *m.begin() && el.matched) {
                                                string s = boost::algorithm::trim_copy (el.str());

                                                if (s.size() > 8) {
                                                        string resn = s.substr (0, 3);
                                                        char chain_id = s.at (4);
                                                        int resi = stoi (s.substr (5, 4));
                                                        char ins_code = (s.size() == 10 ? s.at (9) : ' ');

                                                        if (help::amino_acids.find (resn) != help::amino_acids.end()) {
                                                                __generate_molecule (mols, found_molecule, "");
                                                                mols.last().add_site (site_name, Residue::res_tuple2 (chain_id, resn, resi, ins_code));
                                                        }
                                                }
                                        }
                                }
                        } else if (line.compare (0, 16, "REMARK   7 ALRES") == 0) {
                                // added ALRES since some PDB's (e.g. 1abw) have remark 7
                                string text (line.substr (17));
                                vector<string> vs = help::ssplit (text, ",");

                                if (vs.empty()) {
                                        throw ("die: remark 7 is empty");
                                }

                                vector<string> token = help::ssplit (vs[0], ":");
                                char lchain_id = token[0].at (0);
                                string lresn = token[1];
                                int lresi = stoi (token[2]);
                                char lins_code = token[3].at (0);
                                Model &model = mols.last().last().last();
                                vector<string>::iterator it = vs.begin();
                                it++;

                                for (; it != vs.end(); it++) {
                                        string &s = *it;
                                        vector<string> tok = help::ssplit (s, "=");
                                        vector<string> tok1 = help::ssplit (tok[0], ":");
                                        vector<string> tok2 = help::ssplit (tok[1], ":");
                                        dbgmsg ("ADDING TO REMARK " << "[" << lchain_id << ":" << lresn << ":" << lresi << ":" << lins_code << "] => ["
                                                << tok1[0].at (0) << ":" << tok1[1] << ":" << stoi (tok1[2]) << ":" << tok1[3].at (0)
                                                << "="
                                                << tok2[0].at (0) << ":" << tok2[1] << ":" << stoi (tok2[2]) << ":" << tok2[3].at (0)
                                                << "]");
                                        model.set_remark (7, Residue::res_tuple2 (lchain_id, lresn, lresi, lins_code),
                                                          make_pair (Residue::res_tuple2 (tok1[0].at (0), tok1[1], stoi (tok1[2]), tok1[3].at (0)),
                                                                     Residue::res_tuple2 (tok2[0].at (0), tok2[1], stoi (tok2[2]), tok2[3].at (0))));
                                }
                        } else if (line.compare (0, 16, "REMARK   8 RIGID") == 0) {
                                Model &model = mols.last().last().last();
                                vector<string> vs = help::ssplit (line.substr (17), " ");
                                Atom::Set core, join;

                                for (size_t i = 1; i < vs.size(); ++i)  {
                                        const int atom_number = stoi (vs[i].substr (1));

                                        if (vs[i][0] == 'c') {
                                                core.insert (atom_number_to_atom[&model][atom_number]);
                                        } else if (vs[i][0] == 'j') {
                                                join.insert (atom_number_to_atom[&model][atom_number]);
                                        }
                                }

                                const int seed_id = stoi (vs[0]);
                                model.get_rigid().push_back (Fragmenter::Fragment (core, join, seed_id));
#ifndef NDEBUG
                                dbgmsg ("seed_id = " << seed_id);

                                for (size_t i = 1; i < vs.size(); i++) {
                                        dbgmsg (vs[i]);
                                }

                                for (auto &a : core) {
                                        dbgmsg ("CORE ATOM : " << *a);
                                }

                                for (auto &a : join) {
                                        dbgmsg ("JOIN ATOM : " << *a);
                                }

                                dbgmsg (model.get_rigid());
#endif
                        } else if (line.compare (0, 15, "REMARK   8 ROTA") == 0) {
                                Model &model = mols.last().last().last();
                                vector<string> vs = help::ssplit (line.substr (16), " ");
                                Atom &a1 = *atom_number_to_atom[&model][stoi (vs[1])];
                                Atom &a2 = *atom_number_to_atom[&model][stoi (vs[2])];
                                a1.connect (a2).set_rotatable (vs[0]);
                        } else if (line.compare (0, 19, "REMARK   8 BONDTYPE") == 0) {

                                if ( (__hm & docked_poses_only) && (! ter_found)) {
                                        continue;
                                }

                                if ( (__hm & protein_poses_only)&& (  ter_found)) {
                                        continue;
                                }

                                vector<string> vs = help::ssplit (line.substr (20), " ");

                                 // to avoid when bond_gaff_type is undefined
                                if (vs.size() == 4) {
                                        Model &model = mols.last().last().last();
                                        Atom &a1 = *atom_number_to_atom[&model][stoi (vs[2])];
                                        Atom &a2 = *atom_number_to_atom[&model][stoi (vs[3])];
                                        auto &bond = a1.connect (a2);
                                        bond.set_bo (stoi (vs[0]));
                                        bond.set_bond_gaff_type (vs[1]);
                                }
                        } else if (line.compare (0, 6, "CONECT") == 0) {
                                const string ln = boost::algorithm::trim_right_copy (line.substr (6));
                                dbgmsg ("--" << ln << "--");
                                vector<int> anum;

                                for (size_t i = 0; i < ln.size(); i+=5) {
                                        const int atom_number = stoi (ln.substr (i, 5));
                                        dbgmsg (atom_number);
                                        anum.push_back (atom_number);
                                }

                                // the CONECT words are valid for the whole molecule
                                Molecule &molecule = mols.last();

                                for (auto &assembly : molecule) {
                                        for (auto &model : assembly) {
                                                vector<int>::iterator it = anum.begin();
                                                auto atom_iterator1 = atom_number_to_atom[&model].find (*anum.begin());
                                                dbgmsg ( (atom_iterator1 == atom_number_to_atom[&model].end()
                                                          ? "warning: atom1 not found (check CONECT words)" : ""));

                                                if (atom_iterator1 != atom_number_to_atom[&model].end()) {
                                                        dbgmsg ("found atom = " << boolalpha
                                                                << (atom_iterator1 != atom_number_to_atom[&model].end()));
                                                        Atom &a1 = *atom_iterator1->second;
                                                        it++;

                                                        while (it != anum.end()) {
                                                                auto atom_iterator2 = atom_number_to_atom[&model].find (*it);
                                                                dbgmsg ( (atom_iterator2 == atom_number_to_atom[&model].end()
                                                                          ? "warning: atom2 not found (check CONECT words)" : ""));

                                                                if (atom_iterator2 != atom_number_to_atom[&model].end()) {
                                                                        dbgmsg ("found atom = " << boolalpha
                                                                                << (atom_iterator2 != atom_number_to_atom[&model].end()));
                                                                        Atom &a2 = *atom_iterator2->second;
                                                                        dbgmsg (a1 << " " << a2);
                                                                        const string &resn1 = a1.br().resn();
                                                                        const string &resn2 = a2.br().resn();

                                                                        // don't allow inter-residue bonds e.g., between MG and ATP, or between FAD and CYS (issue #114)
                                                                        if (help::ions.count (resn1) == 0 
                                                                        &&  help::ions.count (resn2) == 0 ){
                                                                                a1.connect (a2);
                                                                        }

#ifndef NDEBUG
                                                                        else {
                                                                                dbgmsg ("Ignoring interresidue bond between " << resn1 << " and " << resn2 << ". Check if that is what you intended!");
                                                                        }

#endif
                                                                }

                                                                it++;
                                                        }
                                                }
                                        }
                                }

                                connect_bonds (get_bonds_in (molecule.get_atoms()));
                        } else if (line.compare (0, 21, "REMARK  20 non-binder") == 0 && (__hm & docked_poses_only)) {
                                return;
                        }
                }
        }
}
