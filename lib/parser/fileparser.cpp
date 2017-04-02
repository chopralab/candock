#include "fileparser.hpp"

#include <memory>
#include <set>
#include <map>
#include <vector>
#include <string>

#include "geom3d/coordinate.hpp"
#include "helper/inout.hpp"
#include "helper/help.hpp"
#include "helper/debug.hpp"
#include "molib/nrset.hpp"
#include "molib/bond.hpp"
#include "fragmenter/fragmenter.hpp"

using namespace std;
using namespace Molib;

namespace Parser {
        void FileParser::rewind() {
                p->set_pos (0);
        }
        void FileParser::set_flags (unsigned int hm) {
                p->set_hm (hm);
        }
        bool FileParser::parse_molecule (Molecules &mols) {
                p->parse_molecule (mols);
                dbgmsg ("PARSED MOLECULES : " << endl << mols);
                return !mols.empty();
        }

        Molecules FileParser::parse_molecule() {
                Molecules mols;
                p->parse_molecule (mols);
                dbgmsg ("PARSED MOLECULES : " << endl << mols);

                if (mols.empty()) {
                        throw Error ("die : could not read molecules from file");
                }

                return mols;
        }

        FileParser::FileParser (const string &molecule_file, unsigned int hm,
                                const int num_occur) {
                prepare_parser (molecule_file, hm, num_occur);
        }

        void FileParser::prepare_parser (const string &molecule_file, unsigned int hm,
                                         const int num_occur) {
                auto ret = molecule_file.find_last_of (".");

                if (ret == string::npos) {
                        throw Error ("die : could not determine the file type of the input molecule");
                }

                string extension = molecule_file.substr (ret + 1);
                transform (extension.begin(), extension.end(),
                           extension.begin(), ::toupper);

                dbgmsg ("pdb reader options = " << hm << " molecule is a " << boolalpha
                        << (extension == "PDB" || extension == "ENT" ? "PDB file" :
                            (extension == "MOL2" ? "Mol2 file" : "undetermined file type")));

                if (extension == "PDB" || extension == "ENT") {
                        p = new PdbParser (molecule_file, hm, num_occur);
                } else if (extension == "MOL2") {
                        p = new Mol2Parser (molecule_file, hm, num_occur);
                } else {
                        throw Error ("die : could not determine the file type of the input molecule");
                }

                if (Inout::file_size (molecule_file) <= 0) {
                        throw Error (string ("die : file not valid: ") + molecule_file + ". Check to see if it exists and has contents!");
                }
        }
};

