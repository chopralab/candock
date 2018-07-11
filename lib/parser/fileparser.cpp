#include "candock/parser/fileparser.hpp"

#include <map>
#include <memory>
#include <set>
#include <string>
#include <vector>

#include "candock/fragmenter/fragmenter.hpp"
#include "candock/geometry/coordinate.hpp"
#include "candock/helper/debug.hpp"
#include "candock/helper/help.hpp"
#include "candock/helper/inout.hpp"
#include "candock/molib/bond.hpp"
#include "candock/molib/nrset.hpp"

using namespace std;
using namespace candock::molib;

namespace candock {
namespace parser {
void FileParser::set_flags(unsigned int hm) { p->set_hm(hm); }
bool FileParser::parse_molecule(Molecules& mols) {
    p->parse_molecule(mols);
    dbgmsg("PARSED MOLECULES : " << endl << mols);
    return !mols.empty();
}

Molecules FileParser::parse_molecule() {
    Molecules mols;
    p->parse_molecule(mols);
    dbgmsg("PARSED MOLECULES : " << endl << mols);

    if (mols.empty()) {
        throw Error("die : could not read molecules from file");
    }

    return mols;
}

FileParser::FileParser(const string& molecule_file, unsigned int hm,
                       const int num_occur) {
    prepare_parser(molecule_file, hm, num_occur);
}

void FileParser::prepare_parser(const string& molecule_file, unsigned int hm,
                                const int num_occur) {
    if (Inout::file_size(molecule_file) <= 0) {
        throw Error("die : file not valid: " + molecule_file +
                    ". Check to see if it exists and has contents!");
    }

    auto ret = molecule_file.find_last_of(".");

    if (ret == string::npos) {
        throw Error(
            "die : could not determine the file type of the input molecule");
    }

    string extension = molecule_file.substr(ret + 1);
    transform(extension.begin(), extension.end(), extension.begin(), ::toupper);

    std::shared_ptr<istream> temp_molecule_stream =
        std::make_shared<ifstream>(molecule_file, std::ios::in);

    prepare_parser(temp_molecule_stream, extension, hm, num_occur);
}

void FileParser::prepare_parser(std::shared_ptr<std::istream>& stream,
                                const std::string& extension, unsigned int hm,
                                const int num_occur) {
    molecule_stream = stream;

    dbgmsg("pdb reader options = "
           << hm << " molecule is a " << boolalpha
           << (extension == "PDB" || extension == "ENT"
                   ? "PDB file"
                   : (extension == "MOL2" ? "Mol2 file"
                                          : "undetermined file type")));

    if (extension == "PDB" || extension == "ENT") {
        p = std::unique_ptr<Parser>(
            new PdbParser(*molecule_stream, hm, num_occur));
    } else if (extension == "MOL2") {
        p = std::unique_ptr<Parser>(
            new Mol2Parser(*molecule_stream, hm, num_occur));
    } else {
        throw Error(
            "die : could not determine the file type of the input molecule");
    }
}
};
}
