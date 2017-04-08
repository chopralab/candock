#include "interface.hpp"
#include "parser/fileparser.hpp"
#include "molib/molecules.hpp"
#include "molib/molecule.hpp"
#include "parser/fileparser.hpp"
#include "score/score.hpp"

#include <memory>

std::unique_ptr <Molib::Molecules> __receptor;
std::unique_ptr <Molib::Molecules> __ligand;
std::unique_ptr <Molib::Score> __score;
std::unique_ptr <Molib::Atom::Grid> __gridrec;

extern "C" const char* initialize_receptor(const char* filename) {
        Parser::FileParser rpdb(filename, Parser::first_model);
        
        __receptor = std::unique_ptr<Molib::Molecules>(new Molib::Molecules);

        rpdb.parse_molecule(*__receptor);

        __receptor->compute_idatm_type()
                .compute_hydrogen()
                .compute_bond_order()
                .compute_bond_gaff_type()
                .refine_idatm_type()
                .erase_hydrogen()  // needed because refine changes connectivities
                .compute_hydrogen()   // needed because refine changes connectivities
                .compute_ring_type()
                .compute_gaff_type()
                .compute_rotatable_bonds() // relies on hydrogens being assigned
                .erase_hydrogen();

        __gridrec = std::unique_ptr<Molib::Atom::Grid>(new Molib::Atom::Grid(__receptor->get_atoms()));

        std::stringstream ss;
        ss << (*__receptor)[0];

        return ss.str().c_str();
}

extern "C" const char*  initialize_ligand(const char* filename) {
        Parser::FileParser rpdb(filename, Parser::first_model, 1);
        
        __ligand = std::unique_ptr<Molib::Molecules>(new Molib::Molecules);

        rpdb.parse_molecule(*__ligand);

        __ligand->compute_idatm_type()
                .compute_hydrogen()
                .compute_bond_order()
                .compute_bond_gaff_type()
                .refine_idatm_type()
                .erase_hydrogen()  // needed because refine changes connectivities
                .compute_hydrogen()   // needed because refine changes connectivities
                .compute_ring_type()
                .compute_gaff_type()
                .compute_rotatable_bonds() // relies on hydrogens being assigned
                .erase_hydrogen();

        std::stringstream ss;
        ss << (*__ligand)[0];

        return ss.str().c_str();
}

extern "C" void   initialize_scoring(const char* filename) {
        __score = std::unique_ptr<Molib::Score> (
                        new  Molib::Score( "reduced", "mean", "radial", 6, 0.10));

        __score->define_composition(__receptor->get_idatm_types(),
                                    __ligand->get_idatm_types())
                .process_distributions_file("data/csd_complete_distance_distributions.txt")
                .compile_scoring_function()
                .parse_objective_function(filename, 10.0 );
        
}
