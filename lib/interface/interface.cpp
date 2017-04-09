#include "interface.hpp"
#include "parser/fileparser.hpp"
#include "molib/molecules.hpp"
#include "molib/molecule.hpp"
#include "parser/fileparser.hpp"
#include "score/score.hpp"
#include "modeler/forcefield.hpp"

#include <boost/filesystem/path.hpp>
#include <memory>

std::unique_ptr <Molib::Molecules> __receptor;
std::unique_ptr <Molib::Molecules> __ligand;
std::unique_ptr <Molib::Score> __score;
std::unique_ptr <Molib::Atom::Grid> __gridrec;
std::unique_ptr <OMMIface::ForceField> __ffield;

const char* initialize_receptor(const char* filename) {
        Parser::FileParser rpdb(filename, Parser::first_model, 1);
        
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
        ss << __receptor->element(0);

        return ss.str().c_str();
}

const char*  initialize_ligand(const char* filename) {
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
        ss << __ligand->element(0);

        return ss.str().c_str();
}

void   initialize_scoring(const char* filename) {

        __score = std::unique_ptr<Molib::Score> (
                        new  Molib::Score( "mean", "reduced", "radial", 6.0, 0.01));

        __score->define_composition(__receptor->get_idatm_types(),
                                    __ligand->get_idatm_types())
                .process_distributions_file("data/csd_complete_distance_distributions.txt")
                .compile_scoring_function()
                .parse_objective_function(filename, 10.0 );
}

void   initialize_ffield(const char* data_dir) {
        boost::filesystem::path p(data_dir);

        __ffield = std::unique_ptr<OMMIface::ForceField>(new OMMIface::ForceField);

        __ffield->parse_gaff_dat_file( (p / "gaff.dat").string() )
              .add_kb_forcefield(*__score, 0.01)
              .parse_forcefield_file( (p / "amber10.xml").string() )
              .parse_forcefield_file( (p / "tip3p.xml").string() );

        __receptor->element(0).prepare_for_mm(*__ffield, *__gridrec);
}

float  calculate_score() {
        return __score->non_bonded_energy(*__gridrec, (*__ligand)[0]);
}

void   set_positions_ligand(const unsigned long* atoms, const float* positions, unsigned long size) {
        Molib::Residue* residue = __ligand->element(0).get_residues().at(0);

        for ( unsigned long i = 0; i < size; ++i) {
                residue->element(atoms[i]).set_crd( Geom3D::Point(positions[i * 3 + 0],
                                                                  positions[i * 3 + 1],
                                                                  positions[i * 3 + 2]
                                                                 ) );
        }
}

void  set_positions_receptor(const unsigned long* atoms, const float* positions, unsigned long size) {
        unsigned long resi = 0;
        Molib::Residue* residue = __receptor->element(0).get_residues().at(resi++);

        for ( unsigned long i = 0; i < size; ++i) {
                
                while ( ! residue->has_element(atoms[i]) ) {
                        residue = __receptor->element(0).get_residues().at(resi++);
                }

                if ( ! residue->has_element(atoms[i]) ) {
                        cout << "Error: could not atom = " << atoms[i] << endl;
                        break;
                }

                residue->element(atoms[i]).set_crd( Geom3D::Point(positions[i * 3 + 0],
                                                                  positions[i * 3 + 1],
                                                                  positions[i * 3 + 2]
                                                                 ) );
        }

        __gridrec.release();
        __gridrec = std::unique_ptr<Molib::Atom::Grid>(new Molib::Atom::Grid(__receptor->get_atoms()));
}