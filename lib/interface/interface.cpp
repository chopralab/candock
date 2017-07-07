#include "interface.hpp"
#include "parser/fileparser.hpp"
#include "molib/molecules.hpp"
#include "molib/molecule.hpp"
#include "parser/fileparser.hpp"
#include "score/score.hpp"
#include "modeler/forcefield.hpp"
#include "modeler/modeler.hpp"

#include <boost/filesystem/path.hpp>
#include <memory>

std::unique_ptr <Molib::Molecules> __receptor;
std::unique_ptr <Molib::Molecules> __ligand;
std::unique_ptr <Molib::Score> __score;
std::unique_ptr <Molib::Atom::Grid> __gridrec;
std::unique_ptr <OMMIface::ForceField> __ffield;
std::string __error_string = "";

const char* cd_get_error() {
        return __error_string.c_str();
}

size_t initialize_receptor(const char* filename) {
        try {
                Parser::FileParser rpdb (filename, Parser::first_model, 1);

                __receptor = std::unique_ptr<Molib::Molecules> (new Molib::Molecules);

                rpdb.parse_molecule (*__receptor);

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

                __gridrec = std::unique_ptr<Molib::Atom::Grid> (new Molib::Atom::Grid (__receptor->get_atoms()));

                return 1;
        } catch ( std::exception &e ) {
                __error_string = std::string("Error in loading receptor ") + e.what();
                return 0;
        }
}

size_t receptor_atom_count() {
        if ( __receptor == nullptr ) {
                __error_string = std::string("You must run initialize_receptor first");
                return 0;
        }

        return __receptor->element(0).get_atoms().size();
}

size_t receptor_atoms(size_t* idx, float* pos) {
        if ( __receptor == nullptr ) {
                __error_string = std::string("You must run initialize_receptor first");
                return 0;
        }
        
        try {
                Molib::Atom::Vec atoms = __receptor->element(0).get_atoms();
                for ( size_t i=0; i < atoms.size(); ++i ) {
                        idx[ i ] = atoms[i]->atom_number();
                        
                        Geom3D::Point& crd = atoms[i]->crd();
                        pos[ i * 3 + 0] = crd.x();
                        pos[ i * 3 + 1] = crd.y();
                        pos[ i * 3 + 2] = crd.z();
                }

                return atoms.size();
        } catch( std::exception &e ) {
                __error_string = std::string("Error creating receptor atom arrays: ") + e.what();
                return 0;
        }
}

size_t receptor_bond_count() {
        if ( __receptor == nullptr ) {
                __error_string = std::string("You must run initialize_receptor first");
                return 0;
        }

        try {
                Molib::BondSet bondset = Molib::get_bonds_in(__receptor->get_atoms());
                return bondset.size();
        } catch (std::exception &e) {
                __error_string = std::string("Error creating receptor bond arrays: ") + e.what();
                return 0;
        }
}

size_t receptor_bonds( size_t* bonds ) {
        if ( __receptor == nullptr ) {
                __error_string = std::string("You must run initialize_receptor first");
                return 0;
        }

        try {
                Molib::BondSet bondset = Molib::get_bonds_in(__receptor->get_atoms());

                size_t i = 0;

                for ( const auto a : bondset) {
                        bonds[ i * 3 + 0 ] = a->atom1().atom_number();
                        bonds[ i * 3 + 1 ] = a->atom2().atom_number();

                        bonds[ i * 3 + 2 ]  = 0;
                        bonds[ i * 3 + 2 ] |= a->is_single()   ? SINGLE_BOND : 0;
                        bonds[ i * 3 + 2 ] |= a->is_double()   ? DOUBLE_BOND : 0;
                        bonds[ i * 3 + 2 ] |= a->is_triple()   ? TRIPLE_BOND : 0;
                        bonds[ i * 3 + 2 ] |= a->is_ring()     ? INRING_BOND : 0;
                        bonds[ i * 3 + 2 ] |= a->is_rotatable()? ROTATE_BOND : 0;
                        bonds[ i * 3 + 2 ] |= a->is_aromatic() ? AROMAT_BOND : 0;

                        ++i;
                }

                return i;
        } catch (std::exception &e) {
                __error_string = std::string("Error creating receptor bond arrays: ") + e.what();
                return 0;
        }
}

size_t initialize_ligand(const char* filename) {
        try {
                Parser::FileParser rpdb (filename, Parser::first_model, 1);

                __ligand = std::unique_ptr<Molib::Molecules> (new Molib::Molecules);

                rpdb.parse_molecule (*__ligand);

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

                return 1;

        } catch ( std::exception &e ) {
                __error_string = std::string("Error in loading ligand ") + e.what();
                return 0;
        }
}

size_t ligand_atom_count() {
        if ( __ligand == nullptr ) {
                __error_string = std::string("You must run initialize_ligand first");
                return 0;
        }

        return __ligand->element(0).get_atoms().size();
}

size_t ligand_atoms(size_t* idx, float* pos) {
        if ( __ligand == nullptr ) {
                __error_string = std::string("You must run initialize_ligand first");
                return 0;
        }
        
        try {
                Molib::Atom::Vec atoms = __ligand->element(0).get_atoms();
                for ( size_t i=0; i < atoms.size(); ++i ) {
                        idx[ i ] = atoms[i]->atom_number();
                        
                        Geom3D::Point& crd = atoms[i]->crd();
                        pos[ i * 3 + 0] = crd.x();
                        pos[ i * 3 + 1] = crd.y();
                        pos[ i * 3 + 2] = crd.z();
                }

                return atoms.size();
        } catch( std::exception &e ) {
                __error_string = std::string("Error creating atom arrays");
                return 0;
        }
}

size_t ligand_bond_count() {
        if ( __ligand == nullptr ) {
                __error_string = std::string("You must run initialize_ligand first");
                return 0;
        }

        try {
                Molib::BondSet bondset = Molib::get_bonds_in(__ligand->get_atoms());
                return bondset.size();
        } catch (std::exception &e) {
                __error_string = std::string("Error creating ligand bond arrays: ") + e.what();
                return 0;
        }
}

size_t ligand_bonds( size_t* bonds ) {
        if ( __ligand == nullptr ) {
                __error_string = std::string("You must run initialize_ligand first");
                return 0;
        }

        try {
                Molib::BondSet bondset = Molib::get_bonds_in(__ligand->get_atoms());

                size_t i = 0;

                for ( const auto a : bondset) {
                        bonds[ i * 3 + 0 ] = a->atom1().atom_number();
                        bonds[ i * 3 + 1 ] = a->atom2().atom_number();

                        bonds[ i * 3 + 2 ]  = 0;
                        bonds[ i * 3 + 2 ] |= a->is_single()   ? SINGLE_BOND : 0;
                        bonds[ i * 3 + 2 ] |= a->is_double()   ? DOUBLE_BOND : 0;
                        bonds[ i * 3 + 2 ] |= a->is_triple()   ? TRIPLE_BOND : 0;
                        bonds[ i * 3 + 2 ] |= a->is_ring()     ? INRING_BOND : 0;
                        bonds[ i * 3 + 2 ] |= a->is_rotatable()? ROTATE_BOND : 0;

                        bonds[ i * 3 + 2 ] |= a->is_aromatic() ? AROMAT_BOND : 0;

                        ++i;
                }
        
                return i;
        } catch (std::exception &e) {
                __error_string = std::string("Error creating ligand bond arrays: ") + e.what();
                return 0;
        }
}

size_t initialize_scoring(const char* obj_dir) {

        if ( __ligand == nullptr  || __receptor == nullptr ) {
                __error_string = std::string("You must run initialize_ligand and initialize_receptor first");
                return 0;
        }

        try {
                __score = std::unique_ptr<Molib::Score> (
                                  new  Molib::Score ("mean", "reduced", "radial", 6.0, 0.01));

                __score->define_composition (__receptor->get_idatm_types(),
                                             __ligand->get_idatm_types())
                .process_distributions_file ("data/csd_complete_distance_distributions.txt")
                .compile_scoring_function()
                .parse_objective_function (obj_dir, 10.0, 1501);

                return 1;

        } catch ( std::exception &e) {
                __error_string = std::string("Error in creating scoring function: ") + e.what();
                return 0;
        }
}

size_t initialize_plugins(const char* plugin_dir) {

        try {

                OMMIface::SystemTopology::loadPlugins(plugin_dir);

                return 1;

        } catch( std::exception& e ) {
                __error_string = std::string("Error in creating KBPlugin: ") + e.what();
                return 0;
        }
}

size_t initialize_ffield(const char* data_dir) {

        if ( __score == nullptr ) {
                __error_string = std::string("You must run initialize_score first");
                return 0;
        }

        try {

                boost::filesystem::path p (data_dir);

                __ffield = std::unique_ptr<OMMIface::ForceField> (new OMMIface::ForceField);

                __ffield->parse_gaff_dat_file ( (p / "gaff.dat").string())
                .add_kb_forcefield (*__score, 0.01, 6.0)
                .parse_forcefield_file ( (p / "amber10.xml").string())
                .parse_forcefield_file ( (p / "tip3p.xml").string());

                __receptor->element (0).prepare_for_mm (*__ffield, *__gridrec);

                __ffield->insert_topology(__receptor->element (0));
                __ffield->insert_topology(__ligand->element(0));

                return 1;

        } catch( std::exception& e ) {
                __error_string = std::string("Error in creating forcefield: ") + e.what();
                return 0;
        }
}

float calculate_score() {
        if ( __ffield == nullptr ) {
                __error_string = std::string("You must run initialize_ffield first");
                return 0;
        }

        try {
                return __score->non_bonded_energy(*__gridrec, (*__ligand)[0]);
        } catch( std::exception& e) {
                __error_string = std::string("Error in scoring: ") + e.what();
                return 0;
        }
}

size_t set_positions_ligand(const size_t* atoms, const float* positions, size_t size) {

        if ( __ligand == nullptr ) {
                __error_string = std::string("You must run initialize_ligand first");
                return 0;
        }

        try { 
        
                Molib::Residue *residue = __ligand->element (0).get_residues().at (0);

                for (size_t i = 0; i < size; ++i) {
                        residue->element (atoms[i]).set_crd (Geom3D::Point (positions[i * 3 + 0],
                                                                            positions[i * 3 + 1],
                                                                            positions[i * 3 + 2]
                                                                           ));
                }

                return 1;
        
        } catch ( std::exception &e ) {
                __error_string = std::string("Error in setting ligand coordinates: ") + e.what();
                return 0;
        }
}

size_t set_positions_receptor(const size_t* atoms, const float* positions, size_t size) {

        if ( __receptor == nullptr ) {
                __error_string = std::string("You must run initialize_receptor first");
                return 0;
        }

        try {

                size_t resi = 0;
                Molib::Residue *residue = __receptor->element (0).get_residues().at (resi++);

                for (size_t i = 0; i < size; ++i) {

                        while (! residue->has_element (atoms[i])) {
                                residue = __receptor->element (0).get_residues().at (resi++);
                        }

                        if (! residue->has_element (atoms[i])) {
                                __error_string = std::string("Error: could not atom = ") + std::to_string(atoms[i]);
                                break;
                        }

                        residue->element (atoms[i]).set_crd (Geom3D::Point (positions[i * 3 + 0],
                                                                            positions[i * 3 + 1],
                                                                            positions[i * 3 + 2]
                                                                           ));
                }

                __gridrec.reset (new Molib::Atom::Grid (__receptor->get_atoms()));
                
                return 1;
                
        } catch (std::exception &e) {
                __error_string = std::string("Error in setting receptor coordinates: ") + e.what();
                return 0;
        }
}

size_t minimize_complex (size_t max_iter, size_t update_freq) {

        if ( __ffield == nullptr ) {
                __error_string = std::string("You must run initialize_ffield first");
                return 0;
        }

        try {

                OMMIface::Modeler modeler (*__ffield, "kb", 6,
                                           0.0001, max_iter,
                                           update_freq, 0.00000000001, false, 2.0, 300, 91);

                Molib::Atom::Vec rec_atoms = __receptor->get_atoms();
                Molib::Atom::Vec lig_atoms = __ligand->get_atoms();

                modeler.add_topology (rec_atoms);
                modeler.add_topology (lig_atoms);

                modeler.init_openmm();

                modeler.add_crds (rec_atoms, __receptor->get_crds());
                modeler.add_crds (lig_atoms, __ligand->get_crds());

                modeler.init_openmm_positions();

                modeler.minimize_state ();

                // init with minimized coordinates
                Geom3D::Point::Vec rec_coords = modeler.get_state (rec_atoms);

                for (size_t i = 0; i < rec_atoms.size(); ++i)
                        rec_atoms[i]->set_crd (rec_coords[i]);

                Geom3D::Point::Vec lig_coords = modeler.get_state (lig_atoms);

                for (size_t j = 0; j < lig_atoms.size(); ++j)
                        lig_atoms[j]->set_crd (lig_coords[j]);

                __gridrec.reset (new Molib::Atom::Grid (__receptor->get_atoms()));

                return 1;

        } catch (std::exception &e) {
                __error_string = std::string("Error in setting receptor coordinates: ") + e.what();
                return 0;
        }
}
