#include "forcefield.hpp"
#include "molib/molecule.hpp"
#include "score/score.hpp"
#include "helper/inout.hpp"
#include "helper/error.hpp"
#include "rapidxml/rapidxml.hpp"
#include "rapidxml/rapidxml_print.hpp"
#include <boost/regex.hpp>
#include <boost/algorithm/string.hpp>
#include <boost/filesystem.hpp>
#include <boost/date_time/posix_time/posix_time.hpp>
#include <boost/algorithm/string/split.hpp>
#include <boost/algorithm/string.hpp>
#include <boost/asio/ip/host_name.hpp>
#include <openmm/Units.h>
#include <string>
#include <sstream>
using namespace std;
//
namespace OMMIface {
	ostream& operator<< (ostream& stream, const ForceField::ResidueTopology& r) {
		for (auto &kv : r.atom) {
			stream << "atom " << kv.first << " = " << kv.second << endl;
		}
		for (auto &kv1 : r.bond) {
			for (auto &kv2 : kv1.second) {
				stream << "bond " << kv1.first << " " << kv2.first << " = " << kv2.second << endl;
			}
		}
		for (auto &from : r.external_bond) {
			stream << "external bond " << from << endl;
		}
		return stream;
	}
	ostream& operator<<(ostream& os, const ForceField& ff) {
		for (auto &kv : ff.residue_topology) {
			os << "Residue: " << kv.first << endl << kv.second << endl;
		}
		return os;
	}
	ostream& operator<<(ostream& os, const ForceField::KBForces &kb) {
		for (auto &kv1 : kb) {
			const int aclass1 = kv1.first;
			for (auto &kv2 : kv1.second) {
				const int aclass2 = kv2.first;
				const auto &pot = kv2.second.potential;
				const auto &der = kv2.second.derivative;
				os << "KB forcefield for " << help::idatm_unmask[aclass1] 
					<< " and " << help::idatm_unmask[aclass2] << " : " << endl; 
				for (size_t i = 0; i < pot.size(); ++i)
					os << "i = " << i << " pot = " << pot[i] << " der = "
						<< der[i] << endl;
			}
		}
		return os;
	}

        const ForceField::KBType& ForceField::get_kb_force_type(const Molib::Atom &atom1, 
                const Molib::Atom &atom2) const {
                const int &aclass1 = atom1.idatm_type(); 
                const int &aclass2 = atom2.idatm_type();
#ifndef NDEBUG
                dbgmsg("add kb force"
                        << " aclass1 = " << help::idatm_unmask[aclass1] 
                        << " aclass2 = " << help::idatm_unmask[aclass2]);
#endif
                if ( kb_force_type.count(aclass1) != 0 && kb_force_type.at(aclass1).count(aclass2) != 0) {
                        return kb_force_type.at(aclass1).at(aclass2);
                }

                if ( kb_force_type.count(aclass2) != 0 && kb_force_type.at(aclass2).count(aclass1) != 0) {
                        return kb_force_type.at(aclass2).at(aclass1);
                }

                stringstream ss;
                ss << "warning : missing kb force type " << help::idatm_unmask[aclass1] 
                   << "-" << help::idatm_unmask[aclass2];
                throw ParameterError(ss.str());

        }

        const ForceField::KBType& ForceField::get_kb_force_type(int aclass1, int aclass2) const {
                if ( kb_force_type.count(aclass1) != 0 && kb_force_type.at(aclass1).count(aclass2) != 0) {
                        return kb_force_type.at(aclass1).at(aclass2);
                }

                if ( kb_force_type.count(aclass2) != 0 && kb_force_type.at(aclass2).count(aclass1) != 0) {
                        return kb_force_type.at(aclass2).at(aclass1);
                }

                stringstream ss;
                ss << "warning : missing kb force type " << help::idatm_unmask[aclass1] 
                   << "-" << help::idatm_unmask[aclass2];
                throw ParameterError(ss.str());

        }

        bool ForceField::has_atom_type(const int aclass1) const {
                return (atom_type.count(aclass1) != 0);
        }


        const ForceField::AtomType& ForceField::get_atom_type(const int type) const {
                if (has_atom_type(type) != 0)
                        return atom_type.at(type);
                throw ParameterError("warning : missing atom type " 
                       + std::to_string(type));
        }

        bool ForceField::has_bond_type(const string& aclass1, const string& aclass2) const {
                return ( bond_type.count(aclass1) != 0 &&
                         bond_type.at(aclass1).count(aclass2) != 0 );
        }

        const ForceField::BondType& ForceField::get_bond_type(const int type1, const int type2) const {

                const AtomType &atype1 = atom_type.at(type1);
                dbgmsg(atype1.cl);
                const AtomType &atype2 = atom_type.at(type2);
                dbgmsg(atype2.cl);
                const string &aclass1 = atype1.cl;
                const string &aclass2 = atype2.cl;
                dbgmsg("add bond type1 = " << type1 << " type2 = " << type2 << " aclass1 = " 
                       << aclass1 << " aclass2 = " << aclass2);
                
                if ( has_bond_type(aclass1, aclass2) ) {
                        return bond_type.at(aclass1).at(aclass2);
                }

                if ( has_bond_type(aclass2, aclass1) ) {
                        return bond_type.at(aclass2).at(aclass1);
                }

                // if this fails, we try finding a similar bond parameter
                for (auto &sclass : help::get_replacement({aclass1, aclass2})) {
                    
                        if ( has_bond_type(sclass[0], sclass[1]) ) {
                            return bond_type.at(sclass[0]).at(sclass[1]);
                        }

                        if ( has_bond_type(sclass[1], sclass[0]) ) {
                            return bond_type.at(sclass[1]).at(sclass[0]);
                        }
                }
                throw ParameterError("warning : missing bond type " + aclass1 
                        + "-" + aclass2);
        }

        bool ForceField::has_angle_type( const string &aclass1, const string &aclass2, const string &aclass3 ) const {
                return ( angle_type.count(aclass1) !=0 &&
                         angle_type.at(aclass1).count(aclass2) != 0 &&
                         angle_type.at(aclass1).at(aclass2).count(aclass3) );
        }

        const ForceField::AngleType& ForceField::get_angle_type(const int type1, 
                const int type2, const int type3) const {

                const AtomType &atype1 = atom_type.at(type1);
                const AtomType &atype2 = atom_type.at(type2);
                const AtomType &atype3 = atom_type.at(type3);
                const string &aclass1 = atype1.cl;
                const string &aclass2 = atype2.cl;
                const string &aclass3 = atype3.cl;
                dbgmsg("add angle type1 = " << type1 << " type2 = " << type2 
                        << " type3 = " << type3 << " aclass1 = " << aclass1 
                        << " aclass2 = " << aclass2 << " aclass3 = " << aclass3);
                // See note under bond stretch above regarding the factor of 2 here.
                if ( has_angle_type(aclass1, aclass2, aclass3) ) return angle_type.at(aclass1).at(aclass2).at(aclass3);
                if ( has_angle_type(aclass3, aclass2, aclass1) ) return angle_type.at(aclass3).at(aclass2).at(aclass1);
                // if this fails, we try finding a similar angle parameter
                for (auto &sclass : help::get_replacement({aclass1, aclass2, aclass3})) {
                        if ( has_angle_type(sclass[0], sclass[1], sclass[2]) )
                            return angle_type.at(sclass[0]).at(sclass[1]).at(sclass[2]);
                        if ( has_angle_type(sclass[2], sclass[1], sclass[0]) )
                            return angle_type.at(sclass[2]).at(sclass[1]).at(sclass[0]);
                }
                throw ParameterError("warning : missing angle type " + aclass1 
                                     + "-" + aclass2 + "-" + aclass3);
        }

        bool ForceField::has_dihedral_type( const string &aclass1, const string &aclass2,
                                            const string &aclass3, const string &aclass4 ) const {
                return ( torsion_type.count(aclass1) !=0 &&
                         torsion_type.at(aclass1).count(aclass2) != 0 &&
                         torsion_type.at(aclass1).at(aclass2).count(aclass3) !=0 &&
                         torsion_type.at(aclass1).at(aclass2).at(aclass3).count(aclass4) != 0 );
        }

        const ForceField::TorsionTypeVec& ForceField::get_dihedral_type(const int type1, 
                const int type2, const int type3, const int type4) const {

                const AtomType &atype1 = atom_type.at(type1);
                const AtomType &atype2 = atom_type.at(type2);
                const AtomType &atype3 = atom_type.at(type3);
                const AtomType &atype4 = atom_type.at(type4);
                const string &aclass1 = atype1.cl;
                const string &aclass2 = atype2.cl;
                const string &aclass3 = atype3.cl;
                const string &aclass4 = atype4.cl;

                dbgmsg("add dihedral type1 = " << type1 << " type2 = " << type2 
                        << " type3 = " << type3 << " type4 = " << type4 << " aclass1 = " << aclass1 
                        << " aclass2 = " << aclass2 << " aclass3 = " << aclass3 << " aclass4 = " << aclass4);

                if ( has_dihedral_type(aclass1, aclass2, aclass3, aclass4) )
                        return torsion_type.at(aclass1).at(aclass2).at(aclass3).at(aclass4);

                if ( has_dihedral_type(aclass4, aclass3, aclass2, aclass1) )
                        return torsion_type.at(aclass4).at(aclass3).at(aclass2).at(aclass1);

                if ( has_dihedral_type("X",aclass2,aclass3,"X") )
                        return torsion_type.at("X").at(aclass2).at(aclass3).at("X");

                if ( has_dihedral_type("X",aclass3,aclass2,"X") )
                        return torsion_type.at("X").at(aclass3).at(aclass2).at("X");

                // if this fails, we try finding a similar dihedral parameter
                for (auto &sclass : help::get_replacement({aclass1, aclass2, aclass3, aclass4})) {

                        if ( has_dihedral_type(sclass[0], sclass[1], sclass[2], sclass[3]) )
                                return torsion_type.at(sclass[0]).at(sclass[1]).at(sclass[2]).at(sclass[3]);

                        if ( has_dihedral_type(sclass[3], sclass[2], sclass[1], sclass[0]) )
                                return torsion_type.at(sclass[3]).at(sclass[2]).at(sclass[1]).at(sclass[0]);

                        if ( has_dihedral_type("X", sclass[1], sclass[2], "X") )
                                return torsion_type.at("X").at(sclass[1]).at(sclass[2]).at("X");

                        if ( has_dihedral_type("X", sclass[2], sclass[1], "X") )
                                return torsion_type.at("X").at(sclass[2]).at(sclass[1]).at("X");
                }
                throw ParameterError("warning : missing dihedral type " + aclass1 
                                     + "-" + aclass2 + "-" + aclass3 + "-" + aclass4);
        }

        bool ForceField::has_improper_type( const string &aclass1, const string &aclass2,
                                            const string &aclass3, const string &aclass4 ) const {
                return ( improper_type.count(aclass1) !=0 &&
                         improper_type.at(aclass1).count(aclass2) != 0 &&
                         improper_type.at(aclass1).at(aclass2).count(aclass3) !=0 &&
                         improper_type.at(aclass1).at(aclass2).at(aclass3).count(aclass4) != 0 );
        }

        const ForceField::TorsionTypeVec& ForceField::get_improper_type(const int type1, 
                const int type2, const int type3, const int type4) const {

                const AtomType &atype1 = atom_type.at(type1);
                const AtomType &atype2 = atom_type.at(type2);
                const AtomType &atype3 = atom_type.at(type3);
                const AtomType &atype4 = atom_type.at(type4);
                const string &aclass1 = atype1.cl;
                const string &aclass2 = atype2.cl;
                const string &aclass3 = atype3.cl;
                const string &aclass4 = atype4.cl;

                dbgmsg("add improper (third atom is central) type1 = " << type1 << " type2 = " << type2 
                        << " type3 = " << type3 << " type4 = " << type4 << " aclass1 = " << aclass1 
                        << " aclass2 = " << aclass2 << " aclass3 = " << aclass3 << " aclass4 = " << aclass4);

                if ( has_improper_type( aclass1, aclass2, aclass3, aclass4 ) )
                        return improper_type.at(aclass1).at(aclass2).at(aclass3).at(aclass4);

                if ( has_improper_type( aclass2, aclass1, aclass3, aclass4 ) )
                        return improper_type.at(aclass2).at(aclass1).at(aclass3).at(aclass4);

                if ( has_improper_type( aclass2, aclass4, aclass3, aclass1 ) )
                        return improper_type.at(aclass2).at(aclass4).at(aclass3).at(aclass1);

                if ( has_improper_type( aclass4, aclass2, aclass3, aclass1 ) )
                        return improper_type.at(aclass4).at(aclass2).at(aclass3).at(aclass1);

                if ( has_improper_type( aclass4, aclass1, aclass3, aclass2 ) )
                        return improper_type.at(aclass4).at(aclass1).at(aclass3).at(aclass2);

                if ( has_improper_type( aclass1, aclass4, aclass3, aclass2 ) )
                        return improper_type.at(aclass1).at(aclass4).at(aclass3).at(aclass2);

                // Missing one atom
                if ( has_improper_type( "X", aclass1, aclass3, aclass2 ) )
                        return improper_type.at("X").at(aclass1).at(aclass3).at(aclass2);

                if ( has_improper_type( "X", aclass1, aclass3, aclass4 ) )
                        return improper_type.at("X").at(aclass1).at(aclass3).at(aclass4);

                if ( has_improper_type( "X", aclass2, aclass3, aclass1 ) )
                        return improper_type.at("X").at(aclass2).at(aclass3).at(aclass1);

                if ( has_improper_type( "X", aclass2, aclass3, aclass4 ) )
                        return improper_type.at("X").at(aclass2).at(aclass3).at(aclass4);

                if ( has_improper_type( "X", aclass4, aclass3, aclass1 ) )
                        return improper_type.at("X").at(aclass4).at(aclass3).at(aclass1);

                if ( has_improper_type( "X", aclass4, aclass3, aclass2 ) )
                        return improper_type.at("X").at(aclass4).at(aclass3).at(aclass2);

                // Missing two atoms
                if ( has_improper_type( "X", "X", aclass3, aclass1 ) )
                        return improper_type.at("X").at("X").at(aclass3).at(aclass1);

                if ( has_improper_type( "X", "X", aclass3, aclass2 ) )
                        return improper_type.at("X").at("X").at(aclass3).at(aclass2);

                if ( has_improper_type( "X", "X", aclass3, aclass4 ) )
                        return improper_type.at("X").at("X").at(aclass3).at(aclass4);

                throw ParameterError("warning : missing improper type " + aclass1 
                                     + "-" + aclass2 + "-" + aclass3 + "-" + aclass4);
        }

	ForceField& ForceField::insert_topology(const Molib::Molecule &molecule) {
		map<const string, const int> atom_name_to_type;
		for (auto &kv : this->atom_type) {
			const int &type = kv.first;
			AtomType &at = kv.second;
			atom_name_to_type.insert({at.cl, type});
		}
		for (auto &presidue : molecule.get_residues()) {
			auto &residue = *presidue;
			// additionally check if residue not already in the ffield (this is for adding receptor
			// topology where amino acids are already in the ffield, but cofactors are not
			if (!this->residue_exists(residue.resn())) {

				if ( help::ions.count(residue.resn()) != 0 ) {
					this->residue_topology.insert({residue.resn(), ResidueTopology()});
					ResidueTopology &rtop = this->residue_topology.at(residue.resn());
					rtop.atom.insert({residue.first().atom_name(), 1});
					continue;
				}

				this->residue_topology.insert({residue.resn(), ResidueTopology()});
				ResidueTopology &rtop = this->residue_topology.at(residue.resn());
				for (auto &atom : residue) {
					if (!atom_name_to_type.count(atom.gaff_type())) { // see issue #113
						throw Error("die : insert topology failed for atom with gaff type (check gaff.dat) = " 
							+ atom.gaff_type());
					}
					rtop.atom.insert({atom.atom_name(), atom_name_to_type.at(atom.gaff_type())});
				}
				for (auto &atom : residue) {
					for (auto &adj_a : atom) {
						if (atom.atom_number() < adj_a.atom_number())
							rtop.bond[atom.atom_name()][adj_a.atom_name()] = true;
					}
				}
				dbgmsg("inserted topology for residue " << residue.resn()
					<< endl << rtop);
			}
		}
		return *this;
	}

	ForceField& ForceField::erase_topology(const Molib::Molecule &molecule) {
		for (auto &presidue : molecule.get_residues()) {
			auto &residue = *presidue;
			dbgmsg("erasing topology for residue " << residue.resn()
				<< endl << residue_topology.at(residue.resn()));
			this->residue_topology.erase(residue.resn());
		}
		return *this;
	}

	ForceField& ForceField::add_kb_forcefield(const Molib::Score &score, 
		const double step, const double cutoff) {
                this->cutoff = cutoff * OpenMM::NmPerAngstrom;
		this->step   = step * OpenMM::NmPerAngstrom;
		for (auto &kv : score.get_energies()) {
			auto &atom_pair = kv.first;
			auto &energy = kv.second;
			this->kb_force_type[atom_pair.first][atom_pair.second].potential = energy;
		}
		for (auto &kv : score.get_derivatives()) {
			auto &atom_pair = kv.first;
			auto &der = kv.second;
			this->kb_force_type[atom_pair.first][atom_pair.second].derivative = der;
		}
		return *this;
	}

	ForceField& ForceField::parse_gaff_dat_file(const string &gaff_dat_file) {
		boost::smatch m;
		vector<string> gdf;
		dbgmsg("reading gaff file = " << gaff_dat_file);
		Inout::read_file(gaff_dat_file, gdf);
		bool non_bonded = false;
		int semaphore = 0;
		int type = 10000; // start ligand types with 10000 to NOT overlap with receptor
		map<const string, const int> atom_name_to_type;
		for (auto &line : gdf) {
			if (line.compare(0, 5, "AMBER") == 0) {
				continue;
			}
			if (line.empty()) {
				semaphore++;
				continue;
			}
			if (line.compare(0, 4, "MOD4") == 0) {
				non_bonded = true;
				continue;
			}
			if (line.compare(0, 3, "END") == 0)
				break;
			// read atom types
			if (semaphore == 0 && boost::regex_search(line, m, boost::regex("^(\\S{1,2})\\s+(\\S+)\\s+(\\S+)"))) {
				if (m[1].matched && m[2].matched && m[3].matched) {
					AtomType &at = this->atom_type[type];
					at.cl = m[1].str();
					at.mass = stod(m[2].str());
					at.element = "";
					//~ at.polarizability = stod(m[3].str()); // not used ...
					atom_name_to_type.insert({at.cl, type});
					dbgmsg("type = " << type << " class = " << at.cl << " at.element = " 
							<< at.element << " at.mass = " << at.mass);
					type++;
				}
			}
			// parse harmonic bonds
			if (semaphore == 1 && boost::regex_search(line, m, boost::regex("^(\\S{1,2})\\s*-(\\S{1,2})\\s+(\\S+)\\s+(\\S+)"))) {
				if (m[1].matched && m[2].matched && m[3].matched && m[4].matched) {
					const string cl1 = m[1].str();
					const string cl2 = m[2].str();
					const double k = stod(m[3].str()) * 2 * OpenMM::KJPerKcal 
                                    * OpenMM::AngstromsPerNm * OpenMM::AngstromsPerNm;
					const double length = stod(m[4].str()) * OpenMM::NmPerAngstrom;
					const bool can_constrain = (cl1.at(0) == 'h' || cl2.at(0) == 'h' ? true : false); // X-H bonds can be constrained
					this->bond_type[cl1][cl2] = BondType{length, k, can_constrain};
					dbgmsg("harmonic bond force between " << cl1 << " and " << cl2 << " is "
							<< " length = " << length
							<< " k = " << k
							<< " can_constrain = " << can_constrain);
				}
			}
			// parse harmonic angles
			if (semaphore == 2 && boost::regex_search(line, m, boost::regex("^(\\S{1,2})\\s*-(\\S{1,2})\\s*-(\\S{1,2})\\s+(\\S+)\\s+(\\S+)"))) {
				if (m[1].matched && m[2].matched && m[3].matched && m[4].matched && m[5].matched) {
					const string cl1 = m[1].str();
					const string cl2 = m[2].str();
					const string cl3 = m[3].str();
					const double k = stod(m[4].str()) * 2 * OpenMM::KJPerKcal;
					const double angle = stod(m[5].str()) * OpenMM::RadiansPerDegree;
					this->angle_type[cl1][cl2][cl3] = AngleType{angle, k};
					dbgmsg("harmonic angle force between " << cl1 << " and " << cl2 << " and " << cl3 <<" is "
							<< " angle = " << angle
							<< " k = " << k);
				}
			}
			// parse proper torsions
			if (semaphore == 3 && boost::regex_search(line, m, boost::regex("^(\\S{1,2})\\s*-(\\S{1,2})\\s*-(\\S{1,2})\\s*-(\\S{1,2})\\s+(\\S+)\\s+(\\S+)\\s+(\\S+)\\s+(\\S+)"))) {
				if (m[1].matched && m[2].matched && m[3].matched && m[4].matched 
					&& m[5].matched && m[6].matched && m[7].matched && m[8].matched) {
					const string cl1 = m[1].str();
					const string cl2 = m[2].str();
					const string cl3 = m[3].str();
					const string cl4 = m[4].str();
					const double idivf = stod(m[5].str());
					const double k = (1 / idivf) * stod(m[6].str()) * OpenMM::KJPerKcal;
					const double phase = stod(m[7].str())  * OpenMM::RadiansPerDegree;
					const int periodicity = abs(stoi(m[8].str()));
					TorsionTypeVec &dihedral = this->torsion_type[cl1][cl2][cl3][cl4];
					dihedral.push_back(TorsionType{periodicity, phase, k});
					dbgmsg("periodic (proper) torsion force between " 
							<< cl1 << " and " << cl2 << " and " << cl3 << " and " << cl4 << " is "
							<< " periodicity" << dihedral.size() << " = " << periodicity
							<< " phase" << dihedral.size() << " = " << phase
							<< " k" << dihedral.size() << " = " << k);
				}
			}
			// parse improper torsions (there is no IDIVF)
			if (semaphore == 4 && boost::regex_search(line, m, boost::regex("^(\\S{1,2})\\s*-(\\S{1,2})\\s*-(\\S{1,2})\\s*-(\\S{1,2})\\s+(\\S+)\\s+(\\S+)\\s+(\\S+)"))) {
				if (m[1].matched && m[2].matched && m[3].matched && m[4].matched 
					&& m[5].matched && m[6].matched && m[7].matched) {
					const string cl1 = m[1].str();
					const string cl2 = m[2].str();
					const string cl3 = m[3].str();
					const string cl4 = m[4].str();
					const double k = stod(m[5].str()) * OpenMM::KJPerKcal;
					const double phase = stod(m[6].str())  * OpenMM::RadiansPerDegree;
					const int periodicity = abs(stoi(m[7].str()));
					TorsionTypeVec &improper = this->improper_type[cl1][cl2][cl3][cl4];
					improper.push_back(TorsionType{periodicity, phase, k});
					dbgmsg("periodic (improper) torsion force between (third atom is central) " 
							<< cl1 << " and " << cl2 << " and " << cl3 << " and " << cl4 << " is "
							<< " periodicity" << improper.size() << " = " << periodicity
							<< " phase" << improper.size() << " = " << phase
							<< " k" << improper.size() << " = " << k);
				}
			}
			if (non_bonded && boost::regex_search(line, m, boost::regex("^\\s{2}(\\S{1,2})\\s+(\\S+)\\s+(\\S+)"))) {
				if (m[1].matched && m[2].matched && m[3].matched) {
					const int type = atom_name_to_type.at(m[1].str());
					AtomType &at = this->atom_type[type];
					at.charge = 0; // for now set charge to 0
					at.sigma = stod(m[2].str()) * OpenMM::NmPerAngstrom 
												* OpenMM::SigmaPerVdwRadius;
					at.epsilon = stod(m[3].str()) * OpenMM::KJPerKcal;
					dbgmsg("nonbonded force for atom type " << type
							<< " charge = " << at.charge
							<< " sigma = " << at.sigma
							<< " epsilon = " << at.epsilon);
				}
			}
		}
		// parse nonbonded force
		this->coulomb14scale = 0.833333;
		this->lj14scale = 0.5;
		return *this;
	}

	template<typename T>
	const char* str(rapidxml::xml_document<> &doc, const T &i) {
		std::ostringstream s;
		s << i;
		return doc.allocate_string(s.str().c_str());
	}

	void ForceField::output_forcefield_file(const string &fn) { // output XML forcefield & topology file
		using namespace rapidxml;
		xml_document<> doc; // character type defaults to char
		xml_node<> *ff_node = doc.allocate_node(node_element, "ForceField");
		doc.append_node(ff_node);
		// atom types to XML
		xml_node<> *atom_type_node = doc.allocate_node(node_element, "AtomTypes");
		ff_node->append_node(atom_type_node);
		for (auto &kv : this->atom_type) {
			const int &name = kv.first;
			const AtomType &at = kv.second;
			xml_node<> *type_node = doc.allocate_node(node_element, "Type");
			atom_type_node->append_node(type_node);
			type_node->append_attribute(doc.allocate_attribute("name", str(doc, name)));
			type_node->append_attribute(doc.allocate_attribute("class", at.cl.c_str()));
			type_node->append_attribute(doc.allocate_attribute("element", at.element.c_str()));
			type_node->append_attribute(doc.allocate_attribute("mass", str(doc, at.mass)));
		}
		// residue topology to XML
		xml_node<> *residues_node = doc.allocate_node(node_element, "Residues");
		ff_node->append_node(residues_node);
		for (auto &kv : this->residue_topology) {
			const string &name = kv.first;
			const ResidueTopology &rtop = kv.second;
			xml_node<> *residue_node = doc.allocate_node(node_element, "Residue");
			residues_node->append_node(residue_node);
			residue_node->append_attribute(doc.allocate_attribute("name", name.c_str()));
			map<const string, const int> atom_name_to_index;
			int i = 0;
			for (auto &kv2 : rtop.atom) {
				const string &name = kv2.first;
				const int &type = kv2.second;
				atom_name_to_index.insert({name, i++});
				xml_node<> *atom_node = doc.allocate_node(node_element, "Atom");
				residue_node->append_node(atom_node);				
				atom_node->append_attribute(doc.allocate_attribute("name", name.c_str()));
				atom_node->append_attribute(doc.allocate_attribute("type", str(doc, type)));
			}
			for (auto &kv2 : rtop.bond) {
				const int from = atom_name_to_index.at(kv2.first);
				for (auto &kv3 : kv2.second) {
					const int to = atom_name_to_index.at(kv3.first);
					xml_node<> *bond_node = doc.allocate_node(node_element, "Bond");
					residue_node->append_node(bond_node);				
					bond_node->append_attribute(doc.allocate_attribute("from", str(doc, from)));
					bond_node->append_attribute(doc.allocate_attribute("to", str(doc, to)));
				}
			}
		}
		// harmonic bonds to XML
		xml_node<> *hb_node = doc.allocate_node(node_element, "HarmonicBondForce");
		ff_node->append_node(hb_node);
		for (auto &kv : this->bond_type) {
			const string &cl1 = kv.first;
			for (auto &kv2 : kv.second) {
				const string &cl2 = kv2.first;
				const BondType &bt = kv2.second;
				xml_node<> *bond_node = doc.allocate_node(node_element, "Bond");
				hb_node->append_node(bond_node);
				bond_node->append_attribute(doc.allocate_attribute("class1", cl1.c_str()));
				bond_node->append_attribute(doc.allocate_attribute("class2", cl2.c_str()));
				bond_node->append_attribute(doc.allocate_attribute("length", str(doc, bt.length)));
				bond_node->append_attribute(doc.allocate_attribute("k", str(doc, bt.k)));
			}
		}
		// harmonic angles to XML
		xml_node<> *ha_node = doc.allocate_node(node_element, "HarmonicAngleForce");
		ff_node->append_node(ha_node);
		for (auto &kv : this->angle_type) {
			const string &cl1 = kv.first;
			for (auto &kv2 : kv.second) {
				const string &cl2 = kv2.first;
				for (auto &kv3 : kv2.second) {
					const string &cl3 = kv3.first;
					const AngleType &at = kv3.second;
					xml_node<> *angle_node = doc.allocate_node(node_element, "Angle");
					ha_node->append_node(angle_node);
					angle_node->append_attribute(doc.allocate_attribute("class1", cl1.c_str()));
					angle_node->append_attribute(doc.allocate_attribute("class2", cl2.c_str()));
					angle_node->append_attribute(doc.allocate_attribute("class3", cl3.c_str()));
					angle_node->append_attribute(doc.allocate_attribute("angle", str(doc, at.angle)));
					angle_node->append_attribute(doc.allocate_attribute("k", str(doc, at.k)));
				}
			}
		}
		// periodic proper torsions to XML
		xml_node<> *pt_node = doc.allocate_node(node_element, "PeriodicTorsionForce");
		ff_node->append_node(pt_node);
		for (auto &kv : this->torsion_type) {
			const string &cl1 = kv.first;
			for (auto &kv2 : kv.second) {
				const string &cl2 = kv2.first;
				for (auto &kv3 : kv2.second) {
					const string &cl3 = kv3.first;
					for (auto &kv4 : kv3.second) {
						const string &cl4 = kv4.first;
						xml_node<> *torsion_node = doc.allocate_node(node_element, "Proper");
						pt_node->append_node(torsion_node);
						torsion_node->append_attribute(doc.allocate_attribute("class1", cl1.c_str()));
						torsion_node->append_attribute(doc.allocate_attribute("class2", cl2.c_str()));
						torsion_node->append_attribute(doc.allocate_attribute("class3", cl3.c_str()));
						torsion_node->append_attribute(doc.allocate_attribute("class4", cl4.c_str()));
						int i = 0;
						for (auto &tt : kv4.second) {
							const string &num = std::to_string(++i);
							torsion_node->append_attribute(doc.allocate_attribute(str(doc, ("periodicity" + num)), str(doc, tt.periodicity)));
							torsion_node->append_attribute(doc.allocate_attribute(str(doc, "phase" + num), str(doc, tt.phase)));
							torsion_node->append_attribute(doc.allocate_attribute(str(doc, "k" + num), str(doc, tt.k)));
						}
					}
				}
			}
		}
		// periodic improper torsions to XML
		for (auto &kv : this->improper_type) {
			const string &cl1 = kv.first;
			for (auto &kv2 : kv.second) {
				const string &cl2 = kv2.first;
				for (auto &kv3 : kv2.second) {
					const string &cl3 = kv3.first;
					for (auto &kv4 : kv3.second) {
						const string &cl4 = kv4.first;
						xml_node<> *torsion_node = doc.allocate_node(node_element, "Improper");
						pt_node->append_node(torsion_node);
						// be careful because in xml file FIRST atom is central and in gaff dat type
						// (and in improper_type array) the THIRD atom is central!
						torsion_node->append_attribute(doc.allocate_attribute("class1", cl3.c_str()));
						torsion_node->append_attribute(doc.allocate_attribute("class2", cl1.c_str()));
						torsion_node->append_attribute(doc.allocate_attribute("class3", cl2.c_str()));
						torsion_node->append_attribute(doc.allocate_attribute("class4", cl4.c_str()));
						int i = 0;
						for (auto &tt : kv4.second) {
							const string &num = std::to_string(++i);
							torsion_node->append_attribute(doc.allocate_attribute(str(doc, ("periodicity" + num)), str(doc, tt.periodicity)));
							torsion_node->append_attribute(doc.allocate_attribute(str(doc, "phase" + num), str(doc, tt.phase)));
							torsion_node->append_attribute(doc.allocate_attribute(str(doc, "k" + num), str(doc, tt.k)));
						}
					}
				}
			}
		}
		// nonbonded force to XML
		xml_node<> *nb_force_node = doc.allocate_node(node_element, "NonbondedForce");
		ff_node->append_node(nb_force_node);
		nb_force_node->append_attribute(doc.allocate_attribute("coulomb14scale", str(doc, this->coulomb14scale)));
		nb_force_node->append_attribute(doc.allocate_attribute("lj14scale", str(doc, this->lj14scale)));
		for (auto &kv : this->atom_type) {
			const int &type = kv.first;
			const AtomType &at = kv.second;
			xml_node<> *atom_node = doc.allocate_node(node_element, "Atom");
			nb_force_node->append_node(atom_node);
			atom_node->append_attribute(doc.allocate_attribute("type", str(doc, type)));
			atom_node->append_attribute(doc.allocate_attribute("charge", str(doc, at.charge)));
			atom_node->append_attribute(doc.allocate_attribute("sigma", str(doc, at.sigma)));
			atom_node->append_attribute(doc.allocate_attribute("epsilon", str(doc, at.epsilon)));
		}
		// knowledge-based force to XML
		xml_node<> *kb_node = doc.allocate_node(node_element, "KBForce");
		ff_node->append_node(kb_node);
		//~ kb_node->append_attribute(doc.allocate_attribute("step", str(doc, this->step)));
		//~ kb_node->append_attribute(doc.allocate_attribute("scale", str(doc, this->scale)));
		for (auto &kv : this->kb_force_type) {
			const auto &cl1 = kv.first;
			for (auto &kv2 : kv.second) {
				const auto &cl2 = kv2.first;
				const KBType &kbt = kv2.second;
				xml_node<> *kb_bond_node = doc.allocate_node(node_element, "Bond");
				kb_node->append_node(kb_bond_node);
				kb_bond_node->append_attribute(doc.allocate_attribute("class1", help::idatm_unmask[cl1]));
				kb_bond_node->append_attribute(doc.allocate_attribute("class2", help::idatm_unmask[cl2]));
				stringstream ss;
				copy(kbt.potential.begin(), kbt.potential.end(), ostream_iterator<double>(ss, " "));
				kb_bond_node->append_attribute(doc.allocate_attribute("potential", str(doc, ss.str())));
				ss.str("");ss.clear();
				copy(kbt.derivative.begin(), kbt.derivative.end(), ostream_iterator<double>(ss, " "));
				kb_bond_node->append_attribute(doc.allocate_attribute("derivative", str(doc, ss.str())));
			}
		}
		// print XML to file
		stringstream ss;
		ss << doc;
		Inout::file_open_put_stream(fn, ss);
	}
	ForceField& ForceField::parse_forcefield_file(const string &fn) { // read XML forcefield & topology file
		using namespace rapidxml;
		string ff_file;
		Inout::read_file(fn, ff_file);
		char *c_ff_file = new char[ff_file.size() + 1];
		strcpy(c_ff_file, ff_file.c_str());
		xml_document<> doc; // character type defaults to char
		doc.parse<0>(c_ff_file); // 0 means default parse flags
		// parse atom types
		for (xml_node<> *node = doc.first_node("ForceField")->first_node("AtomTypes")->first_node(); node; node = node->next_sibling()) {
			const int name = stoi(node->first_attribute("name")->value());
			if (this->atom_type.count(name))
				throw Error("die : ligands's types (names) overlap with receptor's (see xml file)");
			AtomType &at = this->atom_type[name];
			at.cl = node->first_attribute("class")->value();
			at.element = node->first_attribute("element")->value();
			at.mass = stod(node->first_attribute("mass")->value());
			dbgmsg("name = " << name << " class = " << at.cl << " at.element = " 
					<< at.element << " at.mass = " << at.mass);
		}
		// parse residue topology
		for (xml_node<> *node = doc.first_node("ForceField")->first_node("Residues")->first_node();	node; node = node->next_sibling()) {
			const string name = node->first_attribute("name")->value();
			dbgmsg("parsing residue topology name = " << name);
			ResidueTopology &rtop = this->residue_topology[name];
			vector<string> names;
			dbgmsg("before reading atoms first_node = " << node->first_node("Atom")->name());
			for (xml_node<> *node2 = node->first_node("Atom"); node2; node2 = node2->next_sibling("Atom")) {
				dbgmsg("parsing residue topology atom name = " << node2->first_attribute("name")->value()
						<< " type = " << stoi(node2->first_attribute("type")->value()));
				rtop.atom[node2->first_attribute("name")->value()] = stoi(node2->first_attribute("type")->value());
				names.push_back(node2->first_attribute("name")->value());
			}
			dbgmsg("after reading atoms");
			for (xml_node<> *node2 = node->first_node("Bond"); node2; node2 = node2->next_sibling("Bond")) {
				const string atom_name1 = names[stoi(node2->first_attribute("from")->value())];
				const string atom_name2 = names[stoi(node2->first_attribute("to")->value())];
				dbgmsg("parsing residue topology bond atom_name1 = " << atom_name1
						<< " atom_name2 = " << atom_name2);
				rtop.bond[atom_name1][atom_name2] = true;
			}
			for (xml_node<> *node2 = node->first_node("ExternalBond"); node2; node2 = node2->next_sibling("ExternalBond")) {
				dbgmsg("parsing residue topology external_bond from = " 
						<< names[stoi(node2->first_attribute("from")->value())]);
				rtop.external_bond.insert( names[stoi(node2->first_attribute("from")->value())] );
			}
			dbgmsg("name = " << node->name());
		}
		// parse harmonic bonds
		for (xml_node<> *node = doc.first_node("ForceField")->first_node("HarmonicBondForce")->first_node(); node; node = node->next_sibling()) {
			const string cl1 = node->first_attribute("class1")->value();
			const string cl2 = node->first_attribute("class2")->value();
			const bool can_constrain = (cl1.at(0) == 'H' || cl2.at(0) == 'H' ? true : false); // X-H bonds can be constrained
			this->bond_type[cl1][cl2] = BondType{stod(node->first_attribute("length")->value()),
											stod(node->first_attribute("k")->value()),
											can_constrain};
			dbgmsg("harmonic bond force between " << cl1 << " and " << cl2 << " is "
					<< " length = " << stod(node->first_attribute("length")->value()) 
					<< " k = " << stod(node->first_attribute("k")->value())
					<< " can_constrain = " << can_constrain);
		}
		// parse harmonic angles
		for (xml_node<> *node = doc.first_node("ForceField")->first_node("HarmonicAngleForce")->first_node(); node; node = node->next_sibling()) {
			const string cl1 = node->first_attribute("class1")->value();
			const string cl2 = node->first_attribute("class2")->value();
			const string cl3 = node->first_attribute("class3")->value();
			this->angle_type[cl1][cl2][cl3] = AngleType{stod(node->first_attribute("angle")->value()),
												stod(node->first_attribute("k")->value())};
			dbgmsg("harmonic angle force between " << cl1 << " and " << cl2 << " and " << cl3 <<" is "
					<< " angle = " << stod(node->first_attribute("angle")->value())
					<< " k = " << stod(node->first_attribute("k")->value()));
		}
		// parse torsions
		xml_node<> *torsions_node = doc.first_node("ForceField")->first_node("PeriodicTorsionForce");
		if (torsions_node) {
			for (xml_node<> *node = torsions_node->first_node(); node; node = node->next_sibling()) {
				const string cl1 = (node->first_attribute("class1")->value_size() == 0 ? "X" : node->first_attribute("class1")->value());
				const string cl2 = (node->first_attribute("class2")->value_size() == 0 ? "X" : node->first_attribute("class2")->value());
				const string cl3 = (node->first_attribute("class3")->value_size() == 0 ? "X" : node->first_attribute("class3")->value());
				const string cl4 = (node->first_attribute("class4")->value_size() == 0 ? "X" : node->first_attribute("class4")->value());
				const bool is_proper = (string("Proper").compare(node->name()) == 0);
				dbgmsg("node_name = " << node->name() << " is_proper = " << is_proper);
				xml_attribute<> *attr = node->first_attribute("periodicity1");
				while(attr) {
					const int periodicity = stoi(attr->value());
					attr = attr->next_attribute();
					const double phase = stod(attr->value());
					attr = attr->next_attribute();
					const double k = stod(attr->value());
					attr = attr->next_attribute();
					if (is_proper) {
						TorsionTypeVec &dihedral = this->torsion_type[cl1][cl2][cl3][cl4];
						dihedral.push_back(TorsionType{periodicity, phase, k});
						dbgmsg("periodic torsion (dihedral) force between " 
							<< cl1 << " and " << cl2 << " and " << cl3 << " and " << cl4 << " is "
							<< " periodicity" << dihedral.size() << " = " << periodicity
							<< " phase" << dihedral.size() << " = " << phase
							<< " k" << dihedral.size() << " = " << k);
					} else {
						TorsionTypeVec &improper = this->improper_type[cl2][cl3][cl1][cl4];
						improper.push_back(TorsionType{periodicity, phase, k});
						dbgmsg("periodic torsion (improper) force between " 
							<< cl1 << " and " << cl2 << " and " << cl3 << " and " << cl4 << " is "
							<< " periodicity" << improper.size() << " = " << periodicity
							<< " phase" << improper.size() << " = " << phase
							<< " k" << improper.size() << " = " << k
							<< " (here the FIRST atom is central)");
					}
				}
			}
		}
		dbgmsg("before parsing nonbonded force");
		// parse nonbonded force
		xml_node<> *nb_node = doc.first_node("ForceField")->first_node("NonbondedForce");
		if (nb_node) {
			this->coulomb14scale = stod(nb_node->first_attribute("coulomb14scale")->value());
			this->lj14scale = stod(nb_node->first_attribute("lj14scale")->value());
			for (xml_node<> *node = nb_node->first_node(); node; node = node->next_sibling()) {
				const int type = stoi(node->first_attribute("type")->value());
				AtomType &at = this->atom_type[type];
				at.charge = stod(node->first_attribute("charge")->value());
				at.sigma = stod(node->first_attribute("sigma")->value());
				at.epsilon = stod(node->first_attribute("epsilon")->value());
				dbgmsg("nonbonded force for atom type " << type
						<< " charge = " << at.charge
						<< " sigma = " << at.sigma
						<< " epsilon = " << at.epsilon);
			}
		}
		// parse knowledge-based force
		xml_node<> *kb_node = doc.first_node("ForceField")->first_node("KBForce");
		if (kb_node) {
			//~ this->step = stod(kb_node->first_attribute("step")->value());
			//~ this->scale = stod(kb_node->first_attribute("scale")->value());
			for (xml_node<> *node = kb_node->first_node(); node; node = node->next_sibling()) {
				const int cl1 = help::idatm_mask.at(node->first_attribute("class1")->value());
				const int cl2 = help::idatm_mask.at(node->first_attribute("class2")->value());
				vector<string> pot = help::ssplit(node->first_attribute("potential")->value());
				vector<string> der = help::ssplit(node->first_attribute("derivative")->value());
				vector<double> potential(pot.size()), derivative(der.size());
				transform(pot.begin(), pot.end(), potential.begin(), [this](const std::string& val) { return stod(val); });
				transform(der.begin(), der.end(), derivative.begin(), [this](const std::string& val) { return stod(val); });
				this->kb_force_type[cl1][cl2] = KBType{potential, derivative};
				dbgmsg("knowledge-based force between " << help::idatm_unmask[cl1] 
					<< " and " << help::idatm_unmask[cl2] << " is "
					<< " potential = " << stod(node->first_attribute("potential")->value()) 
					<< " derivative = " << stod(node->first_attribute("derivative")->value()));
			}
		}
		delete [] c_ff_file;
		return *this;
	}
}
