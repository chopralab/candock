#ifndef TARGET_H
#define TARGET_H

#include "candock/design/design.hpp"
#include "candock/program/fragmentligands.hpp"

#include "candock/docker/gpoints.hpp"
#include "candock/molib/molecules.hpp"

#include "candock/modeler/forcefield.hpp"
#include "candock/program/dockfragments.hpp"
#include "candock/program/findcentroids.hpp"
#include "candock/program/linkfragments.hpp"
#include "candock/score/kbff.hpp"

#include <string>

namespace candock {

namespace Program {

class CANDOCK_EXPORT Target {
    std::unique_ptr<molib::Molecule> __protein;
    std::unique_ptr<score::Score> __score;
    std::unique_ptr<score::KBFF> __ff_score;
    std::unique_ptr<OMMIface::ForceField> __ffield;
    std::unique_ptr<molib::Atom::Grid> __gridrec;
    std::unique_ptr<Program::FindCentroids> __centroids;
    std::unique_ptr<Program::DockFragments> __prepseeds;
    std::unique_ptr<Program::LinkFragments> __dockedlig;

    void __initialize_score(const FragmentLigands& ligand_fragments);
    void __initialize_ffield();
    void __initialize_kbforce(const FragmentLigands& ligand_fragments);

   public:
    Target(const std::string& input_name);
    // Target(const Target &) {};
    // Target &operator= (const Target &) {};

    std::set<int> get_idatm_types(
        const std::set<int>& previous = std::set<int>()) const;

    void find_centroids();
    void make_gridhcp(const FragmentLigands& ligand_fragments);
    void dock_fragments(const FragmentLigands& ligand_fragments);
    void link_fragments(const FragmentLigands& ligand_fragments);
    void link_fragments(const molib::Molecules& ligand_fragments);
    void make_scaffolds(const std::set<std::string>& seeds_to_add,
                        molib::Molecules& all_designs_out);
    void design_ligands(const std::set<std::string>& seeds_to_add,
                        molib::Molecules& all_designs_out);

    const DockFragments& get_seeds() const { return *__prepseeds; }

    const std::string& name() const { return __protein->name(); }
};
}
}

#endif  // TARGET_H
