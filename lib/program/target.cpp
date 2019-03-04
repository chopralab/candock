#include "candock/program/target.hpp"

#include <boost/filesystem.hpp>

#include "statchem/fragmenter/unique.hpp"
#include "statchem/helper/path.hpp"
#include "statchem/helper/logger.hpp"
#include "statchem/molib/molecule.hpp"
#include "statchem/parser/fileparser.hpp"

#include "statchem/score/kbff.hpp"

#include "statchem/modeler/modeler.hpp"
#include "statchem/modeler/systemtopology.hpp"

#include "candock/program/dockfragments.hpp"
#include "candock/program/findcentroids.hpp"
#include "candock/program/linkfragments.hpp"

#include "statchem/fileio/fileout.hpp"
#include "statchem/fileio/inout.hpp"

namespace candock {
namespace Program {

Target::Target(const std::string& input_name) {
    /* Initialize parsers for receptor and read
     * the receptor molecule(s)
     *
     */
    if (fileio::file_size(input_name) <= 0) {
        throw Error("Provided file or directory does not exist: " + input_name);
    }

    // If the option given is a regular file, act like previous versions
    parser::FileParser rpdb(input_name, parser::first_model);
    molib::Molecules receptors = rpdb.parse_molecule();

    /* Compute atom types for receptor and cofactor(s): gaff types for protein,
     * Mg ions, and water are read from the forcefield xml file later on while
     * gaff types for cofactors (ADP, POB, etc.) are calculated de-novo here
     *
     */
    receptors.compute_idatm_type()
        .compute_hydrogen()
        .compute_bond_order()
        .compute_bond_gaff_type()
        .refine_idatm_type()
        .erase_hydrogen()    // needed because refine changes connectivities
        .compute_hydrogen()  // needed because refine changes connectivities
        .compute_ring_type()
        .compute_gaff_type()
        .compute_rotatable_bonds()  // relies on hydrogens being assigned
        .erase_hydrogen();

    __protein = std::unique_ptr<molib::Molecule>(
        new molib::Molecule(std::move(receptors[0])));

    // Emulate the original version of candock
    __protein->set_name(input_name.substr(0, input_name.length() - 4));
    boost::filesystem::create_directory(__protein->name());

    /* Create receptor grid
     *
     */
    molib::Atom::Grid* gridrec = new molib::Atom::Grid(__protein->get_atoms());
    __gridrec = std::unique_ptr<molib::Atom::Grid>(gridrec);

    // Just initialize, here (copies strings only) and do not run the
    // calculation
    Program::FindCentroids* centroids = new FindCentroids(
        input_name, __protein->get_chain_ids(molib::Residue::protein | molib::Residue::nucleic),
        __protein->name());
    __centroids = std::unique_ptr<Program::FindCentroids>(centroids);
}

/* Read distributions file and initialize scores
 *
 */

void Target::__initialize_score(const FragmentLigands& ligand_fragments) {
    if (__score != nullptr) {
        return;
    }

    score::Score* score = new score::Score(
        cmdl.get_string_option("ref"), cmdl.get_string_option("comp"),
        cmdl.get_string_option("func"), cmdl.get_int_option("cutoff"));
    __score = std::unique_ptr<score::Score>(score);

    __score
        ->define_composition(__protein->get_idatm_types(),
                             ligand_fragments.ligand_idatm_types())
        .process_distributions(cmdl.get_string_option("dist"))
        .compile_scoring_function();
}

void Target::__initialize_ffield() {
    if (__ffield != nullptr) {
        return;
    }

    // Prepare the receptor for docking to
    OMMIface::ForceField* ffield = new OMMIface::ForceField;
    __ffield = std::unique_ptr<OMMIface::ForceField>(ffield);

    __ffield->parse_gaff_dat_file(cmdl.get_string_option("gaff_dat"))
        .parse_forcefield_file(cmdl.get_string_option("amber_xml"))
        .parse_forcefield_file(cmdl.get_string_option("water_xml"));

    if (!cmdl.get_string_option("gaff_heme").empty()) {
        dbgmsg("Adding " << cmdl.get_string_option("gaff_heme") << endl);
        __ffield->parse_gaff_dat_file(cmdl.get_string_option("gaff_heme"));
    }

    __protein->prepare_for_mm(*__ffield, *__gridrec);
}

void Target::__initialize_kbforce(const FragmentLigands& ligand_fragments) {
    OMMIface::SystemTopology::loadPlugins();

    score::KBFF* score = new score::KBFF(
        cmdl.get_string_option("ff_ref"), cmdl.get_string_option("ff_comp"),
        cmdl.get_string_option("ff_func"), cmdl.get_int_option("ff_cutoff"),
        cmdl.get_double_option("step"));

    __ff_score = std::unique_ptr<score::KBFF>(score);

    __ff_score
        ->define_composition(__protein->get_idatm_types(),
                             ligand_fragments.ligand_idatm_types())
        .process_distributions(cmdl.get_string_option("dist"))
        .compile_scoring_function();

    if (cmdl.get_string_option("obj_dir").empty()) {
        __ff_score->compile_objective_function();
    } else {
        __ff_score->parse_objective_function(
            cmdl.get_string_option("obj"), cmdl.get_double_option("scale"), 15);
    }

    __ffield->add_kb_forcefield(*__ff_score,
                                cmdl.get_double_option("dist_cutoff"));
}

std::set<int> Target::get_idatm_types(const std::set<int>& previous) const {
    return __protein->get_idatm_types(previous);
}

void Target::find_centroids() {
    /* Run section of Candock designed to find binding site1s
     * Currently, this runs ProBIS and does not require any
     * previous step to be competed.
     *
     */

    __centroids->run_step();
}

void Target::make_gridhcp(const FragmentLigands& ligand_fragments) {
    __initialize_score(ligand_fragments);
    __initialize_ffield();

    find_centroids();

    if (__prepseeds == nullptr) {
        Program::DockFragments* prepseeds =
            new Program::DockFragments(*__centroids, ligand_fragments, *__score,
                                       *__gridrec, __protein->name());
        __prepseeds = std::unique_ptr<Program::DockFragments>(prepseeds);
    }

    docker::Gpoints gpoints = __prepseeds->get_gridhcp();
    fileio::output_file(
        gpoints,
        Path::join(__protein->name(), cmdl.get_string_option("gridpdb_hcp")));
}

void Target::dock_fragments(const FragmentLigands& ligand_fragments) {
    __initialize_score(ligand_fragments);
    __initialize_ffield();

    find_centroids();

    if (__prepseeds == nullptr) {
        Program::DockFragments* prepseeds =
            new Program::DockFragments(*__centroids, ligand_fragments, *__score,
                                       *__gridrec, __protein->name());
        __prepseeds = std::unique_ptr<Program::DockFragments>(prepseeds);
    }

    __prepseeds->run_step();
}

void Target::link_fragments(const FragmentLigands& ligand_fragments) {
    dock_fragments(ligand_fragments);

    __initialize_kbforce(ligand_fragments);

    if (__dockedlig == nullptr) {
        __ffield->insert_topology(*__protein);
        Program::LinkFragments* dockedlig = new LinkFragments(
            *__protein, *__score, *__ffield, *__prepseeds, *__gridrec);
        __dockedlig = std::unique_ptr<Program::LinkFragments>(dockedlig);
    }

    __dockedlig->run_step();
}

void Target::link_fragments(const molib::Molecules& ligands) {
    if (__prepseeds == nullptr) {
        throw Error("You must run dock_fragments first!");
    }

    if (__dockedlig == nullptr) {
        __ffield->insert_topology(*__protein);
        Program::LinkFragments* dockedlig = new LinkFragments(
            *__protein, *__score, *__ffield, *__prepseeds, *__gridrec);
        __dockedlig = std::unique_ptr<Program::LinkFragments>(dockedlig);
    }

    __dockedlig->link_ligands(ligands);
}

void Target::make_scaffolds(const std::set<std::string>& seeds_to_add,
                            molib::Molecules& all_designs_out) {
    molib::Unique created_design("designed.txt");

    molib::NRset nr = __prepseeds->get_top_seeds(
        seeds_to_add, cmdl.get_double_option("top_percent"));
    for (auto& molecules : nr) {
        design::Design designer(molecules.first(), created_design);
        designer.change_original_name(molecules.name());
        designer.functionalize_hydrogens_with_fragments(
            nr, cmdl.get_double_option("tol_seed_dist"),
            cmdl.get_double_option("clash_coeff"),
            std::make_tuple(cmdl.get_double_option("lipinski_mass"),
                            cmdl.get_int_option("lipinski_hbd"),
                            cmdl.get_int_option("lipinski_nhs"),
                            cmdl.get_int_option("lipinski_ohs")));
#ifndef NDEBUG
        fileio::output_file(designer.get_internal_designs(),
                           "internal_designs.pdb", std::ios_base::app);
#endif
        all_designs_out.add(designer.prepare_designs());
    }

    if (all_designs_out.size() == 0) {
        log_step << "No new designs, exiting" << std::endl;
        return;
    }

    created_design.write_out();
}

void Target::design_ligands(const std::set<std::string>& seeds_to_add,
                            molib::Molecules& all_designs_out) {
    molib::Unique created_design("designed.txt");

    for (auto& molecule : __dockedlig->top_poses()) {
        design::Design designer(molecule, created_design);
        if (!seeds_to_add.empty())
            designer.functionalize_hydrogens_with_fragments(
                __prepseeds->get_top_seeds(
                    seeds_to_add, cmdl.get_double_option("top_percent")),
                cmdl.get_double_option("tol_seed_dist"),
                cmdl.get_double_option("clash_coeff"),
                std::make_tuple(cmdl.get_double_option("lipinski_mass"),
                                cmdl.get_int_option("lipinski_hbd"),
                                cmdl.get_int_option("lipinski_nhs"),
                                cmdl.get_int_option("lipinski_ohs")));

        // const vector<string>& h_single_atoms =
        // cmdl.get_string_vector("add_single_atoms");
        const std::vector<std::string>& a_single_atoms =
            cmdl.get_string_vector("change_terminal_atom");

        if (!a_single_atoms.empty())
            designer.functionalize_extremes_with_single_atoms(a_single_atoms);
// if (!h_single_atoms.empty())
//        designer.functionalize_hydrogens_with_single_atoms(h_single_atoms);
#ifndef NDEBUG
        fileio::output_file(designer.get_internal_designs(),
                           "internal_designs.pdb", std::ios_base::app);
#endif
        all_designs_out.add(designer.prepare_designs());
    }

    created_design.write_out();
}
}
}
