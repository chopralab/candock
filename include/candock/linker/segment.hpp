#ifndef SEGMENT_H
#define SEGMENT_H
#include <functional>
#include <tuple>

#include "statchem/fragmenter/fragmenter.hpp"
#include "statchem/geometry/coordinate.hpp"
#include "statchem/graph/graph.hpp"
#include "statchem/helper/debug.hpp"
#include "statchem/molib/internal.hpp"
#include "statchem/molib/it.hpp"
#include "statchem/molib/molecule.hpp"

#include "candock/linker/state.hpp"

namespace candock {

namespace linker {

class Segment : public molib::template_vector_container<Segment*, Segment> {
   public:
    typedef std::set<Segment*> Set;
    typedef std::set<const Segment*> ConstSet;
    typedef std::vector<Segment*> Vec;
    typedef std::pair<const Segment*, const Segment*> ConstPair;
    typedef graph::Graph<Segment> Graph;
    typedef std::map<ConstPair, Graph::Path> Paths;
    typedef int Id;

   private:
    std::map<const molib::Atom*, int> __amap;
    const molib::Atom::Vec __atoms;
    int __seed_id;  // != -1 if segment is a seed
    Id __id;
    std::vector<std::unique_ptr<State>> __state;  // only seed states here!
    ConstSet __adjacent_seed_segments;
    std::map<const Segment*, double> __max_linker_length;
    std::map<const Segment*, molib::Bond> __bond;
    std::map<const Segment*, Segment*> __next;
    std::vector<bool> __join_atom, __common_atom;

    Segment(const molib::Atom::Vec atoms, const int& seed_id,
            const Segment::Id idx);
    static Paths __find_paths(
        const std::vector<std::unique_ptr<Segment>>& segments);
    static void __set_branching_rules(const Paths& paths);
    static bool __link_adjacent(const Graph::Path& path);
    static void __init_max_linker_length(const Paths& paths);
    static void __compute_max_linker_length(Segment::Graph::Path& path);

   public:
    int get_seed_id() const { return __seed_id; }
    bool has_next(const Segment& goal) const {
        return __next.count(&goal) != 0;
    }
    Segment& get_next(const Segment& goal) const {
        return *__next.at(&goal);
    }  // get next seg in the direction of goal
    void set_next(Segment& goal, Segment& next) {
        __next.insert({&goal, &next});
    }
    const ConstSet& get_adjacent_seed_segments() const {
        return __adjacent_seed_segments;
    }
    void set_adjacent_seed_segments(Segment& seed_seg) {
        __adjacent_seed_segments.insert(&seed_seg);
    }
    bool is_seed_adjacent(const Segment& other) const {
        return __adjacent_seed_segments.count(&other) != 0;
    }
    const molib::Atom::Vec& get_atoms() const { return __atoms; }
    const molib::Atom& get_atom(const int i) const { return *__atoms[i]; }
    bool has_atom(const molib::Atom& atom) const {
        return __amap.count(&atom) != 0;
    }
    int get_idx(const molib::Atom& atom) const { return __amap.at(&atom); }
    std::string get_label() const {
        std::stringstream ss;
        ss << *this;
        return ss.str();
    }                                 // graph ostream operator
    int weight() const { return 0; }  // dummy for graph ostream operator
    void add_state(std::unique_ptr<State> s) {
        __state.push_back(std::move(s));
    }
    bool is_seed() const { return __seed_id != -1; }
    bool is_leaf() const { return size() == 1; }
    bool is_branch() const { return size() > 2; }
    const std::vector<std::unique_ptr<State>>& get_states() const {
        return __state;
    }
    State& get_first_state() const { return *__state[0]; }
    State& get_last_state() const { return *__state.back(); }
    bool is_adjacent(const Segment& other) const {
        for (auto& adj : *this)
            if (&adj == &other) return true;
        return false;
    }
    double get_max_linker_length(const Segment& other) const {
        return __max_linker_length.at(&other);
    }
    void set_max_linker_length(const Segment& other, const double d) {
        if (d > __max_linker_length[&other]) __max_linker_length[&other] = d;
    }
    const molib::Bond& get_bond(const Segment& other) const {
        return __bond.at(&other);
    }
    void set_bond(const Segment& other, molib::Atom& a1, molib::Atom& a2);
    int adjacent_in_segment(const molib::Atom& atom,
                            const molib::Atom& forbidden) const;
    Id get_id() const { return __id; }
    void set_join_atom(const molib::Atom& atom) {
        __join_atom[get_idx(atom)] = true;
    }
    bool is_join_atom(const int i) const { return __join_atom[i]; }
    void set_common_atom(const molib::Atom& atom) {
        __common_atom[get_idx(atom)] = true;
    }
    bool is_common_atom(const int i) const { return __common_atom[i]; }
    friend std::ostream& operator<<(std::ostream& stream, const Segment& s);

    static Graph create_graph(const molib::Molecule& molecule);
};
};
}

#endif
