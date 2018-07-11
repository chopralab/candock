#ifndef MODEL_H
#define MODEL_H
#include "candock/fragmenter/fragmenter.hpp"
#include "candock/geometry/geometry.hpp"
#include "candock/geometry/matrix.hpp"
#include "candock/molib/atom.hpp"
#include "candock/molib/chain.hpp"
#include "candock/molib/element.hpp"
#include "candock/molib/grid.hpp"
#include "candock/molib/it.hpp"
#include "candock/molib/residue.hpp"

namespace candock {

namespace molib {
class Chain;
class Residue;
class Atom;
class Assembly;

class Model : public template_map_container<Chain, Model, Assembly, char> {
    int __number;
    std::map<std::pair<int, Residue::res_tuple2>,
             std::map<Residue::res_tuple2, Residue::res_tuple2>>
        __remarks;
    Fragmenter::Fragment::Vec __rigid;

   public:
    Model(int number) : __number(number) {}
    Model(const Model& rhs) : __number(rhs.__number), __remarks(rhs.__remarks) {
        for (auto& chain : rhs) {
            dbgmsg("Copy constructor : model");
            add(new Chain(chain));
        }
    }
    Fragmenter::Fragment::Vec& get_rigid() { return __rigid; }
    const Fragmenter::Fragment::Vec& get_rigid() const { return __rigid; }
    void set_rigid(const Fragmenter::Fragment::Vec& rigid) { __rigid = rigid; }
    void init_bio(const Model& model_asym, const geometry::Matrix& matrix,
                  const std::set<char>& chains);
    void rotate(const geometry::Matrix& rota, const bool inverse = false);
    Chain& add(Chain* chain) {
        return this->aadd(chain->chain_id(), chain, this);
    }
    void set_remark(
        const int remark_number, const Residue::res_tuple2& ligand,
        std::pair<const Residue::res_tuple2&, const Residue::res_tuple2&>
            rpair) {
        this->__remarks[make_pair(remark_number, ligand)].insert(
            make_pair(rpair.first, rpair.second));
    }
    Model& regenerate_bonds(const Model&);
    void set_number(const int& i) { __number = i; }
    Chain& chain(const char chain_id) const { return this->element(chain_id); }
    bool has_chain(const char chain_id) const {
        return this->has_element(chain_id);
    }
    int number() const { return __number; }

    bool remarks(const int remark_number, const Residue::res_tuple2 ligand,
                 std::map<Residue::res_tuple2, Residue::res_tuple2>& r) {
        auto i = __remarks.find(make_pair(remark_number, ligand));
        if (i == __remarks.end()) return false;
        r = i->second;
        return true;
    }

    Atom::Vec get_atoms(
        const std::string& chain_ids = "",
        const Residue::res_type& rest = Residue::res_type::notassigned) const;
    Model& erase_properties() {
        for (auto& chain : *this) chain.erase_properties();
        return *this;
    }
    friend std::ostream& operator<<(std::ostream& stream, const Model& m);
};

}  // molib
}

#endif
