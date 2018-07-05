#ifndef XSCORE_H
#define XSCORE_H

#include <vector>
#include "candock/molib/molecule.hpp"

namespace candock {

namespace Score {

        std::array<double, 5> vina_xscore(const Molib::Atom::Grid &gridrec, const Molib::Atom::Vec &atoms);

}

}

#endif
