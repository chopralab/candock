#ifndef XSCORE_H
#define XSCORE_H

#include <vector>
#include "molib/molecule.hpp"

namespace Score {

        std::array<double, 5> vina_xscore(const Molib::Atom::Grid &gridrec, const Molib::Atom::Vec &atoms);

}

#endif
