#ifndef XSCORE_H
#define XSCORE_H

#include <vector>
#include "candock/molib/molecule.hpp"

namespace candock {

namespace score {

std::array<double, 5> vina_xscore(const molib::Atom::Grid& gridrec,
                                  const molib::Atom::Vec& atoms);
}
}

#endif
