#ifndef MODELER_H
#define MODELER_H
#include <string>
#include <map>
#include <set>
#include <vector>

#include "candock/geometry/coordinate.hpp"
#include "candock/molib/molecule.hpp"
#include "candock/helper/help.hpp"
#include "candock/helper/benchmark.hpp"
#include "candock/modeler/topology.hpp"
#include "candock/modeler/systemtopology.hpp"

using namespace std;

namespace Score
{
class Score;
};

namespace OMMIface
{

class Modeler
{
public:
  class MinimizationError : public Error
  {
  public:
    MinimizationError(const string &msg) : Error(msg) {}
  };

private:
  const ForceField *__ffield;
  string __fftype;
  double __tolerance;
  int __max_iterations;
  bool __use_constraints;
  double __step_size_in_ps;
  double __temperature;
  double __friction;

  geometry::Point::Vec __positions;
  Topology __topology;
  SystemTopology __system_topology;

  int __dynamics_steps;

public:
  Modeler(const ForceField &ffield,
          const string &fftype = "none",
          double tolerance = 0.0001,
          int max_iterations = 100,
          bool use_constraints = false,
          double step_size_in_fs = 2.0,
          double temperature = 300.0,
          double friction = 91.0
  );

  void mask(const Molib::Atom::Vec &atoms);
  void unmask(const Molib::Atom::Vec &atoms);

  void add_topology(const Molib::Atom::Vec &atoms);
  void add_crds(const Molib::Atom::Vec &atoms, const geometry::Point::Vec &crds);
  void add_random_crds(const Molib::Atom::Vec &atoms);

  geometry::Point::Vec get_state(const Molib::Atom::Vec &atoms);

#ifndef NDEBUG
  void minimize_knowledge_based(Molib::Molecule &ligand, Molib::Molecule &receptor, Score::Score &score);
#endif

  void minimize_state();
  void dynamics();

  void init_openmm_positions();
  void init_openmm(const std::string& platform,
    const std::string& precision, const std::string& accelerators,
    SystemTopology::integrator_type type = SystemTopology::integrator_type::none);

  void set_max_iterations(const int max_iterations) { __max_iterations = max_iterations; }
  void set_num_steps_to_run(const int num_steps_to_run) { __dynamics_steps = num_steps_to_run; }

  double potential_energy();

  const ForceField& get_forcefield() {
    return *__ffield;
  }
};
}

#endif
