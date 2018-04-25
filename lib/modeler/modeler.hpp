#ifndef MODELER_H
#define MODELER_H
#include <string>
#include <map>
#include <set>
#include <vector>

#include "geom3d/coordinate.hpp"
#include "molib/molecule.hpp"
#include "helper/help.hpp"
#include "helper/benchmark.hpp"
#include "topology.hpp"
#include "systemtopology.hpp"

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
  double __dist_cutoff_in_nm;
  double __tolerance;
  int __max_iterations;
  int __update_freq;
  double __position_tolerance_in_nm;
  bool __use_constraints;
  double __step_size_in_ps;
  double __temperature;
  double __friction;

  Geom3D::Point::Vec __positions;
  Topology __topology;
  SystemTopology __system_topology;

  bool __run_dyanmics;
  int __dynamics_steps;

public:
  Modeler(const ForceField &ffield, const string &fftype, double dist_cutoff,
          double tolerance, int max_iterations, int update_freq, double position_tolerance,
          bool use_constraints, double step_size_in_fs, double temperature, double friction);

  void mask(const Molib::Atom::Vec &atoms);
  void unmask(const Molib::Atom::Vec &atoms);

  void add_topology(const Molib::Atom::Vec &atoms);
  void add_crds(const Molib::Atom::Vec &atoms, const Geom3D::Point::Vec &crds);
  void add_random_crds(const Molib::Atom::Vec &atoms);

  Geom3D::Point::Vec get_state(const Molib::Atom::Vec &atoms);

#ifndef NDEBUG
  void minimize_knowledge_based(Molib::Molecule &ligand, Molib::Molecule &receptor, Score::Score &score);
#endif

  void minimize_state();
  void dynamics();

  void knowledge_based_calculation();
  void physical_calculation();

  void init_openmm_positions();
  void init_openmm(SystemTopology::integrator_type type = SystemTopology::integrator_type::none, string platform = "CPU");

  void set_max_iterations(const int max_iterations) { __max_iterations = max_iterations; }
  void set_num_steps_to_run(const int num_steps_to_run) { __dynamics_steps = num_steps_to_run; }

  double potential_energy();
};
}

#endif
