#ifndef GENLIG_H
#define GENLIG_H
#include <ostream>
#include "statchem/molib/molecules.hpp"

namespace candock {
namespace genlig {

typedef std::map<int, statchem::molib::Molecules> BindingSiteClusters;
typedef std::map<int, double> BindingSiteScores;
void generate_ligands(const std::string& receptor_file,
                      const std::string& receptor_chain_id,
                      const std::string& json_file, const std::string& bio_dir,
                      const std::string& lig_code, const std::string& lig_file,
                      const std::string& bsite_file);
std::pair<BindingSiteClusters, BindingSiteScores>
generate_binding_site_prediction(const std::string& json_with_ligs_file,
                                 const std::string& bio_dir,
                                 const int num_bsites);

// Operator overloading for typedef types can lead to issues.
// See http://blog.mezeske.com/?p=170 for details.

}

std::ostream& operator<<(std::ostream& os,
                         const std::map<int, statchem::molib::Molecules>& bsites);

std::ostream& operator<<(std::ostream& os,
                         const std::map<int, double>& bscores);

}

#endif
