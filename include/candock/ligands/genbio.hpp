#ifndef GENBIO_HPP
#define GENBIO_HPP

#include <string>

namespace candock {

namespace genbio {
	void generate_biological_assemblies(const std::string &models, const bool hydrogens,
		const std::string &pdb_dirname, const std::string &qpdb_file, const std::string &qcid,
		const bool neighb, const bool rnolig, const std::string &bio, const bool ralch,
		const std::string &json_file, const std::string &bio_file, const bool noalch,
		const std::string &mols_name, const bool rasym);
}

}

#endif
