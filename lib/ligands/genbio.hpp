#ifndef GENBIO_HPP
#define GENBIO_HPP

namespace genbio {
	void generate_biological_assemblies(const string &models, const bool hydrogens,
		const string &pdb_dirname, const string &qpdb_file, const string &qcid,
		const bool neighb, const bool rnolig, const string &bio, const bool ralch,
		const string &json_file, const string &bio_file, const bool noalch,
		const string &mols_name, const bool rasym);
}

#endif
