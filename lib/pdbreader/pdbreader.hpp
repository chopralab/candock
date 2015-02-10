#ifndef PDBREADER_H
#define PDBREADER_H
#include <memory>
#include <set>
#include <map>
#include <vector>
#include <string>
#include <boost/regex.hpp>
#include <boost/algorithm/string.hpp>
#include <boost/filesystem.hpp>
#include <boost/date_time/posix_time/posix_time.hpp>
#include <boost/algorithm/string/split.hpp>
#include <boost/algorithm/string.hpp>
#include <boost/asio/ip/host_name.hpp>
#include "helper/debug.hpp"
using namespace std;

namespace Molib {
	class Molecules;
	class PDBreader {
	public:
		typedef enum {first_model=1, all_models=2, hydrogens=4, skip_hetatm=8} pdb_read_options;
	private:
		class Parser {
		protected:
			const string __molecule_file;
			unsigned int __hm;
			const int __num_occur;
			streampos __pos;
			bool __giant_molecule;
			void __generate_molecule(Molecules&, bool&, const string&);
			void __generate_assembly(Molecules&, bool&, int, const string&);
			void __generate_model(Molecules&, bool&, int);
		public:
			Parser(const string &molecule_file, unsigned int hm=PDBreader::all_models, 
				const int num_occur=-1)	: __molecule_file(molecule_file), __hm(hm), 
				__num_occur(num_occur), __pos(0), __giant_molecule(false) {} 
			virtual void parse_molecule(Molecules&) = 0;
			virtual void set_pos(streampos pos) { __pos = pos; };
			virtual void set_hm(unsigned int hm) { __hm = hm; };
		};
		class PdbParser : public Parser {
		public:
			using Parser::Parser;
			void parse_molecule(Molecules&);
		};
		class Mol2Parser : public Parser {
		public:
			using Parser::Parser;
			void parse_molecule(Molecules&);
		};
		Parser *p;
	public:
		PDBreader(const string &molecule_file, unsigned int hm=all_models, 
			const int num_occur=-1);
		~PDBreader() { delete p; }
		void rewind();
		void set_flags(unsigned int hm);
		void parse_molecule(Molecules &mols);
		Molecules parse_molecule();
	};
};
#endif
