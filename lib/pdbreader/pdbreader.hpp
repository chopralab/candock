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

#include "parser.hpp"
#include "helper/debug.hpp"

using namespace std;

namespace Parser {
	class PDBreader {
	private:
		class PdbParser : public Parser {
		public:
			using Parser::Parser;
			void parse_molecule(Molib::Molecules&);
		};
		class Mol2Parser : public Parser {
		public:
			using Parser::Parser;
			void parse_molecule(Molib::Molecules&);
		};
		Parser *p;
	public:
		PDBreader() : p(nullptr){};
		PDBreader(const string &molecule_file, unsigned int hm=all_models, 
			const int num_occur=-1);
		~PDBreader() { delete p; }
		void prepare_parser(const string &molecule_file, unsigned int hm=all_models, 
			const int num_occur=-1);
		void rewind();
		void set_flags(unsigned int hm);
		bool parse_molecule(Molib::Molecules &mols);
		Molib::Molecules parse_molecule();
	};
}

#endif
