#ifndef PDBREADER_H
#define PDBREADER_H
#include <memory>
#include <set>
#include <map>
#include <vector>
#include <string>

#include "parser.hpp"
#include "helper/debug.hpp"

using namespace std;

namespace Parser {
        class FileParser {
        private:
                class PdbParser : public Parser {
                public:
                        using Parser::Parser;
                        void parse_molecule (Molib::Molecules &);
                };
                class Mol2Parser : public Parser {
                public:
                        using Parser::Parser;
                        void parse_molecule (Molib::Molecules &);
                };
                Parser *p;
        public:
                FileParser() : p (nullptr) {};
                FileParser (const string &molecule_file, unsigned int hm=all_models,
                            const int num_occur=-1);
                ~FileParser() {
                        delete p;
                }
                void prepare_parser (const string &molecule_file, unsigned int hm=all_models,
                                     const int num_occur=-1);
                void rewind();
                void set_flags (unsigned int hm);
                bool parse_molecule (Molib::Molecules &mols);
                Molib::Molecules parse_molecule();
        };
}

#endif
