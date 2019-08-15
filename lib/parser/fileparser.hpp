/* Copyright (c) 2016-2019 Chopra Lab at Purdue University, 2013-2016 Janez Konc at National Institute of Chemistry and Samudrala Group at University of Washington
 *
 * This program is free for educational and academic use
 *
 * This program is free software; you can redistribute it and/or modify
 * it under the terms of the GNU General Public License as published by
 * the Free Software Foundation version 3 of the License.
 *
 * This program is distributed in the hope that it will be useful,
 * but WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the
 * GNU General Public License for more details.
 */

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
