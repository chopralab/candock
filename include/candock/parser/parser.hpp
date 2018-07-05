#ifndef PARSER_H
#define PARSER_H

#include <string>
#include "candock/molib/molecules.hpp"

namespace candock {

namespace Parser {
        enum pdb_read_options {
                first_model=        1 << 0,
                all_models=         1 << 1,
                hydrogens=          1 << 2,
                skip_hetatm=        1 << 3,
                skip_atom=          1 << 4,
                sparse_macromol=    1 << 5,
                docked_poses_only=  1 << 6,
                protein_poses_only= 1 << 7,
        };

        class Parser {
        protected:
                std::istream &__stream;
                std::mutex __concurrent_read_mtx;
                unsigned int __hm;
                const int __num_occur;
                bool __giant_molecule;
                void __generate_molecule(Molib::Molecules&, bool&, const std::string&);
                void __generate_assembly(Molib::Molecules&, bool&, int, const std::string&);
                void __generate_model(Molib::Molecules&, bool&, int);
        public:
                Parser( std::istream &molecule_file, unsigned int hm= all_models, 
                       const int num_occur=-1)
                       : __stream(molecule_file), __hm(hm), 
                         __num_occur(num_occur), __giant_molecule(false) {}
                virtual ~Parser() {}
                virtual void parse_molecule(Molib::Molecules&) = 0;
                virtual void set_pos(std::streampos pos);
                virtual void set_hm(unsigned int hm);
        };
}

}

#endif // PARSER_H
