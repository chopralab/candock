include_directories($CMAKE_CURRENT_SOURCE_DIR)

set( pdbreader_SRC
        assembly.cpp
        assembly.hpp
        atom.cpp
        atom.hpp
        atomtype.cpp
        atomtype.hpp
        bond.cpp
        bond.hpp
        bondtype.cpp
        bondtype.hpp
        chain.cpp
        chain.hpp
        element.cpp
        element.hpp
        grid.hpp
        hydrogens.cpp
        hydrogens.hpp
        internal.cpp
        internal.hpp
        it.hpp
        model.cpp
        model.hpp
        molecule.cpp
        molecule.hpp
        molecules.cpp
        molecules.hpp
        nrset.cpp
        nrset.hpp
        pdbreader.cpp
        pdbreader.hpp
        residue.cpp
        residue.hpp
)

include_directories( ../modeler )

add_library(pdbreader OBJECT ${pdbreader_SRC})
