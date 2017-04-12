#include "interface.hpp"

#include <stdlib.h>
#include <stdio.h>
#include <string.h>
#include <time.h>

int main( int argc, char **argv) {

        if (argc < 5) {
               puts("Error, please give at least 4 arguments\n");
               return 1;
        }

        clock_t start = clock();

        initialize_receptor(argv[1]);
        char* recetor = (char*)malloc( receptor_string_size() * sizeof(char) );
        copy_receptor_string(recetor);

        initialize_ligand(argv[2]);
        char* ligand  = (char*)malloc( receptor_string_size() * sizeof(char) );
        copy_ligand_string(ligand);

        initialize_scoring(argv[3]);
        initialize_ffield(argv[4]);

        printf("Time to initialize: %f\n", (double) ( clock() - start) / CLOCKS_PER_SEC );

        start = clock();
        printf("Score before setting: %f\n", calculate_score());
        printf("Time to do one scoring: %f\n", (double) ( clock() - start) / CLOCKS_PER_SEC );

        size_t  lig_atom_count = ligand_atom_count();
        size_t* idx = (size_t*)malloc( lig_atom_count * sizeof(size_t) );
        float * positions = (float*)malloc( lig_atom_count * sizeof(float) * 3 );
        ligand_atoms(idx, positions);

        for (int i = 0; i < lig_atom_count; ++i) {
                positions[ i * 3 + 0 ] = positions[ i * 3 + 0 ] + 0.5f;
                positions[ i * 3 + 1 ] = positions[ i * 3 + 1 ] + 0.5f;
                positions[ i * 3 + 2 ] = positions[ i * 3 + 2 ] + 0.5f;
        }

        start = clock();

        set_positions_ligand( idx, positions, lig_atom_count );
        
        printf("Time to set ligand: %f\n", (double) ( clock() - start) / CLOCKS_PER_SEC );

        start = clock();
        printf("Score after setting ligand: %f\n", calculate_score());
        printf("Time to score again: %f\n", (double) ( clock() - start) / CLOCKS_PER_SEC );

        size_t rec_atom_count = receptor_atom_count();
        size_t* idx_rec = (size_t*)malloc(rec_atom_count * sizeof(size_t));
        float * positions_rec = (float*)malloc(rec_atom_count * sizeof(float) * 3);
        receptor_atoms(idx_rec, positions_rec);

        for (int j = 0; j < rec_atom_count; ++j) {
                positions_rec[ j * 3 + 0 ] = positions_rec[ j * 3 + 0 ] + 0.5f;
                positions_rec[ j * 3 + 1 ] = positions_rec[ j * 3 + 1 ] + 0.5f;
                positions_rec[ j * 3 + 2 ] = positions_rec[ j * 3 + 2 ] + 0.5f;
        }

        start = clock();
        set_positions_receptor(idx_rec, positions_rec, rec_atom_count);
        printf("Time to set receptor: %f\n", (double) ( clock() - start) / CLOCKS_PER_SEC );

        start = clock();
        printf("Score after setting receptor: %f\n", calculate_score());
        printf("Time to score again: %f\n", (double) ( clock() - start) / CLOCKS_PER_SEC );

        return 0;
}

