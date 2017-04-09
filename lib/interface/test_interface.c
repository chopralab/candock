#include "interface.hpp"

#include <stdio.h>
#include <string.h>
#include <time.h>

int main( int argc, char **argv) {

        if (argc < 5) {
               puts("Error, please give at least 4 arguments\n");
               return 1;
        }

        clock_t start = clock();
        const char* my_recetor = initialize_receptor(argv[1]);
        const char* ligand = initialize_ligand(argv[2]);

        initialize_scoring(argv[3]);
        initialize_ffield(argv[4]);

        printf("Time to initialize: %f\n", (double) ( clock() - start) / CLOCKS_PER_SEC );

        start = clock();
        printf("Score before setting: %f\n", calculate_score());
        printf("Time to do one scoring: %f\n", (double) ( clock() - start) / CLOCKS_PER_SEC );

        unsigned long idx[24];
        float positions[24 * 3];

        for (int i = 0; i <= 23; ++i) {
                idx[ i ] = i + 2153;
                positions[ i * 3 + 0 ] = 0.0f;
                positions[ i * 3 + 1 ] = 0.0f;
                positions[ i * 3 + 2 ] = 0.0f;
        }

        start = clock();
        set_positions_ligand( idx, positions, 24 );
        printf("Time to set ligand: %f\n", (double) ( clock() - start) / CLOCKS_PER_SEC );

        start = clock();
        printf("Score after setting ligand: %f\n", calculate_score());
        printf("Time to score again: %f\n", (double) ( clock() - start) / CLOCKS_PER_SEC );

        unsigned long idx_rec[2152];
        float positions_rec[2152 * 3];

        for (int j = 0; j <= 2151; ++j) {
                idx_rec[ j ] = j + 1;
                positions_rec[ j * 3 + 0 ] = 0.0f;
                positions_rec[ j * 3 + 1 ] = 0.0f;
                positions_rec[ j * 3 + 2 ] = 0.0f;
        }

        start = clock();
        set_positions_receptor(idx_rec, positions_rec, 2152);
        printf("Time to set receptor: %f\n", (double) ( clock() - start) / CLOCKS_PER_SEC );

        start = clock();
        printf("Score after setting receptor: %f\n", calculate_score());
        printf("Time to score again: %f\n", (double) ( clock() - start) / CLOCKS_PER_SEC );

        return 0;
}

