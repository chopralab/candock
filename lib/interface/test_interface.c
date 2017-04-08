#include "interface.hpp"

#include <stdio.h>
#include <string.h>

int main( int argc, char **argv) {
        
        if (argc < 5) {
               puts("Error, please give at least 4 arguments\n");
               return 1;
        }

        const char* my_recetor = initialize_receptor(argv[1]);
        const char* ligand = initialize_ligand(argv[2]);
        
        initialize_scoring(argv[3]);
        initialize_ffield(argv[4]);

        printf("Score before setting: %f\n", calculate_score());

        unsigned long idx[24];
        float positions[24 * 3];

        for (int i = 0; i <= 23; ++i) {
                idx[ i ] = i + 2153;
                positions[ i * 3 + 0 ] = 0.0f;
                positions[ i * 3 + 1 ] = 0.0f;
                positions[ i * 3 + 2 ] = 0.0f;
        }

        set_positions_ligand( idx, positions, 24 );

        printf("Score after setting: %f\n", calculate_score());

        printf("Done!\n");

        return 0;
}

