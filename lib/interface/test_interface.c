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

        printf("%f\n", calculate_score());

        printf("Done!\n");

        return 0;
}

