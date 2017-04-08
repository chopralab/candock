#include "interface.hpp"

#include <stdio.h>
#include <string.h>

int main( int argc, char **argv) {
        
        if (argc < 3) {
               puts("Error, please give at least 3 arguments\n");
               return 1;
        }

        printf("%s\n",initialize_receptor(argv[1]));
        const char * ligand = initialize_ligand(argv[2]);
        
        for ( int i = 0; i < strlen(ligand); ++i) {
                putchar(ligand[i]);
        }

        printf("Done!\n");

        return 0;
}

