#include <candock/interface/interface.hpp>

#include <stdlib.h>
#include <stdio.h>
#include <string.h>
#include <time.h>

int main( int argc, char **argv) {

        if (argc < 2) {
               puts("Error, please give at least 1 arguments\n");
               return 1;
        }

        clock_t start = clock();

        if (! initialize_complex(argv[1])) {
                printf("%s\n", cd_get_error());
                return 1;
        }

        printf("Time to initialize: %f\n", (double) ( clock() - start) / CLOCKS_PER_SEC );

        size_t rec_atom_count = receptor_atom_count();
        char*   chain= (char*)  malloc(rec_atom_count * sizeof(char) );
        size_t* resi = (size_t*)malloc(rec_atom_count * sizeof(size_t));
        size_t* rest = (size_t*)malloc(rec_atom_count * sizeof(size_t));
        char*   resn = (char*)  malloc(rec_atom_count * sizeof(char) );
        size_t* elem = (size_t*)malloc(rec_atom_count * sizeof(size_t));
        receptor_atom_details(chain, resi, rest, resn, elem);

        printf("Number of atoms: %zu\n", rec_atom_count);

        char current_chain = 0;
        size_t current_resi = 0;
        size_t print_count = 0;
        for ( size_t i = 0; i < rec_atom_count; ++i) {
                if (current_chain != chain[i]) {
                        current_chain = chain[i];
                        printf("\nChain %c:\n", current_chain);
                        print_count = 0;
                }
                if (current_resi != resi[i]) {
                        current_resi = resi[i];
                        putchar(resn[i]);
                        ++print_count;
                }
                if (print_count >= 80) {
                        printf("\n");
                        print_count = 0;
                }
        }
        
        printf("\n");
        
        size_t lig_atom_count = ligand_atom_count();
        char*   lchain= (char*)  malloc(lig_atom_count * sizeof(char) );
        size_t* lresi = (size_t*)malloc(lig_atom_count * sizeof(size_t));
        size_t* lrest = (size_t*)malloc(lig_atom_count * sizeof(size_t));
        size_t* lelem = (size_t*)malloc(lig_atom_count * sizeof(size_t));
        ligand_atom_details(lchain, lresi, lrest, lelem);

        current_chain = 0;
        current_resi = 0;
        print_count = 0;
        for ( size_t i = 0; i < lig_atom_count; ++i) {
                if (current_chain != lchain[i]) {
                        current_chain = lchain[i];
                        printf("\nChain %c:\n", current_chain);
                        print_count = 0;
                }
                if (current_resi != lresi[i]) {
                        printf("%zu\n", lresi[i]);
                        current_resi = lresi[i];                        
                        ++print_count;
                }
        }

        printf("\n");

        size_t* neighs = (size_t*) malloc(4 * sizeof(size_t));
        
        size_t neigh_count = ligand_get_neighbors(5,neighs);
        printf("%zu %zu %zu %zu %zu\n", neighs[0], neighs[1], neighs[2], neighs[3], neigh_count);
        
        size_t* neighs2 = (size_t*) malloc(4 * sizeof(size_t));
        ligand_get_neighbors(neighs[1], neighs2);
        printf("%zu %zu %zu %zu\n", neighs2[0], neighs2[1], neighs2[2], neighs2[3]);

        size_t* neighs3 = (size_t*) malloc(4 * sizeof(size_t));
        ligand_get_neighbors(neighs2[1], neighs3);
        printf("%zu %zu %zu %zu\n", neighs3[0], neighs3[1], neighs3[2], neighs3[3]);
        
        free(chain);
        free(resi);
        free(resn);
        free(rest);
        free(elem);
        free(neighs);
        free(neighs2);
        free(neighs3);
        
        return 0;
}

