/*extern "C" __global__ void calcNeighborList(
    real4 *__restrict__ posq,
    int2 *dev_neighbors,
    bool *dev_bonded_exclusions_matrix)
{
    int block = blockIdx.x;
    int thread = threadIdx.x;
    printf("block %d thread %d \n", block, thread);
    size_t counter = 0;

    for (size_t interaction = block + 1; interaction < num_atoms; interaction++)
    {

        int atom1 = block;
        int atom2 = interaction;

        if (atom1 > atom2)
        {
            real d_sq = pow(posq[atom1].x - posq[atom2].x, 2) + pow(posq[atom1].y - posq[atom2].y, 2) + pow(posq[atom1].z - posq[atom2].z, 2);

            printf("atom1 %d atom2 %d\n", atom1, atom2);
            if (d_sq < dist_sq && dev_bonded_exclusions_matrix[atom1 * num_atoms + atom2] == false)
            {
                if (counter >= 100)
                    printf("counter is greater than 100!!!\n");

                dev_neighbors[atom1 * 100 + counter++] = make_int2(atom1, atom2);
            }
        }
    }
*/
    /* 
     * We want to make sure any spots not used will not be confused for atoms
     * therefore we set them to -1, because that is not a valid atom number.
     */

  /*  while(counter < num_atoms)
        dev_neighbors[block * 100 + counter++] = make_int2(-1, -1);
    
}
*/