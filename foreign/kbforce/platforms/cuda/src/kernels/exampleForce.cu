extern "C" __global__ void calcNonBonded(real4* __restrict__ posq, unsigned long long* __restrict__ forceBuffers, int *particle1, int *particle2, int *dev_type1, int *dev_type2, double *dev_potential_lookup, double *dev_derivative_lookup, double *dev_energy) {

    for (int atom = blockIdx.x*blockDim.x+threadIdx.x; atom < num_bonds; atom += blockDim.x*gridDim.x) {
       
      
        
        int p1 = particle1[atom];
        int p2 = particle2[atom];
        
        real4 pos1 = posq[p1];
        real4 pos2 = posq[p2];
        
        
   
        real3 delta = make_real3(pos2.x-pos1.x, pos2.y-pos1.y, pos2.z-pos1.z);
        real r = SQRT(delta.x*delta.x + delta.y*delta.y + delta.z*delta.z);
        size_t dist = static_cast<size_t>(floor(r / step));
      
        
        const int internal_min = (dev_type1[atom] <= dev_type2[atom]) ? dev_type1[atom] : dev_type2[atom];
        const int internal_max = (dev_type1[atom] <= dev_type2[atom]) ? dev_type2[atom] : dev_type1[atom];
    
        const int offset = (internal_max + internal_min*(number_of_types-1)-internal_min*(internal_min-1)/2) * number_of_steps;

        if (dist >= number_of_steps) 
            continue; // effectively add zero to energy and force
            
        atomicAdd(&dev_energy[0], dev_potential_lookup[offset + dist]);
         
         
        real dEdR = dev_derivative_lookup[offset + dist];
        dEdR = (r > 0) ? (dEdR/r) : 0.0;
        
        real3 force1 = make_real3(delta.x * dEdR, delta.y * dEdR, delta.z * dEdR);
        real3 force2 = make_real3(-delta.x * dEdR, -delta.y * dEdR, -delta.z * dEdR);
        
        
        atomicAdd(&forceBuffers[p1], static_cast<unsigned long long>((long long) (force1.x *0x100000000)));
        atomicAdd(&forceBuffers[p2], static_cast<unsigned long long>((long long) (force2.x *0x100000000)));
        
        atomicAdd(&forceBuffers[p1 + PADDED_NUM_ATOMS], static_cast<unsigned long long>((long long) (force1.y *0x100000000)));
        atomicAdd(&forceBuffers[p2 + PADDED_NUM_ATOMS], static_cast<unsigned long long>((long long) (force2.y *0x100000000)));
        
        atomicAdd(&forceBuffers[p1 + PADDED_NUM_ATOMS*2], static_cast<unsigned long long>((long long) (force1.z *0x100000000)));
        atomicAdd(&forceBuffers[p2 + PADDED_NUM_ATOMS*2], static_cast<unsigned long long>((long long) (force2.z *0x100000000)));
    
         
    }
}
    
