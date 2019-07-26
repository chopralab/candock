{
    real3 delta = posq2 - posq1;
    //real r = sqrt(delta.x * delta.x + delta.y * delta.y + delta.z * delta.z);
    int dist = static_cast<int>(floor(r / step));

    //const int internal_min = dev_atoms[atom1] <= dev_atoms[atom2] ? dev_atoms[atom1] : dev_atoms[atom2];
    //const int internal_max = dev_atoms[atom1] <= dev_atoms[atom2] ? dev_atoms[atom2] : dev_atoms[atom1];
    const int internal_min = atom1 <= atom2 ? atom1 : atom2;
    const int internal_max = atom1 <= atom2 ? atom2 : atom1;


    const int offset = (internal_max + internal_min * (number_of_types - 1) - internal_min * (internal_min - 1) / 2) * number_of_steps;
            
           

    if (dist >= number_of_steps)
        continue; // effectively add zero to energy and force
            

    tempEnergy += dev_potential_lookup[offset + dist];

    dEdR += (r > 0) ? (dev_derivative_lookup[offset + dist] / r) : 0.0;


    real3 force1 = make_real3(delta.x * dEdR, delta.y * dEdR, delta.z * dEdR);
    real3 force2 = make_real3(-delta.x * dEdR, -delta.y * dEdR, -delta.z * dEdR);

    atomicAdd(&forceBuffers[atom1], static_cast<unsigned long long>((long long)(force1.x * 0x100000000)));
    atomicAdd(&forceBuffers[atom2], static_cast<unsigned long long>((long long)(force2.x * 0x100000000)));

    atomicAdd(&forceBuffers[atom1 + PADDED_NUM_ATOMS], static_cast<unsigned long long>((long long)(force1.y * 0x100000000)));
    atomicAdd(&forceBuffers[atom2 + PADDED_NUM_ATOMS], static_cast<unsigned long long>((long long)(force2.y * 0x100000000)));

    atomicAdd(&forceBuffers[atom1 + PADDED_NUM_ATOMS * 2], static_cast<unsigned long long>((long long)(force1.z * 0x100000000)));
    atomicAdd(&forceBuffers[atom2 + PADDED_NUM_ATOMS * 2], static_cast<unsigned long long>((long long)(force2.z * 0x100000000)));
            
   

}


