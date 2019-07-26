extern "C" __global__ void reduceEnergy(double* __restrict__ dev_energy, double* __restrict__ dev_output_energy) {
    extern __shared__ mixed tempBuffer[];
    const unsigned int thread = threadIdx.x;
    mixed sum = 0;
    for (unsigned int index = thread; index < bufferSize; index += blockDim.x)
        sum += dev_energy[index];
    tempBuffer[thread] = sum;
    for (int i = 1; i < workGroupSize; i *= 2) {
        __syncthreads();
        if (thread%(i*2) == 0 && thread+i < workGroupSize)
            tempBuffer[thread] += tempBuffer[thread+i];
    }
    if (thread == 0)
        *dev_output_energy = tempBuffer[0];
}











