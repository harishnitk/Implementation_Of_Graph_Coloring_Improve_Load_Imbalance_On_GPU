#ifndef UTILITY_CUDA_H_
#define UTILITY_CUDA_H_

#include<bits/stdc++.h>
#include "cuda.h"
using namespace std;

void catchCudaError(cudaError_t, const char *);
void ReadColFile(const char[],long int *,long int ***,long int **,long int *);
void ReadMMFile(const char[],long int *,long int ***,long int **,long int *);
#endif
