#ifndef UTILITY_H_INCLUDED
#define UTILITY_H_INCLUDED
#include<bits/stdc++.h>
using namespace std;

void ReadColFile(const char[],long int *,long int ***,long int **,long int *);
void ReadMMFile(const char[],long int *,long int ***,long int **,long int *);
int CountColors(long int ,long int *);
bool IsValidColoring(long int,long int *,long int *,long int *);
bool IsSafe(long int,bool *,long int,long int *, long int);

#endif // UTILITY_H_INCLUDED
