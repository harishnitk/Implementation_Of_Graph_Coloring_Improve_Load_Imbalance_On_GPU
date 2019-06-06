#include<bits/stdc++.h>
#include "cuda.h"
using namespace std;

void catchCudaError(cudaError_t error, const char *function)
{
    if(error!=cudaSuccess)
    {
        printf("\n====== Cuda Error Code %i ======\n %s in CUDA %s\n", error, cudaGetErrorString(error), function);
        exit(-1);
    }
}

void ReadColFile(const char filename[],long int *V,long int ***st_Column,long int **st_degree,long int *counter){
   string s;
   ifstream infile(filename);
   if(infile.fail()){
      cout<<"Fail to open the file\n";
      exit(0);
   }

   long int n_rows,n_edges;
   //Maintain Hash for NNZ,preSum,colIndex
   //allocate dynamic size

   while(getline(infile,s)){
       istringstream iss(s);
       string str;
       long int u,v;
       iss>>str;
       if(str=="p"){
          iss>>s;
          iss>>n_rows;
          iss>>n_edges;
          *V = n_rows;
          *st_degree = new long int[n_rows];
          *st_Column = new long int*[n_rows];
          memset(*st_degree,0,n_rows*sizeof(long int));
          continue;
       }else if(str!="e"){
          continue; 
       }

       iss>>u>>v;
         if(u!=v){
         long int u_len = (*st_degree)[u-1];
         long int v_len = (*st_degree)[v-1];
         (*st_Column)[u-1] = (long int*)realloc((*st_Column)[u-1],sizeof(long int)*(u_len+1));
         (*st_Column)[v-1] = (long int*)realloc((*st_Column)[v-1],sizeof(long int)*(v_len+1));
         (*st_Column)[u-1][u_len] = v-1; 
         (*st_Column)[v-1][v_len] = u-1;
         (*st_degree)[u-1]++;
         (*st_degree)[v-1]++;
         *counter+=2;
       }
        
   }

   infile.close();
}

void ReadMMFile(const char filename[], long int *V,long int ***st_Column,long int **st_degree,long int *counter){
   string s;
   ifstream infile(filename);
   if(infile.fail()){
      cout<<"Fail to open the file\n";
      return;
   }

   //content
   while(getline(infile,s)){
     istringstream iss(s);
     if(s.find("%")==string::npos){
        break;
     }
   }

   istringstream iss(s);
   //Maintain Hash for NNZ,preSum,colIndex
   //allocate dynamic size

   long int n_rows,n_cols,n_edges;
   iss>>n_rows>>n_cols>>n_edges;
   *st_degree = new long int[n_rows];
   *st_Column = new long int*[n_rows];
   memset(*st_degree,0,n_rows*sizeof(long int));
   *V = n_rows;

   //reading edges

   while(getline(infile,s)){
      istringstream iss(s);
      long int u,v,w;
      iss>>u>>v>>w;
      if(u!=v){
        long int u_len = (*st_degree)[u-1];
        long int v_len = (*st_degree)[v-1];
        (*st_Column)[u-1] = (long int*)realloc((*st_Column)[u-1],sizeof(long int)*(u_len+1));
        (*st_Column)[v-1] = (long int*)realloc((*st_Column)[v-1],sizeof(long int)*(v_len+1));
        (*st_Column)[u-1][u_len] = v-1; 
        (*st_Column)[v-1][v_len] = u-1;
        (*st_degree)[u-1]++;
        (*st_degree)[v-1]++;
        *counter+=2;
      }
   }

   infile.close();
}
