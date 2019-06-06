/*
# Edge Base Approach
#{class} {use the bit manipulations for less memeory requirements}
@Working forbidden modified if condition revome cudaSynchronize in function
@Algo 4-opt Modified on 20-03-2019
@Add Degree Base
*/

#include<bits/stdc++.h>
#include<cuda.h>
#include<thrust/count.h>
#include<thrust/extrema.h>
#include<thrust/device_ptr.h>
#include<curand_kernel.h>
#include "Utility.cuh"

#define MAXBLOCKS 1<<32
#define MOD 32

using namespace std;

__global__ void AssignColors(long int V,long int *delta_Degree,long int *colIndex,long int *colors,
                             long int *CS,unsigned long long int *vforbidden){
  
   long int threadId = blockDim.x*blockIdx.x+threadIdx.x;

   if(threadId<V){
      
      //conflicts colors
      if(colors[threadId]==0){
        
        //First Available Color
        if(vforbidden[threadId]==0){

           colors[threadId] = CS[threadId]+1;
        
        }else{
           CS[threadId] = CS[threadId]+1;
           vforbidden[threadId] = 0;
        }

      }

   }


}

__global__ void DetectConflicts(long int V,long int *preSum,long int *colIndex,long int *colors,
                bool *checkConflict,unsigned long long int *vforbidden,long int *deltaDegree,long int *degree){

    long int threadId = blockIdx.x*blockDim.x+threadIdx.x;

    if(threadId<V){

        if(colors[threadId]==0){
           *checkConflict = true;
        }

        for(long int k=preSum[threadId];k<preSum[threadId+1];k++){

            long int j = colIndex[k];

            if((colors[threadId]!=0)&&(colors[j]!=0)&&(colors[threadId]==colors[j])&&(degree[j]>degree[threadId])){

                colors[threadId] = 0;
                *checkConflict = true;
                return;

            }else if((colors[threadId]!=0)&&(colors[j]!=0)&&(colors[threadId]==colors[j])&&(j>threadId)){

                colors[threadId] = 0;
                *checkConflict = true;
                return;

            }
        }
    
    }
  
}

__global__ void ForbiddenColors(long int V,long int *preSum,long int *colIndex,long int *colors,
                                unsigned long long int *vforbidden,long int *CS){
 
   long int threadId = blockDim.x*blockIdx.x+threadIdx.x;

   if(threadId<V&&colors[threadId]==0){
      
      for(long int k=preSum[threadId];k<preSum[threadId+1];k++){
          
          long int j = colIndex[k];

          if(CS[j]==CS[threadId]){
             
             if(colors[j]!=0&&colors[threadId]==0){
                unsigned long long int value = (vforbidden[threadId]|colors[j])-vforbidden[threadId];
                atomicAdd(&vforbidden[threadId],value);
             }else if(colors[j]==0&&colors[threadId]!=0){
                unsigned long long int value = (vforbidden[j]|colors[threadId])-vforbidden[j];
                atomicAdd(&vforbidden[j],value);
             }
          }        

      }
   
   }
  
}

void preSumLength(int V,long int *d_preSum,long int *degree,long int *delta_Degree){
   
    for(long int i=0;i<V;i++){
       
       d_preSum[i+1] = d_preSum[i]+degree[i];
       
       if(*delta_Degree<degree[i]){
          *delta_Degree = degree[i];
       }
    }

    *delta_Degree = *delta_Degree+1;
}

__global__ void IsValidgraph_Coloring(long int V,long int *colors,long int *preSum,long int *colIndex,bool *flag){
   
   long int threadId = blockDim.x*blockIdx.x+threadIdx.x;
   if(threadId<V){
       for(long int i=preSum[threadId];i<preSum[threadId+1];i++){
          if(colors[threadId]==colors[colIndex[i]]||colors[threadId]==-1){
             *flag = false;
          }
       }
   }
   
}


long int EdgeBased_Algorithm(long int V,long int *preSum,long int *colIndex,long int *colors,long int *degree,long int n_zero_counter,long int *delta_Degree){
    
    long int *CS;
    cudaMallocManaged(&CS,sizeof(long int)*V);
    unsigned long long int *d_vforbidden;
    cudaMallocManaged(&d_vforbidden,sizeof(unsigned long long int)*V);
    thrust::fill(d_vforbidden,d_vforbidden+V,0); 
    thrust::fill(CS,CS+V,0);

    /*
    @ step 2 Initialize the colors to 0
    @ until all are colored
    */
    thrust::fill(colors,colors+V,0);
    
    long int minimumColor;
    long int n_threads =  256;
    long int n_blocks = min((V+n_threads-1)/n_threads,(long)MAXBLOCKS);
    bool *checkConflict;
    cudaMallocManaged(&checkConflict,sizeof(bool));
    clock_t gpu_time = clock();

    do{
       
       *checkConflict = false;

       AssignColors<<<n_blocks,n_threads>>>(V,delta_Degree,colIndex,colors,CS,d_vforbidden);
       DetectConflicts<<<n_blocks,n_threads>>>(V,preSum,colIndex,colors,checkConflict,d_vforbidden,delta_Degree,degree);
       ForbiddenColors<<<n_blocks,n_threads>>>(V,preSum,colIndex,colors,d_vforbidden,CS);
       catchCudaError(cudaMemcpy(checkConflict,checkConflict,sizeof(bool),cudaMemcpyDeviceToHost),"conflict");
       
    }while(*checkConflict);
    gpu_time = clock()-gpu_time;

    //Assigned Colors
    /*
    @ last step to print the assigned colors
    */
    cout<<endl;
    for(long int i=0;i<V;i++){
       printf("vertex --> %i Assigned Color --> %d\n",i,colors[i]);
    }
    cout<<endl;
    printf("gpu_time is %.6lf\n",(double)gpu_time/((double)CLOCKS_PER_SEC/1000));
    
    thrust::device_ptr<long int> d_ptr = thrust::device_pointer_cast(colors);
    minimumColor = *(thrust::max_element(d_ptr, d_ptr+V));

    cudaFree(checkConflict);
    cudaFree(d_vforbidden);
    cudaFree(CS);

    //required colors needed
    return minimumColor;
}

void GraphColoring_GPUAllocation(const char filename[]){
   
   //@difficult to allocate memory for large complete dataset not assume complete graph
   long int V; //No. of verties
   long int n_zero_counter = 0;   
   long int **st_Column;
   long int *st_degree;
   
   if(string(filename).find("col")!=string::npos){
     ReadColFile(filename,&V,&st_Column,&st_degree,&n_zero_counter);
   }else{
     ReadMMFile(filename,&V,&st_Column,&st_degree,&n_zero_counter); 
   }
   
   long int *degree;
   catchCudaError(cudaMallocManaged(&degree,sizeof(long int)*V),"Degree Allocation");
   thrust::copy(st_degree,st_degree+V,degree);
   
   long int *d_preSum;
   catchCudaError(cudaMallocManaged(&d_preSum,sizeof(long int)*(V+1)),"preSum Allocation");
   d_preSum[0] = 0;
   //store all the index of non zero element
   long int *d_colIndex;
   catchCudaError(cudaMallocManaged(&d_colIndex,sizeof(long int)*n_zero_counter),"colIndex Allocation");
   
   //Allocation
   long int *colors;
   catchCudaError(cudaMallocManaged(&colors,sizeof(long int)*V),"Color Allocation");
   long int *delta_Degree;
   catchCudaError(cudaMallocManaged(&delta_Degree,sizeof(long int)),"Delta Degree Allocation");
   *delta_Degree = 0;
   
   preSumLength(V,d_preSum,degree,delta_Degree);

   for(int i=0;i<V;i++){
      //Remove the cudaMemcpy it will take more time
      thrust::copy(st_Column[i],st_Column[i]+degree[i],d_colIndex+d_preSum[i]);
   }
   
   
   /*
   @begin CSR
   */
   long int width=16,height=16;
   long int threads_per_blocks = width*height;
   
   //Call the EdgeBase Algorithm
   long int number_Of_Colors_Needed = EdgeBased_Algorithm(V,d_preSum,d_colIndex,colors,degree,n_zero_counter,delta_Degree);  
   
   cout<<"EdgeBase Algorithm coloring found solution with "<<number_Of_Colors_Needed<<" colors"<<endl;
   cout<<"Valid coloring "; 
   
   bool *d_isValidColors; 
   catchCudaError(cudaMallocManaged(&d_isValidColors,sizeof(bool)),"IsValid Allocation");
   *d_isValidColors = true;

   IsValidgraph_Coloring<<<ceil(V/threads_per_blocks)+1,threads_per_blocks>>>(V,colors,d_preSum,d_colIndex,d_isValidColors); 
   catchCudaError(cudaMemcpy(d_isValidColors,d_isValidColors,sizeof(bool),cudaMemcpyDeviceToHost),"IsValid Host Copy");
   
   if(*d_isValidColors){
     cout<<"yes"<<endl;
   }else{
     cout<<"No"<<endl;
   }
   
   //catchCudaError(cudaDeviceSynchronize(),"GraphColoring DeviceSync");
   
   cudaFree(d_preSum);
   cudaFree(d_colIndex);
   cudaFree(colors);
   cudaFree(degree);
   cudaFree(d_isValidColors);
   
}

/* Reading Argument with command line opetion */
int main(int argc,char *argv[])
{
     if(argc<2){
       cout<<"Invalid Input Parameter"<<endl;
       exit(1);
     }else{
      
     /*
     @Adding the clock
     */
     clock_t time = clock();
     GraphColoring_GPUAllocation(argv[1]);
     time = clock()-time; 
     
     cout<<"Total execution time is "<<(double)time/(double)CLOCKS_PER_SEC<<endl;
    
   }   

   return 0;
}
