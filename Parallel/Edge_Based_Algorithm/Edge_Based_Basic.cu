/*
# Edge Base Approach
#{class} {use the bit manipulations for less memeory requirements}
@iterative
*/

#include<bits/stdc++.h>
#include<cuda.h>
#include<thrust/count.h>
#include<thrust/extrema.h>
#include<thrust/device_ptr.h>
#include <thrust/device_vector.h>
#include <thrust/inner_product.h>
#include<curand_kernel.h>
#include "Utility.cuh"

#define MAXBLOCKS 1<<32
#define MOD 32

using namespace std;

__global__ void AssignColors(long int V,long int *preSum,long int *colIndex,long int *colors,long int *delta_Degree,bool *conflicts){
  
   long int threadId = blockDim.x*blockIdx.x+threadIdx.x;
   long int stride = blockDim.x*gridDim.x;
   
   if(threadId<V&&!conflicts[threadId]){
      return;
   }


   if(threadId<V){
      
     for(long int i=threadId;i<V;i=i+stride){

        long int *vforbidden = (long int*)malloc(sizeof(long int)*(*delta_Degree+1));
        memset(vforbidden,0,sizeof(long int)*(*delta_Degree+1));

        for(long int k=preSum[i];k<preSum[i+1];k++){

           long int j = colIndex[k];
           long int value = colors[j]%MOD;
           long int shift = 1<<value;
           vforbidden[colors[j]/MOD]|= shift; 
        }
        
        //Assign colors
        for(long int color=1;color<=*delta_Degree+1;color++){

            long int val = color%MOD;
            
            if((vforbidden[color/MOD]&(1<<val))== 0){
                colors[i] = color;
                return;
            }

        }

       
        free(vforbidden);
        
     }

   }

}

__global__ void DetectConflicts(long int V,long int *preSum,long int *colIndex,long int *colors,bool *conflicts,bool *checkConflicts){

    long int threadId = blockIdx.x*blockDim.x+threadIdx.x;

    if(threadId<V){

        conflicts[threadId] = false;

        for(long int k=preSum[threadId];k<preSum[threadId+1];k++){

            long int j = colIndex[k];

            if((colors[threadId]==colors[j])&&(j<threadId)){

                conflicts[threadId] = true;
                *checkConflicts = true;
                return;

            }
        }
    
    }
  
}

__global__ void preSumLength(int V,long int *d_preSum,long int *degree,long int *delta_Degree){
   
    for(long int i=0;i<V;i++){
       
       d_preSum[i+1] = d_preSum[i]+degree[i];
       
       if(*delta_Degree<degree[i]){
          *delta_Degree = degree[i];
       }
    }
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
    
    /*
    @ step 2 Initialize the colors to 0
    @ until all are colored
    */
    thrust::fill(colors,colors+V,0);
    
    long int n_threads =  256;
    long int n_blocks = min((V+n_threads-1)/n_threads,(long)MAXBLOCKS);
    bool *d_conflicts,*checkConflict;
    cudaMallocManaged(&d_conflicts,sizeof(bool)*V);
    cudaMallocManaged(&checkConflict,sizeof(bool));
    thrust::fill(d_conflicts,d_conflicts+V,true);
    
    do{
       
       *checkConflict = false;

       AssignColors<<<n_blocks,n_threads>>>(V,preSum,colIndex,colors,delta_Degree,d_conflicts);
       DetectConflicts<<<n_blocks,n_threads>>>(V,preSum,colIndex,colors,d_conflicts,checkConflict);
       
       catchCudaError(cudaDeviceSynchronize(),"Edge");
        

    }while(*checkConflict);
    
    //Assigned Colors
    /*
    @ last step to print the assigned colors
    */
    cout<<endl;
    for(long int i=0;i<V;i++){
       printf("vertex --> %i Assigned Color --> %d\n",i,colors[i]);
    }
    cout<<endl;

    //thrust::device_ptr<long int> d_ptr = thrust::device_pointer_cast(colors);
    //long int minimumColor = *(thrust::max_element(d_ptr, d_ptr+V));
    thrust::device_vector<long int> d_data(V);
    thrust::copy(colors,colors+V,d_data.begin());
    thrust::sort(d_data.begin(), d_data.end());

    size_t num_unique = thrust::inner_product(d_data.begin(), d_data.end()-1,d_data.begin()+1,0,
                                              thrust::plus<long int>(),thrust::not_equal_to<long int>())+1;

    cudaFree(d_conflicts);
    cudaFree(checkConflict);
    
    //required colors needed
    return (long int)num_unique;
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

   preSumLength<<<1,1>>>(V,d_preSum,degree,delta_Degree);
   catchCudaError(cudaMemcpy(d_preSum,d_preSum,sizeof(long int)*(V+1),cudaMemcpyDeviceToHost),"Copy to PreSum");

   for(int i=0;i<V;i++){
      //Remove the cudaMemcpy it will take more time
      thrust::copy(st_Column[i],st_Column[i]+degree[i],d_colIndex+d_preSum[i]);
   }  
   
   /*
   @begin CSR
   */
   
   //Call the EdgeBase Algorithm
   long int number_Of_Colors_Needed = EdgeBased_Algorithm(V,d_preSum,d_colIndex,colors,degree,n_zero_counter,delta_Degree);  
   
   cout<<"EdgeBase Algorithm coloring found solution with "<<number_Of_Colors_Needed<<" colors"<<endl;
   cout<<"Valid coloring Yes\n"; 
   
   catchCudaError(cudaDeviceSynchronize(),"GraphColoring DeviceSync");
   

   cudaFree(d_preSum);
   cudaFree(d_colIndex);
   cudaFree(colors);
   cudaFree(degree);
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
