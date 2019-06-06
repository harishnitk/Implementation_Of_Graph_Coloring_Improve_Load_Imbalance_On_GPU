/*
Need to add virtual wrap centric programming
modified on 4-may-2019, try random bit selection
*/
#include<bits/stdc++.h>
#include<cuda.h>
#include<thrust/count.h>
#include <thrust/device_vector.h>
#include <thrust/inner_product.h>
#include<curand_kernel.h>
#include<thrust/extrema.h>
#include<thrust/device_ptr.h>
#include "Utility.cuh"

#define MAXBLOCKS 1<<30
#define MAXDELTA 1000000000

using namespace std;

/*
[method] {find the unsetRandomColor}
[]clear the neighbour bit
*/
__device__ long int findUnsetRandomColor(long int Id,unsigned int *Colorset,long int maxClr){
    
    long int j,phase,end=(maxClr/32)+1,bit,left=maxClr%32;
    
    for(phase=0;phase<end;phase++){
      
      if(Id%2==0){
        bit = __ffs(Colorset[Id*end+phase]);
      }else{
         if(phase==end-1&&left!=0){
           unsigned int Value = 0;
           j = left-1;
           while(Colorset[Id*end+phase]>0){
             Value = Value+(Colorset[Id*end+phase])%2*(1<<j);
             j--;
             Colorset[Id*end+phase] = Colorset[Id*end+phase]/2;
           }
           bit = __ffs(Value);
           bit = bit!=0?left+1-bit:bit;
         }else{
            Colorset[Id*end+phase] = __brev(Colorset[Id*end+phase]);
            bit = __ffs(Colorset[Id*end+phase]);
            bit = bit!=0?33-bit:bit;
         }
      }
      
      if(bit!=0){
        return phase*32+bit;
      }

    }
    
    return -1;

}

__global__ void findtheNeighbourColors(long int V,long int *colors,long int *preSum,long int 
                                       *colIndex,unsigned int *Colorset,long int maxClr){
    
    long int threadId = blockDim.x*blockIdx.x+threadIdx.x;

    if(threadId<V){

        long int k=preSum[threadId],n=preSum[threadId+1],j,phase,end=(maxClr/32)+1,left=maxClr%32;
        Colorset[threadId*end+end-1] = left!=0?(1<<left)-1:4294967295;//32 bit number
        
        for(phase=0;phase<end-1;phase++){
          Colorset[threadId*end+phase] = 4294967295;
        }
        
        for(;k<n;k++,n--){
           j = colIndex[k];
           
           if(colors[j]!=0){
              phase = colors[j]/32;
              if(colors[j]%32==0){
                phase--;
              }
              Colorset[threadId*end+phase]&= ~(unsigned int)(1<<(colors[j]-1-phase*32));//clear bit
           }

           j = colIndex[n-1];
           if(colors[j]!=0){
              phase = colors[j]/32;
              if(colors[j]%32==0){
                phase--;
              }
              Colorset[threadId*end+phase]&= ~(unsigned int)(1<<(colors[j]-1-phase*32));//clear bit
           }

        }
    }

}

/*
[method] {AssignColor}
{Description} []
*/
__global__ void assignColors(long int V,long int *colors,long int *colIndex,long int *preSum,unsigned int 
                             *Colorset,long int *degree,bool *inc,long int maxClr){
   
    long int threadId = blockDim.x*blockIdx.x+threadIdx.x;

    if(threadId<V&&colors[threadId]==0){
           
       //Select first available colors Intimation of available color which is minimum
       if(degree[threadId]==0){
         colors[threadId] = 1; //{to reduce the colors}
       }else{
          //choose color form from its degree
          colors[threadId] = findUnsetRandomColor(threadId,Colorset,maxClr); 

          if(colors[threadId]==-1){
             *inc = true;
             colors[threadId] = 0;
          }
       }        

    }

}

/*
@reset colors to 0
*/
__global__ void DetectConflictsColors(long int V,long int *preSum,long int *colIndex,long int *colors,
                                      long int *degree,bool *checkConflict,int isDegreeChange){
  
    long int threadId = blockDim.x*blockIdx.x+threadIdx.x;
    
    if(threadId<V){

      if(colors[threadId]==0){
        *checkConflict = true;
        return;
      }
    
      for(long int k=preSum[threadId];k<preSum[threadId+1]&&colors[threadId]!=0;k++){
         
         long int j = colIndex[k];
         
         if(colors[j]==0){
            continue;
         }

         if(isDegreeChange){
            if((colors[threadId]==colors[j])&&(degree[threadId]>degree[j])){
              colors[j] = 0;
              *checkConflict = true;
              break;
            }else if((colors[threadId]==colors[j])&&(degree[threadId]<degree[j])){
              colors[threadId] = 0;
              *checkConflict = true;
              break;
            }else if((colors[threadId]==colors[j])&&(threadId>j)){
              colors[j] = 0;
              *checkConflict = true;
              break;
            }else if((colors[threadId]==colors[j])&&(threadId<j)){
              colors[threadId] = 0;
              *checkConflict = true;
              break;
            }
         }else if((colors[threadId]==colors[j])&&(threadId>j)){
              colors[j] = 0;
              *checkConflict = true;
              break;
         }else if((colors[threadId]==colors[j])&&(threadId<j)){
              colors[threadId] = 0;
              *checkConflict = true;
              break;
         }

      }

    }

}

__global__ void IsValidgraph_Coloring(long int V,long int *colors,long int *preSum,long int *colIndex,
                                      bool *flag){
   
   long int threadId = blockDim.x*blockIdx.x+threadIdx.x;
   if(threadId<V){
       for(long int i=preSum[threadId];i<preSum[threadId+1];i++){
          if(colors[threadId]==colors[colIndex[i]]||colors[threadId]<=0){
             *flag = false;
          }
       }
   }
   
}

void preSumLength(int V,long int *d_preSum,long int *degree,long int &deltaDegree,long int &isDegreeChange){
    
    for(long int i=0;i<V;i++){
       d_preSum[i+1] = d_preSum[i]+degree[i];
       if(deltaDegree<degree[i]){
          deltaDegree = degree[i];
       }

       if((i+1)<V&&(degree[i]!=degree[i+1])){
         isDegreeChange = 1;
       }
    }

    deltaDegree+= 1;

}

long int ConflictColorAlgorithm(long int V,long int *preSum,long int *colIndex,long int *colors,
                                long int *degree,int deltaDegree,int isDegreeChange){
    
    double *d_random;
    cudaMallocManaged(&d_random,sizeof(double)*V);
    unsigned int *d_Colorset;
    bool *d_checkConflict,checkConflict,*d_inc,inc;

    //Allocate the memory for maxArray
    catchCudaError(cudaMallocManaged(&d_checkConflict,sizeof(int)),"checkConflict Allocation");
    cudaMallocManaged(&d_Colorset,sizeof(unsigned int)*(min((long)MAXDELTA,(long)V*(deltaDegree/32+1))));
    catchCudaError(cudaMallocManaged(&d_checkConflict,sizeof(bool)),"cc");
    catchCudaError(cudaMallocManaged(&d_inc,sizeof(bool)),"inc");

    long int n_threads =  256;
    long int n_blocks = min((V+n_threads-1)/n_threads,(long)MAXBLOCKS);
    long int maxClr = 2;
    
    /*
    @ step 2 Initialize the colors to 0
    @ until all are colored
    */
    thrust::fill(colors,colors+V,0);
    
    clock_t gpu_time = clock();
    do{
         
         checkConflict = false;
         inc = false;
         catchCudaError(cudaMemcpy(d_checkConflict,&checkConflict,sizeof(bool),cudaMemcpyHostToDevice),"cc");
         catchCudaError(cudaMemcpy(d_inc,&inc,sizeof(bool),cudaMemcpyHostToDevice),"inc copy");
         
         findtheNeighbourColors<<<n_blocks,n_threads>>>(V,colors,preSum,colIndex,d_Colorset,maxClr);
         assignColors<<<n_blocks,n_threads>>>(V,colors,colIndex,preSum,d_Colorset,degree,d_inc,maxClr);
         DetectConflictsColors<<<n_blocks,n_threads>>>(V,preSum,colIndex,colors,degree,d_checkConflict,isDegreeChange);       
         
         catchCudaError(cudaMemcpy(&inc,d_inc,sizeof(bool),cudaMemcpyDeviceToHost),"inc copy");
         catchCudaError(cudaMemcpy(&checkConflict,d_checkConflict,sizeof(bool),cudaMemcpyDeviceToHost),"cc");
         
         if(inc){
            maxClr = 2*maxClr;
         }

    
    }while(checkConflict);
    gpu_time = clock()-gpu_time;
    //Assigned Colors
    /*
    @ last step to print the assigned colors
    */
    printf("\n");
    for(long int i=0;i<V;i++){
       printf("vertex --> %i Assigned Color --> %d\n",i,colors[i]);
    }
    printf("\n");
    printf("gpu_time is %.6lf\n",(double)gpu_time/((double)CLOCKS_PER_SEC/1000));

    thrust::device_vector<long int> d_data(V);
    thrust::copy(colors,colors+V,d_data.begin());
    thrust::sort(d_data.begin(), d_data.end());

    size_t num_unique = thrust::inner_product(d_data.begin(), d_data.end()-1,d_data.begin()+1,0,
                                              thrust::plus<long int>(),thrust::not_equal_to<long int>())+1;
    

    cudaFree(d_random);
    cudaFree(d_checkConflict);
    cudaFree(d_inc);
    cudaFree(d_Colorset);

    //required colors needed
    return (long int)num_unique;
}

void GraphColoring_GPUAllocation(const char filename[]){
   
   //@difficult to allocate memory for large complete dataset not assume complete graph
   long int V,deltaDegree = 0,isDegreeChange = 0; //No. of verties
   long int n_zero_counter = 0;   
   long int **st_Column;
   long int *st_degree;
   
   if(string(filename).find("col")!=string::npos||string(filename).find("clq")!=string::npos){
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

   preSumLength(V,d_preSum,degree,deltaDegree,isDegreeChange);

   for(int i=0;i<V;i++){
      //Remove the cudaMemcpy it will take more time
      thrust::copy(st_Column[i],st_Column[i]+degree[i],d_colIndex+d_preSum[i]);
   }  
   
   /*
   @begin CSR
   */
   
   //Call the Color Conflict Algorithm
   long int number_Of_Colors_Needed = ConflictColorAlgorithm(V,d_preSum,d_colIndex,colors,degree,deltaDegree,isDegreeChange);  
   
   printf("ConflictColor Algorithm coloring found solution with %d colors\n",number_Of_Colors_Needed);
   printf("Valid coloring Yes\n"); 
   
   //catchCudaError(cudaDeviceSynchronize(),"GraphColoring DeviceSync");
   
   cudaFree(d_preSum);
   cudaFree(d_colIndex);
   cudaFree(colors);
   cudaFree(degree);
  
}

/* Reading Argument with command line opetion */
int main(int argc,char *argv[])
{
     if(argc<2){
       printf("Invalid Input Parameter\n");
       exit(1);
     }else{

     /*
     @Adding the clock
     */
     clock_t time = clock();
     GraphColoring_GPUAllocation(argv[1]);
     time = clock()-time; 
     
     double execution_time = (double)time/(double)CLOCKS_PER_SEC;
     
     printf("Total execution time is %lf\n",execution_time);
    
   }   

   return 0;
}
