/*
# Edge Base Approach
#{class} {use the bit manipulations for less memeory requirements}
@working Latest more fast but give more color than previous
@Recolor function Added 21-03-2019, resolve Degree Based
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


__global__ void AssignColors(long int V,long int *delta_Degree,long int *colIndex,long int *colors,
                             unsigned int *CS,unsigned long long int *vforbidden){
  
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
                                bool *checkConflict,long int *degree,long int *modify){

    long int threadId = blockIdx.x*blockDim.x+threadIdx.x;

    if(threadId<V){
        
        modify[threadId] = 0;

        if(colors[threadId]==0){
           *checkConflict = true;
           return;
        }

        for(long int k=preSum[threadId];k<preSum[threadId+1];k++){

            long int j = colIndex[k];

            if(colors[j]==0){
               continue;
            }

            if((colors[threadId]==colors[j])&&degree[j]<degree[threadId]){

                colors[j] = 0;
                *checkConflict = true;

            }else if((colors[threadId]==colors[j])&&degree[j]>degree[threadId]){

                colors[threadId] = 0;
                *checkConflict = true;
                return;

            }else if((colors[threadId]==colors[j])&&j<threadId){

                colors[j] = 0;
                *checkConflict = true;

            }else if((colors[threadId]==colors[j])&&j>threadId){

                colors[threadId] = 0;
                *checkConflict = true;
                return;

            }
        }
    
    }
  
}

__global__ void ForbiddenColors(long int V,long int *preSum,long int *colIndex,long int *colors,
                                unsigned long long int *vforbidden,unsigned int *CS,int flag){
 
   long int threadId = blockDim.x*blockIdx.x+threadIdx.x;

   if(threadId<V){
      
      if(colors[threadId]){
         return;
      }  
      
      long int maxColor = LONG_MIN;
      
      for(long int k=preSum[threadId];k<preSum[threadId+1];k++){
          
          long int j = colIndex[k];
          
          if(colors[j]==0){
            continue;
          }
          
          if(CS[j]==CS[threadId]){
             
             if(colors[j]!=0){
                unsigned long long int value = (vforbidden[threadId]|colors[j])-vforbidden[threadId];
                atomicAdd(&vforbidden[threadId],value);
             }
          }

          if(colors[j]!=0&&maxColor<colors[j]){
             maxColor = colors[j];
          }        

      }
      
      if(flag==0){
        if((colors[threadId]==0)){//&&(vforbidden[threadId]==0)){
           CS[threadId] = maxColor;
        }
      }
   
   }
  
}

/*
@Larger degree first approch
*/
__global__ void resolveAdj_LargeIndex(long int V,long int *preSum,long int *colIndex,long int *colors,
                              unsigned long long int *vforbidden,unsigned int *CS,bool *change,int flag,long int *degree,long int *modify){

     long int node = blockDim.x*blockIdx.x+threadIdx.x;

     if(node<V){

        if(vforbidden[node]==0||colors[node]||modify[node]==1){
           return;
        }
        
        if(flag){
          
          long int larger=node,seclarger=node,start=preSum[node],end=preSum[node+1]-1,maxColor=LONG_MIN;

          for(;start<=end;){
              
              long int j = colIndex[start];
              
              if((CS[node]==CS[j])&&colors[j]==0&&degree[j]>degree[node]&&modify[j]==0){
                 
                 if(larger<j&&degree[larger]<degree[j]){
                    seclarger = larger;
                    larger = j;
                    *change = true;
                    CS[larger] = CS[larger]+1;
                    if(modify[seclarger]==0){
                      CS[seclarger] = CS[seclarger]+1;  
                    }
                    modify[seclarger] = 1;
                    modify[larger] = 1;
                 }else if(degree[seclarger]<degree[j]){
                    seclarger = j;
                    modify[seclarger] = 1;
                    CS[seclarger] = CS[seclarger]+1;
                 }
                 

              }else if(colors[j]!=0&&maxColor<colors[j]){
                 maxColor = colors[j];
              }

              j = colIndex[end];
              if((CS[node]==CS[j])&&colors[j]==0&&degree[j]>degree[node]&&modify[j]==0){
                 
                 if(larger<j&&degree[larger]<degree[j]){
                    seclarger = larger;
                    larger = j;
                    *change = true;
                    CS[larger] = CS[larger]+1;
                    if(modify[seclarger]==0){
                      CS[seclarger] = CS[seclarger]+1;  
                    }
                    modify[seclarger] = 1;
                    modify[larger] = 1;
                 }else if(degree[seclarger]<degree[j]){
                    seclarger = j;
                    modify[seclarger] = 1;
                    CS[seclarger] = CS[seclarger]+1;
                 }
                 

              }else if(colors[j]!=0&&maxColor<colors[j]){
                 maxColor = colors[j];
              }

              start++;
              end--;

          }

          if(larger!=seclarger){
             colors[larger] = maxColor+1;
             colors[seclarger] = colors[larger]+1;
          }else if(modify[node]==0){
             modify[node] = 1;
             CS[node] = CS[node]+1;
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
    
    unsigned int *CS;
    cudaMallocManaged(&CS,sizeof(unsigned int)*V);
    unsigned long long int *d_vforbidden;
    cudaMallocManaged(&d_vforbidden,sizeof(unsigned long long int)*V);
    thrust::fill(d_vforbidden,d_vforbidden+V,0); 
    thrust::fill(CS,CS+V,0);
    bool *d_change;
    int flag=1;
    cudaMallocManaged(&d_change,sizeof(bool));
    long int *d_modify;
    cudaMallocManaged(&d_modify,sizeof(long int)*V);
    
    /*
    @ step 2 Initialize the colors to 0
    @ until all are colored
    */
    thrust::fill(colors,colors+V,0);
    
    long int n_threads =  256;
    long int n_blocks = min((V+n_threads-1)/n_threads,(long)MAXBLOCKS);
    bool *checkConflict;
    cudaMallocManaged(&checkConflict,sizeof(bool));  
        
    clock_t gpu_time = clock();
    do{
       
       *checkConflict = false;
       *d_change = false;

       AssignColors<<<n_blocks,n_threads>>>(V,delta_Degree,colIndex,colors,CS,d_vforbidden);
       DetectConflicts<<<n_blocks,n_threads>>>(V,preSum,colIndex,colors,checkConflict,degree,d_modify);
       ForbiddenColors<<<n_blocks,n_threads>>>(V,preSum,colIndex,colors,d_vforbidden,CS,flag);
       catchCudaError(cudaMemcpy(checkConflict,checkConflict,sizeof(bool),cudaMemcpyDeviceToHost),"Conflict Copy");
  
       if(flag){
          
          resolveAdj_LargeIndex<<<n_blocks,n_threads>>>(V,preSum,colIndex,colors,d_vforbidden,CS,d_change,flag,degree,d_modify);
          catchCudaError(cudaMemcpy(d_change,d_change,sizeof(bool),cudaMemcpyDeviceToHost),"d_change Copy");
          
          if(*d_change==false){
             flag = 0;
          }

       }


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

    //thrust::device_ptr<long int> d_ptr = thrust::device_pointer_cast(colors);
    //long int minimumColor = *(thrust::max_element(d_ptr, d_ptr+V));
    thrust::device_vector<long int> d_data(V);
    thrust::copy(colors,colors+V,d_data.begin());
    thrust::sort(d_data.begin(), d_data.end());

    size_t num_unique = thrust::inner_product(d_data.begin(), d_data.end()-1,d_data.begin()+1,0,
                                              thrust::plus<long int>(),thrust::not_equal_to<long int>())+1;

    cudaFree(checkConflict);
    cudaFree(d_vforbidden);
    cudaFree(CS);

    //required colors needed
    return (long int)num_unique;//minimumColor;
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
   clock_t t = clock();
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
