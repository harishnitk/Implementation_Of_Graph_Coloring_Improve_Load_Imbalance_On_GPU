/*
@Optimization of Hybrid-6-opt, remove variables
Need to optimize
@Add Sorting,heapSort
*/
#include<bits/stdc++.h>
#include<cuda.h>
#include <thrust/device_vector.h>
#include <thrust/inner_product.h>
#include<thrust/count.h>
#include<curand_kernel.h>
#include<thrust/extrema.h>
#include<thrust/device_ptr.h>
#include "Utility.cuh"

#define MAXBLOCKS 1<<30

using namespace std;


//Heap Sort
__global__ void heapsortBasedOnDegree(long int V,long int *preSum,long int *colIndex,long int *degree){
  
    long int threadId = blockDim.x*blockIdx.x+threadIdx.x;
    
    if(threadId<V){

        long int n = preSum[threadId+1]-preSum[threadId];
        long int largest,previous,temp;
        
        for (int i = preSum[threadId]+n/2-1;i>=preSum[threadId];i--){ 
            
            largest = i,previous=-1;     
            while(largest!=previous){
              int l = 2*largest-preSum[threadId]+1;
              int r = 2*largest-preSum[threadId]+2;
              previous = largest;
              if(l<preSum[threadId+1]&&degree[colIndex[l]]<degree[colIndex[largest]]) 
                  largest = l; 
     
              if(r<preSum[threadId+1]&&degree[colIndex[r]]<degree[colIndex[largest]]) 
                largest = r; 
     
              if(largest!=previous) 
              { 
                temp = colIndex[previous];
                colIndex[previous] = colIndex[largest];
                colIndex[largest] = temp;
              }
           }
     
        }
 
        // One by one extract an element from heap 
        for (int i=preSum[threadId+1]-1;i>=preSum[threadId];i--) 
        { 
            temp = colIndex[preSum[threadId]];
            colIndex[preSum[threadId]] = colIndex[i];
            colIndex[i] = temp;
            largest = preSum[threadId],previous=-1;
            while(largest!=previous){
              int l = 2*largest-preSum[threadId]+1;
              int r = 2*largest-preSum[threadId]+2;
              previous = largest;
              if(l<i&&degree[colIndex[l]]<degree[colIndex[largest]]) 
                  largest = l; 
     
              if(r<i&&degree[colIndex[r]]<degree[colIndex[largest]]) 
                largest = r; 
     
              if(largest!=previous) 
              { 
                temp = colIndex[previous];
                colIndex[previous] = colIndex[largest];
                colIndex[largest] = temp;
             
              }
           }
     
        } 
    }

}

//Heap Sort
__global__ void heapsortBasedOnRandom(long int V,long int *preSum,long int *colIndex,double *random){
  
    long int threadId = blockDim.x*blockIdx.x+threadIdx.x;
    
    if(threadId<V){

        long int n = preSum[threadId+1]-preSum[threadId];
        long int largest,previous,temp;
        
        for (int i = preSum[threadId]+n/2-1;i>=preSum[threadId];i--){ 
            
            largest = i,previous=-1;     
            while(largest!=previous){
              int l = 2*largest-preSum[threadId]+1;
              int r = 2*largest-preSum[threadId]+2;
              previous = largest;
              if(l<preSum[threadId+1]&&random[colIndex[l]]<random[colIndex[largest]]) 
                  largest = l; 
     
              if(r<preSum[threadId+1]&&random[colIndex[r]]<random[colIndex[largest]]) 
                largest = r; 
     
              if(largest!=previous) 
              { 
                temp = colIndex[previous];
                colIndex[previous] = colIndex[largest];
                colIndex[largest] = temp;
              }
           }
     
        }
 
        // One by one extract an element from heap 
        for (int i=preSum[threadId+1]-1;i>=preSum[threadId];i--) 
        { 
            temp = colIndex[preSum[threadId]];
            colIndex[preSum[threadId]] = colIndex[i];
            colIndex[i] = temp;
            largest = preSum[threadId],previous=-1;
            while(largest!=previous){
              int l = 2*largest-preSum[threadId]+1;
              int r = 2*largest-preSum[threadId]+2;
              previous = largest;
              if(l<i&&random[colIndex[l]]<random[colIndex[largest]]) 
                  largest = l; 
     
              if(r<i&&random[colIndex[r]]<random[colIndex[largest]]) 
                largest = r; 
     
              if(largest!=previous) 
              { 
                temp = colIndex[previous];
                colIndex[previous] = colIndex[largest];
                colIndex[largest] = temp;
             
              }
           }
     
        } 
    }

}


__global__ void generateRandom(curandState* globalState, double* randomArray,long int V,unsigned long seed) 
{
    long int idx = blockIdx.x*blockDim.x+threadIdx.x;
    if(idx<V){
       curandState localState = globalState[idx];
       curand_init(seed,idx,0,&localState);
       double RANDOM = curand_uniform(&localState);
       randomArray[idx] = RANDOM;
       globalState[idx] = localState;
    }
}


__global__ void assignColors(long int V,long int c,double *random,long int *colors,double *max_Array,
                             int *left){
   
    long int threadId = blockDim.x*blockIdx.x+threadIdx.x;

    if(threadId<V){
        if(colors[threadId]==-1){
            if(random[threadId]>max_Array[threadId]){
                colors[threadId] = c;
            }else{
                *left = 1;
            }
        }
    }
}

__global__ void assign_DegreeColors(long int V,long int c,long int *degree,bool *flag,long int *colors,
                                    long int *max_Degree,bool *isDegree,int *left){
    
    long int threadId = blockDim.x*blockIdx.x+threadIdx.x;
    
    if(threadId<V){
        if(colors[threadId]==-1){
           if(flag[threadId]==false){
              if(degree[threadId]>max_Degree[threadId]){
                  colors[threadId] = c;
                  *isDegree = true;
                  //reAssign Values
              }else{
                *left = 1;
              } 
           }else{
              flag[threadId] = false;
              *left = 1;
           }
        }
    }
}

//Degree Coloring
__global__ void degree_Coloring(long int V,long int c,long int *preSum,long int *colIndex,
                                long int *degree,long int *colors,long int *max_Degree,bool *d_flag){
  
   long int threadId = blockDim.x*blockIdx.x+threadIdx.x;
   long int maxDegree = -1;

   if(threadId<V&&colors[threadId]!=-1){
      return;
   }
   else if(threadId<V){
            
       if(colors[threadId]==-1){
          
           maxDegree = -1;
          
           for(long int k=preSum[threadId];k<preSum[threadId+1];k++){
               long int j=colIndex[k];
               long int jc = colors[j];
               
               if(colors[j]!=-1){
                 continue;
               }
              
              //check all neighborhood which all are uncolored
              if(jc==-1){
                 if(degree[threadId]==degree[j]){
                    d_flag[threadId] = true;
                    break;  
                 }else if(degree[j]>maxDegree){
                   maxDegree = degree[j];
                   break;
                 }
              }
           }              
           
           max_Degree[threadId] = maxDegree;
       }
        
        
   }

}

// How to arrange the program to assign the different color to adjacent
__global__ void randomized_Based_Label(long int V,long int c,long int *preSum,long int *colIndex,
                                       double *random,long int *colors,double *max_Array){
     
     long int threadId = blockDim.x*blockIdx.x+threadIdx.x;
     double max_Value = -1;
     
     if(threadId<V&&colors[threadId]!=-1){
        return;
     }
     else if(threadId<V){
        
         if(colors[threadId]==-1){
          
           max_Value = -1;

           for(long int k=preSum[threadId];k<preSum[threadId+1];k++){
             long int j = colIndex[k];
             long int jc = colors[j];
             
             if(colors[j]!=-1){
               continue;
             }

             if((threadId!=j&&colors[threadId]==colors[j])||(jc==-1)){        
               //check with random value
                if(random[j]>max_Value){
                    max_Value = random[j];
                    break;
                }
             }
             
           }
       
           max_Array[threadId] = max_Value;

          }     
         
     }
} 

void preSumLength(int V,long int *d_preSum,long int *degree){
   
    for(long int i=0;i<V;i++){
       d_preSum[i+1] = d_preSum[i]+degree[i];
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


long int Hybrid_Algorithm(long int V,long int *preSum,long int *colIndex,long int *colors,long int *degree){
    
    curandState* devStates;
    double *d_random;
    cudaMallocManaged(&d_random,sizeof(double)*V);
    double *d_maxArray;
    long int *d_maxDegree;
    bool *d_flag;
    bool *d_hybrid,hybrid;
    
    //Allocate the memory for maxArray
    catchCudaError(cudaMallocManaged(&d_maxArray,sizeof(double)*V),"MaxArray Allocation");
    catchCudaError(cudaMallocManaged(&d_maxDegree,sizeof(long int)*V),"MaxDegree Allocation");
    catchCudaError(cudaMallocManaged(&d_flag,sizeof(bool)*V),"Flag Allocation");
    catchCudaError(cudaMallocManaged(&d_hybrid,sizeof(bool)),"Hybrid Allocation");
    cudaMallocManaged(&devStates,V*sizeof(curandState));
    
    int *d_left,left;
    cudaMallocManaged(&d_left,sizeof(int));
    
    long int n_threads =  256;
    long int n_blocks = min((V+n_threads-1)/n_threads,(long)MAXBLOCKS);
    //Default case;
    hybrid = true;
    
     /*
    @ step 2 Initialize the colors to -1
    @ until all are colored
    */
    thrust::fill(colors,colors+V,-1);
    long int cr = 1;

    generateRandom<<<n_blocks,n_threads>>>(devStates,d_random,V,time(NULL));
    heapsortBasedOnDegree<<<n_blocks,n_threads>>>(V,preSum,colIndex,degree);

    /*
      Call until all vetex are colored 
    */ 
    clock_t gpu_time = clock();
    do{      
        
       //For Degree Coloring
       if(hybrid){
          left = 0;
          hybrid = false;
          catchCudaError(cudaMemcpy(d_hybrid,&hybrid,sizeof(bool),cudaMemcpyHostToDevice),"ishybrid Copy");
          cudaMemcpy(d_left,&left,sizeof(int),cudaMemcpyHostToDevice);
          degree_Coloring<<<n_blocks,n_threads>>>(V,cr,preSum,colIndex,degree,colors,d_maxDegree,d_flag);       
          assign_DegreeColors<<<n_blocks,n_threads>>>(V,cr,degree,d_flag,colors,d_maxDegree,d_hybrid,d_left);
          catchCudaError(cudaMemcpy(&hybrid,d_hybrid,sizeof(bool),cudaMemcpyDeviceToHost),"ishybrid Copy");
          cudaMemcpy(&left,d_left,sizeof(int),cudaMemcpyDeviceToHost);
           
          if(hybrid==true){
            cr++;
          }else{
            heapsortBasedOnRandom<<<n_blocks,n_threads>>>(V,preSum,colIndex,d_random);
          }   
       }else{
      
         //For Randomized Based Labelling
         left = 0;
         cudaMemcpy(d_left,&left,sizeof(int),cudaMemcpyHostToDevice);
         randomized_Based_Label<<<n_blocks,n_threads>>>(V,cr,preSum,colIndex,d_random,colors,d_maxArray);       
         assignColors<<<n_blocks,n_threads>>>(V,cr,d_random,colors,d_maxArray,d_left);
         cudaMemcpy(&left,d_left,sizeof(int),cudaMemcpyDeviceToHost);
         cr++;
       }             

    
    }while(left);
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

    //thrust::device_ptr<long int> d_ptr = thrust::device_pointer_cast(colors);
    //long int minimumColor = *(thrust::max_element(d_ptr, d_ptr+V));
    thrust::device_vector<long int> d_data(V);
    thrust::copy(colors,colors+V,d_data.begin());
    thrust::sort(d_data.begin(), d_data.end());

    size_t num_unique = thrust::inner_product(d_data.begin(), d_data.end()-1,d_data.begin()+1,0,
                                              thrust::plus<long int>(),thrust::not_equal_to<long int>())+1;

    cudaFree(d_random);
    cudaFree(d_flag);
    cudaFree(d_maxArray);
    cudaFree(d_maxDegree);
    cudaFree(d_hybrid);
    cudaFree(d_left);

    return (long int)num_unique;
}

void GraphColoring_GPUAllocation(const char filename[]){
   
   //@difficult to allocate memory for large complete dataset not assume complete graph
   long int V; //No. of verties
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
   //Allocatio
   long int *colors;
   catchCudaError(cudaMallocManaged(&colors,sizeof(long int)*V),"Color Allocation");

   preSumLength(V,d_preSum,degree);

   for(int i=0;i<V;i++){
      //Remove the cudaMemcpy it will take more time
      thrust::copy(st_Column[i],st_Column[i]+degree[i],d_colIndex+d_preSum[i]);
   }  
   
   /*
   @begin CSR
   */
   long int width=16,height=16;
   long int threads_per_blocks = width*height;
   
   //Call the Hybrid Algorithm
   long int number_Of_Colors_Needed = Hybrid_Algorithm(V,d_preSum,d_colIndex,colors,degree);  
   
   printf("Hybrid Algorithm coloring found solution with %ld colors\n",number_Of_Colors_Needed);
   printf("Valid coloring "); 
   
   bool *d_isValidColors; 
   catchCudaError(cudaMallocManaged(&d_isValidColors,sizeof(bool)),"IsValid Allocation");
   *d_isValidColors = true;

   IsValidgraph_Coloring<<<ceil(V/(width*height))+1,threads_per_blocks>>>(V,colors,d_preSum,d_colIndex,d_isValidColors); 
   catchCudaError(cudaMemcpy(d_isValidColors,d_isValidColors,sizeof(bool),cudaMemcpyDeviceToHost),"Copy isValid Host");
   
   if(*d_isValidColors){
     printf("yes\n");
   }else{
     printf("No\n");
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
       printf("Invalid Input Parameter\n");;
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
