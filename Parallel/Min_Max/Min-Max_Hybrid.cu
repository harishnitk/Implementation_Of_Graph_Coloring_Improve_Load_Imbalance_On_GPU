/*
Need optimization add sorting
opt of 4-opt
bottlneck of sorting kernel, Modified 30-3-2019 heapSort added
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


__global__ void assignColors(long int V,long int c,double *random,long int *colors,double *randomValue,
                             int *maxe,int *mine,int *left){
   
    long int threadId = blockDim.x*blockIdx.x+threadIdx.x;

    if(threadId<V&&colors[threadId]!=-1){
       
       return;
    
    }else if(threadId<V){
        
        if(colors[threadId]==-1){
            
            if(randomValue[threadId]==1){
                
                colors[threadId] = c;
                atomicAdd(maxe,1);
            
            }else if(randomValue[threadId]==-1){
                
                colors[threadId] = c+1;
                atomicAdd(mine,1);
            
            }else{
             
                *left = 1;
            
            }

        }

    }

}

__global__ void assign_DegreeColors(long int V,long int c,long int *degree,bool *flag,
                                    long int *colors,long int *degreeValue,bool *isDegree,int *maxe,int *mine,int *left){
    
    long int threadId = blockDim.x*blockIdx.x+threadIdx.x;
    
    if(threadId<V&&colors[threadId]!=-1){
       
       return;

    }else if(threadId<V){
        
        if(colors[threadId]==-1){
           
           if(flag[threadId]==false){
              
              if(degreeValue[threadId]==1){
                  
                  colors[threadId] = c;
                  *isDegree = true;
                  atomicAdd(maxe,1);
                  //reAssign Values
              
              }else if(degreeValue[threadId]==-1){
                  
                  colors[threadId] = c+1;
                  *isDegree = true;
                  atomicAdd(mine,1);
              
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
                                long int *degree,long int *colors,long int *degreeValue,bool *d_flag){
  
   long int threadId = blockDim.x*blockIdx.x+threadIdx.x;
   long int minmaxDegree;

   if(threadId<V&&colors[threadId]!=-1){
      
      return;

   }
   else if(threadId<V){
        
       if(colors[threadId]==-1){
          
           minmaxDegree = -1;
           long int start = preSum[threadId],end = preSum[threadId+1]-1;

           for(;start<=end;){
               long int j=colIndex[start];
               long int jc = colors[j];
               
               if(colors[j]!=-1){
                 start++;
                 continue;
               }

              //check all neighborhood which all are uncolored
              if(jc==-1){
                 if(degree[threadId]==degree[j]){
                    d_flag[threadId] = true;
                    break;  
                 }else if(degree[j]>minmaxDegree){
                   minmaxDegree = degree[j];
                   break;
                 }

                 start++;

              }

           }

           degreeValue[threadId] = degree[threadId]>minmaxDegree?1:0;
           
           if(degreeValue[threadId]==0&&d_flag[threadId]==false){
              
              minmaxDegree = LONG_MAX;
              for(;start<=end;){
                 long int j=colIndex[end];
                 long int jc = colors[j];
                 
                 if(colors[j]!=-1){
                   end--;
                   continue;
                 }

                //check all neighborhood which all are uncolored
                if(jc==-1){
                   if(degree[threadId]==degree[j]){
                      d_flag[threadId] = true;
                      break;  
                   }else if(degree[j]<minmaxDegree){
                      minmaxDegree = degree[j];
                      break;
                   }

                   end--;

                }

             }               
             
             degreeValue[threadId] = degree[threadId]<minmaxDegree?-1:0;
          }

       }  
        
   }

}

// How to arrange the program to assign the different color to adjacent
__global__ void minmax_Based_Label(long int V,long int c,long int *preSum,long int *colIndex,
                                   double *random,long int *colors,double *randomValue){
     
     long int threadId = blockDim.x*blockIdx.x+threadIdx.x;
     double minMaxValue;
     
     if(threadId<V&&colors[threadId]!=-1){
        return;
     }
     else if(threadId<V){
              
         if(colors[threadId]==-1){
          
            minMaxValue = -1;
            long int start = preSum[threadId],end = preSum[threadId+1]-1,jcs,jce,js,je;
            
            for(;start<=end;){

                js = colIndex[start];
                jcs = colors[js];
                
                if(jcs!=-1){
                  start++;
                  continue;
                }
                else if(jcs==-1){       
                  //check with random value
                   if(random[js]>minMaxValue){
                      minMaxValue = random[js];
                      break;
                   }
                   start++;
                }
             
            }
            
            randomValue[threadId] = random[threadId]>minMaxValue?1:0;
            
            if(randomValue[threadId]==0){
              
              minMaxValue = LONG_MAX;

              for(;start<=end;){

                  je = colIndex[end];
                  jce = colors[je];
                  
                  if(jce!=-1){
                    end--;
                    continue;
                  }else if(jce==-1){

                    //check with random value
                    if(random[je]<minMaxValue){
                      minMaxValue = random[je];
                      break;
                    }

                    end--;

                   }
               
              }
              
              randomValue[threadId] = random[threadId]<minMaxValue?-1:0;
            
            }

     
          }

     }
} 

int preSumLength(int V,long int *d_preSum,long int *degree){
    
    long int prev = degree[0];
    int flag = 0;

    for(long int i=0;i<V;i++){
       d_preSum[i+1] = d_preSum[i]+degree[i];
       if(prev!=degree[i]){
         flag = 1;
       }
    }
    return flag;
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


long int Min_Max_Algorithm(long int V,long int *preSum,long int *colIndex,long int *colors,
                           long int *degree,int isSwitching){
    curandState* devStates;
    double *d_random;
    cudaMallocManaged(&d_random,sizeof(double)*V);
    long int *d_degreeValue;
    double *d_randomValue;
    bool *d_flag;
    bool *d_isDegree,isDegree;
    int *d_mine,*d_maxe,*d_left,left,maxe,mine;

    //Allocate the memory for maxArray
    catchCudaError(cudaMallocManaged(&d_randomValue,sizeof(double)*V),"d_randomValue Allocation");
    catchCudaError(cudaMallocManaged(&d_degreeValue,sizeof(long int)*V),"d_degreeValue Allocation");
    catchCudaError(cudaMallocManaged(&d_flag,sizeof(bool)*V),"Flag Allocation");
    catchCudaError(cudaMallocManaged(&d_isDegree,sizeof(bool)),"Check Allocation");
    catchCudaError(cudaMallocManaged(&d_maxe,sizeof(int)),"maxe Allocation");
    catchCudaError(cudaMallocManaged(&d_mine,sizeof(int)),"mine Allocation");
    catchCudaError(cudaMallocManaged(&d_left,sizeof(int)),"left Allocation");
    cudaMallocManaged(&devStates,V*sizeof(curandState));

    long int n_threads =  256;
    long int n_blocks = min((V+n_threads-1)/n_threads,(long)MAXBLOCKS);
    /*
    @ step 2 Initialize the colors to -1
    @ until all are colored
    */
    thrust::fill(colors,colors+V,-1);
    long int cr = 1;

    //Default case;
    isDegree = false;
    clock_t gpu_time = clock();

    //Step 1 assign the random value to all the vertices
    if(isSwitching){
       generateRandom<<<n_blocks,n_threads>>>(devStates,d_random,V,time(NULL));
       heapsortBasedOnDegree<<<n_blocks,n_threads>>>(V,preSum,colIndex,degree);
    }else{
       generateRandom<<<n_blocks,n_threads>>>(devStates,d_random,V,time(NULL));
       heapsortBasedOnRandom<<<n_blocks,n_threads>>>(V,preSum,colIndex,d_random);
    }
    
    do{
      
       //For Degree Coloring
       if(isDegree||isSwitching){
          
          isDegree = false;
          mine = 0;
          maxe = 0;
          left = 0;
          catchCudaError(cudaMemcpy(d_isDegree,&isDegree,sizeof(bool),cudaMemcpyHostToDevice),"isDegreeD Copy");
          catchCudaError(cudaMemcpy(d_maxe,&maxe,sizeof(int),cudaMemcpyHostToDevice),"maxeD copy");
          catchCudaError(cudaMemcpy(d_maxe,&mine,sizeof(int),cudaMemcpyHostToDevice),"mineD copy");
          catchCudaError(cudaMemcpy(d_left,&left,sizeof(int),cudaMemcpyHostToDevice),"leftD copy");
          degree_Coloring<<<n_blocks,n_threads>>>(V,cr,preSum,colIndex,degree,colors,d_degreeValue,d_flag);       
          assign_DegreeColors<<<n_blocks,n_threads>>>(V,cr,degree,d_flag,colors,d_degreeValue,d_isDegree,d_maxe,d_mine,d_left);
          catchCudaError(cudaMemcpy(&isDegree,d_isDegree,sizeof(bool),cudaMemcpyDeviceToHost),"isDegree Copy");
          catchCudaError(cudaMemcpy(&maxe,d_maxe,sizeof(int),cudaMemcpyDeviceToHost),"maxe copy");
          catchCudaError(cudaMemcpy(&mine,d_mine,sizeof(int),cudaMemcpyDeviceToHost),"mine copy");;
          catchCudaError(cudaMemcpy(&left,d_left,sizeof(int),cudaMemcpyDeviceToHost),"left copy");
          
          cr = (maxe&&mine)?cr+2:(maxe||mine)?cr+1:cr;

          if(isDegree==false){
            isSwitching = 0;
            heapsortBasedOnRandom<<<n_blocks,n_threads>>>(V,preSum,colIndex,d_random);
          }

       }else{
         
         //For Randomized Based Labelling
         maxe = 0;
         mine = 0;
         left = 0;
         catchCudaError(cudaMemcpy(d_left,&left,sizeof(int),cudaMemcpyHostToDevice),"left copy");
         catchCudaError(cudaMemcpy(d_maxe,&maxe,sizeof(int),cudaMemcpyHostToDevice),"maxe copy");
         catchCudaError(cudaMemcpy(d_maxe,&mine,sizeof(int),cudaMemcpyHostToDevice),"mine copy");
         minmax_Based_Label<<<n_blocks,n_threads>>>(V,cr,preSum,colIndex,d_random,colors,d_randomValue);       
         assignColors<<<n_blocks,n_threads>>>(V,cr,d_random,colors,d_randomValue,d_maxe,d_mine,d_left);
         catchCudaError(cudaMemcpy(&left,d_left,sizeof(int),cudaMemcpyDeviceToHost),"left copy");
         catchCudaError(cudaMemcpy(&maxe,d_maxe,sizeof(int),cudaMemcpyDeviceToHost),"maxe copy");
         catchCudaError(cudaMemcpy(&mine,d_mine,sizeof(int),cudaMemcpyDeviceToHost),"mine copy");
         
         cr = (maxe&&mine)?cr+2:(maxe||mine)?cr+1:cr;

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
    cudaFree(d_randomValue);
    cudaFree(d_degreeValue);
    cudaFree(d_isDegree);
    cudaFree(d_left);
    cudaFree(d_maxe);
    cudaFree(d_mine);
    
    //required colors needed
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
   
   //Allocation
   long int *colors;
   catchCudaError(cudaMallocManaged(&colors,sizeof(long int)*V),"Color Allocation");

   int isSwitching = preSumLength(V,d_preSum,degree);

   for(int i=0;i<V;i++){
      //Remove the cudaMemcpy it will take more time
      thrust::copy(st_Column[i],st_Column[i]+degree[i],d_colIndex+d_preSum[i]);
   }  
   
   /*
   @begin CSR
   */
   long int width=16,height=16;
   long int threads_per_blocks = width*height;
   
   //Call the Min-Max Algorithm
   long int number_Of_Colors_Needed = Min_Max_Algorithm(V,d_preSum,d_colIndex,colors,degree,isSwitching);  
   
   printf("Min-Max Algorithm coloring found solution with %d colors\n",number_Of_Colors_Needed);
   printf("Valid coloring "); 
   
   bool *d_isValidColors; 
   catchCudaError(cudaMallocManaged(&d_isValidColors,sizeof(bool)),"IsValid Allocation");
   *d_isValidColors = true;
   
   IsValidgraph_Coloring<<<ceil(V/threads_per_blocks)+1,threads_per_blocks>>>(V,colors,d_preSum,d_colIndex,d_isValidColors); 
   catchCudaError(cudaMemcpy(d_isValidColors,d_isValidColors,sizeof(bool),cudaMemcpyDeviceToHost),"IsValid Host Copy");
   
   if(*d_isValidColors){
     printf("yes\n");
   }else{
     printf("No");
   }
   
   //catchCudaError(cudaDeviceSynchronize(),"GraphColoring DeviceSync");
   
   cudaFree(d_preSum);
   cudaFree(d_colIndex);
   cudaFree(colors);
   cudaFree(degree);
   //cudaFree(d_isValidColors);
  
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
