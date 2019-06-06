#include<bits/stdc++.h>
#include "Utility.h"
using namespace std;

//Heap Sort
void heapsortBasedOnDegree(long int V,long int *preSum,long int *colIndex,long int *degree){
            
    for(long int Id = 0;Id<V;Id++){ 

        long int n = preSum[Id+1]-preSum[Id];
        long int largest,previous,temp;
        
        for(long int i = preSum[Id]+n/2-1;i>=preSum[Id];i--){ 
            
            largest = i,previous=-1;     
            while(largest!=previous){
              long int l = 2*largest-preSum[Id]+1;
              long int r = 2*largest-preSum[Id]+2;
              previous = largest;
              if(l<preSum[Id+1]&&degree[colIndex[l]]<degree[colIndex[largest]]) 
                  largest = l; 
     
              if(r<preSum[Id+1]&&degree[colIndex[r]]<degree[colIndex[largest]]) 
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
        for(long int i=preSum[Id+1]-1;i>=preSum[Id];i--) 
        { 
            temp = colIndex[preSum[Id]];
            colIndex[preSum[Id]] = colIndex[i];
            colIndex[i] = temp;
            largest = preSum[Id],previous=-1;
            while(largest!=previous){
              int l = 2*largest-preSum[Id]+1;
              int r = 2*largest-preSum[Id]+2;
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
void heapsortBasedOnRandom(long int V,long int *preSum,long int *colIndex,long int *random){
            
    for(long int Id = 0;Id<V;Id++){ 

        long int n = preSum[Id+1]-preSum[Id];
        long int largest,previous,temp;
        
        for(long int i = preSum[Id]+n/2-1;i>=preSum[Id];i--){ 
            
            largest = i,previous=-1;     
            while(largest!=previous){
              long int l = 2*largest-preSum[Id]+1;
              long int r = 2*largest-preSum[Id]+2;
              previous = largest;
              if(l<preSum[Id+1]&&random[colIndex[l]]<random[colIndex[largest]]) 
                  largest = l; 
     
              if(r<preSum[Id+1]&&random[colIndex[r]]<random[colIndex[largest]]) 
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
        for(long int i=preSum[Id+1]-1;i>=preSum[Id];i--) 
        { 
            temp = colIndex[preSum[Id]];
            colIndex[preSum[Id]] = colIndex[i];
            colIndex[i] = temp;
            largest = preSum[Id],previous=-1;
            while(largest!=previous){
              int l = 2*largest-preSum[Id]+1;
              int r = 2*largest-preSum[Id]+2;
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

void degreeBased_Approach(long int V,long int c,long int *NNZ,long int *preSum,long int *colIndex,long int *degree,long int *colors,
                          bool &ind_hybrid){

    /*
    @check all the vertex
    */
    long int max_Degree = INT_MIN;
    bool *flag = (bool*)malloc(sizeof(bool)*V);
    long int *maxDegree_Array = (long int*)malloc(sizeof(long int)*V);
    memset(maxDegree_Array,-1,sizeof(long int)*V);
    memset(flag,false,sizeof(bool)*V);
    ind_hybrid = false;

    for(long int i=0;i<V;i++){

        if(colors[i]==-1){
           max_Degree = -1;

           for(long int k=preSum[i];k<preSum[i+1];k++){
              long int j = colIndex[k];
              long int jc = colors[j];

              if(jc==-1){
                 if(degree[j]==degree[i]){
                    flag[i] = true;
                    break;
                }
                if(degree[j]>max_Degree){
                    max_Degree = degree[j];
                    break;
                }
              }
           }

          maxDegree_Array[i] = max_Degree;
        }
    }

    for(long int i=0;i<V;i++){
        if(colors[i]==-1){
            if(flag[i]==false){
                if(degree[i]>maxDegree_Array[i]){
                    colors[i] = c;
                    ind_hybrid = true;
                }
            }
        }
    }

}

void randomizedBased_Approach(long int V,long int c,long int *NNZ,long int *preSum,long int *colIndex,long int *random,long int *colors){

    /*
    @check all the vertex
    */
    long int max_Value = INT_MIN;
    long int *max_Array = (long int*)malloc(sizeof(long int)*V);
    memset(max_Array,-1,sizeof(int)*V);
    
    for(long int i=0;i<V;i++){

        if(colors[i]!=-1){
            continue;
        }else{
           max_Value = -1;
        }

        for(long int k=preSum[i];k<preSum[i+1];k++){
            long int j = colIndex[k];
            long int jc = colors[j];

            if(((jc!=-1)&&(jc!=c))||(i==j)){
               continue;
            }else{
                if(random[j]>max_Value){
                    max_Value = random[j];
                    break;
                }
            }
        }
        max_Array[i] = max_Value;
    }

    for(long int i=0;i<V;i++){
        if(colors[i]==-1){
            if(random[i]>max_Array[i]){
                colors[i] = c;
            }
        }
    }

}

long int Hybrid_Algorithm(long int V,long int *NNZ,long int *preSum,long int *colIndex,long int *colors,
                          long int *degree,long int isDegreeChange){

   long int *random = (long int*)malloc(sizeof(long int)*V);
   long int minimumColor = 0;
   bool hybrid = false;

   srand(time(NULL));
   //Allocation Baseline Algorithm
   for(long int i=0;i<V;i++){
      random[i] = rand();
   }

   memset(colors,-1,sizeof(long int)*V);
   long int cr = 1;
   
   clock_t cpu_time = clock();
   heapsortBasedOnDegree(V,preSum,colIndex,degree);

   for(;cr<=V;cr++){
      
      long int left = 0;

      if(hybrid||isDegreeChange){
        
        degreeBased_Approach(V,cr,NNZ,preSum,colIndex,degree,colors,hybrid);

        if(hybrid==false){
           heapsortBasedOnRandom(V,preSum,colIndex,random);
           isDegreeChange = 0;
        }

      }else{
        
        randomizedBased_Approach(V,cr,NNZ,preSum,colIndex,random,colors);

      }

      for(long int k=0;k<V;k++){
        if(colors[k]==-1){
           left++;
        }
      }

      if(left==0){
         break;
      }

   }
   cpu_time = clock()-cpu_time;

   cout<<endl;
   for(long int i=0;i<V;i++){
      cout<<"Vertex --> "<<i<<" Assigned Color --> "<<colors[i]<<endl;
   }
   cout<<endl;
   
   unordered_set<long int> us(colors,colors+V);
   minimumColor = us.size();
   printf("Cpu_time is %.6lf\n",(double)cpu_time/((double)CLOCKS_PER_SEC/1000));

   free(random);
   //required colors
   return minimumColor;
   
}

void CSRConvert(long int rows,long int **NNZ,long int *preSum,long int **colIndex,long int *counter,long int **st_Column,
                long int &isDegreeChange,long int *degree){
   long int cols = rows;
   *counter = 0;
   preSum[0] = 0;

   for(long int i=0;i<rows;i++){
      for(long int j=0;j<degree[i];j++){
         
          *NNZ = (long int*)realloc(*NNZ,sizeof(long int)*(*counter+1));
          (*NNZ)[*counter] = 1;
          *colIndex = (long int*)realloc(*colIndex,sizeof(long int)*(*counter+1));
          (*colIndex)[*counter] = st_Column[i][j];
          *counter+=1;
         
      }
   }

   for(long int i=0;i<rows;i++){
      preSum[i+1] = preSum[i]+degree[i];
      if(i+1<rows&&(degree[i]!=degree[i+1])){
         isDegreeChange = 1;
      }
   }

}

void Hybrid_GraphColoring(long int V,long int **st_Column,long int *degree,long int n_zero_counter){

   //convert the sparse array into compact sparse row format
    long int *NNZ = (long int*)malloc(sizeof(long int));
    long int *preSum = (long int*)malloc(sizeof(long int)*(V + 1));
    long int *colIndex= (long int*)malloc(sizeof(long int));
    long int counter = 0,isDegreeChange=0;
    long int *colors = (long int*)malloc(sizeof(long int)*V);
    CSRConvert(V, &NNZ, preSum, &colIndex, &counter,st_Column,isDegreeChange,degree);
    
    //calling the baseline algorithm
    long int number_Of_Colors_Needed = Hybrid_Algorithm(V,NNZ,preSum,colIndex,colors,degree,isDegreeChange);

    cout<<"Hybrid Algorithm coloring found solution with "<<number_Of_Colors_Needed<<" colors"<<endl;
    cout<<"Is Valid coloring ";
    if(IsValidColoring(V,colors,colIndex,preSum)){
        cout<<"Yes"<<endl;
    }else{
        cout<<"No"<<endl;
    }

    free(NNZ);
    free(preSum);
    free(colIndex);
    free(colors);
}

int main(int argc,char *argv[])
{
   /*
   @Assume like a complete graph [V-vertex]
   @colors Array
   */
   long int V; //No. of verties
   long int n_zero_counter = 0;   
   long int **st_Column;
   long int *st_degree;

   /*
   @Adding the clock
   */
   clock_t time = clock();
   

   if(argc<2){
      cout<<"Invalid CommandLine Argument-->  ./out filename"<<endl; 
   }else{
      if(string(argv[1]).find("col")!=string::npos){
         ReadColFile(argv[1],&V,&st_Column,&st_degree,&n_zero_counter);
      }else{
         ReadMMFile(argv[1],&V,&st_Column,&st_degree,&n_zero_counter); 
      }
   }
   
   Hybrid_GraphColoring(V,st_Column,st_degree,n_zero_counter);

   time = clock() - time;
   cout<<"Total execution time is "<<(double)time/(double)CLOCKS_PER_SEC<<endl;

   free(*st_Column);
   free(st_degree);

   return 0;
}
