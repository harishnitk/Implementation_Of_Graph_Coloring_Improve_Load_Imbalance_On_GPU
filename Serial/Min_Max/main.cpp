#include<bits/stdc++.h>
#include "Utility.h"
using namespace std;

void degreeBased_Approach(long int V,long int c,long int *NNZ,long int *preSum,long int *colIndex,
                          long int *degree,long int *colors,long int &both,bool &hybrid){

    /*
    @check all the vertex
    */
    long int max_Degree = INT_MIN,min_Degree=LONG_MAX;
    bool *flag = (bool*)malloc(sizeof(bool)*V);
    long int *maxDegree_Array = (long int*)malloc(sizeof(long int)*V);
    long int *minDegree_Array = (long int*)malloc(sizeof(long int)*V);
    memset(maxDegree_Array,-1,sizeof(long int)*V);
    memset(flag,false,sizeof(bool)*V);
    hybrid = false;

    for(long int i=0;i<V;i++){

        if(colors[i]==-1){
           max_Degree = -1;
           min_Degree = LONG_MAX;

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
                }
                if(degree[j]<min_Degree){
                    min_Degree = degree[j];
                }
              }
           }

          maxDegree_Array[i] = degree[i]>max_Degree?1:0;
          minDegree_Array[i] = degree[i]<min_Degree?1:0;

        }
    }

    for(long int i=0;i<V;i++){
        if(colors[i]==-1){
            if(flag[i]==false){
                if(maxDegree_Array[i]){
                    colors[i] = c;
                    hybrid = true;
                }else if(minDegree_Array[i]){
                    colors[i] = c+1;
                    hybrid = true;
                }
            }
        }
    }

}

void minMaxBased_Approach(long int V,long int c,long int *NNZ,long int *preSum,long int *colIndex,long int *random,long int *colors,long int &both){

    /*
    @check all the vertex
    */
    long int max_Value = INT_MIN,min_Value=LONG_MAX;
    long int *max_Array = (long int*)malloc(sizeof(long int)*V);
    long int *min_Array = (long int*)malloc(sizeof(long int)*V);
    memset(max_Array,-1,sizeof(long int)*V);
    
    for(long int i=0;i<V;i++){

        if(colors[i]!=-1){
            continue;
        }else{
           max_Value = -1;
           min_Value = LONG_MAX;
        }

        for(long int k=preSum[i];k<preSum[i+1];k++){
            long int j = colIndex[k];
            long int jc = colors[j];

            if(((jc!=-1)&&(jc!=c))||(i==j)){
               continue;
            }else{
                if(random[j]>max_Value){
                    max_Value = random[j];
                }
                if(random[j]<min_Value){
                    min_Value = random[j];
                }
            }
        }
        max_Array[i] = random[i]>max_Value?1:0;
        min_Array[i] = random[i]<min_Value?1:0;
    }

    for(long int i=0;i<V;i++){
        if(colors[i]==-1){
            if(max_Array[i]){
                colors[i] = c;
                both++;
            }else if(min_Array[i]){
                colors[i] = c+1;
                both++;
            }
        }
    }

}

long int MinMax_Algorithm(long int V,long int *NNZ,long int *preSum,long int *colIndex,long int *colors,
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

   long int cr = 1,both=0;
   clock_t cpu_time = clock();

   do{
      
       long int left = 0,both=0;
       if(hybrid||isDegreeChange){

          degreeBased_Approach(V,cr,NNZ,preSum,colIndex,degree,colors,both,hybrid);

          if(hybrid==false){
             isDegreeChange = 0;
          }

       }else{

          both = 0;
          minMaxBased_Approach(V,cr,NNZ,preSum,colIndex,random,colors,both);

       }
      
      for(long int k=0;k<V;k++){
         if(colors[k]==-1){
            left++;
         }
      }

      if(left==0){
         break;
      }else{
          cr = cr+2;
      }

   }while(left);
   cpu_time = clock()-cpu_time;

   cout<<endl;
   for(long int i=0;i<V;i++){
      cout<<"Vertex --> "<<i<<" Assigned Color --> "<<colors[i]<<endl;
   }
   cout<<endl;
   
   unordered_set<long int> us(colors,colors+V);
   minimumColor = us.size();
   printf("Cpu_Time is %.6lf\n",(double)cpu_time/((double)CLOCKS_PER_SEC/1000));

   free(random);
   //required colors
   return minimumColor;
   
}

void CSRConvert(long int rows,long int **NNZ,long int *preSum,long int **colIndex,long int *counter,
                long int **st_Column,long int *degree,long int &isDegreeChange){
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

void MinMax_GraphColoring(long int V,long int **st_Column,long int *degree,long int n_zero_counter){

   //convert the sparse array into compact sparse row format
    long int *NNZ = (long int*)malloc(sizeof(long int));
    long int *preSum = (long int*)malloc(sizeof(long int)*(V + 1));
    long int *colIndex= (long int*)malloc(sizeof(long int));
    long int counter = 0,isDegreeChange=0;
    long int *colors = (long int*)malloc(sizeof(long int)*V);
    CSRConvert(V, &NNZ, preSum, &colIndex, &counter,st_Column,degree,isDegreeChange);
    
    //calling the baseline algorithm
    long int number_Of_Colors_Needed = MinMax_Algorithm(V,NNZ,preSum,colIndex,colors,degree,isDegreeChange);

    cout<<"MinMax Algorithm coloring found solution with "<<number_Of_Colors_Needed<<" colors"<<endl;
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
   

   MinMax_GraphColoring(V,st_Column,st_degree,n_zero_counter);

   time = clock() - time;
   cout<<"Total execution time is "<<(double)time/(double)CLOCKS_PER_SEC<<endl;


   return 0;
}
