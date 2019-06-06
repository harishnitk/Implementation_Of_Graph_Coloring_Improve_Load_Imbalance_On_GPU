#include<bits/stdc++.h>
#include "Utility.h"
#define MOD 32

using namespace std;

void AssignColors(long int V,long int deltaDegree,long int *colIndex,long int *colors,
                  long int *CS,long int *vforbidden){
   
   

   for(long int i=0;i<V;i++){

      //conflicts colors
      if(colors[i]==0){
        
        //First Available Color
        if(vforbidden[i]==0){

           colors[i] = CS[i]+1;
        
        }else{
           CS[i] = CS[i]+1;
           vforbidden[i] = 0;
        }

     }

   }

}

/*
@smaller index be conflicts
*/
void DetectConflicts(long int V,long int *preSum,long int *colIndex,long int *colors,
                     bool &checkConflict,long int *vforbidden,long int *degree){
   
   for(long int node=0;node<V;node++){
  
        if(colors[node]==0){
           checkConflict = true;
        }

        for(long int k=preSum[node];k<preSum[node+1];k++){

            long int j = colIndex[k];

            if((colors[node]!=0)&&(colors[j]!=0)&&(colors[node]==colors[j])&&(j>node)){

                colors[node] = 0;
                checkConflict = true;
                break;

            }
        }
    
    }

}

/*
@ store the conflicts colors
*/
void ForbiddenColors(long int V,long int *preSum,long int *colIndex,long int *colors,long int *vforbidden,
                     long int *CS){
 

   for(long int node=0;node<V;node++){
      
      for(long int k=preSum[node];k<preSum[node+1];k++){
          
          long int j = colIndex[k];

          if(CS[j]==CS[node]){
             
             if(colors[j]!=0&&colors[node]==0){
                vforbidden[node]|= colors[j];
             }else if(colors[j]==0&&colors[node]!=0){
                vforbidden[j]|= colors[node];
             }
          }        

      }
      

   }

  
}


long int EdgeBased_Algorithm(long int V,long int *NNZ,long int *preSum,long int *colIndex,
                             long int *colors,long int *degree,long int n_zero_counter,long int deltaDegree){
     
    long int *CS = (long int *)malloc(sizeof(long int)*V);
    long int *vforbidden = (long int*)malloc(sizeof(long int)*V);
    memset(vforbidden,0,sizeof(long int)*V);
    memset(CS,0,sizeof(long int)*V); 
    long int minimumColor = 0;
    
    /*
    @ step 2 Initialize the colors to 0
    @ until all are colored
    */
    memset(colors,0,sizeof(long int)*V);
    
    bool checkConflict;
    clock_t cpu_time = clock();
    do{

      
       checkConflict = false;
       
       AssignColors(V,deltaDegree,colIndex,colors,CS,vforbidden);
       DetectConflicts(V,preSum,colIndex,colors,checkConflict,vforbidden,degree);
       ForbiddenColors(V,preSum,colIndex,colors,vforbidden,CS);
          
    }while(checkConflict);
    
    cpu_time = clock()-cpu_time;
    
    //Assigned Colors
    /*
    @ last step to print the assigned colors
    */
    cout<<endl;
    for(long int i=0;i<V;i++){
       printf("vertex --> %i Assigned Color --> %d\n",i,colors[i]);
    }
    cout<<endl;
    
    printf("Cpu_time is %.6lf\n",(double)cpu_time/((double)CLOCKS_PER_SEC/1000));

    unordered_set<long int> us(colors,colors+V);

    free(CS);
    free(vforbidden);
    
    minimumColor = us.size();
    //required colors needed
    return minimumColor;

}

void CSRConvert(long int rows,long int **NNZ,long int *preSum,long int **colIndex,long int *counter,
                long int **st_Column,long int *degree,long int &deltaDegree){
   
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

     if(deltaDegree<degree[i]){
         deltaDegree = degree[i];
     }
   }

   deltaDegree = deltaDegree+1;
}

void EdgeBased_GraphColoring(long int V,long int **st_Column,long int *degree,long int n_zero_counter){

   //convert the sparse array into compact sparse row format
    long int *NNZ = (long int*)malloc(sizeof(long int));
    long int *preSum = (long int*)malloc(sizeof(long int)*(V + 1));
    long int *colIndex= (long int*)malloc(sizeof(long int));
    long int counter = 0;
    long int *colors = (long int*)malloc(sizeof(long int)*V);
    long int deltaDegree = 0;
    CSRConvert(V, &NNZ, preSum, &colIndex, &counter,st_Column,degree,deltaDegree);
    
    //calling the baseline algorithm
    long int number_Of_Colors_Needed = EdgeBased_Algorithm(V,NNZ,preSum,colIndex,colors,degree,n_zero_counter,deltaDegree);

    cout<<"Algorithm coloring found solution with "<<number_Of_Colors_Needed<<" colors"<<endl;
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
   

   EdgeBased_GraphColoring(V,st_Column,st_degree,n_zero_counter);

   time = clock() - time;
   cout<<"Total execution time is "<<(double)time/(double)CLOCKS_PER_SEC<<endl;


   return 0;
}
