#include<bits/stdc++.h>
#include "Utility.h"
#define MOD 32

using namespace std;

void AssignColors(long int V,long int *preSum,long int *colIndex,long int *colors,
                  long int deltaDegree,bool *d_conflicts){
   
   
   for(long int i=0;i<V;i++){

        
        if(!d_conflicts[i]){
          continue;
        }
        
        long int *vforbidden = (long int*)malloc(sizeof(long int)*(deltaDegree+1));
        memset(vforbidden,0,sizeof(long int)*(deltaDegree+1));

        for(long int k=preSum[i];k<preSum[i+1];k++){

           long int j = colIndex[k];
           long int value = colors[j]%MOD;
           long int shift = 1<<value;
           vforbidden[colors[j]/MOD]|= shift; 
        }
        
        //Assign colors
        for(long int color=1;color<=deltaDegree+1;color++){

            long int val = color%MOD;
            
            if((vforbidden[color/MOD]&(1<<val))== 0){
                colors[i] = color;
                break;
            }

        }
    }


}

void DetectConflicts(long int V,long int *preSum,long int *colIndex,long int *colors,
                     bool *conflicts,bool &checkConflict){
   
   for(long int node=0;node<V;node++){
  
        conflicts[node] = false;

        for(long int k=preSum[node];k<preSum[node+1];k++){

            long int j = colIndex[k];

            if((colors[node]==colors[j])&&(j<node)){

                conflicts[node] = true;
                checkConflict = true;
                break;

            }
        }
    
    }

}

long int EdgeBased_Algorithm(long int V,long int *NNZ,long int *preSum,long int *colIndex,
                             long int *colors,long int *degree,long int deltaDegree){
     
      /*
    @ step 2 Initialize the colors to 0
    @ until all are colored
    */
    memset(colors,0,sizeof(long int)*V);

    long int minimumColor = 0;
    long int n_threads =  256;
    bool *d_conflicts;
    bool checkConflict;
    d_conflicts = (bool*)malloc(sizeof(bool)*V);
    memset(d_conflicts,true,sizeof(bool)*V);
    
    do{
       
       checkConflict = false;

       AssignColors(V,preSum,colIndex,colors,deltaDegree,d_conflicts);
       DetectConflicts(V,preSum,colIndex,colors,d_conflicts,checkConflict);
        

    }while(checkConflict);
    
    //Assigned Colors
    /*
    @ last step to print the assigned colors
    */
    cout<<endl;
    for(long int i=0;i<V;i++){
       printf("vertex --> %i Assigned Color --> %d\n",i,colors[i]);
    }
    cout<<endl;

    long int *req_clr = max_element(colors, colors+V);

    free(d_conflicts);
    
    minimumColor = *req_clr;
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
    long int number_Of_Colors_Needed = EdgeBased_Algorithm(V,NNZ,preSum,colIndex,colors,degree,deltaDegree);

    cout<<"EdgeBase Algorithm coloring found solution with "<<number_Of_Colors_Needed<<" colors"<<endl;
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
