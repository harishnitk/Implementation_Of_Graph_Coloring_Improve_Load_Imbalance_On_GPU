#include<bits/stdc++.h>
#include "Utility.h"
#define MOD 32

using namespace std;

int __popc(unsigned int x){
    
    int count = 0;
    unsigned int temp = x;

    while(temp>0){
      if(temp%2==1){
        count++;
      }
      temp = temp/2;
    }
    
    return count;
}


long int find_FirstUnsetBit(long int Id,long int *preSum,long int *colIndex,long int *colors,unsigned int
                            *Colorset,long int maxClr){

    long int k=preSum[Id],j,phase,end=(maxClr/32)+1,bit=0,left=maxClr%32,count=0;
    
    Colorset[end-1] = left!=0?(1<<left)-1:4294967295;//32 bit number
    
    for(phase=0;phase<end-1;phase++){
      Colorset[phase] = 4294967295;
    }

    for(;k<preSum[Id+1];k++){
       j = colIndex[k];
       
       if(colors[j]!=0){
          phase = colors[j]/32;
          if(colors[j]%32==0){
            phase--;
          }
          Colorset[phase]&= ~(unsigned int)(1<<(colors[j]-1-phase*32));//clear bit
       }
    }
    
    for(phase=0;phase<end;phase++){
      
      
      count = __popc(Colorset[phase]);
      if(count!=0){
        
        count = Id%count;
        j = 0;
        while(count>=0&&Colorset[phase]>0){
            if(Colorset[phase]%2==1){
               bit = j+1;
               count--;
            }
            j++;
            Colorset[phase] = Colorset[phase]/2;
        } 

        if(bit!=0){
          return phase*32+bit;
        }
      
      }

    }
    
    return -1;

}

void AssignColors(long int V,long int *preSum,long int *colIndex,long int *colors,unsigned int *Colorset,
                  long int maxClr,int &inc,long int *degree){

     for(long int u=0;u<V;u++){
         
         //need to recolor
         if(colors[u]==0){
            
            if(degree[u]==0){
               
               colors[u] = 1;

            }else{
               
               colors[u] = find_FirstUnsetBit(u,preSum,colIndex,colors,Colorset,maxClr);
               if(colors[u]==-1){
                  colors[u] = 0;
                  inc = 1;
               }

            }
         } 
     }

}

void DetectConflicts(long int V,long int *preSum,long int *colIndex,long int *colors,
                     bool &checkConflict,long int *degree){
   
   for(long int node=0;node<V;node++){
        
        if(colors[node]==0){
           checkConflict = true;
           continue;
        }

        for(long int k=preSum[node];k<preSum[node+1];k++){

            long int j = colIndex[k];

            if(colors[j]==0){
              continue;
            }
            
            if((colors[node]==colors[j])&&(degree[j]<degree[node])){
                colors[j] = 0;
                checkConflict = true;
                break;
            }else if((colors[node]==colors[j])&&(degree[j]>degree[node])){
                colors[node] = 0;
                checkConflict = true;
                break;
            }if((colors[node]==colors[j])&&(j>node)){
                colors[node] = 0;
                checkConflict = true;
                break;
            }else if((colors[node]==colors[j])&&(j<node)){
                colors[j] = 0;
                checkConflict = true;
                break;
            }
        }
    
    }

}

long int CConflict_Algorithm(long int V,long int *NNZ,long int *preSum,long int *colIndex,
                             long int *colors,long int *degree,long int deltaDegree){
     
    /*
    @ step 2 Initialize the colors to 0
    @ until all are colored
    */
    memset(colors,0,sizeof(long int)*V);

    long int minimumColor = 0,maxClr=2;
    long int n_threads =  256;
    bool checkConflict;
    int inc;
    unsigned int *d_Colorset = (unsigned int*)malloc(sizeof(unsigned int)*deltaDegree);
    /*
    Default Case
    */
    clock_t cpu_time = clock();
    do{
       
       inc = 0;
       checkConflict = false;
       AssignColors(V,preSum,colIndex,colors,d_Colorset,maxClr,inc,degree);
       DetectConflicts(V,preSum,colIndex,colors,checkConflict,degree);

       if(inc==1){
          maxClr = 2*maxClr;
       }

    }while(checkConflict);
    cpu_time = clock()-cpu_time;

    //Assigned Colors
    /*
    @ last step to print the assigned colors
    */
    cout<<endl;
    //for(long int i=0;i<V;i++){
      // printf("vertex --> %i Assigned Color --> %d\n",i,colors[i]);
    //}
    cout<<endl;
    
    printf("cpu_time is %.6lf\n",(double)cpu_time/((double)CLOCKS_PER_SEC/1000));

    //long int *req_clr = max_element(colors, colors+V);
    unordered_set<long int> us(colors,colors+V);

    free(d_Colorset);

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

   deltaDegree+= 1;
}

void Conflict_GraphColoring(long int V,long int **st_Column,long int *degree,long int n_zero_counter){

   //convert the sparse array into compact sparse row format
    long int *NNZ = (long int*)malloc(sizeof(long int));
    long int *preSum = (long int*)malloc(sizeof(long int)*(V + 1));
    long int *colIndex= (long int*)malloc(sizeof(long int));
    long int counter = 0;
    long int *colors = (long int*)malloc(sizeof(long int)*V);
    long int deltaDegree = 0;
    
    CSRConvert(V, &NNZ, preSum, &colIndex, &counter,st_Column,degree,deltaDegree);
    
    //calling the baseline algorithm
    long int number_Of_Colors_Needed = CConflict_Algorithm(V,NNZ,preSum,colIndex,colors,degree,deltaDegree);

    cout<<"CConflict Algorithm coloring found solution with "<<number_Of_Colors_Needed<<" colors"<<endl;
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
   

   Conflict_GraphColoring(V,st_Column,st_degree,n_zero_counter);

   time = clock() - time;
   cout<<"Total execution time is "<<(double)time/(double)CLOCKS_PER_SEC<<endl;


   return 0;
}
