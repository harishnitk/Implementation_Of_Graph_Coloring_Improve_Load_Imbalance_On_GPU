#include<bits/stdc++.h>
using namespace std;

//Read MatrixMarket graphs
//Assumes input nodes are numbered starting from 1
void ReadMMFile(const char filename[], long int *V,long int ***st_Column,long int **st_degree,long int *counter){
   string s;
   ifstream infile(filename);
   if(infile.fail()){
      cout<<"Fail to open the file\n";
      return;
   }

   //content
   while(getline(infile,s)){
     istringstream iss(s);
     if(s.find("%")==string::npos){
        break;
     }
   }

   istringstream iss(s);
   //Maintain Hash for NNZ,preSum,colIndex
   //allocate dynamic size

   long int n_rows,n_cols,n_edges;
   iss>>n_rows>>n_cols>>n_edges;
   *st_degree = new long int[n_rows];
   *st_Column = new long int*[n_rows];
   memset(*st_degree,0,n_rows*sizeof(long int));
   *V = n_rows;

   //reading edges

   while(getline(infile,s)){
      istringstream iss(s);
      long int u,v,w;
      iss>>u>>v>>w;
      if(u!=v){
        long int u_len = (*st_degree)[u-1];
        long int v_len = (*st_degree)[v-1];
        (*st_Column)[u-1] = (long int*)realloc((*st_Column)[u-1],sizeof(long int)*(u_len+1));
        (*st_Column)[v-1] = (long int*)realloc((*st_Column)[v-1],sizeof(long int)*(v_len+1));
        (*st_Column)[u-1][u_len] = v-1; 
        (*st_Column)[v-1][v_len] = u-1;
        (*st_degree)[u-1]++;
        (*st_degree)[v-1]++;
        *counter+=2;
      }
   }

   infile.close();
}


//Read DIMACS graphs
//Assumes input nodes are numbered starting from 1 in the column file
void ReadColFile(const char filename[],long int *V,long int ***st_Column,long int **st_degree,long int *counter){
   string s;
   ifstream infile(filename);
   if(infile.fail()){
      cout<<"Fail to open the file\n";
      exit(0);
   }

   long int n_rows,n_edges;
   //Maintain Hash for NNZ,preSum,colIndex
   //allocate dynamic size
   while(getline(infile,s)){
       istringstream iss(s);
       string str;
       long int u,v;
       iss>>str;
       if(str=="p"){
          iss>>s;
          iss>>n_rows;
          iss>>n_edges;
          *V = n_rows;
          *st_degree = new long int[n_rows];
          *st_Column = new long int*[n_rows];
          memset(*st_degree,0,n_rows*sizeof(long int));
          continue;
       }else if(str!="e"){
          continue; 
       }

       iss>>u>>v;
       if(u!=v){
         long int u_len = (*st_degree)[u-1];
         long int v_len = (*st_degree)[v-1];
         (*st_Column)[u-1] = (long int*)realloc((*st_Column)[u-1],sizeof(long int)*(u_len+1));
         (*st_Column)[v-1] = (long int*)realloc((*st_Column)[v-1],sizeof(long int)*(v_len+1));
         (*st_Column)[u-1][u_len] = v-1; 
         (*st_Column)[v-1][v_len] = u-1;
         (*st_degree)[u-1]++;
         (*st_degree)[v-1]++;
         *counter+=2;
       }
   }

   infile.close();
}

/*
@[description] Counts the number of unique colors in a solution
*/
int CountColors(long int V,long int *colors)
{
   long int n_colors=0;
   set<long int>seen_colors;

   for(long int i=0;i<V;i++) {
      if(seen_colors.find(colors[i])==seen_colors.end()) {
         seen_colors.insert(colors[i]);
         n_colors++;
      }
   }

   return n_colors;
}


/*
@[description]Returns true if the color assignment is valid for the graph
*/
bool IsValidColoring(long int V,long int *colors,long int *colIndex,long int *preSum)
{
   for(long int i=0;i<V;i++) {
      for(long int j=preSum[i];j<preSum[i+1];j++) {
         
          if(colors[i]==colors[colIndex[j]]) {
              printf("Vertex %d and Vertex %d are connected and have the same color %d\n", i, colIndex[j], colors[i]);
              return false;
          }
          if(colors[i]<1) {
              printf("Vertex %d has invalid color %d\n", i, colors[i]);
             return false;
          }
         
      }
   }

   return true;
}

/*
@[description] A utility function to check if the current color assignment
               is safe for vertex v
*/
bool IsSafe(long int v,bool *graph,long int V,long int *colors, long int c)
{
   for(long int i=0;i<V;i++)
      if(graph[v*V+i]&&c==colors[i])
         return false;

   return true;
}
