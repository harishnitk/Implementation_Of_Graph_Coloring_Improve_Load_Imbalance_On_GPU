#include<bits/stdc++.h>
using namespace std;

int main(int argc,char *argv[])
{
	
	string str[] = {"./main ../../../../DataSet/amazon.col","./main ../../../../DataSet/DSJC250.9.col",
                  "./main ../../../../DataSet/DSJC500.1.col","./main ../../../../DataSet/DSJC500.9.col",
                  "./main ../../../../DataSet/DSJC1000.9.col","./main ../../../../DataSet/DSJR500.1.col",
                  "./main ../../../../DataSet/DSJR500.5.col","./main ../../../../DataSet/flat300_20_0.col",
                  "./main ../../../../DataSet/flat300_26_0.col","./main ../../../../DataSet/flat1000_50_0.col",
                  "./main ../../../../DataSet/flat1000_76_0.col","./main ../../../../DataSet/fpsol2.i.2.col",
                  "./main ../../../../DataSet/frb100-40.col","./main ../../../../DataSet/frb100-40-1.col",
                  "./main ../../../../DataSet/inithx.i.3.col","./main ../../../../DataSet/le450_15a.col",
                  "./main ../../../../DataSet/le450_15d.col","./main ../../../../DataSet/le450_25d.col",
                  "./main ../../../../DataSet/miles500.col","./main ../../../../DataSet/miles1000.col",
                  "./main ../../../../DataSet/movie_taste_net.col","./main ../../../../DataSet/mulsol.i.5.col",
                  "./main ../../../../DataSet/myciel4.col","./main ../../../../DataSet/queen16_16.col"};
	
	for(int itr=0;itr<4;itr++){
	  
	  const char *command = str[itr].c_str();
	  for(int i=0;i<5;i++){
         system(command);
      }

      cout<<"--------------------------------------------------------------------"<<endl;
    
    }
    
    /*string str2 = "./main ../../../../DataSet/audikw_1.mtx";

    const char *command2 = str2.c_str();
	  for(int i=0;i<5;i++){
         system(command2);
    }

    cout<<"--------------------------------------------------------------------"<<endl;
    

    
    string str3 = "./main ../../../../DataSet/Bump_2911.mtx";

    const char *command3 = str3.c_str();
	  for(int i=0;i<5;i++){
         system(command3);
    }

    cout<<"--------------------------------------------------------------------"<<endl;
    

    string str4 = "./main ../../../../DataSet/cage12.mtx";

    const char *command4 = str4.c_str();
	  for(int i=0;i<5;i++){
         system(command4);
    }

    cout<<"--------------------------------------------------------------------"<<endl;

    string str5 = "./main ../../../../DataSet/circuit5M.mtx";

    const char *command5 = str5.c_str();
	  for(int i=0;i<5;i++){
         system(command5);
    }

    cout<<"--------------------------------------------------------------------"<<endl;
    
    string str6 = "./main ../../../../DataSet/G3_circuit.mtx";

    const char *command6 = str6.c_str();
	  for(int i=0;i<5;i++){
         system(command6);
    }

    cout<<"--------------------------------------------------------------------"<<endl;
    
    string str7 = "./main ../../../../DataSet/hollywood-2009.mtx";

    const char *command7 = str7.c_str();
	  for(int i=0;i<1;i++){
         system(command7);
    }

    cout<<"--------------------------------------------------------------------"<<endl;
    
    string str8 = "./main ../../../../DataSet/poisson3Da.mtx";

    const char *command8 = str8.c_str();
	  for(int i=0;i<5;i++){
         system(command8);
    }

    cout<<"--------------------------------------------------------------------"<<endl;
    

    string str9 = "./main ../../../../DataSet/rgg_n_2_24_s0.mtx";

    const char *command9 = str9.c_str();
	  for(int i=0;i<1;i++){
         system(command9);
    }

    cout<<"--------------------------------------------------------------------"<<endl;*/
    

	return 0;
}