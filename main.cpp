#include <iostream>
#include <vector>
#include <math.h>
#include <fstream>
#include "Eigen/Dense"
#include "cheb.h"

std::fstream plotfile1;
std::fstream plotfile2;
int main(void){
	
	//if N􏰚􏰚􏰃􏰖 ==0 return D =0􏰚􏰃􏰘 x=1􏰚􏰄􏰘

	 
	 
		std::vector <double> xx ;//This the an array with x used for construct the exact solution 
	
	 	for (signed int i =-100; i<100; i++){
	 		double temp = double(i)/100.0;
	 		xx.push_back(temp);

	 	}
	
	 	std::vector <double> uu ;//This the an array for the exact solution 
	
	 	for (int i =0; i<xx.size();i++){
		
	 		uu.push_back(exp(xx[i])*sin(5*xx[i]));
	 	}
		
	
	 	//for (int N =10; N<20 ;N++){
		int N=20;
		Eigen::MatrixXd D = Pre_chebD(N);
		Eigen::MatrixXd Xi = Pre_chebXi(N);
	
		cheb(N, D, Xi);
	 		//}
	 	std::vector <double> u; 	//Soluton of approximation 	
	 	//<<x->numDimensions<<endl;
	
	 	for(int i = 0;i<N+1;i++){
	 		
	 		u.push_back(exp(Xi(i,0))*sin(5*Xi(i,0)));
	 	} 
	 	
	
	 	//plotting the exact solution 
	 	plotfile1.open("smoothfuncExactSol.dat",  std::fstream::out);
	
	         if(plotfile1.fail()){

	               std::cout<<"Error opening file"<< std::endl;
	               exit(0);
	           }
	
	           for (int i=0;i<xx.size();i++){
		
	 		 plotfile1<<xx[i]<<"\t"<<uu[i]<<"\t"<< std::endl;
                
	             }
	    
	 	plotfile1.close();
	
	 	
	 	plotfile2.open("smoothfunc.dat",  std::fstream::out);
	
	         if(plotfile2.fail()){

	 		std::cout<<"Error opening file"<< std::endl;
	               exit(0);
	           }

	         for (int i=0;i<N+1;i++){
		
	 		 plotfile2<<Xi(i,0)<<"\t"<<u[i]<<"\t"<< std::endl;
                
	 	 }
	 	 plotfile2.close();
	 	 
	        return 0;


}