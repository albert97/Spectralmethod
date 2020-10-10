#include <iostream>
#include <vector>
#include <math.h>
#include "cheb.h"

//differentiation matrix to work out Chebyshev grid


//Let us staet with a 2D matrix

Eigen::MatrixXd Pre_chebD(int N){
	
	Eigen::MatrixXd D = Eigen::MatrixXd::Constant(N+1,N+1, 0);
	
	return D;
}

Eigen::MatrixXd Pre_chebXi(int N){
	
    	
	Eigen::MatrixXd Xi = Eigen::MatrixXd::Constant(N+1,1, 0);
	
	return Xi;
	
}


void cheb(int N, Eigen::MatrixXd  &D, Eigen::MatrixXd &Xi){
	
	if (N == 0){
		
		D =  Eigen::MatrixXd::Constant(N+1,N+1, 1);
		Xi = Eigen::MatrixXd::Constant(N+1,1, 1);
		
	}
	else{
	    	//Declare the transpose of X

	   	Eigen::ArrayXXf Xj = Eigen::ArrayXXf:: Zero(N+1,1);
    	
		// To generate Chebyshev points, we let xj= transpose(cos(􏰒pi􏰔􏰒􏰃􏰍 j /N)􏰓􏰑􏰘)
	
		for (int j =0; j<N+1; j++){
     	
		     Xj(j,0) = cos((M_PI * j)/N);
		     Xi(j,0) = Xj(j,0);
	     
		}
	
	
	   	 //The off diagonal elements of Chebyshev matrix has a "c" coefficient 􏰚 􏰜􏰅􏰘
	
		 Eigen:: MatrixXd C= Eigen:: MatrixXd::Constant(N+1, 1, 1);
		 C(0,0)=2;
		 C(N,0)=2;
 	
		 Eigen::ArrayXXf Crighthand (N+1,1);
		 Crighthand.col(0) = Eigen::ArrayXf::LinSpaced(N+1, 0, N);
	 
		 for (int i =0; i< N+1;i++){
	 	
			 C(i,0) = C(i,0)* pow(-1, Crighthand(i,0));
		
		 }
	
  	
	  	  //Work out dX through Xi and Xj where Xj is a transpose of Xi
	 
		 Eigen::ArrayXXf X= Eigen::ArrayXXf:: Zero(N+1,N+1);
	 
		 for (int i =0; i<N+1;i++){
	 	
			 X.col(i) = Xj.col(0);
		
		 }
   	
	  	
		Eigen::ArrayXXf XTranpose= Eigen::ArrayXXf:: Zero(N+1,N+1);
		XTranpose= X.transpose() ;
   	

		Eigen::ArrayXXf dX= X-XTranpose;

	
	
	   	 //Find the off-diagonal elements
	
	
		Eigen::MatrixXd  Dnumerator = Eigen:: MatrixXd::Constant(N+1,N+1, 0);
		Eigen::MatrixXd  Ddenomenator = Eigen:: MatrixXd::Identity(N+1,N+1);
	
		Eigen:: MatrixXd Ctemp;
		Eigen:: MatrixXd Ctemp2;
		Ctemp.array()= 1 / C.array();
		Ctemp2 =Ctemp.transpose();
	
   	
	 	Dnumerator = C * Ctemp2 ;
	
   	
	
		Eigen::MatrixXd Ddenomenatortemp (N+1,N+1);
	
		for (int i =0; i< N+1;i++){
			for (int j =0 ; j<N+1;j++){
			
				Ddenomenatortemp(i,j)=dX(i,j)+Ddenomenator(i,j);
			
			}
		}
	
   
		for (int i =0; i< N+1;i++){
			for (int j =0 ; j<N+1;j++){
			
				D(i,j)=Dnumerator(i,j)/Ddenomenatortemp(i,j) ;
			
			}
		}
	
   	

	   	 //Find the diagonal entries
	
		Eigen::MatrixXd DiagonalEntry = Eigen:: MatrixXd::Identity(N+1,N+1);
	
		Eigen::MatrixXd Dtranpose = D.transpose();
	
   	
	
		Eigen::MatrixXd SumColDtranspose(1,N+1);
		for (int i =0; i<N+1;i++){
		
			SumColDtranspose(0, i)=(Dtranpose.col(i)).sum();
		
		}

	
		for (int i=0;i<N+1 ;i++){
			for (int j =0;j<N+1;j++){
				if (i==j){
				
					DiagonalEntry (i,j) = SumColDtranspose(0,i);
				
				}
			
			}
		}

		D = D - DiagonalEntry;
	
		
	}

    
    
}
