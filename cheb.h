#ifndef CHEB_H
#define CHEB_H

#include "Eigen/Dense"

Eigen::MatrixXd Pre_chebD(int );
Eigen::MatrixXd Pre_chebXi(int );


void cheb(int , Eigen::MatrixXd &, Eigen::MatrixXd &);

#endif