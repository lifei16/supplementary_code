#ifndef FUSED_LASSO_E_H
#define FUSED_LASSO_E_H



#include "Matrix.h"
#include "LADMM.h"
#include <string>
#include <gsl/gsl_rng.h>
#include <gsl/gsl_randist.h>
#include <stdio.h>      /* printf */
#include <time.h>
#include <fstream>
#include <algorithm>
#include <iomanip>
#include <ctime>
#include <math.h>

//This class solves problem: min_x f(x)+ g(x)+ h(Mx) by LADMM;
// where g(x)=\frac{lambda2}{2}\|x\|_2^2 +lambda1\|x\|_1;
// f(x)= 1/2*\|Ax- b\|_2^2= \sum_i phi_i(A_ix);
// f_i(x)= 0.5*(x- b_i)^2, h_i(x)= \|x\|, Mx= (x^1- x^2, ..., x_n- x_1).





template<typename L, typename D>
class Fused_lasso_e: public LADMM<L, D>
{
private:


D lambda1;

D lambda2;

D lambda3;

Matrix<L,D> my_M;

Matrix<L,D> my_A;




protected:

public:

  Fused_lasso_e(const char* Matrix_file,const char* Matrix_file2, D val_lambda1, D val_lambda2, D val_lambda3)
  :LADMM<L,D>(),my_A(Matrix_file), my_M(Matrix_file2)
  {
    lambda1=val_lambda1;
    lambda2=val_lambda2;
    lambda3=val_lambda3;
    this->lambda_f= 1;
  }
  
Fused_lasso_e(const char* Matrix_file,D val_lambda1, D val_lambda2, D val_lambda3)
  :LADMM<L,D>(),my_A(Matrix_file)
  {
    lambda1=val_lambda1;
    lambda2=val_lambda2;
    lambda3=val_lambda3;
    this->lambda_f= 1;
    my_M.construct_fused_matrix(my_A);
  }

  L get_n(){return my_M.nfeatures;}
  L get_m(){return my_M.nsamples;}

  inline D gradient_of_f_j(D x1, L i){
      return x1- my_A.b[i];
  }



  inline D value_of_g_j(D x, L i){
        return lambda2*x*x/2+lambda1*fabs(x);
      }

  inline D value_of_f_j(D x1, L i){
  		return 0.5*(x1- my_A.b[i])*(x1- my_A.b[i]);
	}

  inline D value_of_h_j(D x, L i){
  	return lambda3*fabs(x);
  }


  inline D prox_of_h_j(D x1,D x2, L i){
  	D newx;
  	if (x1- lambda3/x2> my_M.b[i]){
  		newx= x1- lambda3/x2;
	  }
	  else if(x1+ lambda3/x2< my_M.b[i]){
	  	newx= x1+ lambda3/x2;
	  }
	  else{
	  	newx= my_M.b[i];
	  }
	  return newx;
  }

  inline D prox_of_g_j(D x1,D x2, L i){
        D new_x;
        if(x1*x2> lambda1)
         new_x=(x1*x2- lambda1)/(x2+lambda2);
        else if(x1*x2< -lambda1)
         new_x=(x1*x2+ lambda1)/(x2+ lambda2);
        else
         new_x=0;
        return new_x;
      }

  inline void set_matrix_M(){
    this->data_M=my_M;
  }
  
  inline void set_matrix_A(){
    this->data_A=my_A;
  }
/*
  inline void set_matrix_A(){
    this->data_A.nsamples=0;
    this->data_A.nfeatures=this->data_M.nfeatures;
    this->data_A.nnz=0;
    this->data_A.ptr.resize(1,0);
    this->data_A.ptr_t.resize(this->data_M.nfeatures+1,0);
  }
*/
  
  
  
  

  void LADMM_solver(D beta_0,  D rho,vector<D> & x0,vector<D> & y0, vector<D> & lambda0, L max_nb_outer, L p_N_1, string filename1, D time){
    this->ADMM_solve_with_Linear(beta_0, rho,x0,y0, lambda0, max_nb_outer,  p_N_1, filename1, time);
}




};

#endif
