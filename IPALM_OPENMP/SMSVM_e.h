#ifndef SMSVM_E_H
#define SMSVM_E_H



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

//This class solves problem: min_{x,omega} \sum_i max(1- b_i^T(A_ix+ omega), 0)+ P(x) by LADMM
// where P(x)=\frac{lambda2}{2}\|x\|_2^2 +lambda1\|x\|_1.

// Let X= (x; omega), M= (b.*A b), b_i^T(A_ix+ omega)= M_iX
// phi_i(x)= 0, h_i(x)= max(1- x,0), g(X)= P(x).






template<typename L, typename D>
class SMSVM_e: public LADMM<L, D>
{
private:


D lambda1;

D lambda2;

Matrix<L,D> my_M;

D val_lambda_f;



protected:

public:

  SMSVM_e(const char* Matrix_file,D val_lambda1, D val_lambda2)
  :LADMM<L,D>(),my_M(Matrix_file)
  {
    lambda1=val_lambda1;
    lambda2=val_lambda2;
    this->lambda_f=0;
    my_M.matrix_modified_to_svm(my_M); 
      
  }

  L get_n(){return my_M.nfeatures;}
  L get_m(){return my_M.nsamples;}

  inline D gradient_of_f_j(D x1, L i){
      return 0;
  }



  inline D value_of_g_j(D x, L i){
  	    if(i== this->m_1-1){
  	    	return 0;
		  } 
		else{
           return lambda2*x*x/2+lambda1*fabs(x);
       }
   }

  inline D value_of_f_j(D x1, L i){return 0;}

  inline D value_of_h_j(D x, L i){
  	if(x<= 1.0){
  		return 1.0- x;
	  }
	  else{
	  	return 0;
	  }
  }
  

  inline D prox_of_g_j(D x1,D x2, L i){
  	    if(i== this->m_1 - 1){
  	    	return x1;
		  }
		else{
            D new_x;
        if(x1*x2> lambda1)
         new_x=(x1*x2- lambda1)/(x2+lambda2);
        else if(x1*x2< -lambda1)
         new_x=(x1*x2+ lambda1)/(x2+ lambda2);
        else
         new_x=0;
        return new_x;
        }
      }

      
// compute argmin_x h(x)+ x2/2*(x- x1)*(x- x1)
  inline D prox_of_h_j(D x1, D x2, L j){
  	if (x1+ 1.0/x2< 1.0){
  		return x1+ 1.0/x2;
   }
	  else if (x1> 1.0){
	  	return x1;
	  }
	  else{
	  	return 1.0;;
	  }
	}


  inline void set_matrix_M(){
    this->data_M=my_M;
  }

  inline void set_matrix_A(){
    this->data_A.nsamples=0;
    this->data_A.nfeatures=this->data_M.nfeatures;
    this->data_A.nnz=0;
    this->data_A.ptr.resize(1,0);
    this->data_A.ptr_t.resize(this->data_M.nfeatures+1,0);
  }

  
  
  

  void LADMM_solver(D beta_0,  D rho,vector<D> & x0,vector<D> & y0, vector<D> & lambda0, L max_nb_outer, L p_N_1, string filename1, D time){
    this->ADMM_solve_with_Linear(beta_0, rho,x0,y0, lambda0, max_nb_outer,  p_N_1, filename1, time);
}




};

#endif
