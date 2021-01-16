#ifndef BASIS_PURSUIT_C_H
#define BASIS_PURSUIT_C_H



#include "Matrix.h"
#include "PDCD_I.h"
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

//This class solves problem of the form f(x)+g(x) under the constraint Mx=c by SMART_CD;
// where f(x)=0
//and g(x)=frac{lambda2}{2}|x|_2+lambda1|x|_1.





template<typename L, typename D>
class Basis_pursuit_c: public PDCD_I<L, D>
{
private:


  D lambda1;

  D lambda2;

  Matrix<L,D> my_M;

  D val_lambda_f;



protected:

public:

  Basis_pursuit_c(const char* Matrix_file,D val_lambda1, D val_lambda2)
  :PDCD_I<L,D>(),my_M(Matrix_file)
  {
    lambda1=val_lambda1;
    lambda2=val_lambda2;
    val_lambda_f=0;
    this->mu_f=0;
    this->mu_g=lambda2;
  }

  L get_n(){return my_M.nfeatures;}
  L get_m(){return my_M.nsamples;}

  inline D gradient_of_phi_j(D x1, L i){
    return 0;
  }



  inline D value_of_g_i(D x, L i){
    return lambda2*x*x/2+lambda1*fabs(x);
  }

  inline D value_of_phi_j(D x1, L i){return 0;}


  inline D prox_of_g_i(D x1,D x2,D x3, L i){
    return compute_one_step(1.0/x2, x1, x3)-x3;
  }

  D compute_one_step(D tau, D u, D x){
    D new_x;
    if(x>tau*(lambda1+u))
    new_x=(x-tau*(lambda1+u))/(1+lambda2*tau);
    else if(x<tau*(u-lambda1))
    new_x=(x-tau*(u-lambda1))/(1+lambda2*tau);
    else
    new_x=0;
    return new_x;
  }


  inline D value_of_phistar_i(D x,L i) {return 0;}


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


  inline D distance_to_subgradient_of_g2(vector<D> & gra, vector<D> & x){
    D res=0;
    D tmp;
    D xv;
    for(L i=0;i<this->n;i++){
      tmp=-gra[i];
      xv=x[i];
      if(xv>0) res+=(tmp-lambda2*xv-lambda1)*(tmp-lambda2*xv-lambda1);
      else if(xv<0) res+=(tmp-lambda2*xv+lambda1)*(tmp-lambda2*xv+lambda1);
      else if(tmp-lambda2*xv>lambda1) res+=(tmp-lambda2*xv-lambda1)*(tmp-lambda2*xv-lambda1);
      else if(tmp-lambda2*xv<-lambda1) res+=(tmp-lambda2*xv+lambda1)*(tmp-lambda2*xv+lambda1);
      else res+=0;
    }
    return sqrt(res);
  }





  void SMART_CD_solver(D beta_0, vector<D> & x0,vector<D> & y0,L val_tau, L max_nb_outer, L p_N, string filename1, D time){
    this->PDCD_I_solver(beta_0, x0,y0, val_tau, max_nb_outer,  p_N, val_lambda_f, filename1,time);
  }




};

#endif
