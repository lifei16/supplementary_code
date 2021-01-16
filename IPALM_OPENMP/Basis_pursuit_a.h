#ifndef BASIS_PURSUIT_A_H
#define BASIS_PURSUIT_A_H



#include "Matrix.h"
#include "APPROX2.h"
#include "ALM_I_APPROX.h"
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

/*This class solves problem of the form f(x)+g(x) under the constraint Mx=c by IPALM_APPROX;
where f(x)=0 and g(x)=frac{lambda2}{2}|x|_2+lambda1 |x|_1.
*/




template<typename L, typename D>
class Basis_pursuit_a: public ALM_I_APPROX<L, D>
{
private:


  D lambda1;

  D lambda2;

  Matrix<L,D> my_M;

  Matrix<L,D> my_A;

  D val_lambda_f;



protected:

public:

  Basis_pursuit_a(const char* Matrix_file,D val_lambda1, D val_lambda2)
  :ALM_I_APPROX<L,D>(), my_M(Matrix_file)
  {
    lambda1=val_lambda1;
    lambda2=val_lambda2;
    val_lambda_f=0;
    this->val_mu_f=0;
    this->val_mu_g=lambda2;
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

  inline D feasible_dual(vector<D> & x, vector<D> & y){

    if(lambda2+this->beta_s>0)
    {
      return 1;
    }
    else
    {
      D scal=1;
      L l=x.size();
      for(L i=0;i<l;i++)
      if(fabs(x[i]+y[i])>lambda1)
      scal=min(lambda1/fabs(x[i]+y[i]),scal);
      return scal;
    }
  }

  inline D rescale(vector<D> & x, vector<D> & y){
    D scal=1;
    L l=x.size();
    for(L i=0;i<l;i++){
      if(fabs(x[i]+y[i])>lambda1){
        scal=min(lambda1/fabs(x[i]+y[i]),scal);
      }
    }
    return scal;
  }

  inline D value_of_phistar_i(D x,L i) {return 0;}


  inline D value_of_g_tilde_star(D scal,vector<D> & x, vector<D> & y){

    if(lambda2+this->tau_s*this->beta_s>0)
    {
      D res=0;
      D ip=0;
      L l=x.size();
      for(L i=0;i<l;i++)
      {
        ip=max(fabs(x[i]+y[i]+this->tau_s*this->beta_s*this->x_s[i])-lambda1,0.0);
        res+=ip*ip;
      }
      D xsnorm=0;
      for(L i=0;i<l;i++)
      xsnorm+=this->tau_s*this->beta_s/2*this->x_s[i]*this->x_s[i];
      return res/2/(lambda2+this->tau_s*this->beta_s)-xsnorm;
    }
    else
    {
      return 0;
    }
  }

  inline D value_of_g_star(D scal,vector<D> & x, vector<D> & y){
    return 0;
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

  inline D distance_to_subgradient_of_g(){
    D res=0;
    D tmp;
    D xv;
    for(L i=0;i<this->n;i++){
      tmp=-this->gradient_of_f[i]-this->M_tlambda_s[i];
      xv=this->x_s[i];
      if(xv>0) res+=(tmp-lambda2*xv-lambda1)*(tmp-lambda2*xv-lambda1);
      else if(xv<0) res+=(tmp-lambda2*xv+lambda1)*(tmp-lambda2*xv+lambda1);
      else if(tmp-lambda2*xv>lambda1) res+=(tmp-lambda2*xv-lambda1)*(tmp-lambda2*xv-lambda1);
      else if(tmp-lambda2*xv<-lambda1) res+=(tmp-lambda2*xv+lambda1)*(tmp-lambda2*xv+lambda1);
      else res+=0;
    }
    return sqrt(res);
  }

  inline D subgradient_of_Tx(){
    D res=0;
    D tmp;
    D xv;
    for(L i=0;i<this->n;i++){
      xv=this->Tx[i];
      tmp=this->gradient_of_Tx[i]+ lambda2*xv+this->tau_s*this->beta_s*(xv- this->x_s[i]);
      if(xv>0) res+=(tmp+ lambda1)*(tmp+ lambda1);
      else if(xv<0) res+=(tmp-lambda1)*(tmp-lambda1);
      else if(tmp>lambda1) res+=(tmp-lambda1)*(tmp-lambda1);
      else if(tmp<-lambda1) res+=(tmp+lambda1)*(tmp+lambda1);
      else res+=0;
    }
    return sqrt(res);
  }





  void ALM_I_APPROX_solver(D beta_0, D epsilon_0,  D eta, D rho,vector<D> & x0,vector<D> & y0,L val_tau, L max_nb_outer, L p_N_1, L p_N_2,string filename1, string filename2, D time){
    this->ALM_I_APPROX_solve_with_APPROX(beta_0, epsilon_0,  eta, rho,x0,y0, val_tau, max_nb_outer,  p_N_1,  p_N_2,  val_lambda_f, filename1,  filename2, time);
  }




};

#endif
