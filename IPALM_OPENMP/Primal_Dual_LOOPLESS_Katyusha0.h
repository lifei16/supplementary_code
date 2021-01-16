#ifndef PRIMAL_DUAL_LOOPLESS_KATYUSHA0_H
#define PRIMAL_DUAL_LOOPLESS_KATYUSHA0_H

#include "Primal_Dual_LOOPLESS.h"


#include <string>
#include <gsl/gsl_rng.h>
#include <gsl/gsl_randist.h>
#include <stdio.h>      /* printf */
#include <time.h>
#include <fstream>
#include <algorithm>
#include <iomanip>
#include <ctime>
#include <sstream>


//This class implements the method loopless Katyusha with arbitrary sampling, in the paper L-SVRG and L-Katyusha with arbitrary sampling. arXiv:1906.01481


/*
The optimization problem to solve is:

\sum_{i=1}^n \lambda_i\phi_i(A_i^{\top} w)+ g(w)
Assumption 1: For each i, \phi_i is 1-smooth
By default, lambda_i=1/n for all i.
*
* The dual problem is
*         g*(-sum_{i=1}^n \lambda_i \alpha_i A_i)+\sum_{i=1}^n \lambda_i \phi*_i(\alpha_i)
*/

template<typename L, typename D>
class Primal_Dual_LOOPLESS_Katyusha0: public Primal_Dual_LOOPLESS<L, D>
{


protected:
  D theta1;
  D theta2;
  D theta3;
  D eta; //the stepsize parameter alpha in the paper



  std::vector<D> primal_y;     // the y variable

  std::vector<D> primal_z;     // the z variable

  std::vector<D> primal_w;     // the w variable


public:


  virtual inline D gradient_of_phi_i(D, L){return D(NULL);}
  virtual inline D gradient_of_gstar_j(D, L){return D(NULL);}
  virtual inline D value_of_phi_i(D, L) {return D(NULL);}
  virtual inline D value_of_g_j(D, L){return D(NULL);}
  virtual inline D value_of_phistar_i(D,L) {return D(NULL);}
  virtual inline D value_of_gstar(vector<D> &){return D(NULL);}
  virtual inline D value_of_gstar_minus(vector<D> &, D){return D(NULL);}
  virtual inline D feasible_dual(vector<D> &){return D(NULL);}
  virtual inline D compute_delta_alpha(D,D,L){return D(NULL);}
  virtual inline void set_auxiliary_v(){}
  virtual inline D get_lambda1(){return D(NULL);}
  virtual inline D get_lambda2(){return D(NULL);}
  virtual inline void compute_just_in_time_prox_grad(D, D, D &, L,L, D, D &, L){}

  Primal_Dual_LOOPLESS_Katyusha0()
  : Primal_Dual_LOOPLESS<L,D>()
  {

  }


  Primal_Dual_LOOPLESS_Katyusha0(const char* matrix_file, const char* vector_file)
  : Primal_Dual_LOOPLESS<L,D>(matrix_file, vector_file)
  {

  }

  Primal_Dual_LOOPLESS_Katyusha0(const char* matrix_file)
  : Primal_Dual_LOOPLESS<L,D>(matrix_file)
  {

  }


  void set_theta(){
    theta2=this->L2/2/max(this->Lf,this->L2);
    if(this->Lf<=this->L2/this->p){
      D tmp=this->mu/this->L2/this->p;
      cout<<"this->L2="<<this->L2<<endl;
      cout<<"tmp="<<tmp<<endl;
      if(tmp>=1){
        theta1=theta2;
      }
      else{
        theta1=sqrt(tmp)*theta2;
      }
    }else
    theta1=min(sqrt(this->mu/this->Lf),this->p/2);
    theta3=1-theta1-theta2;
    eta=1./(theta1*(this->Lf+2*max(this->L2,this->Lf)));
    cout<<"this->Lf="<<this->Lf<<"; this->L2="<<this->L2<<"; theta1="<<theta1<<"; theta2="<<theta2<<"; theta3="<<theta3<<endl;
    cout<<"eta="<<eta<<endl;
  }

  inline D compute_current_xj_value_without_update(D ui, L j, L t0 , L t1)
  {
    D w=primal_w[j];
    D y=primal_y[j];
    D z=primal_z[j];
    compute_just_in_time_prox_grad(eta, ui , z, t0, t1, w, y, j);
    D x=theta1*z+theta2*w+theta3*y;
    return x;
  }

  inline D compute_current_xj_value(D ui, L j, L t0 , L t1)
  {
    D w=primal_w[j];
    D y=primal_y[j];
    D z=primal_z[j];
    compute_just_in_time_prox_grad(eta, ui , z, t0, t1, w, y, j);
    D x=theta1*z+theta2*w+theta3*y;
    this->primal_x[j]=x;
    primal_z[j]=z;
    primal_y[j]=y;
    return x;
  }

  inline void set_stepsize(){
    primal_y.clear();
    primal_z.clear();
    primal_w.clear();
    primal_y.resize(this->nfeatures,0);
    primal_z.resize(this->nfeatures,0);
    primal_w.resize(this->nfeatures,0);
    for(L j=0; j<this->nfeatures; j++){
      primal_y[j]=this->primal_x[j];
      primal_z[j]=this->primal_x[j];
      primal_w[j]=this->primal_x[j];
    }
    cout<<"w size="<<primal_w.size()<<endl;
    set_theta();
    cout<<"w 2="<<primal_w.size()<<endl;
  }
/*
  inline void update_x(D xj,L j){
    this->primal_x[j]=xj;
    primal_z[j]=z;
    primal_y[j]=y;
  }
  */


  inline void update_baralpha()
  {
    D yi=gsl_rng_uniform(this->rng);
    if(yi<=this->p){
      this->nb_iters+=this->nsamples/this->tau;
      #pragma omp parallel for
      for(L i=0;i<this->nsamples;i++)
      {
        D aitg=this->compute_AiTxk(i);
        D deltaalphai=gradient_of_phi_i(aitg,i)-this->dual_alpha[i];
        this->dual_alpha[i]+=deltaalphai;
        for (L k = this->ptr[i]; k < this->ptr[i + 1];k++)
        {
          L j=this->row_idx[k];
          this->baralpha[j]+=this->lambda_f[i]*deltaalphai*this->A[k];
        }
      }

      for(L j=0;j<this->nfeatures;j++)
      primal_w[j]=this->primal_x[j];
    }
  }



  void loopless_Katyusha(vector<D> & x0, vector<D> & w0, string filename, vector<D> & L_phi, D val_mu, D val_epsilon, L max_nb, L nb_tau, L nb_c, L u, L p_mod, D scal_p)
  {
    filename="Katyusha"+filename;
    this->loopless( x0,  w0,  filename,  L_phi,  val_mu,  val_epsilon, max_nb, nb_tau, nb_c, u, p_mod,scal_p);
  }

  void loopless_Katyusha2(vector<D> & x0, vector<D> & w0, string filename, vector<D> & L_phi, D val_mu, D val_epsilon, L max_nb, L nb_tau, L nb_c, L u, L p_mod, D scal_p)
  {
    filename="Katyusha"+filename;
    this->loopless2( x0,  w0,  filename,  L_phi,  val_mu,  val_epsilon, max_nb, nb_tau, nb_c, u, p_mod,scal_p);
  }












};

#endif /* MIN_SMOOTH_CONVEX_H */
