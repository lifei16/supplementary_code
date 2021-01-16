#ifndef ALM_KATYUSHA_H
#define ALM_KATYUSHA_H

#include "Primal_Dual_LOOPLESS_Katyusha0.h"


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
//#include "cmd_line.h"


//This class implements the method IPALM_Katyusha

/*
The optimization problem to solve is:

min \sum_i f^i(A_ix)+ \sum_i h^i(M_ix)+ P(x)
Assumption 1: For each i, f_i is 1-smooth, P(x) is seperable

*/

// Each subproblem solves problem of the form sum_i f^i(A_ix)+ sum_i h^i_{beta_s}(M_ix;lambda_s^i) +P(x)+ beta_s/2|x-x_s|^2 by L_Katyuhsa


template<typename L, typename D>
class ALM_Katyusha: public Primal_Dual_LOOPLESS_Katyusha0<L, D>
{
private:

  std::vector<D> Ax_s;
  std::vector<D> old_Ax_s;

  std::vector<D> Mx_s;
  std::vector<D> old_Mx_s;

protected:
  Matrix<L, D>	data_A;

  Matrix<L, D> data_M;

  std::vector<D> x_s;

  std::vector<D> old_x_s;

  std::vector<D> lambda_s;

  std::vector<D> old_lambda_s;


  D beta_s;

  D epsilon_s;

  D m_s;

  D m_0;

  L m_1;

  L m_2;

  L m_3;

  D eta;

  D rho;

  D tau_s;

  std::vector<D> L_phi;

  D function_value;

  L print_every_N_ALM;

  D running_time_ALM;

  L nb_outer_iters;

  ofstream samp_ALM;


public:

  D mu_g;

  D lambda1;

  D lambda2;

  D L_h;


  virtual inline D value_of_f_j(D, L){return D(NULL);}
  virtual inline D value_of_h_j(D, L){return D(NULL);}
  virtual inline D gradient_of_f_j(D, L){return D(NULL);}
  virtual inline D value_of_f_star_j(D, L){return D(NULL);}
  virtual inline D value_of_h_star_j(D, L){return D(NULL);}
  virtual inline D prox_of_h_star_j(D,D, L){return D(NULL);}
  virtual inline D value_of_P_j(D, L){return D(NULL);}
  virtual inline D prox_of_P_j(D, D, L){return D(NULL);}
  virtual inline D feasible_dual(vector<D> &){return D(NULL);}
  virtual inline void compute_just_in_time_prox_grad_without_x_s(D, D, D &, L,L, D, D &, L){}
  virtual inline void rescale(){}
  virtual inline void set_matrix_M(){}
  virtual inline void set_matrix_A(){}

  ALM_Katyusha()
  : Primal_Dual_LOOPLESS_Katyusha0<L,D>(),data_A(),data_M()
  {
  	this->gamma= 1;
  }

  ALM_Katyusha(const char* matrix_file, const char* vector_file)
  : Primal_Dual_LOOPLESS_Katyusha0<L,D>(matrix_file, vector_file),data_A(), data_M()
  {
    this->gamma=1;
  }

  ALM_Katyusha(const char* matrix_file)
  : Primal_Dual_LOOPLESS_Katyusha0<L,D>(matrix_file),data_A(), data_M()
  {
    this->gamma=1;
  }
/*
  ALM_Katyusha(const char* matrix_file, const char* matrix_file2)
  : Primal_Dual_LOOPLESS_Katyusha0<L,D>(),data_A(matrix_file), data_M(matrix_file2)
  {
  	this->matrix_merge(data_A,data_M);
    this->gamma=1;
  }
 */
  inline void set_L_phi(){
  	for(int i=m_3; i<m_2+ m_3; i++){
  		L_phi[i]= this->nsamples;
	  }
  	for(int i= 0; i<m_3; i++){
    	L_phi[i]= this->nsamples/beta_s;
	}
  }

  L get_nb_features(){
    return data_A.get_d();
  }

  inline D value_of_phi_i(D x, L i) {
    if ( i>= data_M.nsamples){
      return value_of_f_j(x, i- data_M.nsamples);
    }
    else{
      D tmp= prox_of_h_star_j(1.0/beta_s*x+ lambda_s[i],beta_s,i);
      return beta_s*(x*tmp- value_of_h_star_j(tmp,i)- beta_s/2*(tmp- lambda_s[i])*(tmp- lambda_s[i]));
    }
  }



  inline D gradient_of_phi_i(D x, L i){
    if (i>= data_M.nsamples){
      return gradient_of_f_j(x,i- data_M.nsamples);
    }
    else{
      D tmp= prox_of_h_star_j(1.0/beta_s*x+ lambda_s[i],beta_s,i);
      return beta_s*tmp;
    }

  }


  inline D value_of_phistar_i(D x, L i) {
    if ( i>= data_M.nsamples){
      return value_of_f_star_j(x,i- data_M.nsamples);
    }
    else{
      return beta_s*value_of_h_star_j(x/beta_s,i)+ 0.5*(x- beta_s*lambda_s[i])*(x- beta_s*lambda_s[i]);
    }
  }

  inline D value_of_g_j(D x, L j){
  	return value_of_P_j(x,j)+ beta_s/2*(x- x_s[j])*(x- x_s[j]);
  }

  inline D value_of_gstar(vector<D> &x){
  	D res= 0;
  	L l= x.size();
  	D tmp= 0;
  	for (L i= 0; i< l; i++){
  		tmp= prox_of_P_j(1.0/beta_s*x[i]+ x_s[i],beta_s,i);
        res+= x[i]*tmp- value_of_P_j(tmp,i)- beta_s/2*(tmp- x_s[i])*(tmp- x_s[i]);
	  }
  	return res;
  }

  inline D gradient_of_gstar_j(D x, L j){
  	return prox_of_P_j(1.0/beta_s*x+ x_s[j],beta_s,j);
  }

  inline D value_of_gstar_minus(vector<D> &x, D scal){
  	L l= x.size();
  	vector<D> tmp(l,0);
  	for (L j= 0;j< l; j++){
  		tmp[j]= -scal*x[j];
	  }
  	return value_of_gstar(tmp);
  }
 /* inline D get_lambda1(){return lambda1;}
  inline D get_lambda2(){return lambda2;}
*/

  void compute_x(){
    for (L i=0; i< this->nfeatures; i++){
      x_s[i]= this->primal_x[i];
    }
    for (L j=m_3; j< m_2+ m_3; j++){
      Ax_s[j- m_3]= compute_AiTx_s(j);
    }
    for (L j=0; j< m_3; j++){
      Mx_s[j]= compute_AiTx_s(j);
	}
  }

  D compute_AiTx_s(L i){
    D res=0;
    for (L k = this->ptr[i]; k < this->ptr[i + 1];k++)
    {
      L j=this->row_idx[k];
      res+= this->A[k]*x_s[j];
    }
    return res;
  }

  inline void compute_just_in_time_prox_grad(D tau, D u, D &x, L t0,L t1, D w, D &y, L j){
    D u_s=u-tau_s*beta_s*x_s[j];
    compute_just_in_time_prox_grad_without_x_s(tau, u_s,  x,  t0,  t1,  w, y,j);
  }

  void compute_function_value(){
    D res= 0;
    for (L i= 0; i< this->nfeatures; i++){
      res+= value_of_P_j(x_s[i],i);
    }
    for (L i= 0; i< data_M.nsamples; i++){
      res+= value_of_h_j(Mx_s[i],i);
    }
    for (L i= 0; i<data_A.nsamples; i++){
      res+= value_of_f_j(Ax_s[i],i);
    }
    function_value= res;
  }

  /*void compute_infeasibility(){
  D tmp1= 0;
  D tmp2= 0;
  for (L j=0; j< this->nsamples; j++){
  if (1- Ax_s[j]> tmp1){
  tmp1= 1- Ax_s[j];
}
tmp2+= max(0.0,1- Ax_s[j]);
}
residual1= tmp1;
residual2= tmp2/this->nsamples;
}*/



inline void compute_m0(D beta0, D val_eta, D val_rho, L val_tau){
  cout<<"beta0="<<beta0<<"val_eta="<<val_eta<<"; val_rho="<<val_rho<<"; val_tau="<<val_tau<<endl;
  m_s= this->nsamples/val_tau*1e+6;

  //D tmp7=1-val_tau/this->n*sqrt(beta0/(max_Lf_s+max_M_s/beta0+beta0));
  //cout<<"tmp7: "<<tmp7<<"; "<<val_tau/(this->n+0.)*sqrt(beta0/(max_Lf_s+max_M_s/beta0+beta0))<<endl;
  //m_0=(2*log(val_rho)+log(val_eta)+log(2))/log(tmp7);

  cout<<" m_s="<<m_s<<endl;
}

void update_y(){
  for(L j=0;j<m_3;j++)
  {
    D aiTx= Mx_s[j];
    lambda_s[j]= prox_of_h_star_j(1.0/beta_s*aiTx+ lambda_s[j],beta_s,j);
  }
}

void Initialize(D beta_0, D epsilon_0,  D val_eta, D val_rho,L val_tau, vector<D> & x0,vector<D> & y0){
  cout<<"start initializing"<<endl;
  set_matrix_M();
  set_matrix_A();


  this->tau=val_tau;
  m_1=data_A.get_d();
  m_2=data_A.get_n();
  m_3=data_M.get_n();
  cout<<"m_1="<<m_1<<endl;
  cout<<"m_2="<<m_2<<endl;
  cout<<"m_3="<<m_3<<endl;

  beta_s=beta_0;
  //tau_s= this->nsamples;
  tau_s= 1;

  epsilon_s=epsilon_0;
  compute_m0(beta_s,val_eta,val_rho,val_tau);

  eta=val_eta;
  rho=val_rho;

  x_s.resize(m_1,0);
  old_x_s.resize(m_1,0);
  lambda_s.resize(m_3,0);
  old_lambda_s.resize(m_3,0);
  for(L i=0;i<m_1;i++){
    x_s[i]=x0[i];
    old_x_s[i]= x0[i];
  }
  for(L j=0;j<m_3;j++){
    lambda_s[j]=y0[j];
    old_lambda_s[j]= y0[j];
  }

  Ax_s.clear();
  Ax_s.resize(m_2,0);
  old_Ax_s.clear();
  old_Ax_s.resize(m_2,0);
  Mx_s.clear();
  Mx_s.resize(m_3,0);
  old_Mx_s.clear();
  old_Mx_s.resize(m_3,0);
  L_phi.clear();
  L_phi.resize(m_2+ m_3,0);
  for(int i=0; i< m_3; i++){
  	L_phi[i]= this->nsamples/beta_0;
  }
  for(int i=m_3; i< m_2+ m_3; i++){
  	L_phi[i]= this->nsamples;
  }
}

void reset_everything(){
  epsilon_s*=rho;
  beta_s*=eta;

}

void update_m_s(L val_tau){
  D tmp= 0;
  D tmpx= 0;
  D tmpy= 0;
  D tmpy2= 0;
  for(L j=0;j<m_3;j++)
  {
    D tmp_p= prox_of_h_star_j(Mx_s[j]/beta_s/eta+ lambda_s[j], beta_s/eta, j);
    tmpy+=(tmp_p- lambda_s[j])*(tmp_p- lambda_s[j]);
    tmp+= (lambda_s[j]- old_lambda_s[j])*(lambda_s[j]- old_lambda_s[j]);
    tmpy2+= (beta_s*eta*lambda_s[j]- beta_s*old_lambda_s[j])*(beta_s*eta*lambda_s[j]- beta_s*old_lambda_s[j]);
  }
  for (L i= 0; i< m_1;i++){
    tmpx= (x_s[i]- old_x_s[i])*(x_s[i]- old_x_s[i]);
  }
  D tmp4= 2*epsilon_s+ beta_s*tmp;
  D tmp5= (1- eta)*beta_s/2*tmpy;
  D tmp6= sqrt(tmp)*beta_s*(1+ eta)*L_h;
  D tmp7= sqrt(tmp)*sqrt(tmpy2);
  D tmp1= log(tmp4+ tmp5+ tmp6+ tmp7+ 2/(2*eta- 1)*tau_s*beta_s/2*tmpx);
  D tmp2= log(epsilon_s)+ log(rho)- log(2);
  D tmp3= sqrt(lambda2/this->sumLi);
  m_s= ceil((tmp1- tmp2)/log(2)*4*max(1.0*this->nsamples,1/tmp3)/val_tau);
}

inline void compute_and_record_res(){
  if(nb_outer_iters%print_every_N_ALM==0){
    compute_function_value();
    cout<<setprecision(9)<<"Iteration: "<<nb_outer_iters<<"; time="<<running_time_ALM<< "; function value="<<function_value<<endl;
    samp_ALM<<setprecision(9)<<nb_outer_iters<<" "<<running_time_ALM<<" "<< function_value<<" "<<endl;
  }
}

void ALM_solve_with_L_Katyusha(D beta_0, D epsilon_0,  D eta, D rho,vector<D> & x0,vector<D> & y0, L val_tau, L max_nb_outer, L p_N_1, L p_N_2, string filename1, string filename2, D time){
  Initialize(beta_0, epsilon_0, eta, rho,val_tau, x0, y0);
  nb_outer_iters=0;
  //string sampname2= ALGparam.data_dir +"/results/L_Katyusha_"+filename2;
  string sampname2= "results/L_Katyusha_"+filename2;
  this->samp.open(sampname2.c_str());
  //filename1= ALGparam.data_dir +"/results/ALM_"+filename1;
  filename1= "results/ALM_"+filename1;
  samp_ALM.open(filename1.c_str());
  running_time_ALM=0;
  print_every_N_ALM=p_N_1;
  this->set_print_every_N(p_N_2);
  compute_and_record_res();
  D start;
  //start = std::clock();
  start=omp_get_wtime();
  //cout<<"m0="<<ceil(m_s/this->n*val_tau)<<endl;
  cout<<"m_s= "<< ceil(m_s/this->nsamples*val_tau)<<"; beta_s="<<beta_s<<"; epsilon_s="<<epsilon_s<<endl;
  //this->L_Katyusga_MU(x_s, val_tau, val_mu_f, val_mu_g, 3, p_N_2, ceil(m_s/this->n*val_tau), epsilon_s, filename2,1);
  rescale();
  set_L_phi();
  this->loopless_Katyusha(y0, x_s, filename2, L_phi, mu_g+ tau_s*beta_s, epsilon_s, ceil(m_s/this->nsamples*val_tau), val_tau, 1, 1, 0, 1);
  for(L i=0;i<m_1;i++){
    old_x_s[i]=x_s[i];
  }
  for(L i=0;i<m_2;i++){
    old_Ax_s[i]=Ax_s[i];
  }
  compute_x();
  for(L i=0;i<m_3;i++){
    old_lambda_s[i]=lambda_s[i];
  }
  update_y();
  update_m_s(val_tau);
  nb_outer_iters++;
  //running_time_ALM+=( std::clock() - start ) / (double) CLOCKS_PER_SEC;
  running_time_ALM+=omp_get_wtime()-start;
  compute_and_record_res();
  //start = std::clock();
  start=omp_get_wtime();
  reset_everything();
  //running_time_ALM+=( std::clock() - start ) / (double) CLOCKS_PER_SEC;
  running_time_ALM+=omp_get_wtime()-start;
  while(nb_outer_iters<max_nb_outer){
    //start = std::clock();
    start=omp_get_wtime();
    cout<<"ms="<<ceil(m_s/this->nsamples*val_tau)<<"; beta_s="<<beta_s<<"; epsilon_s"<<epsilon_s<<endl;
    rescale();
    set_L_phi();
    this->loopless_Katyusha(y0, x_s, filename2, L_phi, mu_g+ tau_s*beta_s, epsilon_s, ceil(m_s/this->nsamples*val_tau), val_tau, 1, 1, 0, 1);
    for(L i=0;i<m_1;i++){
      old_x_s[i]=x_s[i];
    }
    for(L i=0;i<m_2;i++){
      old_Ax_s[i]=Ax_s[i];
    }
    for(L i=0;i<m_3;i++){
      old_lambda_s[i]=lambda_s[i];
    }
    compute_x();
    update_y();
    update_m_s(val_tau);
    nb_outer_iters++;
    //running_time_ALM+=( std::clock() - start ) / (double) CLOCKS_PER_SEC;
    running_time_ALM+=omp_get_wtime()-start;
    compute_and_record_res();
    //start = std::clock();
    start=omp_get_wtime();
    reset_everything();
    //running_time_ALM+=( std::clock() - start ) / (double) CLOCKS_PER_SEC;
    running_time_ALM+=omp_get_wtime()-start;
    if (running_time_ALM> time){
      break;
    }
  }

}
};

#endif /* MIN_SMOOTH_CONVEX_H */
