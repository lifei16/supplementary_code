#ifndef ALM_APG_H
#define ALM_APG_H



#include "Matrix.h"
#include "APPROX2.h"
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


//This class implements the method ASGARD_DL.

/*
The optimization problem to solve is of the form min f(x)+g(x)+h(Mx) ;
where f(x)=\sum_{j=1}^m lambda_f[j] \phi_j(<A_j,x>)
and g(x)=sum_{i=1}^n g_i(x_i). We all assume that each \phi_j is 1-smooth.
*/

// Each subproblem solves problem of the form f(x)+ h_{beta_s}(Mx; lambda_s) +g(x) by APG;







template<typename L, typename D>
class ALM_APG: public APPROX2<L, D>
{
private:

  std::vector<D> Au;
  std::vector<D> Az;

  std::vector<D> Mu;
  std::vector<D> Mz;


  std::vector<D> Mx_s;
  std::vector<D> Ax_s;
  std::vector<D> lambda_f;








protected:
  Matrix<L,D> data_A;

  Matrix<L, D> data_M;

  std::vector<D> x_s;

  std::vector<D> lambda_s;

  std::vector<D> M_tlambda_s;

  D beta_s;

  D epsilon_s;

  D m_s;

  D m_0;

  L m_1;

  L m_2;


  D omega;

  D max_Lf_s;

  D max_M_s;

  D Lf;

  D L_m;

  L val_mu_f;
  L val_mu_g;

  std::vector<D> dual_alpha1;   //dual_alpha1[j]=gradient_of_phi_j(<A_j,x>)

  std::vector<D> dual_alpha2;  //dual_alpha2[j]=(<M_j,x>-c_j)

  std::vector<D> baralpha1;  //baralpha1=-sum_{j=1}^m_1 lambda_f[j]*dual_alpha1[j]*A_j

  std::vector<D> baralpha2;   //baralpha2=-sum_{j=1}^m_2 1/beta_s*dual_alpha2[j]*M_j


  D function_value;

  L print_every_N_ALM_APG;

  D running_time_ALM_APG;

  L nb_outer_iters;

  std::vector<D> gradient_of_f;

  ofstream samp_ALM_APG;
public:



  virtual inline D gradient_of_phi_j(D, L){return D(NULL);}

  virtual inline D value_of_g_i(D, L){return D(NULL);}
  virtual inline D value_of_phi_j(D, L){return D(NULL);}
  virtual inline D value_of_h_j(D, L){return D(NULL);}
  virtual inline D value_of_h_star_j(D, L){return D(NULL);}

  virtual inline D prox_of_g_i(D,D,D, L){return D(NULL);}
  virtual inline D prox_of_h_star_j(D,D, L){return D(NULL);}

  virtual inline D feasible_dual(std::vector<D> &,std::vector<D> &) {return D(NULL);}

  virtual inline D value_of_phistar_i(D,L) {return D(NULL);}
  virtual inline D value_of_g_tilde_star(D, vector<D> &,vector<D> &){return D(NULL);}

  virtual inline void set_matrix_M(){}

  virtual inline void set_matrix_A(){}



  virtual inline D distance_to_subgradient_of_g(){return D(NULL);}
  //virtual inline D distance_to_subgradient_of_g(std::vector<D> &, std::vector<D> & ){return D(NULL);}

  ALM_APG(const char* Matrix_file, D val_lambda_f)
  :data_A(), data_M()
  {

  }

  ALM_APG()
  :data_A(), data_M()
  {

  }


  L get_nb_features(){
    return data_A.get_d();
  }


  //This function computes \nabla_i f(gamma u +z)+ \nabla_i h_beta(gamma Mu+ Mz)
  inline D partial_i_of_f(L i){
    D res=0;
    for (L k = data_A.ptr_t[i]; k < data_A.ptr_t[i + 1];k++)
    {
      L j=data_A.col_idx[k];
      D tmp=lambda_f[j]*gradient_of_phi_j(this->gamma*Au[j]+Az[j], j);
      res+=data_A.A_t[k]*tmp;
    }
    for (L k = data_M.ptr_t[i]; k < data_M.ptr_t[i + 1];k++)
    {
      L j=data_M.col_idx[k];
      D tmp=prox_of_h_star_j((this->gamma*Mu[j]+Mz[j])/beta_s+ lambda_s[j],beta_s,j);
      res+=data_M.A_t[k]*tmp;
    }
    return res;
  }

  //This function computes argmin{x1 t+x2 t^2/2+g_i(x3+t)}
  //inline D compute_prox(D x1, D x2, D x3, L i){
  //return prox_operA_tor(x1,x2,x3,i);
  //}


  void compute_x(){
    for(L i=0;i<this->n;i++){
      this->x[i]=this->gamma*this->u[i]+this->z[i];
      x_s[i]=this->x[i];
    }
    for(L j=0;j<m_2;j++)
    Mx_s[j]=this->gamma*this->Mu[j]+this->Mz[j];
    for(L j=0;j<m_1;j++)
    Ax_s[j]=this->gamma*this->Au[j]+this->Az[j];
  }

  inline void compute_primal_value() {

    D res=0;
    for(L i=0;i<this->n;i++){
      this->x[i]=this->gamma*this->u[i]+this->z[i];
      res+=value_of_g_i(this->x[i],i);
    }

    for(L j=0;j<m_1;j++){
      res+=lambda_f[j]*value_of_phi_j(this->gamma*Au[j]+Az[j],j);
    }
    for(L j=0;j<m_2;j++){
      D tmp= prox_of_h_star_j((this->gamma*Mu[j]+Mz[j])/beta_s+ lambda_s[j],beta_s,j);
      res+=(this->gamma*Mu[j]+Mz[j])*tmp- value_of_h_star_j(tmp,j)- beta_s/2*(tmp- lambda_s[j])*(tmp- lambda_s[j]);
    }
    this->primal_value=res;
  }

  inline void compute_dual_value(){
    for(L j=0;j<m_1;j++){
      dual_alpha1[j]=gradient_of_phi_j(this->gamma*Au[j]+Az[j], j);
    }
    for(L j=0;j<m_2;j++){
      dual_alpha2[j]=prox_of_h_star_j((this->gamma*Mu[j]+Mz[j])/beta_s+ lambda_s[j],beta_s,j);
    }
    for(L i=0;i<this->n;i++){
      baralpha1[i]=0;
      for(L k = data_A.ptr_t[i]; k < data_A.ptr_t[i + 1];k++){
        L j=data_A.col_idx[k];
        baralpha1[i]-=dual_alpha1[j]*data_A.A_t[k]*lambda_f[j];
      }
      baralpha2[i]=0;
      for(L k = data_M.ptr_t[i]; k < data_M.ptr_t[i + 1];k++){
        L j=data_M.col_idx[k];
        baralpha2[i]-=dual_alpha2[j]*data_M.A_t[k];
      }
    }
    D res=0;
    D scal=feasible_dual(baralpha1, baralpha2);
    for(L j=0;j<m_1;j++)
    {
      res-=value_of_phistar_i(scal*dual_alpha1[j],j)*lambda_f[j];
    }
    for(L j=0;j<m_2;j++)
    {
      res-=value_of_h_star_j(scal*dual_alpha2[j],j)+ beta_s/2*(scal*dual_alpha2[j]- lambda_s[j])*(scal*dual_alpha2[j]- lambda_s[j]);
    }
    res-=value_of_g_tilde_star(scal,baralpha1,baralpha2);
    this->dual_value=res;
  }

  inline void compute_function_value() {
    D res=0;
    for(L i=0;i<this->n;i++){
      res+=value_of_g_i(x_s[i],i);
    }
    for(L j=0;j<m_1;j++){
      res+=lambda_f[j]*value_of_phi_j(Ax_s[j],j);
    }
    for(L j=0;j<m_2;j++){
      res+=value_of_h_j(Mx_s[j],j);
    }
    function_value=res;
  }


  D compute_prox(D x1,D x2,D x3, L i){
    D res=prox_of_g_i(x1+beta_s*(x3-x_s[i]), beta_s+x2,x3,i);
    if(fabs(res)>1e10)
    {cout<<"larger than 1e10:"<<res<<endl;
    cout<<beta_s<<"; "<<x1<<"; "<<x3<<"; "<<x_s[i]<<"; "<<i<<endl;
    D res2=0;
    for (L k = data_A.ptr_t[i]; k < data_A.ptr_t[i + 1];k++)
    {
      L j=data_A.col_idx[k];
      D tmp=lambda_f[j]*gradient_of_phi_j(this->gamma*Au[j]+Az[j], j);
      res2+=data_A.A_t[k]*tmp;
      cout<<"Ax["<<j<<"]="<<this->gamma*Au[j]+Az[j]<<endl;
    }
    for (L k = data_M.ptr_t[i]; k < data_M.ptr_t[i + 1];k++)
    {
      L j=data_M.col_idx[k];
      D tmp=(this->gamma*Mu[j]+Mz[j]-data_M.b[j])/beta_s+lambda_s[j];
      cout<<"Mu["<<j<<"]="<<Mu[j]<<"; gamma="<<this->gamma<<endl;
      res2+=data_M.A_t[k]*tmp;
    }
  }
  return res;
}

inline void compute_gradient_norm(){
  std::vector<D> Tx(this->n,0);
  this->do_single_step_prox(Tx);
  D res=0;
  for(L i=0;i<this->n;i++)
  res+=(Tx[i]-this->x[i])*(Tx[i]-this->x[i]);
  this->gradient_norm=sqrt(res);
}

inline void set_v()
{
  this->v.resize(this->n,0);
  D maxv=0;
  D minv=std::numeric_limits<double>::max();
  D sumv=0;
  D sumvi1=0;
  L sumw=0;
  L maxw=0;
  L minw=this->n;
  this->sumofLi=0;
  for(L j=0;j<m_1;j++)
  {
    sumw+=data_A.w_t[j];
    maxw=max(maxw,data_A.w_t[j]);
    minw=min(minw,data_A.w_t[j]);
  }
  cout<<"sumw="<<sumw<<";  maxw="<<maxw<<"; minw="<<minw<<endl;
  for(L i=0;i<this->n;i++)
  {
    D vi=0;
    D vi1=0;
    for (L k = data_A.ptr_t[i]; k < data_A.ptr_t[i + 1];k++)
    {
      L j=data_A.col_idx[k];
      vi+=(1.+(data_A.w_t[j]-1.)*(this->tau-1.)/max(this->n-1.,1.))*data_A.A_t[k]*data_A.A_t[k]*lambda_f[j];
      vi1+=data_A.A_t[k]*data_A.A_t[k]*lambda_f[j];
    }
    for (L k = data_M.ptr_t[i]; k < data_M.ptr_t[i + 1];k++)
    {
      L j=data_M.col_idx[k];
      vi+=(1.+(data_M.w_t[j]-1.)*(this->tau-1.)/max(this->n-1.,1.))*data_M.A_t[k]*data_M.A_t[k]/beta_s;
      vi1+=data_M.A_t[k]*data_M.A_t[k]/beta_s;
    }
    this->v[i]=vi;
    sumv+=vi;
    sumvi1+=vi1;
    if(maxv<vi) maxv=vi;
    if(minv>vi) minv=vi;
  }
  if(this->tau==this->n){
    for(L i=0;i<this->n;i++)
    this->v[i]=sumvi1;
  }
  this->sumofLi=sumvi1;
  cout<<"  max of v: "<<maxv<<" ;  min of v: "<<minv<<" ;  sumofv: "<<sumv<<" sumofLi="<<this->sumofLi<<endl;
}

inline void compute_m0(D beta0, D val_omega,  L val_tau){
  cout<<"beta0="<<beta0<<"val_omega="<<val_omega<<"; val_tau="<<val_tau<<endl;
  D tmp1=sqrt(2*max_M_s)/beta0;
  cout<<tmp1<<" ;"<<val_omega<<endl;
  m_s=this->n/val_tau*30;
  cout<<" m_s="<<m_s<<endl;
}

inline void rescale_Matrix(){
  Lf=0;
  max_Lf_s=0;
  max_M_s=0;
  L_m=0;
  cout<<"this->n="<<this->n<<endl;
  for(L i=0;i<this->n;i++)
  {
    D vi=0;
    for (L k = data_A.ptr_t[i]; k < data_A.ptr_t[i + 1];k++)
    {
      L j=data_A.col_idx[k];
      Lf+=data_A.A_t[k]*data_A.A_t[k]*lambda_f[j];
      vi+=(1.+(data_A.w_t[j]-1.)*(this->tau-1.)/max(this->n-1.,1.))*data_A.A_t[k]*data_A.A_t[k]*lambda_f[j];
    }
    if(vi>max_Lf_s) max_Lf_s=vi;
    D vi2=0;
    for (L k = data_M.ptr_t[i]; k < data_M.ptr_t[i + 1];k++)
    {
      L j=data_M.col_idx[k];
      L_m+=data_M.A_t[k]*data_M.A_t[k];
      vi2+=(1.+(data_M.w_t[j]-1.)*(this->tau-1.)/max(this->n-1.,1.))*data_M.A_t[k]*data_M.A_t[k];
    }
    if(vi2>max_M_s) {max_M_s=vi2;cout<<vi2<<" ";}
  }
cout<<"max_M_s="<<max_M_s<<"; L_m="<<L_m<<"; max_Lf_s="<<max_Lf_s<<"; L_f="<<Lf<<endl;
}

inline void set_p(){
  this->proba_vector.resize(this->n,(0.0+this->tau)/this->n);
  this->max_p=(0.0+this->tau)/this->n;
}



void testing(){
  std::vector<D> nablaf(this->n);
  D res1=0;
  compute_Mty();
  for(L i=0;i<this->n;i++)
  {
    nablaf[i]=0;
    for (L k = data_M.ptr_t[i]; k < data_M.ptr_t[i + 1];k++)
    {
      L j=data_M.col_idx[k];
      D tmp=lambda_s[j];
      nablaf[i]+=data_M.A_t[k]*tmp;
    }
    res1+=(nablaf[i]-M_tlambda_s[i])*(nablaf[i]-M_tlambda_s[i]);
  }
}



inline void update_z_coordinate( L i, D dz){
  this->z[i]+=dz;
  for (L k = data_A.ptr_t[i]; k < data_A.ptr_t[i + 1];k++)
  {
    L j=data_A.col_idx[k];
    Az[j]+=dz*data_A.A_t[k];
  }
  for (L k = data_M.ptr_t[i]; k < data_M.ptr_t[i + 1];k++)
  {
    L j=data_M.col_idx[k];
    Mz[j]+=dz*data_M.A_t[k];
  }
}

inline void update_x_coordinate( L i, D dx){
  this->x[i]+=dx;
  x_s[i]=this->x[i];
  L j;
  for (L k = data_M.ptr_t[i]; k < data_M.ptr_t[i + 1];k++)
  {
    j=data_M.col_idx[k];
    Mx_s[j]+=dx*data_M.A_t[k];
  }
  for (L k = data_A.ptr_t[i]; k < data_A.ptr_t[i + 1];k++)
  {
    j=data_A.col_idx[k];
    Ax_s[j]+=dx*data_A.A_t[k];
  }
}

inline void update_u_coordinate( L i, D du){
  this->u[i]+=du;
  for (L k = data_A.ptr_t[i]; k < data_A.ptr_t[i + 1];k++)
  {
    L j=data_A.col_idx[k];
    Au[j]+=du*data_A.A_t[k];
  }
  for (L k = data_M.ptr_t[i]; k < data_M.ptr_t[i + 1];k++)
  {
    L j=data_M.col_idx[k];
    Mu[j]+=du*data_M.A_t[k];
  }
}

void compute_Au(){
  for(L j=0;j<m_1;j++)
  for(L k = data_A.ptr[j]; k < data_A.ptr[j + 1];k++){
    L kj=data_A.row_idx[k];
    Au[j]+=this->u[kj]*data_A.A[k];
  }
}

void compute_Mu(){
  for(L j=0;j<m_2;j++)
  for(L k = data_M.ptr[j]; k < data_M.ptr[j + 1];k++){
    L kj=data_M.row_idx[k];
    Mu[j]+=this->u[kj]*data_M.A[k];
  }
}

void compute_Az(){
  for(L j=0;j<m_1;j++)
  for(L k = data_A.ptr[j]; k < data_A.ptr[j + 1];k++){
    L kj=data_A.row_idx[k];
    Az[j]+=this->z[kj]*data_A.A[k];
  }
}
void compute_Mz(){
  for(L j=0;j<m_2;j++)
  for(L k = data_M.ptr[j]; k < data_M.ptr[j + 1];k++){
    L kj=data_M.row_idx[k];
    Mz[j]+=this->z[kj]*data_M.A[k];
  }
}

void compute_Au(vector<D> & x0){
  for(L j=0;j<m_1;j++)
  for(L k = data_A.ptr[j]; k < data_A.ptr[j + 1];k++){
    L kj=data_A.row_idx[k];
    Au[j]+=x0[kj]*data_A.A[k];
  }

}

void compute_Az(vector<D> & x0){
  for(L j=0;j<m_1;j++)
  for(L k = data_A.ptr[j]; k < data_A.ptr[j + 1];k++){
    L kj=data_A.row_idx[k];
    Az[j]+=x0[kj]*data_A.A[k];
  }
}
void compute_Mu(vector<D> & x0){
  for(L j=0;j<m_2;j++)
  for(L k = data_M.ptr[j]; k < data_M.ptr[j + 1];k++){
    L kj=data_M.row_idx[k];
    Mu[j]+=x0[kj]*data_M.A[k];
  }

}


void compute_Mz(vector<D> & x0){
  for(L j=0;j<m_2;j++)
  for(L k = data_M.ptr[j]; k < data_M.ptr[j + 1];k++){
    L kj=data_M.row_idx[k];
    Mz[j]+=x0[kj]*data_M.A[k];
  }
}

void compute_Mx(vector<D> & x0){
  for(L j=0;j<m_2;j++)
  for(L k = data_M.ptr[j]; k < data_M.ptr[j + 1];k++){
    L kj=data_M.row_idx[k];
    Mx_s[j]+=x0[kj]*data_M.A[k];
  }
}

void compute_Ax(vector<D> & x0){
  for(L j=0;j<m_1;j++)
  for(L k = data_A.ptr[j]; k < data_A.ptr[j + 1];k++){
    L kj=data_A.row_idx[k];
    Ax_s[j]+=x0[kj]*data_A.A[k];
  }
}


void compute_Mty(){
  M_tlambda_s.clear();
  M_tlambda_s.resize(this->n,0);
  for(L i=0;i<this->n;i++)
  {
    for(L k = data_M.ptr_t[i]; k < data_M.ptr_t[i + 1];k++){
      L j=data_M.col_idx[k];
      M_tlambda_s[i]+=lambda_s[j]*data_M.A_t[k];
    }
  }

}

void update_y(){
  D res=0;
  for(L j=0;j<m_2;j++)
  {
    D tmp=prox_of_h_star_j(Mx_s[j]/beta_s+ lambda_s[j],beta_s,j);
    res+=(tmp- lambda_s[j])*(tmp- lambda_s[j]);
    L i;
    for (L k = data_M.ptr[j]; k < data_M.ptr[j + 1];k++)
    {
      i=data_M.row_idx[k];
      M_tlambda_s[i]=tmp*data_M.A[k];
    }
    lambda_s[j]= tmp;
  }
}



void Initialize(D beta_0, D epsilon_0,  D val_omega, L val_tau, vector<D> & x0,vector<D> & y0, D val_lambda_f){
  cout<<"start initializing"<<endl;
  set_matrix_M();
  set_matrix_A();

  this->tau=val_tau;
  m_1=data_A.get_n();
  m_2=data_M.get_n();
  cout<<"m_1="<<m_1<<endl;
  cout<<"m_2="<<m_2<<endl;
  lambda_f.resize(m_1,val_lambda_f);
  dual_alpha1.resize(m_1,0);
  dual_alpha2.resize(m_2,0);
  this->n=data_A.nfeatures;
  baralpha1.resize(this->n,0);
  baralpha2.resize(this->n,0);
  gradient_of_f.resize(this->n,0);

  rescale_Matrix();

  beta_s=beta_0;


  epsilon_s=epsilon_0;
  compute_m0(beta_s,val_omega,val_tau);

  omega=val_omega;

  x_s.resize(this->n,0);
  lambda_s.resize(m_2,0);
  for(L i=0;i<this->n;i++)
  x_s[i]=x0[i];

  for(L j=0;j<m_2;j++)
  lambda_s[j]=y0[j];

  M_tlambda_s.resize(this->n,0);

  compute_Mty();

  Au.clear();
  Au.resize(m_1,0);
  Az.clear();
  Az.resize(m_1,0);
  Ax_s.clear();
  Ax_s.resize(m_1,0);
  Mu.clear();
  Mu.resize(m_2,0);
  Mz.clear();
  Mz.resize(m_2,0);
  Mx_s.clear();
  Mx_s.resize(m_2,0);

  compute_Az(x0);
  compute_Mz(x0);
  compute_Ax(x0);
  compute_Mx(x0);
}

void reset_everything(){
  epsilon_s/=omega;
  m_s*=omega;
  beta_s/=omega;

  compute_Mty();

  Au.clear();
  Au.resize(m_1,0);
  Az.clear();
  Az.resize(m_1,0);
  Mu.clear();
  Mu.resize(m_2,0);
  Mz.clear();
  Mz.resize(m_2,0);
  compute_Az(x_s);
  compute_Mz(x_s);
}



void APPROX_restart(L K, vector<D> & x0, L val_tau, L eval, L p_N,  L max_nb, D eps, string filename){
  std::vector<D> tmp_x(this->n);
  for(L i=0;i<this->n;i++){
    tmp_x[i]=x0[i];
  }
  L k=0;
  while(k<ceil(max_nb/(K+0.0))){
    k++;
    Au.clear();
    Au.resize(m_1,0);
    Az.clear();
    Az.resize(m_1,0);
    Mu.clear();
    Mu.resize(m_2,0);
    Mz.clear();
    Mz.resize(m_2,0);
    compute_Az(tmp_x);
    compute_Mz(tmp_x);
    cout<<"p_N="<<tmp_x[0]<<" "<<tmp_x[1]<<" "<<tmp_x[2]<<" "<<p_N<<endl;
    this->APPROX_MU(tmp_x, val_tau, 0, 0, eval, p_N, K, eps, filename,1);

    for(L i=0;i<this->n;i++){
      tmp_x[i]=this->x[i];
    }
    cout<<"restarting..."<<k<<"; "<<(max_nb/(K+0.0))<<endl;
    if(this->delta<eps) k=ceil(max_nb/(K+0.0))+1;
  }

}


inline void compute_and_record_res(){
  if(nb_outer_iters%print_every_N_ALM_APG==0){
    compute_function_value();
    cout<<setprecision(9)<<"Iteration: "<<nb_outer_iters<<"; time="<<running_time_ALM_APG<<"; function value="<<function_value<<endl;
    samp_ALM_APG<<setprecision(9)<< nb_outer_iters<<" "<<running_time_ALM_APG<<" "<<function_value<<" "<<endl;
  }
}




void ALM_APG_solve_with_APPROX(D beta_0, D epsilon_0,  D omega, vector<D> & x0,vector<D> & y0, L val_tau, L max_nb_outer, L p_N_1, L p_N_2, D val_lambda_f,string filename1, string filename2, D time){
  Initialize(beta_0, epsilon_0, omega,val_tau, x0, y0, val_lambda_f);

  nb_outer_iters=0;
  string sampname2="results/APPROXMU_"+filename2;
  this->samp.open(sampname2.c_str());
  filename1="results/DLRCSGR2_"+filename1;
  samp_ALM_APG.open(filename1.c_str());
  running_time_ALM_APG=0;
  print_every_N_ALM_APG=p_N_1;
  compute_and_record_res();
  D start;
  //start = std::clock();
  start=omp_get_wtime();
  std::vector<D> record_x_s(this->n);
  std::vector<D> record_x_sp1(this->n);
  cout<<"m0="<<ceil((m_s)/this->n*val_tau)<<endl;
  cout<<"ms="<<ceil(m_s/this->n*val_tau)<<"; beta_s="<<beta_s<<"; epsilon_s="<<epsilon_s<<endl;
  this->APPROX_MU(x_s, val_tau, val_mu_f, val_mu_g, 3, p_N_2, ceil((m_s)/this->n*val_tau), epsilon_s, filename2,1);
  compute_x();
  update_y();
  nb_outer_iters++;
  //running_time_ALM_APG+=( std::clock() - start ) / (double) CLOCKS_PER_SEC;
  running_time_ALM_APG+=omp_get_wtime()-start;
  compute_and_record_res();
  //start = std::clock();
  start=omp_get_wtime();
  reset_everything();
  //running_time_ALM_APG+=( std::clock() - start ) / (double) CLOCKS_PER_SEC;
  running_time_ALM_APG+=omp_get_wtime()-start;
  while(nb_outer_iters<max_nb_outer){
    //start = std::clock();
    start=omp_get_wtime();
    cout<<"ms="<<ceil(m_s/this->n*val_tau)<<"; beta_s="<<beta_s<<"; epsilon_s"<<epsilon_s<<endl;
    this->APPROX_MU(x_s, val_tau, val_mu_f, val_mu_g, 3, p_N_2, ceil(m_s/this->n*val_tau), epsilon_s, filename2,1);
    compute_x();
    update_y();
    nb_outer_iters++;
    //running_time_ALM_APG+=( std::clock() - start ) / (double) CLOCKS_PER_SEC;
    running_time_ALM_APG+=omp_get_wtime()-start;
    testing();
    compute_and_record_res();
    //start = std::clock();
    start=omp_get_wtime();
    reset_everything();
    //running_time_ALM_APG+=( std::clock() - start ) / (double) CLOCKS_PER_SEC;
    running_time_ALM_APG+=omp_get_wtime()-start;
    if (running_time_ALM_APG> time){
      break;
    }
  }

}






};

#endif
