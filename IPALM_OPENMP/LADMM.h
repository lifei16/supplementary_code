#ifndef LADMM_H
#define LADMM_H

#include "Matrix.h"


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
#include <omp.h>
//#include "cmd_line.h"


//This class implements the method LADMM

/*
The optimization problem to solve is:

min \sum_i f^i(A_ix)+ \sum_i h^i(y^i)+ g(x)
s.t. Mx= y
Assumption 1: For each i, f_i is smooth, g(x) is seperable
*/

template<typename L, typename D>
class LADMM
{
private:

  std::vector<D> Ax;
  std::vector<D> old_Ax;

  std::vector<D> Mx;
  std::vector<D> old_Mx;

  std::vector<D> Mty;
  std::vector<D> old_Mty;

  std::vector<D> Mtlambda;
  std::vector<D> old_Mtlambda;

  std::vector<D> gradient;

  std::vector<D> MtMx;

protected:
  Matrix<L, D> data_A;

  Matrix<L, D> data_M;

  std::vector<D> x;

  std::vector<D> old_x;

  std::vector<D> y;

  std::vector<D> old_y;

  std::vector<D> lambda;

  std::vector<D> old_lambda;


  D beta;

  L m_1;

  L m_2;

  L m_3;

  D rho;

  D tau;

  D sigma;

  D L_phi;

  D function_value;

  D infeas;

  L print_every_N_ADMM;

  D running_time_ADMM;

  L nb_outer_iters;

  ofstream samp_ADMM;


public:
  D lambda_f;

  D mu_g;

  D lambda1;

  D lambda2;

  D L_h;


  virtual inline D value_of_f_j(D, L){return D(NULL);}
  virtual inline D value_of_h_j(D, L){return D(NULL);}
  virtual inline D gradient_of_f_j(D, L){return D(NULL);}
  virtual inline D prox_of_h_j(D,D, L){return D(NULL);}
  virtual inline D value_of_g_j(D, L){return D(NULL);}
  virtual inline D prox_of_g_j(D, D, L){return D(NULL);}
  virtual inline void set_matrix_M(){}
  virtual inline void set_matrix_A(){}


/*
  LADMM(const char* matrix_file, const char* matrix_file2)
  : Primal_Dual_LOOPLESS_Katyusha0<L,D>(),data_A(matrix_file), data_M(matrix_file2)
  {
  	this->matrix_merge(data_A,data_M);
    this->gamma=1;
  }
 */
  inline void set_L_phi(){
  	if (data_A.nsamples== 0){
  		L_phi= 0;
	  }
	else{
  		L_phi= compute_lambda_max_A(10);
  	}
  }

  D compute_lambda_max_A(L K){

    std::vector<D> bk(data_A.nfeatures);
    for (L j=0;j<data_A.nfeatures;j++)
    {
      bk[j]=1;
    }
    std::vector<D> yk(data_A.nsamples);
    D normk;
    D tmp;
    for(L kk=0;kk<K;kk++){
      for (L i=0;i<data_A.nsamples;i++){
        tmp=0;
        for (L k = data_A.ptr[i]; k < data_A.ptr[i + 1];	k++)
        {
          L j=data_A.row_idx[k];
          tmp+=data_A.A[k]*bk[j];
        }
        yk[i]=tmp;
      }
      normk=0;
      for (L j=0;j<data_A.nfeatures;j++){
        bk[j]=0;
        for (L k = data_A.ptr_t[j]; k < data_A.ptr_t[j + 1];	k++)
        {
          L i=data_A.col_idx[k];
          bk[j]+=data_A.A_t[k]*yk[i]*lambda_f;
        }
        normk+=bk[j]*bk[j];
      }
      normk=sqrt(normk);
      for (L j=0;j<data_A.nfeatures;j++)
      {bk[j]=bk[j]/normk; }
    }
    cout<<endl;
    D res=0;
    normk=0;
    for (L i=0;i<data_A.nsamples;i++){
      tmp=0;
      for (L k = data_A.ptr[i]; k < data_A.ptr[i + 1];	k++)
      {
        L j=data_A.row_idx[k];
        tmp+=data_A.A[k]*bk[j];
      }
      yk[i]=tmp;
      normk+=yk[i]*yk[i];
    }
    std::vector<D> bk2(data_A.nfeatures);
    for (L j=0;j<data_A.nfeatures;j++){
      bk2[j]=0;
      for (L k = data_A.ptr_t[j]; k < data_A.ptr_t[j + 1];	k++)
      {
        L i=data_A.col_idx[k];
        bk2[j]+=data_A.A_t[k]*yk[i]*lambda_f;
      }
    }
    for (L j=0;j<data_A.nfeatures;j++)
    res+=bk2[j]*bk[j];
    return res;
  }

D compute_lambda_max_M(L K){

    std::vector<D> bk(data_M.nfeatures);
    for (L j=0;j<data_M.nfeatures- 1;j++)
    {
      bk[j]=1;
    }
    bk[data_M.nfeatures- 1]= 2;
    std::vector<D> yk(data_M.nsamples);
    D normk;
    D tmp;
    for(L kk=0;kk<K;kk++){
      for (L i=0;i<data_M.nsamples;i++){
        tmp=0;
        for (L k = data_M.ptr[i]; k < data_M.ptr[i + 1];	k++)
        {
          L j=data_M.row_idx[k];
          tmp+=data_M.A[k]*bk[j];
        }
        yk[i]=tmp;
      }
      normk=0;
      for (L j=0;j<data_M.nfeatures;j++){
        bk[j]=0;
        for (L k = data_M.ptr_t[j]; k < data_M.ptr_t[j + 1];	k++)
        {
          L i=data_M.col_idx[k];
          bk[j]+=data_M.A_t[k]*yk[i];
        }
        normk+=bk[j]*bk[j];
      }
      normk=sqrt(normk);
      for (L j=0;j<data_M.nfeatures;j++)
      {bk[j]=bk[j]/normk; }
    }
    D res=0;
    normk=0;
    for (L i=0;i<data_M.nsamples;i++){
      tmp=0;
      for (L k = data_M.ptr[i]; k < data_M.ptr[i + 1];	k++)
      {
        L j=data_M.row_idx[k];
        tmp+=data_M.A[k]*bk[j];
      }
      yk[i]=tmp;
      normk+=yk[i]*yk[i];
    }
    std::vector<D> bk2(data_M.nfeatures);
    for (L j=0;j<data_M.nfeatures;j++){
      bk2[j]=0;
      for (L k = data_M.ptr_t[j]; k < data_M.ptr_t[j + 1];	k++)
      {
        L i=data_M.col_idx[k];
        bk2[j]+=data_M.A_t[k]*yk[i];
      }
    }
    for (L j=0;j<data_M.nfeatures;j++)
    res+=bk2[j]*bk[j];
    return res;
  }

  L get_nb_features(){
    return data_A.get_d();
  }


 /* inline D get_lambda1(){return lambda1;}
  inline D get_lambda2(){return lambda2;}
*/

  void update_x(){
    #pragma omp parallel for
    for (L i= 0; i< m_1; i++){
    	gradient[i]= 0;
    	for (L k = data_A.ptr_t[i]; k < data_A.ptr_t[i + 1];k++)
    	{
      		L j=data_A.col_idx[k];
      		gradient[i]+= data_A.A_t[k]*gradient_of_f_j(Ax[j],j);
    	}
    	MtMx[i]= 0;
    	for (L k = data_M.ptr_t[i]; k < data_M.ptr_t[i + 1];k++)
    	{
      		L j=data_M.col_idx[k];
      		MtMx[i]+= data_M.A_t[k]*Mx[j];
    	}
    	x[i]= prox_of_g_j(x[i]- (Mtlambda[i]+ gradient[i]+ beta*MtMx[i]- beta*Mty[i])/(beta*tau+ L_phi), beta*tau+ L_phi, i);
	}
  #pragma omp parallel
  {
    #pragma omp  for
	for (L i= 0; i< m_2; i++){
		Ax[i]= compute_AiTx(i);
	}
  #pragma omp  for
	for (L i= 0; i< m_3; i++){
		Mx[i]= compute_MiTx(i);
	}
   }
  }

  void update_y(){
    #pragma omp parallel for
  	for (L i=0; i< m_3; i++){
  		y[i]= prox_of_h_j(y[i]+ (lambda[i]+ beta*Mx[i]- beta*y[i])/(beta*sigma), beta*sigma, i);
	  }
    #pragma omp parallel for
	for (L i=0; i< m_1; i++){
		Mty[i]= compute_MTiTy(i);
	}
  }

  D compute_AiTx(L i){
    D res=0;
    for (L k = data_A.ptr[i]; k < data_A.ptr[i + 1];k++)
    {
      L j=data_A.row_idx[k];
      res+= data_A.A[k]*x[j];
    }
    return res;
  }

  D compute_MiTx(L i){
    D res=0;
    for (L k = data_M.ptr[i]; k < data_M.ptr[i + 1];k++)
    {
      L j=data_M.row_idx[k];
      res+= data_M.A[k]*x[j];
    }
    return res;
  }

  D compute_MTiTy(L i){
    D res=0;
    for (L k = data_M.ptr_t[i]; k < data_M.ptr_t[i + 1];k++)
    {
      L j=data_M.col_idx[k];
      res+= data_M.A_t[k]*y[j];
    }
    return res;
  }

  D compute_MTiTlambda(L i){
    D res=0;
    for (L k = data_M.ptr_t[i]; k < data_M.ptr_t[i + 1];k++)
    {
      L j=data_M.col_idx[k];
      res+= data_M.A_t[k]*lambda[j];
    }
    return res;
  }


  void compute_function_value(){
    D res= 0;
    for (L i= 0; i< data_A.nfeatures; i++){
      res+= value_of_g_j(x[i],i);
    }
    for (L i= 0; i< data_M.nsamples; i++){
      res+= value_of_h_j(Mx[i],i);
    }
    for (L i= 0; i<data_A.nsamples; i++){
      res+= value_of_f_j(Ax[i],i);
    }
    function_value= res;
  }

  void compute_infeasibility(){
  	D res = 0;
  	for (L i= 0; i< data_M.nsamples; i++){
      res+= (Mx[i]- data_M.b[i])*(Mx[i]- data_M.b[i]);
    }
    infeas= sqrt(res);
  }



void update_lambda(){
  #pragma omp  for
  for(L j=0;j<m_3;j++)
  {
    lambda[j]= lambda[j]+ beta*(Mx[j]- y[j]);
  }
  #pragma omp  for
  for (L i=0; i< m_1; i++){
		Mtlambda[i]= compute_MTiTlambda(i);
	}
}

void Initialize(D beta_0, D val_rho, vector<D> & x0,vector<D> & y0, vector<D> & lambda0){
  cout<<"start initializing"<<endl;
  set_matrix_M();
  set_matrix_A();


  m_1=data_A.get_d();
  m_2=data_A.get_n();
  m_3=data_M.get_n();
  cout<<"m_1="<<m_1<<endl;
  cout<<"m_2="<<m_2<<endl;
  cout<<"m_3="<<m_3<<endl;

  beta=beta_0;
  //tau= data_A.nsamples;

  tau= 1.02*compute_lambda_max_M(10);
  sigma= 1;

  rho=val_rho;

  x.resize(m_1,0);
  old_x.resize(m_1,0);
  y.resize(m_3,0);
  old_y.resize(m_3,0);
  lambda.resize(m_3,0);
  old_lambda.resize(m_3,0);
  for(L i=0;i<m_1;i++){
    x[i]=x0[i];
  }
  for(L j=0;j<m_3;j++){
    y[j]=y0[j];
  }
  for(L j=0;j<m_3;j++){
    lambda[j]=lambda0[j];
  }

  Ax.clear();
  Ax.resize(m_2,0);
  Mx.clear();
  Mx.resize(m_3,0);
  Mty.clear();
  Mty.resize(m_1,0);
  Mtlambda.clear();
  Mtlambda.resize(m_1,0);
  gradient.clear();
  gradient.resize(m_1,0);
  MtMx.clear();
  MtMx.resize(m_1,0);
  L_phi= 0;
  set_L_phi();
  cout<< "L_phi= "<< L_phi << " beta= "<< beta<< " rho= "<< rho<< " tau= "<< tau<< " sigma= "<< sigma<<endl;
}

void reset_everything(){
  beta*=rho;

}


inline void compute_and_record_res(){
  if(nb_outer_iters%print_every_N_ADMM==0){
    compute_function_value();
    compute_infeasibility();
    cout<<setprecision(9)<<"Iteration: "<<nb_outer_iters<<"; time="<<running_time_ADMM<< "; function value="<<function_value<< "; infeasibility="<< infeas<< endl;
    samp_ADMM<<setprecision(9)<<nb_outer_iters<<" "<<running_time_ADMM<<" "<< function_value<<" "<< infeas<< endl;
  }
}

void ADMM_solve_with_Linear(D beta_0, D val_rho,vector<D> & x0,vector<D> & y0, vector<D> & lambda0, L max_nb_outer, L p_N_1, string filename1, D time){
  Initialize(beta_0,val_rho, x0, y0, lambda0);
  nb_outer_iters=0;
  //string sampname2= ALGparam.data_dir +"/results/L_Katyusha_"+filename2;
  //filename1= ALGparam.data_dir +"/results/ADMM_"+filename1;
  filename1= "results/ADMM_"+filename1;
  samp_ADMM.open(filename1.c_str());
  running_time_ADMM=0;
  print_every_N_ADMM=p_N_1;
  compute_and_record_res();
  D start;
  D res_x, res_y, res_l;
  /*
  for(L i=0;i<m_3;i++){
    old_lambda[i]=lambda[i];
  }
  */
  while(nb_outer_iters<max_nb_outer){
    //rescale();
    for(L i=0;i<m_1;i++){
      old_x[i]=x[i];
    }
    for(L i=0;i<m_3;i++){
      old_y[i]=y[i];
    }
    for(L i=0;i<m_3;i++){
      old_lambda[i]=lambda[i];
    }
    //start = std::clock();
    start=omp_get_wtime();
    update_x();
    update_y();
    update_lambda();
    nb_outer_iters++;
    //running_time_ADMM+=( std::clock() - start ) / (double) CLOCKS_PER_SEC;
    running_time_ADMM+=omp_get_wtime()-start;
    compute_and_record_res();
    /*
    res_x= 0;
    for(L i=0;i<m_1;i++){
      res_x+= (old_x[i]- x[i])*(old_x[i]- x[i]);
    }
    res_y= 0;
    for(L i=0;i<m_3;i++){
      res_y= (old_y[i]- y[i])*(old_y[i]- y[i]);
    }
    res_l= 0;
    for(L i=0;i<m_3;i++){
      res_l= (old_lambda[i]- lambda[i])*(old_lambda[i]- lambda[i]);
    }
    cout<< "res_x= "<< res_x<< " res_y= "<< res_y<< " res_l= "<< res_l<< endl;
    system("pause");
    */
    //start = std::clock();
    start=omp_get_wtime();
    reset_everything();
    //running_time_ADMM+=( std::clock() - start ) / (double) CLOCKS_PER_SEC;
    running_time_ADMM+=omp_get_wtime()-start;
    if (running_time_ADMM> time){
      break;
    }
  }

}
};

#endif /* MIN_SMOOTH_CONVEX_H */
