#ifndef PRIMAL_DUAL_LOOPLESS_H
#define PRIMAL_DUAL_LOOPLESS_H

#include "problem_data.h"


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


//This class implements the loopless variance reduced type methods with arbitrary sampling

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
class Primal_Dual_LOOPLESS: public problem_data<L, D>
{


private:

  // involved variables







  std::vector<D> proba_vector;

  std::vector<D>  tilde_proba_vector;


  std::vector<D> group_C;

  std::vector<D> index_group_C;

  std::vector<D> maxp_group_C;

  std::vector<D> isolated_I;

  std::vector<D> sump_group_C;
  vector<D> theta_S;   //theta_S in the paper




  std::vector<D> Li;



  std::vector<D> gk;     // the vector g^k in the paper








  // auxiliary variables
  std::vector<D> deltagk;
  std::vector<D> index_to_do_prox;
  std::vector<D> whether_or_not_to_do_prox;

  L nb_of_iters_per_loop;



  D primal_value;

  D dual_value;

  D epsilon;

  L max_nb_loops;

  D max_p;

  D maxv;



  D lambda_1;
  D lambda_2;

  D running_time;

  L nb_loops;



  L print_every_N;

  D delta;


  vector<L> batch_i;
  vector<L> my_batch;
  vector<D> batch_deltaalphai;

  L batch_size;

  L nb_groups;    //number of groups in group sampling




protected:

  std::vector<D> primal_x;     // the x variable


  string uniform;

  D Lf;

  D L2;

  D L1;

  D sumLi;

  D p; // the probability of changing w to x as in the paper
  L tau;  //number of threads on each node/computer

  D scaler;


  std::vector<D> dual_alpha;  //  dual_alpha[i]=\phi_i'(A_i't)

  std::vector<D> lambda_f;    // lambda_f[i]=lambda_i
  std::vector<D> baralpha;   // baralpha=sum_{i=1}^{nsamples} lambda_f[i]* dual_alpha[i]*A_i

  L current_nb_iters;

  L nb_iters;

  std::vector<D> last_update_i;

  ofstream samp;

public:


  gsl_rng * rng;
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
  virtual inline D compute_current_xj_value(D , L, L , L){return D(NULL);}
  virtual inline D compute_current_xj_value_without_update(D , L, L , L){return D(NULL);}
  virtual inline void set_stepsize(){}
  virtual inline void update_x(D,L){}
  virtual inline void update_baralpha(){}

  Primal_Dual_LOOPLESS()
  : problem_data<L,D>()
  {

  }

  Primal_Dual_LOOPLESS(const char* matrix_file, const char* vector_file)
  : problem_data<L,D>(matrix_file, vector_file)
  {

  }

  Primal_Dual_LOOPLESS(const char* matrix_file)
  : problem_data<L,D>(matrix_file)
  {

  }


  void set_rng()
  {
    gsl_rng_env_setup();
    const gsl_rng_type * T;
    T = gsl_rng_default;
    rng = gsl_rng_alloc(T);
    gsl_rng_set(rng,time(NULL));
    //gsl_rng_set(rng, 27432042);

  }



  void set_print_every_N(L i){print_every_N=i;}

  void set_L1(L p){
    if(p==0){ //batch sampling mode  (stochastic process in the paper)
      L1=(1-1./tau)*Lf+L2;
    }else{ //group sampling mode
      L1=Lf+L2;
    }
    //  alpha_l=1/(4.*L1+2*L2);
  }

  void set_Li_Lf()
  {
    sumLi=0;
    maxv=0;
    D minv=std::numeric_limits<double>::max();
    Li.resize(this->nsamples,0);
    for(L i=0;i<this->nsamples;i++)
    {
      D vi=0;
      for (L k = this->ptr[i]; k < this->ptr[i + 1];k++)
      {
        vi+=lambda_f[i]*this->nsamples*this->A[k]*this->A[k];
      }
      Li[i]=vi;
      if(maxv<vi) maxv=vi;
      if(minv>vi) minv=vi;
    }
    for(L i=0;i<this->nsamples;i++){
      if(Li[i]<maxv*1e-15)
      {
        for (L k = this->ptr[i]; k < this->ptr[i + 1];k++)
        { this->A[k]=0;
          L j=this->row_idx[k];
          for (L kj = this->ptr_t[j]; kj < this->ptr_t[j + 1];	kj++)
          {
            L i2=this->col_idx[kj];
            if(i2==i) this->A_t[kj]=0;
          }
        }
        Li[i]=maxv*1e-15;
      }
      sumLi+=Li[i];
    }
    Lf=compute_lambda_max(10);
    cout<<"  max of v: "<<maxv<<"  min of v: "<<minv<<" sumof Li="<<sumLi<<" Lf= "<<Lf<<endl;
  }


  //compute the maximal eigenvalue of sum_{i=1}^m \lambda_f[i]*A_i*A_i' by power iteration
  D compute_lambda_max(L K){

    std::vector<D> bk(this->nfeatures);
    for (L j=0;j<this->nfeatures;j++)
    {
      bk[j]=1;
    }
    std::vector<D> yk(this->nsamples);
    D normk;
    D tmp;
    for(L kk=0;kk<K;kk++){
      for (L i=0;i<this->nsamples;i++){
        tmp=0;
        for (L k = this->ptr[i]; k < this->ptr[i + 1];	k++)
        {
          L j=this->row_idx[k];
          tmp+=this->A[k]*bk[j];
        }
        yk[i]=tmp;
      }
      normk=0;
      for (L j=0;j<this->nfeatures;j++){
        bk[j]=0;
        for (L k = this->ptr_t[j]; k < this->ptr_t[j + 1];	k++)
        {
          L i=this->col_idx[k];
          bk[j]+=this->A_t[k]*yk[i]*lambda_f[i];
        }
        normk+=bk[j]*bk[j];
      }
      normk=sqrt(normk);
      for (L j=0;j<this->nfeatures;j++)
      {bk[j]=bk[j]/normk; }
    }
    cout<<endl;
    D res=0;
    normk=0;
    for (L i=0;i<this->nsamples;i++){
      tmp=0;
      for (L k = this->ptr[i]; k < this->ptr[i + 1];	k++)
      {
        L j=this->row_idx[k];
        tmp+=this->A[k]*bk[j];
      }
      yk[i]=tmp;
      normk+=yk[i]*yk[i];
    }
    std::vector<D> bk2(this->nfeatures);
    for (L j=0;j<this->nfeatures;j++){
      bk2[j]=0;
      for (L k = this->ptr_t[j]; k < this->ptr_t[j + 1];	k++)
      {
        L i=this->col_idx[k];
        bk2[j]+=this->A_t[k]*yk[i]*lambda_f[i];
      }
    }
    for (L j=0;j<this->nfeatures;j++)
    res+=bk2[j]*bk[j];
    return res;
  }

  void set_optimal_probability()
  {
    tilde_proba_vector.clear();
    tilde_proba_vector.resize(this->nsamples,0);
    proba_vector.clear();
    proba_vector.resize(this->nsamples,0);
    D sum=0;
    for(L i=0; i<this->nsamples;i++)
    {
      tilde_proba_vector[i]=Li[i];
      sum+=Li[i];
    }
    max_p=0;
    for(L i=0; i<this->nsamples;i++)
    {
      tilde_proba_vector[i]=tilde_proba_vector[i]/sum;
      proba_vector[i]=1-pow(1-tilde_proba_vector[i],tau);
      if(max_p<tilde_proba_vector[i])
      {
        max_p=tilde_proba_vector[i];
      }
    }
    cout<<"sum="<<sum<<"; proba"<<tilde_proba_vector[0]<<endl;
  }


  void set_group_sampling_probability(){
    proba_vector.clear();
    proba_vector.resize(this->nsamples,0);
    std::vector<D> q(this->nsamples);
    D sumq=0;
    cout<<"start setting group sampling probablity"<<endl;
    for(L i=0; i<this->nsamples;i++)
    {
      q[i]=Li[i];
      sumq+=q[i];
    }
    D maxq=0;
    L nb=0;
    D tmp=0;
    for(L i=0; i<this->nsamples;i++)
    {
      q[i]=q[i]/sumq*tau;
      if(q[i]>maxq) maxq=q[i];
      if(q[i]>1) {
        nb++; //count the number of elements larger than 1
        tmp+=q[i]-1;
      }
    }
    cout<<"maxq="<<maxq<<endl;
    if(maxq<=1){
      for(L i=0; i<this->nsamples;i++)
      {
        proba_vector[i]=q[i];
      }
    }else{
      for(L i=0; i<this->nsamples;i++)
      {
        if(q[i]>1) q[i]=1;
        else {
          D deltaq=min(1-q[i],tmp);
          tmp=tmp-deltaq;
          q[i]+=deltaq;
        }
        proba_vector[i]=q[i];
      }
    }
  }



  void set_L2(L p_mod, L u){
    if(p_mod==0){ //batch sampling mode
      D tmp=0;
      D st;
      for(L i=0;i<this->nsamples;i++){
        st=Li[i]/tilde_proba_vector[i];
        if (st>tmp) tmp=st;
      }
      L2=tmp/this->nsamples/tau;
      cout<<"L2="<<L2<<"; "<<tilde_proba_vector[0]<<endl;
    }
    else{ //group sampling mode
      D tmp=0;
      D st;
      for(L i=0;i<this->nsamples;i++){
        st=Li[i]/proba_vector[i];
        if(isolated_I[i]==1) st=st-Li[i];
        if(st>tmp) tmp=st;
      }
      L2=tmp/this->nsamples;
    }
    cout<<"set_L2="<<L2<<endl;
  }

  void set_uniform_probability(L p_mod)
  {
    if(p_mod==0) {
      //batch sampling (stochastic process as defined in the paper)
      proba_vector.clear();
      tilde_proba_vector.clear();
      tilde_proba_vector.resize(this->nsamples,1./this->nsamples);
      proba_vector.resize(this->nsamples,1-pow(1-1.0/this->nsamples,tau));
      max_p=1.0/this->nsamples;
      cout<<"uniform="<<tilde_proba_vector[0]<<endl;
    }
    else{ //group sampling
      D pi=(tau*this->c+0.0)/this->nsamples;
      proba_vector.clear();
      proba_vector.resize(this->nsamples,pi);
    }
  }




  void sort_p(){
    vector<pair<D,L> >a;
    for (L i = 0 ;i < this->nsamples ; i++) {
      a.push_back(make_pair(proba_vector[i],i)); // k = value, i = original index
    }
    sort(a.begin(),a.end());
    group_C.clear();
    group_C.resize(this->nsamples);
    isolated_I.clear();
    isolated_I.resize(this->nsamples,0);
    D tmp=0;
    D maxpi=0;
    D previous_i=0;
    index_group_C.clear();
    index_group_C.push_back(0);
    maxp_group_C.clear();
    sump_group_C.clear();

    for (L i = 0 ;i < this->nsamples ; i++){
      group_C[i]=a[this->nsamples-1-i].second;
      D pi=a[this->nsamples-1-i].first;
      if(tmp+pi<=1){
        if(pi>maxpi) maxpi=pi;
        tmp+=pi;
      }
      else{
        index_group_C.push_back(i);
        maxp_group_C.push_back(maxpi);
        sump_group_C.push_back(tmp);
        tmp=pi;
        maxpi=pi;
        if(i-previous_i==1) isolated_I[group_C[previous_i]]=1;
        previous_i=i;
      }
    }
    index_group_C.push_back(this->nsamples);
    maxp_group_C.push_back(maxpi);
    sump_group_C.push_back(tmp);
    nb_groups=sump_group_C.size();
    if(previous_i==this->nsamples-1) isolated_I[group_C[previous_i]]=1;
    cout<<"size of group="<<nb_groups<<endl;
  }


  inline L sampling(L n)
  {
    L i=(floor)(gsl_rng_uniform(rng)*n);
    D y=gsl_rng_uniform(rng);
    while(y*max_p>tilde_proba_vector[i])
    {
      i=(floor)(gsl_rng_uniform(rng)*n);
      y=gsl_rng_uniform(rng);
    }
    return i;
  }






  void set_initial_dual(){
    for(L i=0;i<this->nsamples;i++)
    {
      D res=0;
      for (L k = this->ptr[i]; k < this->ptr[i + 1];k++)
      {
        L j=this->row_idx[k];
        res+=this->A[k]*primal_x[j];
      }
      dual_alpha[i]=gradient_of_phi_i(res,i);
    }
  }




  void update_primal()
  {
    L nb_indices=0;
    //#pragma omp parallel for
    for(L i_t=0;i_t<batch_size;i_t++)
    {
      L i=batch_i[i_t];
      D deltaalphai=batch_deltaalphai[i_t];
      for (L k = this->ptr[i]; k < this->ptr[i + 1];k++)
      {
        L j=this->row_idx[k];
        if (whether_or_not_to_do_prox[j]==0) {
          //#pragma omp critical
          {
          whether_or_not_to_do_prox[j]=1;
          index_to_do_prox[nb_indices]=j;
          nb_indices++;
          deltagk[j]=0;
          }
        }
        deltagk[j]+=lambda_f[i]*deltaalphai*this->A[k]*theta_S[i];
      }
    }
    #pragma omp parallel for
    for(L i_d=0;i_d<nb_indices;i_d++){
      L j=index_to_do_prox[i_d];
      D xj=compute_current_xj_value(baralpha[j]+deltagk[j], j,last_update_i[j],current_nb_iters);
      //update_x(xj,j);
      last_update_i[j]=current_nb_iters;
      whether_or_not_to_do_prox[j]=0;
      index_to_do_prox[i_d]=-1;
    }

  }

  void update_baralpha(L  i, D deltaalphai)
  {
    for (L k = this->ptr[i]; k < this->ptr[i + 1];	k++)
    {
      L j=this->row_idx[k];
      baralpha[j]+=lambda_f[i]*deltaalphai*this->A[k];
    }
  }



  D compute_AiTxk(L i) //This function compute A_i^\top x^k
  {
    D res=0;
    D xj;
    for (L k = this->ptr[i]; k < this->ptr[i + 1];k++)
    {
      L j=this->row_idx[k];
      if(current_nb_iters!=last_update_i[j])
      {
        xj=compute_current_xj_value(baralpha[j], j,last_update_i[j],current_nb_iters);
        //update_x(xj,j);
      }
      res+=this->A[k]*primal_x[j];
      last_update_i[j]=current_nb_iters;
    }
    return res;
  }


  D compute_AiTxk_without_update(L i) //This function compute A_i^\top x^k  without updating x^k
  {
    D res=0;
    D xj;
    //#pragma omp parallel for reduction(+:res)
    for (L k = this->ptr[i]; k < this->ptr[i + 1];k++)
    {
      L j=this->row_idx[k];
      if(current_nb_iters!=last_update_i[j])
      {
        xj=compute_current_xj_value_without_update(baralpha[j], j,last_update_i[j],current_nb_iters);
      }
      else {
        xj=primal_x[j];
      }
      res+=this->A[k]*xj;
    }
    return res;
  }




  void compute_dual_value()
  {
    D res=0;
    //D alphaAx= 0;
    std::vector<D> alpha_tmp(this->nsamples,0);
    std::vector<D> baralpha_tmp(this->nfeatures,0);
    //#pragma omp parallel for schedule(dynamic)
    for(L i=0;i<this->nsamples;i++)
    {
      D aitx=compute_AiTxk_without_update(i);
      alpha_tmp[i]=gradient_of_phi_i(aitx,i);
      //alphaAx+= lambda_f[i]*alpha_tmp[i]*aitx;
      for (L k = this->ptr[i]; k < this->ptr[i + 1];	k++)
      {
        L j=this->row_idx[k];
        baralpha_tmp[j]+=lambda_f[i]*alpha_tmp[i]*this->A[k];
      }
    }

    D scal=feasible_dual(baralpha_tmp);
    //#pragma omp parallel for reduction(-:res)
    for(L i=0;i<this->nsamples;i++)
    {
      res-=value_of_phistar_i(scal*alpha_tmp[i],i)*lambda_f[i];
    }
    res-=value_of_gstar_minus(baralpha_tmp,scal); // should have used value_of_gstar(-baralpha)
    dual_value=res;
  }

  void compute_primal_value()
  {
    D res=0;
    //#pragma omp parallel for reduction(+:res)
    for(L i=0;i<this->nsamples;i++)
    {
      D aitx=compute_AiTxk_without_update(i);
      res+=lambda_f[i]*value_of_phi_i(aitx,i);
    }

    D res2=0;
    //#pragma omp parallel for reduction(+:res2)
    for(L j=0; j<this->nfeatures; j++)
    {
      D xj=compute_current_xj_value_without_update(baralpha[j], j,last_update_i[j],current_nb_iters);
      D gj=value_of_g_j(xj,j);
      res2+=gj;
    }
    primal_value= res+res2;
  }




  void compute_initial_dual_value()
  {

    D res=0;
    D scal=feasible_dual(baralpha);
    #pragma omp parallel for reduction(-:res)
    for(L i=0;i<this->nsamples;i++)
    {
      res-=value_of_phistar_i(scal*dual_alpha[i],i)*lambda_f[i];
    }
    //cout<< "d_1(x)= "<< res<< " d_2(x)= "<< value_of_gstar_minus(baralpha,scal)<< endl;
    res-=value_of_gstar_minus(baralpha,scal);  // should have used value_of_gstar(-baralpha)
    dual_value=res;
  }


  void compute_initial_primal_value()
  {
    D res=0;
    #pragma omp parallel for reduction(+:res)
    for(L j=0; j<this->nfeatures; j++)
    {
      res+=value_of_g_j(primal_x[j],j);
    }
    D res2=0;
    #pragma omp parallel for reduction(+:res2)
    for(L i=0;i<this->nsamples;i++)
    {
      D aitx=0;
      for (L k = this->ptr[i]; k < this->ptr[i + 1];	k++)
      {
        L j=this->row_idx[k];
        aitx+=this->A[k]*primal_x[j];
      }
      res2+=lambda_f[i]*value_of_phi_i(aitx,i);
    }
    //cout<< "p_1(x)= "<< res<< " p_2(x)= "<< res2<< endl;
    primal_value= res+res2;
  }



  void set_p(D scal_p){

    if(scal_p==4.6){
      p=sqrt(this->mu/L1*tau/this->nsamples);
    }
    else if (scal_p==5.6){
      p=this->mu/L1;
    }
    else{
      p=scal_p*tau/(0.0+this->nsamples);
    }
    cout<<"changing probablity: "<<p<<endl;
  }




  void initialize(vector<D> & x0, vector<D> & w0, vector<D> & L_phi, D val_mu, D val_epsilon, L max_nb, L nb_tau, L nb_c, L u, L p_mod, D scal_p)
  {
    cout<<"start initializing"<<" u="<<u<<endl;
    this->distributed_set_up(nb_c);
    if(nb_tau>this->noverc) perror("tau should be less than n over c");
    tau=nb_tau;
    nb_of_iters_per_loop=floor(this->nsamples/(this->c*(tau+0.0)));
    batch_i.clear();
    batch_deltaalphai.clear();
    batch_i.resize(this->nsamples,0);
    batch_deltaalphai.resize(this->nsamples,0);

    deltagk.clear();
    deltagk.resize(this->nfeatures,0);
    index_to_do_prox.clear();
    index_to_do_prox.resize(this->nfeatures,-1);
    whether_or_not_to_do_prox.clear();
    whether_or_not_to_do_prox.resize(this->nfeatures,0);
    gk.clear();
    gk.resize(this->nfeatures,0);

    /**setup parameters**/
    epsilon=val_epsilon;
    max_nb_loops=max_nb;
    this->mu=val_mu;
    cout<<"mu="<<this->mu<<endl;
    lambda_f.clear();
    lambda_f.resize(this->nsamples,0); //L_phi is the Lipschitz constant of the function \phi_i.
    for(int i=0;i< this->nsamples; i++){
      lambda_f[i]= L_phi[i]/this->nsamples;
    }

    set_rng();
    set_auxiliary_v();


    lambda_1=get_lambda1();
    lambda_2=get_lambda2();

    set_Li_Lf();

    running_time=0;
    /**setup probability**/
    if(p_mod==0)   // batch sampling (stochastic process as defined in the paper)
    {
      if(u==0)   // uniform sampling
      {
        set_uniform_probability(p_mod);
        this->uniform="uniform";
      }
      else{
        set_optimal_probability();
        this->uniform="nonuniform";
      }
    }
    else{  //group sampling
      if(u==0){
        set_uniform_probability(p_mod);
        this->uniform="uniform";
      }else{
        this->set_group_sampling_probability();
        this->uniform="nonuniform";
      }
      sort_p();
      cout<<"sort p"<<endl;
    }




    set_L2(p_mod,u);
    scaler=scal_p;

    set_L1(p_mod);
    cout<<"set L1="<<L1<<endl;
    set_p(scal_p);
    set_thetaS(p_mod);
    cout<<"setthetaS"<<endl;

    primal_x=w0;
    dual_alpha=x0;
    set_stepsize();
    cout<<"set stepsize"<<endl;


    set_initial_dual();

    baralpha.clear();
    baralpha.resize(this->nfeatures,0);
    last_update_i.clear();
    last_update_i.resize(this->nfeatures,0);
    for(L i=0;i<this->nsamples;i++)
    {
      update_baralpha(i, dual_alpha[i]);
    }

    current_nb_iters=0;
    D start_time_initial=std::clock();
    compute_initial_primal_value();
    compute_initial_dual_value();
    D end_time_initial=std::clock();
    cout<<"computation of initial primal dual value time="<<(end_time_initial-start_time_initial) / (double) CLOCKS_PER_SEC<<endl;

    cout<<"primal value="<<primal_value<<endl;
    cout<<"dual_value="<<dual_value<<endl;


    cout<<"Initialization is finished!"<<endl;
  }


  void result_record()
  {
    cout<<primal_value<<endl;
    cout<<dual_value<<endl;
    ofstream result;
    result.open("results/x0.dat");
    for (L j=0;j<this->nfeatures;j++){
      result <<setprecision(20)<< primal_x[j] << " ";
      //cout<< primal_x[j] << " ";
    }
    result << endl;
    ofstream resulta;
    resulta.open("results/alpha0.dat");
    for (L i=0;i<this->nsamples;i++)
    resulta << setprecision(20)<<dual_alpha[i] << " ";
    resulta << endl;
    result.close();
  }


  void read_w0(){
    ifstream recordw0("results/x0.dat");
    for(L j=0;j<this->nfeatures;j++)
    {
      recordw0>>primal_x[j];
    }
    recordw0.close();

    ifstream recordalpha0("results/alpha0.dat");
    for(L i=0;i<this->nsamples;i++)
    {
      recordalpha0>>dual_alpha[i];
    }
    recordalpha0.close();
  }







  void compute_and_record_result()
  {
    if(nb_loops%print_every_N==0){
      D start = std::clock();
      compute_primal_value();
      D endt=std::clock();
      cout<<"computation of primal value time="<<(endt-start) / (double) CLOCKS_PER_SEC<<endl;
      compute_dual_value();
      delta=primal_value-dual_value;
      cout<<setprecision(16)<<floor(((0.0+nb_iters)*this->c*tau/(this->nsamples)))<<";;;  "<<running_time<<"; primal dual gap="<<delta<<" primal value="<<primal_value<<"; dual value="<<dual_value<<endl;
      samp<<floor(((0.0+nb_iters)*this->c*tau/(this->nsamples)))<<" "<<delta<<" "<<running_time<<" "<<primal_value<<endl;
    }

  }




  void batch_sampling()
  {
    batch_size=0;
    L i;
    for(L it_t=0;it_t<tau;it_t++)
    {
      i=sampling(this->nsamples);
      batch_i[batch_size]=i;
      batch_size++;
    }
  }


  void set_thetaS(L p_mod){
    theta_S.clear();
    theta_S.resize(this->nsamples,0);
    if(p_mod==0){
      for(L i=0;i<this->nsamples;i++)
      theta_S[i]=1.0/(tau*tilde_proba_vector[i]);
    }
    else{
      for(L i=0;i<this->nsamples;i++)
      theta_S[i]=1.0/proba_vector[i];
    }
  }

  void group_sampling(){
    batch_size=0;
    std::vector<L> sampled_groups;
    for(L i=0;i<nb_groups;i++){
      D y=gsl_rng_uniform(rng);
      if(y<sump_group_C[i])
      {
        sampled_groups.push_back(i);
        batch_size++;
      }
    }
    for(L t=0;t<batch_size;t++){
      L group_i=sampled_groups[t];
      L s1=index_group_C[group_i];
      L s2=index_group_C[group_i+1];
      L i=s1+(floor)(gsl_rng_uniform(rng)*(s2-s1));
      D y=gsl_rng_uniform(rng);
      D maxpi=maxp_group_C[group_i];
      while(y*maxpi>proba_vector[group_C[i]])
      {
        i=s1+(floor)(gsl_rng_uniform(rng)*(s2-s1));
        y=gsl_rng_uniform(rng);
      }
      batch_i[t]=group_C[i];
    }
  }




  void loopless(vector<D> & x0, vector<D> & w0, string filename, vector<D> & L_phi, D val_mu, D val_epsilon, L max_nb, L nb_tau, L nb_c, L u, L p_mod, D scal_p)
  {
    initialize(x0,  w0,  L_phi, val_mu, val_epsilon,  max_nb,  nb_tau, nb_c,  u, p_mod, scal_p);
    string sampname="results/L_"+filename+uniform;
    if(p_mod==0) sampname=sampname+"_batch";
    else sampname=sampname+"_group";
    string scal_str;
    stringstream  scal_convert;
    scal_convert<<scal_p;
    scal_str=scal_convert.str();
    sampname+=scal_str;
    cout<<"running Loopless"<<" ; "<<sampname<<endl;
    samp.open(sampname.c_str());


    delta=primal_value-dual_value;
    nb_loops=0;
    nb_iters=0;

    cout<<setprecision(9)<<"initial: "<<" delta: "<<delta<<" primal: "<<primal_value<<"; dual: "<<dual_value<<"epsilon="<<epsilon<<endl;
    samp<<((0.0+nb_iters)*this->c*tau/(this->nsamples))<<" "<<delta<<running_time<<" "<<primal_value<<endl;
    srand48(27432042);
    //srand(time(NULL));
    printf("There are %d threads\n",omp_get_num_threads());

    D start;
    while(delta>=epsilon && nb_loops<max_nb_loops)
    {
      compute_and_record_result();

      if(delta<epsilon) break;
      start = std::clock();

      nb_loops++;

      //cout<<"before the loop time elapsed="<<duration<<endl;
      start = std::clock();

      //cout<<"nb_of_iters_per_loop="<<nb_of_iters_per_loop<<endl;
      for(L it=0;it<nb_of_iters_per_loop;it++)
      {
        if(p_mod==0)
        batch_sampling();
        else
        group_sampling();

        //start = std::clock();
        start=omp_get_wtime();
        #pragma omp parallel
        {
          #pragma omp for
          for(L i_t=0;i_t<batch_size;i_t++)
          {
          L i=batch_i[i_t];
          D aitg=compute_AiTxk(i);
          D deltaalphai=(gradient_of_phi_i(aitg,i)-dual_alpha[i]);
          batch_deltaalphai[i_t]=deltaalphai;
         }
      }
        //cout<<omp_get_wtime()-start<<endl;
        current_nb_iters++;
        update_primal();
        update_baralpha();
        nb_iters++;
        D end_time = omp_get_wtime();
        //running_time+= ( std::clock() - start ) / (double) CLOCKS_PER_SEC;
        running_time+=end_time-start;

      }

      //cout<<"after the loop time elapsed="<<duration<<endl;
    }
    if(delta>=epsilon) {compute_primal_value();
      compute_dual_value();
      delta=primal_value-dual_value;
      cout<<setprecision(9)<<floor(((0.0+nb_iters)*this->c*tau/(this->nsamples)))<<";  "<<running_time<<"; primal dual gap="<<delta<<" primal value="<<primal_value<<"; dual value="<<dual_value<<endl;
      samp<<floor(((0.0+nb_iters)*this->c*tau/(this->nsamples)))<<" "<<delta<<" "<<running_time<<" "<<primal_value<<endl;
    }
    cout<<setprecision(9)<<"nb_loops: "<<nb_loops<<"nb_iters"<<floor(((0.0+nb_iters)*this->c*tau/(this->nsamples)))<<";  delta: "<<delta<<" primal: "<<primal_value<<"; dual: "<<dual_value<<endl;
    samp.close();
  }

  void loopless2(vector<D> & x0, vector<D> & w0, string filename, vector<D> & L_phi, D val_mu, D val_epsilon, L max_nb, L nb_tau, L nb_c, L u, L p_mod, D scal_p)
  {
    initialize(x0,  w0,  L_phi, val_mu, val_epsilon,  max_nb,  nb_tau, nb_c,  u, p_mod, scal_p);
    string sampname="results/L_"+filename+uniform;
    if(p_mod==0) sampname=sampname+"_batch";
    else sampname=sampname+"_group";
    string scal_str;
    stringstream  scal_convert;
    scal_convert<<scal_p;
    scal_str=scal_convert.str();
    sampname+=scal_str;
    cout<<"running Loopless"<<" ; "<<sampname<<endl;
    samp.open(sampname.c_str());


    delta=primal_value-dual_value;
    nb_loops=0;
    nb_iters=0;

    cout<<setprecision(9)<<"initial: "<<" delta: "<<delta<<" primal: "<<primal_value<<"; dual: "<<dual_value<<endl;
    samp<<((0.0+nb_iters)*this->c*tau/(this->nsamples))<<" "<<delta<<running_time<<" "<<primal_value<<endl;
    //srand48(27432042);
    srand(time(NULL));


    D start;
    while(nb_loops<max_nb_loops)
    {
      compute_and_record_result();

      start = std::clock();

      nb_loops++;

      //cout<<"before the loop time elapsed="<<duration<<endl;
      start = std::clock();
      //cout<<"nb_of_iters_per_loop="<<nb_of_iters_per_loop<<endl;
      for(L it=0;it<nb_of_iters_per_loop;it++)
      {
        if(p_mod==0)
        batch_sampling();
        else
        group_sampling();

        start = std::clock();
        for(L i_t=0;i_t<batch_size;i_t++)
        {
          L i=batch_i[i_t];
          D aitg=compute_AiTxk(i);
          D deltaalphai=(gradient_of_phi_i(aitg,i)-dual_alpha[i]);
          batch_deltaalphai[i_t]=deltaalphai;
        }
        current_nb_iters++;
        update_primal();
        update_baralpha();
        nb_iters++;
        running_time+= ( std::clock() - start ) / (double) CLOCKS_PER_SEC;
      }
      //cout<<"after the loop time elapsed="<<duration<<endl;
    }
    cout<<setprecision(9)<<"nb_loops: "<<floor(((0.0+nb_iters)*this->c*tau/(this->nsamples)))<<";  delta: "<<delta<<" primal: "<<primal_value<<"; dual: "<<dual_value<<endl;
    samp.close();
  }









};

#endif /* MIN_SMOOTH_CONVEX_H */
