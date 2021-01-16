#ifndef APPROX2_H
#define APPROX2_H




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
#include <omp.h>


//This class implements the method APPROX presented in the paper Accelerated, parallel and proximal coordinate descent, by Fercoq & Richtarik.

/*
The optimization problem to solve is:

min f(x)+ P(x)

*/

template<typename L, typename D>
class APPROX2
{


protected:

  // parameters

  D mu_f;

  D mu_psi;

  L n;   // x\in \R^n

  L tau;  //number of threads on each node/computer

  D sumofLi;

  // variables
  std::vector<D> u;

  std::vector<D> z;

  std::vector<D> x;

  std::vector<D> v;

  std::vector<D> t;

  D gamma;

  D theta;

  // sampling variables

  std::vector<D> proba_vector;

  std::vector<D> S;

  std::vector<D> all_n;

  std::vector<D> sampled;

  D max_p;

  // auxiliary variables

  D tauovern;
  D novertau;
  D novertau2;
  D novertau3;
  D novertau4;
  D taumuovern;

  // variables for print
  D primal_value;

  D dual_value;

  D gradient_norm;

  D subgradient;

  D epsilon;

  L max_nb_loops;

  L evaluation;

  D delta;   //primal dual gap

  D delta2; //distance to subgradient

  ofstream samp;

  L nb_iters;

  L nb_of_iters_per_loop;

  D running_time;

  L print_every_N;

  L mod;





public:


  gsl_rng * rng;
  virtual inline D partial_i_of_f(L i){return D(NULL);}
  virtual inline D compute_prox(D, D, D, L){return D(NULL);}
  virtual inline void compute_primal_value() {}
  virtual inline void compute_dual_value(){}
  virtual inline void compute_gradient_norm(){}
  virtual inline void set_v(){}
  virtual inline void set_p(){}
  virtual inline void update_z_coordinate( L, D){}
  virtual inline void update_u_coordinate( L, D){}
  virtual inline void update_x_coordinate( L, D){}
  APPROX2()
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



  // sample i with probability pi=proba_vector[i]
  L sampling()
  {
    //L i=(floor)(gsl_rng_uniform(rng)*n);
    L i=gsl_rng_uniform_int(rng, n);
    if(tau==1)
    {
      D y=gsl_rng_uniform(rng);
      while(y*max_p>proba_vector[i])
      {
        i=(floor)(gsl_rng_uniform(rng)*n);
        y=gsl_rng_uniform(rng);
      }
    }
    return i;
  }


  // sample S
  void batch_sampling()
  {
    if(tau<n)
    {
      L i=sampling();
      for(L k=0;k<tau;k++)
      {
        while(sampled[i]==1)
        {
          i=sampling();
        }
        sampled[i]=1;
        S[k]=i;
      }
      for(L k=0;k<tau;k++)
      {
        sampled[S[k]]=0;
      }
    }
    else {
      S=all_n;
      //cout<<"s=all_n"<<endl;
    }
  }


  void compute_x()
  {
    for(L i=0;i<n;i++)
    x[i]=gamma*u[i]+z[i];
  }

  void update_theta()
  {
    if(mod==1){
      D theta2=theta*theta;
      D tmp=mu_f+mu_psi-theta2*novertau2-mu_psi*(theta+1)*novertau;
      D tmp2=theta2*novertau4+theta*mu_psi*novertau3;
      D tmp3=mu_f+mu_psi-novertau2*theta2-novertau*mu_psi*(theta+1);
      theta=0.5/novertau2*(sqrt(tmp*tmp+4*tmp2)+tmp3);
    }

  }

  void initialize(vector<D> & x0, L val_tau, D val_mu_f, D val_mu_psi, L eval, L p_N, L val_mod)
  {
    tau=val_tau;
    theta=tau/(n+0.0);
    tauovern=tau/(n+0.0);
    novertau=n/(tau+0.0);
    novertau2=novertau*novertau;
    novertau3=novertau2*novertau;
    novertau4=novertau2*novertau2;
    all_n.resize(n,0);
    for(L i=0;i<n;i++)
    all_n[i]=i;
    gamma=1;
    mu_f=val_mu_f;
    mu_psi=val_mu_psi;
    evaluation=eval;
    nb_of_iters_per_loop=floor(max(1.,n/(tau+0.0)));
    print_every_N=p_N;
    mod=val_mod;

    taumuovern=tau*mu_f/(n+0.0);
    delta=std::numeric_limits<double>::max();
    gradient_norm=std::numeric_limits<double>::max();


    u.clear();
    u.resize(n,0);
    z.clear();
    z.resize(n,0);
    x.clear();
    x.resize(n,0);
    for(L i=0;i<n;i++)
    {
      z[i]=x0[i];
      x[i]=x0[i];
    }
    sampled.clear();
    sampled.resize(n,0);
    S.clear();
    S.resize(tau,0);
    t.clear();
    t.resize(n,0);

    set_v();
    set_p();
    set_rng();

    if(evaluation==1)
    {
      compute_primal_value();
      compute_dual_value();
      cout<<"initial primal value="<<primal_value<<endl;
      cout<<"initial dual value="<<dual_value<<endl;
    }
    else if(evaluation==2)
    {
      compute_primal_value();
      compute_gradient_norm();
      cout<<"initial primal value="<<primal_value<<endl;
      cout<<"initial gradient norm="<<gradient_norm<<endl;
    }
    running_time=0;
    cout<<"finished APPROX initializing"<<endl;
  }





  void compute_and_record_result()
  {
    if(evaluation==1&&nb_iters%print_every_N==0)
    {
      compute_primal_value();
      compute_dual_value();
      compute_gradient_norm();
      delta=primal_value-dual_value;
      cout<<setprecision(9)<<(0.0+nb_iters)<<";  "<<running_time<<"; primal dual gap="<<delta<<" primal value="<<primal_value<<"; dual value="<<dual_value<<endl;
      samp<<(0.0+nb_iters)<<" "<<delta<<" "<<running_time<<" "<<primal_value<<" "<<dual_value<<" "<<gradient_norm<<endl;
      gradient_norm=std::numeric_limits<double>::max();
    }
    else if(evaluation==2&&nb_iters%print_every_N==0)
    {
      compute_primal_value();
      compute_gradient_norm();
      cout<<setprecision(9)<<(0.0+nb_iters)<<";  "<<running_time<<"; gradient norm="<<gradient_norm<<" primal value="<<primal_value<<endl;
      samp<<(0.0+nb_iters)<<" "<<gradient_norm<<" "<<running_time<<" "<<primal_value<<endl;
    }
    else if(evaluation==3&&nb_iters%print_every_N==0) //running APPROX for solving the subproblem in ALM
    {
      compute_primal_value();
      compute_dual_value();
      compute_gradient_norm();
      delta=primal_value-dual_value;
      cout<<setprecision(9)<<"      "<<(0.0+nb_iters)<<";  "<<running_time<<"; primal dual gap="<<delta<<" primal value="<<primal_value<<"; dual value="<<dual_value<<"; gradient_norm="<<gradient_norm<<endl;
      samp<<(0.0+nb_iters)<<" "<<delta<<" "<<running_time<<" "<<primal_value<<" "<<dual_value<<" "<<gradient_norm<<endl;
      gradient_norm=std::numeric_limits<double>::max();
    }

  }

  void compute_and_record_result_always(){
    if(evaluation==1)
    {
      compute_primal_value();
      compute_dual_value();
      compute_gradient_norm();
      delta=primal_value-dual_value;
      cout<<setprecision(9)<<(0.0+nb_iters)<<";  "<<running_time<<"; primal dual gap="<<delta<<" primal value="<<primal_value<<"; dual value="<<dual_value<<endl;
      samp<<(0.0+nb_iters)<<" "<<delta<<" "<<running_time<<" "<<primal_value<<" "<<dual_value<<" "<<gradient_norm<<endl;
      gradient_norm=std::numeric_limits<double>::max();
    }
    else if(evaluation==2)
    {
      compute_primal_value();
      compute_gradient_norm();
      cout<<setprecision(9)<<(0.0+nb_iters)<<";  "<<running_time<<"; gradient norm="<<gradient_norm<<" primal value="<<primal_value<<endl;
      samp<<(0.0+nb_iters)<<" "<<gradient_norm<<" "<<running_time<<" "<<primal_value<<endl;
    }
    else if(evaluation==3) //running APPROX for solving the subproblem in ALM
    {
      compute_primal_value();
      compute_dual_value();
      compute_gradient_norm();
      delta=primal_value-dual_value;
      cout<<setprecision(9)<<"      "<<(0.0+nb_iters)<<";  "<<running_time<<"; primal dual gap="<<delta<<" primal value="<<primal_value<<"; dual value="<<dual_value<<"; gradient_norm="<<gradient_norm<<endl;
      samp<<(0.0+nb_iters)<<" "<<delta<<" "<<running_time<<" "<<primal_value<<" "<<dual_value<<" "<<gradient_norm<<endl;
      gradient_norm=std::numeric_limits<double>::max();
    }


  }






  void APPROX_MU(vector<D> & x0, L val_tau, D val_mu_f, D val_mu_psi, L eval, L p_N,  L max_nb_epoch, D eps, string filename, L val_mod)
  {
    initialize(x0,  val_tau,  val_mu_f,  val_mu_psi,  eval, p_N, val_mod);
    cout<<"running APPROX MU"<<" ; "<<filename<<" max_nb_epoch "<<max_nb_epoch<<"; eps="<<eps<<endl;

    nb_iters=0;
    compute_and_record_result();
    //srand48(27432042);
    srand(time(NULL));
    D start;
    while(delta>eps && nb_iters<max_nb_epoch)
    {
      //start = std::clock();
      start=omp_get_wtime();

      for(L it=0;it<nb_of_iters_per_loop;it++)
      {
        gamma*=(1-theta)/(1-taumuovern);
        if(theta==1) gamma=1;
        batch_sampling();

        //#pragma omp parallel
        {
          //#pragma omp for
          for(L it_S=0;it_S<tau;it_S++)
          {
            L i=S[it_S];
            D si=z[i]+taumuovern*gamma/theta*u[i];
            D gi=partial_i_of_f(i);
            D Li=v[i]*novertau*theta;
            t[i]=compute_prox(gi, Li, si,i);
          }
        }
        //#pragma omp parallel for
        for(L it_S=0;it_S<tau;it_S++)
        {
          L i=S[it_S];
          D ti=t[i];
          update_z_coordinate(i, ti+taumuovern*gamma/theta*u[i]);
          update_u_coordinate(i, (-1+novertau*theta)/gamma*ti-taumuovern/theta*u[i]);
        }
        update_theta();
      }
      //running_time+=( std::clock() - start ) / (double) CLOCKS_PER_SEC;
      running_time+=omp_get_wtime()-start;
      nb_iters++;
      compute_and_record_result();
    }

    compute_and_record_result_always();
  }

  void APPROX_MU2(vector<D> & x0, L val_tau, D val_mu_f, D val_mu_psi, L eval, L p_N,  L max_nb_epoch, D eps, string filename, L val_mod)
  {
    initialize(x0,  val_tau,  val_mu_f,  val_mu_psi,  eval, p_N, val_mod);
    cout<<"running APPROX MU"<<" ; "<<filename<<" max_nb_epoch "<<max_nb_epoch<<"; eps="<<eps<<endl;

    nb_iters=0;
    compute_and_record_result();
    //srand48(27432042);
    srand(time(NULL));
    D start;
    while(nb_iters<max_nb_epoch)
    {
      start = std::clock();
      D ti=0;
      D si=0;
      D gi=0;
      D Li=0;
      L i=0;
      for(L it=0;it<nb_of_iters_per_loop;it++)
      {
        gamma*=(1-theta)/(1-taumuovern);
        if(theta==1) gamma=1;
        batch_sampling();
        for(L it_S=0;it_S<tau;it_S++)
        {
          i=S[it_S];
          si=z[i]+taumuovern*gamma/theta*u[i];
          gi=partial_i_of_f(i);
          Li=v[i]*novertau*theta;
          t[i]=compute_prox(gi, Li, si,i);
        }
        for(L it_S=0;it_S<tau;it_S++)
        {
          i=S[it_S];
          ti=t[i];
          update_z_coordinate(i, ti+taumuovern*gamma/theta*u[i]);
          update_u_coordinate(i, (-1+novertau*theta)/gamma*ti-taumuovern/theta*u[i]);
        }
        update_theta();
      }
      running_time+=( std::clock() - start ) / (double) CLOCKS_PER_SEC;
      nb_iters++;
      compute_and_record_result();
    }

    compute_and_record_result_always();
  }


  void prox_grad_step(){
    D gi;
    for(L i=0;i<n;i++){
      gi=partial_i_of_f(i);
      t[i]=compute_prox(gi, sumofLi, x[i],i);
    }
    for(L i=0;i<n;i++){
      update_x_coordinate(i, t[i]);
    }


  }

  void do_single_step_prox(vector<D> & Tx){
    for(L i=0;i<n;i++){
      D gi=partial_i_of_f(i);
      Tx[i]=compute_prox(gi, sumofLi, x[i],i)+x[i];
    }
  }



};

#endif
