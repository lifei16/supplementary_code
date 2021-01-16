#ifndef PDCD_H
#define PDCD_H


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



//This class implements the method SMART_CD.

/*
This class solves problem of the form f(x)+g(x)+h(Mx) ;
where f(x)=sum_{j=1}^m lambda_f[j] phi_j(<A_j,x>) and g(x)=sum_{i=1}^n g_i(x_i).
We all assume that each \phi_j is 1-smooth.
*/



template<typename L, typename D>
class PDCD
{
private:

  std::vector<D> Au;
  std::vector<D> Az;

  std::vector<D> Mu;
  std::vector<D> Mz;


  std::vector<D> Mx;
  std::vector<D> Ax;
  std::vector<D> lambda_f;



protected:

  // parameters

  D mu_f;

  D mu_psi;

  D mu_g;

  L n;   // x\in \R^n

  L tau;  //number of threads on each node/computer

  D sumofLi;

  // variables
  std::vector<D> u;

  std::vector<D> z;

  std::vector<D> x;

  std::vector<D> v;

  std::vector<D> L_M;

  std::vector<D> L_f;

  std::vector<D> t;

  D gamma;

  D theta;

  D theta0;

  L m_1;

  L m_2;

  // sampling variables

  std::vector<D> proba_vector;

  std::vector<D> S;

  std::vector<D> all_n;

  std::vector<D> sampled;

  D max_p;

  D min_p;

  // auxiliary variables

  L max_nb_loops;

  L evaluation;

  ofstream samp;

  L nb_iters;

  L nb_of_iters_per_loop;

  L print_every_N;

  D running_time;

  Matrix<L,D> data_A;

  Matrix<L, D> data_M;

  std::vector<D> lambda;

  std::vector<D> M_tlambda;

  D beta_s;

  D function_value;




public:


  gsl_rng * rng;

  virtual inline D gradient_of_phi_j(D, L){return D(NULL);}

  virtual inline D value_of_g_i(D, L){return D(NULL);}
  virtual inline D value_of_phi_j(D, L){return D(NULL);}
  virtual inline D value_of_h_j(D, L){return D(NULL);}
  virtual inline D value_of_h_star_j(D, L){return D(NULL);}

  virtual inline D prox_of_g_i(D,D,D, L){return D(NULL);}
  virtual inline D prox_of_h_star_j(D,D, L){return D(NULL);}

  virtual inline D value_of_phistar_i(D,L) {return D(NULL);}

  virtual inline void set_matrix_M(){}

  virtual inline void set_matrix_A(){}
  PDCD()
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
    for(L i=0;i<n;i++){
    x[i]=gamma*u[i]+z[i];
	}
	for(L j=0;j<m_2;j++){
        Mx[j]=gamma*Mu[j]+Mz[j];
    }
    for(L j=0;j<m_1;j++){
        Ax[j]=gamma*Au[j]+Az[j];
    }
  }

  inline void compute_function_value() {
     D res=0;
     for(L i=0;i<this->n;i++){
       res+=value_of_g_i(x[i],i);
     }
     for(L j=0;j<m_1;j++){
       res+=lambda_f[j]*value_of_phi_j(Ax[j],j);
     }
     for(L j=0;j<m_2;j++){
       res+=value_of_h_j(Mx[j],j);
     }
     function_value=res;
   }

  void initialize(D beta_0, vector<D> & x0, vector<D> & y0, L val_tau, D p_N, D val_lambda_f)
  {
  	set_matrix_M();
    set_matrix_A();

    this->tau=val_tau;
    m_1=data_A.get_n();
    m_2=data_M.get_n();
    cout<<"m_1="<<m_1<<endl;
    cout<<"m_2="<<m_2<<endl;
    lambda_f.resize(m_1,val_lambda_f);
    n=data_A.nfeatures;
    tau=val_tau;
    all_n.resize(n,0);
    for(L i=0;i<n;i++)
    all_n[i]=i;
    print_every_N= p_N;
    nb_of_iters_per_loop=floor(max(1.,n/(tau+0.0)));
    lambda.resize(m_2,0);
    for(L j=0;j<m_2;j++){
        lambda[j]=y0[j];
    }
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
    M_tlambda.resize(n,0);
    compute_Mty();
    Au.clear();
    Au.resize(m_1,0);
    Az.clear();
    Az.resize(m_1,0);
    Ax.clear();
    Ax.resize(m_1,0);
    Mu.clear();
    Mu.resize(m_2,0);
    Mz.clear();
    Mz.resize(m_2,0);
    Mx.clear();
    Mx.resize(m_2,0);
    compute_Az(x0);
    compute_Mz(x0);
    compute_Ax(x0);
    compute_Mx(x0);
    beta_s= beta_0;
    set_v();
    set_p();
    set_rng();

    theta0= min_p;
    theta= theta0;
    cout<< "theta0= "<< theta0<< endl;
    gamma= 1- theta0;

    cout<<"finished SMART_CD initializing"<<endl;
  }

  inline void compute_and_record_res(){
  	    if (nb_iters%print_every_N== 0){
  	    cout<< "gamma= "<< gamma<< " theta= "<< theta<< " beta_s= "<< beta_s<< endl;
        //compute_KKT_residual();
        compute_function_value();
        cout<<setprecision(9)<<"Iteration: "<<nb_iters<<"; time="<<running_time<<"; function value="<<function_value<<endl;
        samp<<setprecision(9)<< nb_iters<<" "<<running_time<<" "<<function_value<<" "<<endl;
        }
   }


  inline D partial_i_of_f(L i)
  {
  	D res=0;
      for (L k = data_A.ptr_t[i]; k < data_A.ptr_t[i + 1];k++)
      {
        L j=data_A.col_idx[k];
        D tmp=lambda_f[j]*gradient_of_phi_j(gamma*Au[j]+Az[j], j);
        res+=data_A.A_t[k]*tmp;
      }
      for (L k = data_M.ptr_t[i]; k < data_M.ptr_t[i + 1];k++)
      {
        L j=data_M.col_idx[k];
        D tmp=prox_of_h_star_j((gamma*Mu[j]+Mz[j])/beta_s+ lambda[j],beta_s,j);
        res+=data_M.A_t[k]*tmp;
      }
      return res;
  }


  inline void set_v()
   {
     v.resize(n,0);
     L_f.resize(n,0);
     L_M.resize(n,0);
     D maxv=0;
     D minv=std::numeric_limits<double>::max();
     D sumv=0;
     D sumvi1=0;
     L sumw=0;
     L maxw=0;
     L minw=n;
     sumofLi=0;
     for(L j=0;j<m_1;j++)
     {
       sumw+=data_A.w_t[j];
       maxw=max(maxw,data_A.w_t[j]);
       minw=min(minw,data_A.w_t[j]);
     }
     cout<<"sumw="<<sumw<<";  maxw="<<maxw<<"; minw="<<minw<<endl;
     for(L i=0;i<n;i++)
     {
       D lfi=0;
       D lmi=0;
       D lfi1= 0;
       D lmi1= 0;
       for (L k = data_A.ptr_t[i]; k < data_A.ptr_t[i + 1];k++)
       {
         L j=data_A.col_idx[k];
         lfi+=(1.+(data_A.w_t[j]-1.)*(tau-1.)/max(n-1.,1.))*data_A.A_t[k]*data_A.A_t[k]*lambda_f[j];
         lfi1+=data_A.A_t[k]*data_A.A_t[k]*lambda_f[j];
       }
       for (L k = data_M.ptr_t[i]; k < data_M.ptr_t[i + 1];k++)
       {
         L j=data_M.col_idx[k];
         lmi+=(1.+(data_M.w_t[j]-1.)*(tau-1.)/max(n-1.,1.))*data_M.A_t[k]*data_M.A_t[k];
         lmi1+=data_M.A_t[k]*data_M.A_t[k];
       }
       L_f[i]=lfi;
       L_M[i]= lmi;
       sumv+=lfi;
       sumvi1+=lfi1;
       if(maxv<lfi) maxv=lfi;
       if(minv>lfi) minv=lfi;
     }
     if(tau==n){
          for(L i=0;i<n;i++)
              L_f[i]=sumvi1;
     }
     sumofLi=sumvi1;
     cout<<"  max of v: "<<maxv<<" ;  min of v: "<<minv<<" ;  sumofv: "<<sumv<<" sumofLi="<<sumofLi<<endl;
   }

   inline void set_p(){
     proba_vector.resize(n,0.0);
     D res= 0;
     D tmp= 1;
     D tmp2= 0;
     for (L i=0; i< n; i++){
     	res+= sqrt(L_f[i]+ L_M[i]/beta_s);
	 }
	 for (L i= 0; i< n; i++){
	 	proba_vector[i]= sqrt(L_f[i]+ L_M[i]/beta_s)/res;
	 	if (proba_vector[i]< tmp){
	 		tmp= proba_vector[i];
		 }
		if (proba_vector[i]> tmp2){
	 		tmp2= proba_vector[i];
		 }
	 }
	 if (tmp== 0){
	 	for (L i= 0; i< n; i++){
	 	proba_vector[i]= 1.0/n;
	    }
	    min_p= 1.0/n;
	    max_p= 1.0/n;
	 }
	 else{
	    min_p= tmp;
	    max_p= tmp2;
     }
	 cout<< "max_p= "<< max_p<< "; min_p= "<< min_p<< endl;
   }


  inline void update_z_coordinate( L i, D dz){
     z[i]+=dz;
     for (L k = data_A.ptr_t[i]; k < data_A.ptr_t[i + 1];k++)
     {
       L j=data_A.col_idx[k];
       Az[j]+=dz*data_A.A_t[k];
     }
     for (L k = data_M.ptr_t[i]; k < data_M.ptr_t[i + 1];k++)
     {
       L j=data_M.col_idx[k];
       Mz[j]+=dz*data_M.A_t[k];
       //if(j==2869) cout<<"gamma="<<gamma<<"; "<<dz<<"; "<<Mz[j]<<endl;
     }
   }

    inline void update_x_coordinate( L i, D dx){
     x[i]+=dx;
     x[i]=x[i];
     L j;
     for (L k = data_M.ptr_t[i]; k < data_M.ptr_t[i + 1];k++)
     {
       j=data_M.col_idx[k];
       Mx[j]+=dx*data_M.A_t[k];
     }
      for (L k = data_A.ptr_t[i]; k < data_A.ptr_t[i + 1];k++)
     {
       j=data_A.col_idx[k];
       Ax[j]+=dx*data_A.A_t[k];
     }
   }

   inline void update_u_coordinate( L i, D du){
     u[i]+=du;
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
         Au[j]+=u[kj]*data_A.A[k];
       }
   }

   void compute_Mu(){
     for(L j=0;j<m_2;j++)
       for(L k = data_M.ptr[j]; k < data_M.ptr[j + 1];k++){
         L kj=data_M.row_idx[k];
         Mu[j]+=u[kj]*data_M.A[k];
       }
   }

   void compute_Az(){
     for(L j=0;j<m_1;j++)
       for(L k = data_A.ptr[j]; k < data_A.ptr[j + 1];k++){
         L kj=data_A.row_idx[k];
         Az[j]+=z[kj]*data_A.A[k];
       }
   }
   void compute_Mz(){
     for(L j=0;j<m_2;j++)
       for(L k = data_M.ptr[j]; k < data_M.ptr[j + 1];k++){
         L kj=data_M.row_idx[k];
         Mz[j]+=z[kj]*data_M.A[k];
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
         Mx[j]+=x0[kj]*data_M.A[k];
       }
   }

   void compute_Ax(vector<D> & x0){
     for(L j=0;j<m_1;j++)
       for(L k = data_A.ptr[j]; k < data_A.ptr[j + 1];k++){
         L kj=data_A.row_idx[k];
         Ax[j]+=x0[kj]*data_A.A[k];
       }
   }

  void compute_Mty(){
        M_tlambda.clear();
        M_tlambda.resize(this->n,0);
        for(L i=0;i<this->n;i++)
        {
          for(L k = data_M.ptr_t[i]; k < data_M.ptr_t[i + 1];k++){
            L j=data_M.col_idx[k];
            M_tlambda[i]+=lambda[j]*data_M.A_t[k];
           }
        }

    }

    void update_theta()
  {

	 theta= (sqrt(theta*theta*theta*theta+ 4*theta*theta)- theta*theta)/2;
	
  }

  void reset(){
  	beta_s= beta_s/(1+ theta);
  }

  void PDCD_solver(D beta_0, vector<D> & x0,vector<D> & y0, L val_tau, L max_nb_epoch, L p_N, D val_lambda_f,string filename1, D time2)
  {
    initialize(beta_0, x0, y0, val_tau,  p_N, val_lambda_f);
    cout<<"running SMART_CD"<<" ; "<<filename1<<" max_nb_epoch "<<max_nb_epoch<<endl;
    nb_iters=0;
    filename1="results/PDCD_"+filename1;
    samp.open(filename1.c_str());
    //srand48(27432042);
    srand(time(NULL));
    running_time= 0;
    print_every_N=p_N;
    compute_and_record_res();
    D start;
    start = std::clock();
    while(nb_iters<max_nb_epoch)
    {
      start = std::clock();
      D ti=0;
      //D si=0;
      D gi=0;
      D Li=0;
      L i=0;
      for(L it=0;it<nb_of_iters_per_loop;it++)
      {
        gamma*=(1-theta);
        if(theta==1) gamma=1;
        batch_sampling();
        for(L it_S=0;it_S<tau;it_S++)
        {
          i=S[it_S];
          gi=partial_i_of_f(i);
          Li=(L_f[i]+ L_M[i]/beta_s)*theta/theta0;
          t[i]=prox_of_g_i(gi, Li, z[i],i);
        }
        for(L it_S=0;it_S<tau;it_S++)
        {
          i=S[it_S];
          ti=t[i];
          update_z_coordinate(i, ti);
          update_u_coordinate(i, -(1- theta/theta0)/gamma*ti);
        }
        update_theta();
        reset();
      }
      running_time+=( std::clock() - start ) / (double) CLOCKS_PER_SEC;
      nb_iters++;
      compute_x();
      compute_and_record_res();
      if (running_time> time2){
      	break;
	  }
    }
    compute_and_record_res();

  }





};

#endif
