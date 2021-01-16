#ifndef SMSVM_B_H
#define SMSVM_B_H

#include "ALM_Katyusha.h"
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

//This class solves problem: min_{x,omega} sum_i max(1- b_i^T(A_ix+ omega), 0)+ P(x) by IPALM_Katyusha
// where P(x)=frac{sig2}{2}|x|_2^2 +sig1|x|_1.

// Let X= (x; omega), M= (b.*A b), b_i^T(A_ix+ omega)= M_iX
// f^i(x)= 0, h^i(x)= max(1- x,0).





template<typename L, typename D>
class SMSVM_b: public ALM_Katyusha<L, D>
{
private:

	Matrix<L,D> my_A;

	D val_lambda_f;



protected:

public:

	D sig1;
	D sig2;

	SMSVM_b(const char* Matrix_file,D val_lambda1, D val_lambda2)
	:ALM_Katyusha<L,D>(),my_A(Matrix_file)
	{
		sig1=val_lambda1;
		sig2=val_lambda2;
		this->mu_g= sig2;
		this->L_h= 1;
		this->matrix_modified_to_svm(my_A);
	}

	L get_n(){return my_A.nfeatures;}
	L get_m(){return my_A.nsamples;}

	inline void set_matrix_A(){
		this->data_A.nsamples=0;
		this->data_A.nfeatures=this->data_M.nfeatures;
		this->data_A.nnz=0;
		this->data_A.ptr.resize(1,0);
		this->data_A.ptr_t.resize(this->data_M.nfeatures+1,0);
	}

	inline void set_matrix_M(){
		this->data_M.matrix_modified_to_svm(my_A);
	}

	inline D value_of_f_j(D x, L i){
		return 0;
	}

	inline D gradient_of_f_j(D x, L i){
		return 0;
	}

	inline D value_of_f_star_j(D x, L i){
		return 0;
	}

	inline void rescale(){
		this->lambda1= sig1;
		this->lambda2= sig2+ this->tau_s*this->beta_s;
	}


	inline D value_of_P_j(D x, L j){
		if (j== this->nfeatures- 1){
			return 0;
		}
		else{
			return sig1*fabs(x)+ sig2/2*x*x;
		}
	}

	inline D prox_of_P_j(D x1, D x2, L j){
		if (j== this->nfeatures- 1){
			return x1;
		}
		else{
			if (x1*x2- sig1> 0){
				return (x1*x2- sig1)/(sig2+ x2);
			}
			else if(x1*x2+ sig1< 0){
				return (x1*x2+ sig1)/(sig2+ x2);
			}
			else{
				return 0;
			}
		}
	}


	inline D value_of_h_j(D x, L i){
		if(x<= 1.0){
			return 1.0- x;
		}
		else{
			return 0;
		}
	}

	inline D value_of_h_star_j(D x, L j){
		if (x<= 0.0 && x>= -1.0){
			return x;
		}
		else{
			cout<< "error in h*(x)"<< endl;
			return std::numeric_limits<double>::max();
		}
	}

	inline D compute_one_step(D tau, D u, D x, L j){
		D new_x= 0;
		if(j== this->nfeatures- 1){
			new_x= (x-tau*u)/(1+this->lambda2*tau);
		}
		else{
			if(x>tau*(this->lambda1+u))
			new_x=(x-tau*(this->lambda1+u))/(1+this->lambda2*tau);
			else if(x<tau*(u-this->lambda1))
			new_x=(x-tau*(u-this->lambda1))/(1+this->lambda2*tau);
			else
			new_x=0;
		}
		return new_x;
	}

	inline D prox_of_h_star_j(D x1, D x2, L j){
		if (x1- 1.0/x2> 0.0){
			return 0.0;
		}
		else if (x1- 1.0/x2< -1.0){
			return -1.0;
		}
		else{
			return x1- 1.0/x2;
		}
	}


	inline D feasible_dual(vector<D> & x){

		if(this->lambda2>0)
		{
			return 1;
		}
		else
		{
			D scal=1;
			L l=x.size();
			for(L i=0;i<l;i++)
			if(fabs(x[i])>this->lambda1)
			scal=min(this->lambda1/fabs(x[i]),scal);
			return scal;
		}
	}

	/*  inline D distance_to_subgradient_of_g(){
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
*/
//We implement Section 6.2 of SPDC paper for 'just-in-time' update
//t0 is t0+1 in the paper SPDC ; x is the z vector in Katyusha
inline void compute_just_in_time_prox_grad_without_x_s(D tau, D u, D &x, L t0,L t1, D w, D &y, L j){
	if(t0==t1)
	return;
	else{
		if(this->lambda2>0)
		{
			D theta=1./(1+this->lambda2*tau);
			D hplus=(u+this->lambda1)/this->lambda2;
			D hminus=(u-this->lambda1)/this->lambda2;
			D p=pow(theta,t1-t0);
			if(this->lambda1==0||j==this->nfeatures-1){
				hplus=u/this->lambda2;
				y=compute_aid_2(p,theta, hplus, t1-t0, x, y, w);
				x=p*(x)-(1-p)*u/this->lambda2;
			}else{
				if(x==0){
					if(this->lambda1+u<0){
						y=compute_aid_2(p,theta, hplus, t1-t0, x, y, w);
						x=p*x-(1-p)*hplus;
					}
					else if(u-this->lambda1>0){
						y=compute_aid_2(p,theta, hminus, t1-t0, x, y, w);
						x=p*x-(1-p)*hminus;
					}else{
						x=0;
						D p2=pow(this->theta3,t1-t0);
						D tmp2=(1-p2)/(this->theta1+this->theta2);
						y=(this->theta2*w)*tmp2+p2*y;
					}
				}else if(x>0){
					if(this->lambda1+u<=0){
						y=compute_aid_2(p,theta, hplus, t1-t0, x, y, w);
						x=p*x-(1-p)*hplus;
					}else{
						D t=t0-log(1+this->lambda2*x/(this->lambda1+u))/log(theta);
						if(t<t1){
							L t_tmp=floor(t);
							p=pow(theta,t_tmp-t0);
							y=compute_aid_2(p,theta, hplus, t_tmp-t0, x, y, w);
							x=p*x-(1-p)*hplus;
							x=compute_one_step(tau,u,x,j);
							y=this->theta1*x+this->theta2*w+this->theta3*y;
							compute_just_in_time_prox_grad_without_x_s(tau, u,  x,  t_tmp+1, t1, w, y,j);
						}
						else{
							y=compute_aid_2(p,theta, hplus, t1-t0, x, y, w);
							x=p*x-(1-p)*hplus;
						}
					}

				}else if(x<0){
					D minusw=-w;
					D minusx=-x;
					D minusy=-y;
					compute_just_in_time_prox_grad_without_x_s(tau, -u,  minusx,  t0, t1,minusw, minusy,j);
					x=-minusx;
					y=-minusy;
				}

			}

		}
		else{

			if(this->lambda1==0||j==this->nfeatures-1){
				y=compute_aid(tau*u, t1-t0, x, y, w);
				x=x-(t1-t0)*tau*u;
			}
			else{
				if(x==0){
					if(this->lambda1+u<0){
						y=compute_aid(tau*(this->lambda1+u), t1-t0, x, y, w);
						x=x-(t1-t0)*tau*(this->lambda1+u);
					}
					else if(u-this->lambda1>0){
						y=compute_aid(tau*(-this->lambda1+u), t1-t0, x, y, w);
						x=x-(t1-t0)*tau*(u-this->lambda1);
					}else{
						D p2=pow(this->theta3,t1-t0);
						D tmp2=(1-p2)/(this->theta1+this->theta2);
						y=(this->theta2*w)*tmp2+p2*y;
						x=0;
					}
				}else if(x>0){
					if(this->lambda1+u<=0){
						y=compute_aid(tau*(this->lambda1+u), t1-t0, x, y, w);
						x=x-(t1-t0)*tau*(this->lambda1+u);
					}else{
						D t=t0+x/(this->lambda1+u);
						if(t<t1){
							L t_tmp=floor(t);
							y=compute_aid(tau*(this->lambda1+u), t_tmp-t0, x, y, w);
							x=x-(t_tmp-t0)*tau*(this->lambda1+u);
							x=compute_one_step(tau,u,x,j);
							y=this->theta1*x+this->theta2*w+this->theta3*y;
							compute_just_in_time_prox_grad_without_x_s(tau, u,  x,  t_tmp+1, t1, w, y,j);
						}
						else{
							y=compute_aid(tau*(this->lambda1+u), t1-t0, x, y, w);
							x=x-(t1-t0)*tau*(this->lambda1+u);
						}
					}

				}else if(x<0){
					D minusw=-w;
					D minusx=-x;
					D minusy=-y;
					compute_just_in_time_prox_grad_without_x_s(tau, -u,  minusx,  t0, t1,minusw, minusy,j);
					x=-minusx;
					y=-minusy;
				}

			}
		}
	}

};

// compute_aid returns the value of y_{k+s} given by: y_{k+1}=theta1*z_{k+1}+theta2*w+(1-theta1-theta2)*y_k; z_{k+i}=z_k-i*h
D compute_aid(D h, D s, D x, D y, D w){
	D p2=pow(this->theta3,s);
	D tmp2=(1-p2)/(this->theta1+this->theta2);
	D tmp3=(tmp2*this->theta3-p2*s)/(this->theta1+this->theta2);
	y=this->theta1*(tmp2*(x-s*h)+h*tmp3)+this->theta2*w*tmp2+p2*y;
	return y;
}


// compute_aid_2 returns the value of y_{k+s} given by: y_{k+1}=theta1*z_{k+1}+theta2*w+(1-theta1-theta2)*y_k; z_{k+i}=q^i(z_k+h)-h
//p=q^s;
D compute_aid_2(D p, D q, D h, D s, D x, D y, D w){
	D theta3overq=this->theta3/q;
	D p2=pow(this->theta3,s);
	D p3=pow(theta3overq,s);
	D tmp=(1-p3)/(1-theta3overq);
	D tmp2=(1-p2)/(this->theta1+this->theta2);
	y=this->theta1*p*(x+h)*tmp+(this->theta2*w-this->theta1*h)*tmp2+p2*y;
	return y;
}



void Katyusha_solver(D beta_0, D epsilon_0,  D eta, D rho, vector<D> & x0,vector<D> & y0,L val_tau, L max_nb_outer, L p_N_1, L p_N_2,string filename1, string filename2, D time){
	this->ALM_solve_with_L_Katyusha(beta_0, epsilon_0,  eta, rho, x0,y0, val_tau, max_nb_outer,  p_N_1,  p_N_2, filename1,  filename2, time);
}




};

#endif
