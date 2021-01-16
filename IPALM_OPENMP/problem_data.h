#ifndef PROBLEM_DATA_H
#define PROBLEM_DATA_H


using namespace std;
#include "Matrix.h"
#include <math.h>


template<typename L, typename D>
class problem_data: public Matrix<L, D>{
  
public:



protected:





  D lambda;

  D gamma;

  D mu;

  problem_data()
  :Matrix<L,D>()
  {
  }

  problem_data(const char* matrix_file, const char* vector_file)
  :Matrix<L,D>(matrix_file,vector_file)
  {
  }


  problem_data(const char* matrix_file)
  :Matrix<L,D>(matrix_file)
  {
  }



  void set_lambda_and_gamma(D val_lambda, D val_gamma)
  {
    lambda=val_lambda;
    gamma=val_gamma;
  }




  void print_b()
  {
    cout<<this->b[this->nsamples- 3]<<" , "<<this->b[this->nsamples- 1]<<endl;
  }
















};


#endif /*PROBLEM_DATA_H*/
