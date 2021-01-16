#ifndef MATRIX_H
#define	MATRIX_H


#include <cstdlib>
#include <map>
#include <vector>
#include <stdio.h>
#include <stdlib.h>
#include <iostream>
#include <fstream>
#include <limits>

#include "svm_parser.h"

using namespace std;


// Matrix A=(A_1,A_2,...,A_n) where A_i is a vector of length d
/* The matrix A is stored using three vectors: one vector of length nnz containing the values of nonzero elements in A, one vector
of length n+1 containing the number of nonzero elements summed up to each column, one vector of length nnz containing the row index of each nonzero element .
*/
template <typename L, typename D>
class Matrix
{
public:

  typedef vector<L> pointers;
  typedef vector<D> values;
  L nsamples;
  L nfeatures;

  L nnz;
  pointers ptr;
  pointers row_idx;
  values A;
  pointers w;  // number of nonzeros elements in each row of A
  values b;

  values A_t;
  pointers ptr_t;
  pointers col_idx;
  pointers w_t; // number of nonzero elements in each column of A

  L*  index_of_samples;

  vector<L> wprime;

  L c;

  L noverc;

public:
  Matrix()
  {}


    Matrix(const char* filename,const char* vector_file)
    {
      char mystring [100];
      FILE* pFile;
      cout<<filename<<endl;
      pFile=fopen(filename, "r");
      if (pFile == NULL) perror ("Error opening file");
      else {
        if ( fgets (mystring , 100 , pFile) != NULL )
        nsamples=(L)atof(mystring);
        if ( fgets (mystring , 100 , pFile) != NULL )
        nfeatures=(L)atof(mystring);
        if (fgets (mystring , 100 , pFile) != NULL )
        nnz=(L)atof(mystring);
        cout << "Data file " << filename << " contains " << nfeatures << " features, "
        << nsamples << " samples " << "and total "
        << nnz << " non zero elements" << endl;

        ptr.resize(nsamples+1,0);
        A.resize(nnz,0);
        row_idx.resize(nnz,0);
        w.resize(nfeatures,0);

        L k=0;
        while ( k<=nsamples && fgets (mystring , 100 , pFile) != NULL  ){
          ptr[k]=(L)atof(mystring);//cout<<" k: "<<k<<" ptr[k]: "<<ptr[k];
          k++;
        }
        cout<<endl;
        k=0;
        while ( k<nnz && fgets (mystring , 100 , pFile) != NULL  ){
          row_idx[k]=(L)atof(mystring);//cout<<" k: "<<k<<" row_idx[k]: "<<row_idx[k];
          w[row_idx[k]]++;
          k++;
        }
        cout<<endl;
        k=0;
        while (k<nnz && fgets (mystring , 100 , pFile) != NULL  ){
          A[k]=atof(mystring);
          k++;
        }
        cout<<endl;
      }

      b.resize(nsamples);
      char mystring_b [100];
      FILE* pFile_b;
      pFile_b=fopen(vector_file, "r");
      L i=0;
      if (pFile_b == NULL) perror ("Error opening file");
      else
      {
        while ( fgets (mystring_b , 100 , pFile_b) != NULL ){
          if(i>nsamples){
            perror ("wrong size of vector"); break;
          }
          b[i]=atof(mystring_b);
          i++;
        }
      }
      //read_Matrix();
    }




    void write_matrix_to_file() {

      ofstream matrixfile;
      ofstream vectorfile;
      matrixfile.open("datas/matrix_test1");
      matrixfile << nsamples << endl;
      matrixfile << nfeatures << endl;
      matrixfile << nnz << endl;
      for(L i=0;i<nsamples+1;i++)
      matrixfile<<ptr[i]<<endl;
      for(L i=0;i<nnz;i++)
      matrixfile<<row_idx[i]<<endl;
      for(L i=0;i<nnz;i++)
      matrixfile<<A[i]<<endl;
      matrixfile.close();
      vectorfile.open("datas/vector_test1");
      for(L j=0;j<nsamples;j++)
      vectorfile<<b[j]<<endl;
      cout<<"writing to file"<<endl;

    }

    void write_matrix_to_file_t(const char* filename) {
      transpose_matrix();
      ofstream matrixfile;
      ofstream vectorfile;
      matrixfile.open(filename);
      matrixfile << nfeatures << endl;
      matrixfile << nsamples << endl;
      matrixfile << nnz << endl;
      for(L i=0;i<nfeatures+1;i++)
      matrixfile<<ptr_t[i]<<endl;
      for(L i=0;i<nnz;i++)
      matrixfile<<col_idx[i]<<endl;
      for(L i=0;i<nnz;i++)
      matrixfile<<A_t[i]<<endl;
      matrixfile.close();
      vectorfile.open(filename);
      for(L j=0;j<nsamples;j++)
      vectorfile<<b[j]<<endl;
      cout<<"writing to file"<<endl;

    }

    void normalize_matrix(){
      std::vector<D> vi1(nsamples);
      for(L i=0;i<nsamples;i++)
      {
        for (L k = ptr[i]; k < ptr[i + 1];k++)
        {
          vi1[i]+=this->A[k]*this->A[k];
        }
      }
      for(L i=0;i<nsamples;i++)
      {
        for (L k = ptr[i]; k < ptr[i + 1];k++)
        {
          this->A[k]=this->A[k]/sqrt(vi1[i]);
        }
      }
    }

    void transpose_matrix()
    {
      w_t.resize(nsamples,0);
      ptr_t.resize(nfeatures+1,0);
      col_idx.resize(nnz,0);
      A_t.resize(nnz);
      vector<L> nb_nnz(nfeatures);
      for(L i=0;i<nsamples;i++)
      {
        for (L k = ptr[i]; k < ptr[i + 1];k++)
        {
          L j=this->row_idx[k];
          nb_nnz[j]++;
        }
        w_t[i]=ptr[i+1]-ptr[i];
      }
      for(L j=0; j<nfeatures;j++)
      ptr_t[j+1]=ptr_t[j]+nb_nnz[j];
      nb_nnz.clear();
      nb_nnz.resize(nfeatures,0);
      for(L i=0;i<nsamples;i++)
      {
        for (L k = ptr[i]; k < ptr[i + 1];k++)
        {
          L j=this->row_idx[k];
          col_idx[ptr_t[j]+nb_nnz[j]]=i;
          A_t[ptr_t[j]+nb_nnz[j]]=A[k];
          nb_nnz[j]++;
        }
      }
      cout<<"transposation finished"<<endl;
      //read_Matrix();
    }

    void distributed_set_up(L nb_c)
    {
      if(nsamples%nb_c!=0)
      {
        cout<<"nsamples "<<nsamples<<endl;
        cout<<"nb_c "<<nb_c<<endl;
        cout<<"nsamples%nb_c"<<nsamples%nb_c<<endl;
        perror("number of nodes must be divisible by number of samples");
        abort();
      }
      c=nb_c;
      noverc=nsamples/c;
      L* permutation_m=new_order(nsamples,1);
      if(c==1)
      {
        index_of_samples=permutation_m;
        wprime.clear();
        wprime.resize(nfeatures,1);
      }
      else
      {
        wprime.clear();
        wprime.resize(nfeatures,0);
        transpose_matrix();
        permutation(permutation_m, nsamples, 1);
        index_of_samples=new L [nsamples];
        for(L i=0;i<nsamples;i++)
        {
          index_of_samples[permutation_m[i]]=i;
        }
        cout<<"dddddf"<<endl;
        for(L j=0;j<nfeatures;j++)
        {
          vector<int> tmp(c);
          for (L k = ptr_t[j]; k < ptr_t[j + 1];k++)
          {
            L i=col_idx[k];
            L tmpc=(permutation_m[i])/noverc;
            //cout<<"i"<<permutation_m[i]<<"   tmpc"<<tmpc<<" ; c:"<<c<<"noverc: "<<noverc<<endl;
            if(tmp[tmpc]==0)
            {
              wprime[j]++;
              tmp[tmpc]=1;
            }
          }
          if(j%1000==0) cout<<wprime[j]<<" "<<"'"<<w[j]<<"'"<<" ";
        }
        cout<<endl;
      }
    }

    void read_Matrix()
    {
      L t=10;
      if(t>nsamples) t=nsamples;

      for(L i=0;i<t;i++)
      {
        cout<<"column "<<i<<" : ";
        for (L k = ptr[i]; k < ptr[i + 1];k++) {
          L kj=row_idx[k];
          cout<<"A["<<kj<<"]="<<A[k]<<" ";
        }
        cout<<endl;
      }
      for(L j=0;j<10;j++)
      cout<<"w"<<j<<": "<<w[j]<<" ";
      cout<<endl;

    }

    L get_n(){return nsamples;}
    L get_d(){return nfeatures;}
    L get_nnz(){return nnz;}

    void parse_LIB_SVM_data_get_size(const char* filename, L &nsamples,
      L &nfeatures, L &nonzero_elements_of_input_data) {
        nfeatures = 0;
        nsamples = 0;
        nonzero_elements_of_input_data = 0;
        FILE* file
        = fopen(filename, "r");
        if (file == 0) {
          printf("File '%s' not found\n", filename);
          perror ("Error opening file");
        }
        char* stringBuffer = (char*) malloc(65536);
        bool end_of_file = false;
        nsamples = -1;
        //      for (int i = 0; i < nsamples; i++) {
        while (!end_of_file) { //feof(file)
          nsamples++;
          char c;
          int pos = 0;
          char* bufferPointer = stringBuffer;
          do {
            c = fgetc(file);

            if (c == EOF) {
              end_of_file = true;
              break;
            }

            if ((c == ' ') || (c == '\n')) {
              if (pos == 0) {
                //Label found
                *(bufferPointer) = 0;
                int value;
                sscanf(stringBuffer, "%i", &value);
                //printf ("string: %s\n",stringBuffer);

                if (value < 100)
                pos++;

              } else {
                //Feature found
                *(bufferPointer) = 0;
                double value;
                sscanf(stringBuffer, "%lf", &value);
                if (pos > 0) {
                  if (nfeatures < pos)
                  nfeatures = pos;
                  pos = -1;
                  nonzero_elements_of_input_data++;
                }
              }
              bufferPointer = stringBuffer;
            } else if (c == ':') {
              //Position found
              *(bufferPointer) = 0;
              int value;
              sscanf(stringBuffer, "%i", &value);
              pos = value;
              bufferPointer = stringBuffer;
            } else {
              *(bufferPointer) = c;
              bufferPointer++;
            }

          } while (c != '\n');
        }
        free(stringBuffer);
        fclose(file);

      }

      Matrix(const char* filename) {

        int nclasses=0;

        cout << "Going to parse SVM data" << endl;
        parse_LIB_SVM_data_get_size(filename, nsamples, nfeatures,nnz);
        cout << "Data file " << filename << " contains " << nfeatures << " features, "
        << nsamples << " samples " << "and total "
        << nnz << " nonzero elements" << endl;

        FILE* filePtr = fopen(filename, "r");
        if (filePtr == 0) {
          printf("File   '%s' not found\n", filename);
          perror ("no File Found!");
        }

        cout << "Going to process data" << endl;

        ptr.resize(nsamples+1);
        A.resize(nnz);
        row_idx.resize(nnz);
        w.resize(nfeatures);
        b.resize(nsamples);

        L nnzPossition = 0;
        L processedSamples = -1;

        char* stringBuffer = (char*) malloc(65536);
        for (L i = 0; i < nsamples; i++) {
          char c;
          L pos = 0;
          char* bufferPointer = stringBuffer;
          do {
            c = fgetc(filePtr);
            //
            if ((c == ' ') || (c == '\n')) {
              if (pos == 0) {
                //Label found
                *(bufferPointer) = 0;
                double value;
                sscanf(stringBuffer, "%lf", &value);

                D ddval = value;
                if (value < 100) {
                  if (nclasses == 2 && value == 0) {
                    ddval = (float) -1;
                  } else {
                  }

                  processedSamples++;
                  b[processedSamples] = ddval;
                  ptr[processedSamples] = nnzPossition;

                  pos++;
                }
              } else {
                //Feature found
                *(bufferPointer) = 0;
                double value;
                sscanf(stringBuffer, "%lf", &value);

                if (pos > 0) {
                  pos--;  // zero based index in C++
                  row_idx[nnzPossition] = pos;
                  w[pos]++;
                  A[nnzPossition] = value;
                  nnzPossition++;
                  pos = -1;
                }

              }
              bufferPointer = stringBuffer;
            } else if (c == ':') {
              //Position found
              *(bufferPointer) = 0;
              int value;
              sscanf(stringBuffer, "%i", &value);
              pos = value;
              bufferPointer = stringBuffer;
            } else {
              *(bufferPointer) = c;
              bufferPointer++;
            }

          } while (c != '\n');
        }

        processedSamples++;
        ptr[processedSamples] = nnzPossition;

        free(stringBuffer);
        fclose(filePtr);

        L sumnnz=0;
        for(L j=0;j<nfeatures;j++)
        sumnnz+=w[j];

        cout<<"sumnnz"<<sumnnz;
        cout<<endl;
        //write_matrix_to_file();
        //normalize_matrix();
        //read_Matrix();
        transpose_matrix();
      }


      L * new_order(L num, L ic)
      {
        L * order = new L [num*ic];
        for(L j=0;j<ic;j++)
        for (L i=0; i<num; i++)
        {
          order[i+j*num]=i;
        }
        return order;
      }

      void permutation(L * order, L num, L ic)
      {
        L i,j,k;
        for(k=0;k<ic;k++)
        for (i=0; i<num; i++)
        {
          j=i + (int)((num-i)*rand() / (double)RAND_MAX);
          if (j>=num) j=num-1;
          if (j>i)
          {
            int tmp=order[i+k*num];
            order[i+k*num]=order[j+k*num];
            order[j+k*num]=tmp;
          }
        }
      }

      void matrix_merge(Matrix<L,D> M, Matrix<L,D> N){
        L  n1=M.get_d();
        L  n2=M.get_n();
        L nnz_M=M.get_nnz();
        L  n3=N.get_n();
        L nnz_N=N.get_nnz();
        nsamples=n2+ n3;
        nfeatures=n1;
        nnz=nnz_M+nnz_N;
        ptr.clear();
        ptr.resize(nsamples+ 1,0);
        A.clear();
        A.resize(nnz,0);
        row_idx.clear();
        row_idx.resize(nnz,0);
        w.clear();
        w.resize(nfeatures,0);
        L pt_i=0;
        for(L i=0;i<n2;i++){
          ptr[i]= M.ptr[i];
          for(L k=M.ptr[i];k<M.ptr[i+1];k++)
          {
            A[pt_i]=M.A[k];
            row_idx[pt_i]=M.row_idx[k];
            w[row_idx[pt_i]]++;
            pt_i++;
          }
        }
        for(L i=0;i<n3;i++){
          ptr[n2+ i]= nnz_M+ N.ptr[i];
          for(L k=N.ptr[i];k<N.ptr[i+1];k++)
          {
            A[pt_i]=N.A[k];
            row_idx[pt_i]=N.row_idx[k];
            w[row_idx[pt_i]]++;
            pt_i++;
          }
        }
        ptr[nsamples]= nnz;
        b.clear();
        b.resize(nsamples,0);
        for(L i=0; i< n2; i++){
          b[i]= M.b[i];
        }
        for(L i=0; i< n3; i++){
          b[n2+ i]= N.b[i];
        }
        transpose_matrix();
      }

      void matrix_modified_to_svm(Matrix<L,D> M){
        nsamples= M.nsamples;
        nfeatures=M.nfeatures+ 1;
        nnz=M.nnz+ M.nsamples;
        ptr.clear();
        ptr.resize(nsamples+ 1,0);
        A.clear();
        A.resize(nnz,0);
        row_idx.clear();
        row_idx.resize(nnz,0);
        w.clear();
        w.resize(nfeatures,0);
        b.clear();
        b.resize(nsamples,0);
        L pt_i= 0;
        for(L i=0;i<nsamples;i++){
          ptr[i]= M.ptr[i]+ i;
          for(L k=M.ptr[i];k<M.ptr[i+1];k++)
          {
            A[pt_i]=M.A[k];
            row_idx[pt_i]=M.row_idx[k];
            w[row_idx[pt_i]]++;
            pt_i++;
          }
          A[pt_i]= 1;
          row_idx[pt_i]= nfeatures- 1;
          w[row_idx[pt_i]]++;
          pt_i++;
          b[i]= 0;
        }
        ptr[nsamples]= nnz;
        for(L i=0;i<nsamples;i++)
        {
          for (L k = ptr[i]; k < ptr[i + 1];k++) {
            A[k]*=M.b[i];
          }
        }
        transpose_matrix();
      }

      void construct_fused_matrix(Matrix<L,D> M){
        nfeatures= M.nfeatures;
        nsamples= nfeatures;
        nnz= 2*nfeatures;
        ptr.clear();
        ptr.resize(nsamples+ 1,0);
        A.clear();
        A.resize(nnz,0);
        row_idx.clear();
        row_idx.resize(nnz,0);
        w.clear();
        w.resize(nfeatures,2);
        b.clear();
        b.resize(nsamples,0);
        for(L i=0;i<nsamples;i++){
          ptr[i]= 2*i;
          A[2*i]= 1;
          A[2*i+ 1]= -1;
          row_idx[2*i]=i;
          if (i== nsamples- 1){
            row_idx[2*i+ 1]= 0;
          }
          else{
            row_idx[2*i+ 1]= i+ 1;
          }
          b[i]= 0;
        }
        ptr[nsamples]= nnz;
        transpose_matrix();
      }
    };

    #endif /*MATRIX_H*/
