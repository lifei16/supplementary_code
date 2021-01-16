To compile, please do:


>> cd PUT_PATH_TO_ROOT_WHERE_THIS_README_FILE_IS
>> export OMP_NUM_THREADS=4
>> g++-10  -I/usr/local/include -L/usr/local/lib -o main main.cpp -lgsl -lgslcblas -fopenmp


Then run the code following the same commands as in ../readme.txt for each experiment. 

For the following three problems:

1. Basis pursuit: f(x)= 0, g(x)= lambda1*||x||_1, h(x)= I_{x= b}
2. Least absolute deviation: f(x)= 0, g(x)= lambda1*||x||_1, h(x)= lambda3*||x||_1
3. Fused lasso: f(x)= 0.5*||Ax- b||_2^2, g(x)= lambda1*||x||_1, h(x)= lambda3*||x||_1


Parallel implementation is useful for speedup the following three algorithm:

b. IPALM_Katyusha
d. ASGARD_DL
e. LADMM


For the problem 

4. L_1 norm support vector machine: f(x)= 0, g(x,w)= lambda1*||x||_1, h(x)= max(0,1- x)

Parallel implementation is useful for speedup the following two algorithm:

d. ASGARD_DL
e. LADMM
