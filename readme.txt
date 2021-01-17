General instruction

This small package contains the code which was used to implement algorithms and generate plots in the paper

<< F. Li and Z. Qu: An inexact proximal augmented Lagrangian framework with arbitrary linearly convergent inner solver for composite convex optimization. >>


=================================================


The optimization model which can be solved with this package is:
minimize {f(x)+ g(x)+ h(p(x))}
where f(x) is convex and differentiable on an open set containing dom(g), g(x) is proper, convex and closed, h(x) is proper, convex and closed or an indicator function of a convex and closed set, and p(x) is a linear operator (matrix) or quadratic operator.



==============
Problem types:
==============

1. Basis pursuit: f(x)= 0, g(x)= lambda1*||x||_1, h(x)= I_{x= b}
2. Least absolute deviation: f(x)= 0, g(x)= lambda1*||x||_1, h(x)= lambda3*||x||_1
3. Fused lasso: f(x)= 0.5*||Ax- b||_2^2, g(x)= lambda1*||x||_1, h(x)= lambda3*||x||_1
4. L_1 norm support vector machine: f(x)= 0, g(x,w)= lambda1*||x||_1, h(x)= max(0,1- x)
5. Quadratically constrained quadratic program: f(x)= 0.5*x^TQ_0x+ b_0x, g(x)= I_{-b<= x<= b}, h(x)= I_{x<= 0), p(x)= (0.5*x^TQ_ix+ b_ix- 1)_{i=1,...,m}




More detailed description about the problem types and algorithms can be found in the paper.


==================================================


- The folder IPALM/ contains the C++ code implementing the following five algorithms:

a. IPALM-APPROX
b. IPALM-Katyusha
c. SMART_CD
d. ASGARD_DL
e. LADMM

Please open IPALM/readme.txt for more explanation about the structure of the code.


- The folder datas/ contains 20 datasets news20scale2, rcv1, rcv2, rcv1mc, news20binary, news20binary_bp, ijcnn1, w4a, w6a, w7a, w8a, a7a, a8a, a9a, real-sim, covtype, downloaded or modified from https://www.csie.ntu.edu.tw/~cjlin/libsvm/; qcqp1, qcqp2, qcqp3, qcqp4, generated randomly, that are used in our experiments. Since real-sim is too large to upload, readers need to download it by themselves and save it as matrix_realsim. Datasets news20binary, news20binary_bp and covtype are zip files. So readers need unzip them before running the following commands. Since qcqp3 and qcqp4 are too large to upload, readers can generate the datasets randomly by qcqp_generate.m. For detalied instruction, please refer to the readme.txt in folder CVX/.


- The folder MATLAB_code/ contains the Matlab code used to generate plots.



- The folder CVX/ contains the code of directly calling the cvx package to solve the problems. Please open CVX/readme.txt for instruction.



- The folder libsvm/ contains the code of directly calling the LIBSVM package to solve the problems. Please open libsvm/readme.txt for instruction.


==================================================


To compile, please do:

>> cd PUT_PATH_TO_ROOT_WHERE_THIS_README_FILE_IS
>> cd IPALM
>> g++ -o main main.cpp -lgsl -lgslcblas

To run, please type:
>> ./main arg1 arg2 arg3 arg4 arg5 arg6 arg7 arg8 arg9 arg10

arg1: problem_type;
arg2: algorithm;
arg3: filename;
arg4: lambda1;
arg5: lambda3;
arg6: maximum number of outer iterations;
arg7: maximum running time;
arg8: (optional): beta_0, default value= 1; 
arg9: (optional): epsilon_0, default value=0;
arg10: (optional): p_N_2 (print or compute every p_N_2 iterations), default value=100;

===================================================

For example, to solve the least absolute deviation problem with lambda1=0.01, lambda3= 1 in the paper, using IPALM-APPROX, beta_0= 1, epsilon_0= 100, p_N_2= 100, do: 

>>./main 2 a news20scale2 0.01 1 100000 100 1 100 100 

Note if using IPALM, the value epsilon_0 cannot be 0.

If instead you want to solve the problem with SMART_CD, and beta_0= 1000, p_N_2= 50, do:

>>./main 2 c news20scale2 0.01 1 100000 100 1000 50  




The results are saved in the folder results/.

==================================================


==================================================

Reproduction of figures in the paper
==================================================

To reproduce Figure 1(a) in our paper, do:

>> cd PUT_PATH_TO_ROOT_WHERE_THIS_README_FILE_IS
>> cd IPALM
>> g++ -o main main.cpp -lgsl -lgslcblas
>> ./main 2 a news20scale2 0.01 1 100000 800 1 100 100 ; ./main 2 c news20scale2 0.01 1 100000 1200 1000 50 ; ./main 2 d news20scale2 0.01 1 100000 1200 10 ; ./main 2 e news20scale2 0.01 1 1000000 1200 0.1 100


This should take roughly 800+1200+1200+1200=4400 seconds. The results are saved in the folder results/. 

Then run the following commands to generate the plot from the results.

>>cd ..
>>cd MATLAB_code
>>matlab 
>>run plot_news20scale2_lasso





======================

To reproduce Figure 1(b) in our paper, do:

>> cd PUT_PATH_TO_ROOT_WHERE_THIS_README_FILE_IS
>> cd IPALM
>> g++ -o main main.cpp -lgsl -lgslcblas
>> ./main 2 a rcv1 0.01 1 100000 600 1 100 100 ; ./main 2 c rcv1 0.01 1 100000 1200 1000 50 ; ./main 2 d rcv1 0.01 1 100000 1200 10 ; ./main 2 e rcv1 0.01 1 1000000 1200 0.1 100


This should take roughly 600+1200+1200+1200=4200 seconds. The results are saved in the folder results/. 

Then run the following commands to generate the plot from the results.

>>cd ..
>>cd MATLAB_code
>>matlab 
>>run plot_rcv1_lasso





======================
To reproduce Figure 1(c) in our paper, do:

>> cd PUT_PATH_TO_ROOT_WHERE_THIS_README_FILE_IS
>> cd IPALM
>> g++ -o main main.cpp -lgsl -lgslcblas
>> ./main 2 a rcv1mc 0.01 1 100000 600 1 100 100 ; ./main 2 c rcv1mc 0.01 1 100000 1200 1000 50 ; ./main 2 d rcv1mc 0.01 1 100000 1200 10 ; ./main 2 e rcv1mc 0.01 1 1000000 1200 0.1 100


This should take roughly 600+1200+1200+1200=4200 seconds. The results are saved in the folder results/. 

Then run the following commands to generate the plot from the results.


>>cd ..
>>cd MATLAB_code
>>matlab 
>>run plot_rcv1mc_lasso




======================
To reproduce Figure 1(d) in our paper, do:

>> cd PUT_PATH_TO_ROOT_WHERE_THIS_README_FILE_IS
>> cd IPALM
>> g++ -o main main.cpp -lgsl -lgslcblas
>> ./main 2 a news20binary 0.01 1 100000 3600 0.1 100 50 ; ./main 2 c news20binary 0.01 1 100000 3600 10 10 ; ./main 2 d news20binary 0.01 1 100000 3600 0.1 ; ./main 2 e news20binary 0.01 1 1000000 3600 0.1 100


This should take roughly 3600+3600+3600+3600=144000 seconds. The results are saved in the folder results/. 

Then run the following commands to generate the plot from the results.

>>cd ..
>>cd MATLAB_code
>>matlab 
>>run plot_news20binary_lasso





======================
To reproduce Figure 2(a) and 3(a) in our paper, do:

>> cd PUT_PATH_TO_ROOT_WHERE_THIS_README_FILE_IS
>> cd IPALM
>> g++ -o main main.cpp -lgsl -lgslcblas
>> ./main 1 a news20scale2 0.01 1 100000 1000 0.01 1 50 ; ./main 1 c news20scale2 1 1 100000 1800 1 50 ; ./main 1 d news20scale2 1 1 100000 2400 1 ; ./main 1 e news20scale2 0.01 1 1000000 1800 0.1 100


This should take roughly 1000+1800+2400+1800=7000 seconds. The results are saved in the folder results/. 

Then run the following commands to generate the plot from the results.

>>cd ..
>>cd MATLAB_code
>>matlab 
>>run plot_news20scale2_bp





======================
To reproduce Figure 2(b) and 3(b) in our paper, do:

>> cd PUT_PATH_TO_ROOT_WHERE_THIS_README_FILE_IS
>> cd IPALM
>> g++ -o main main.cpp -lgsl -lgslcblas
>> ./main 1 a rcv2 0.01 1 100000 2400 0.01 1 50 ; ./main 1 c rcv2 0.01 1 100000 2400 1 50 ; ./main 1 d rcv2 0.01 1 100000 2400 1 ; ./main 1 e rcv2 0.01 1 1000000 2400 0.1 100


This should take roughly 2400+2400+2400+2400=9600 seconds. The results are saved in the folder results/. 

Then run the following commands to generate the plot from the results.

>>cd ..
>>cd MATLAB_code
>>matlab 
>>run plot_rcv1_bp





======================
To reproduce Figure 2(c) and 3(c) in our paper, do:

>> cd PUT_PATH_TO_ROOT_WHERE_THIS_README_FILE_IS
>> cd IPALM
>> g++ -o main main.cpp -lgsl -lgslcblas
>> ./main 1 a rcv1mc 0.01 1 100000 3600 0.01 1 50 ; ./main 1 c rcv1mc 0.01 1 100000 3600 1 50 ; ./main 1 d rcv1mc 0.01 1 100000 4200 1 ; ./main 1 e rcv1mc 0.01 1 100000 3600 0.1 100


This should take roughly 3600+3600+4200+3600=15000 seconds. The results are saved in the folder results/. 

Then run the following commands to generate the plot from the results.

>>cd ..
>>cd MATLAB_code
>>matlab 
>>run plot_rcv1mc_bp





======================
To reproduce Figure 2(d) and 3(d) in our paper, do:

>> cd PUT_PATH_TO_ROOT_WHERE_THIS_README_FILE_IS
>> cd IPALM
>> g++ -o main main.cpp -lgsl -lgslcblas
>> ./main 1 a news20binary_bp 0.01 1 100000 2400 0.01 1 50 ; ./main 1 c news20binary_bp 0.01 1 100000 2400 0.1 10 ; ./main 1 d news20binary_bp 0.01 1 100000 2400 0.1 ; ./main 1 e news20binary_bp 0.01 1 100000 2400 0.01 100


This should take roughly 2400+2400+2400+2400=9600 seconds. The results are saved in the folder results/. 

Then run the following commands to generate the plot from the results.

>>cd ..
>>cd MATLAB_code
>>matlab 
>>run plot_news20binary_bp





======================
To reproduce Figure 4(a) in our paper, do:

>> cd PUT_PATH_TO_ROOT_WHERE_THIS_README_FILE_IS
>> cd IPALM
>> g++ -o main main.cpp -lgsl -lgslcblas
>> ./main 3 a news20scale2 0.01 0.01 100000 2400 1 100 100 ; ./main 3 b news20scale2 0.01 0.01 100000 600 1 1 40 ; ./main 3 c news20scale2 0.01 0.01 100000 2400 10 50 ; ./main 3 d news20scale2 0.01 0.01 100000 2400 10 ; ./main 3 e news20scale2 0.01 0.01 1000000 2400 0.1 100


This should take roughly 2400+600+2400+2400+2400=10200 seconds. The results are saved in the folder results/. 

Then run the following commands to generate the plot from the results.

>>cd ..
>>cd MATLAB_code
>>matlab 
>>run plot_news20scale2_fl




======================
To reproduce Figure 4(b) in our paper, do:

>> cd PUT_PATH_TO_ROOT_WHERE_THIS_README_FILE_IS
>> cd IPALM
>> g++ -o main main.cpp -lgsl -lgslcblas
>> ./main 3 a rcv1 0.01 0.01 100000 600 1 100 100 ; ./main 3 b rcv1 0.01 0.01 100000 600 1 1 40 ; ./main 3 c rcv1 0.01 0.01 100000 2400 10 50 ; ./main 3 d rcv1 0.01 0.01 100000 2400 10 ; ./main 3 e rcv1 0.01 0.01 1000000 2400 1 100


This should take roughly 600+600+2400+2400+2400=8400 seconds. The results are saved in the folder results/. 

Then run the following commands to generate the plot from the results.

>>cd ..
>>cd MATLAB_code
>>matlab 
>>run plot_rcv1_fl






======================
To reproduce Figure 4(c) in our paper, do:

>> cd PUT_PATH_TO_ROOT_WHERE_THIS_README_FILE_IS
>> cd IPALM
>> g++ -o main main.cpp -lgsl -lgslcblas
>> ./main 3 a rcv1mc 0.01 0.01 100000 1500 1 100 100 ; ./main 3 b rcv1mc 0.01 0.01 100000 600 1 1 40 ; ./main 3 c rcv1mc 0.01 0.01 100000 3600 10 50 ; ./main 3 d rcv1mc 0.01 0.01 100000 3600 10 ; ./main 3 e rcv1mc 0.01 0.01 1000000 3600 1 100


This should take roughly 1500+600+3600+3600+3600=12900 seconds. The results are saved in the folder results/. 

Then run the following commands to generate the plot from the results.

>>cd ..
>>cd MATLAB_code
>>matlab 
>>run plot_rcv1mc_fl





======================
To reproduce Figure 4(d) in our paper, do:

>> cd PUT_PATH_TO_ROOT_WHERE_THIS_README_FILE_IS
>> cd IPALM
>> g++ -o main main.cpp -lgsl -lgslcblas
>> ./main 3 a news20binary 0.01 0.01 100000 3600 1 100 100 ; ./main 3 b news20binary 0.01 0.01 100000 3600 1 10 10 ; ./main 3 c news20binary 0.01 0.01 100000 3600 10 10 ; ./main 3 d news20binary 0.01 0.01 100000 3600 0.1 ; ./main 3 e news20binary 0.01 0.01 1000000 3600 0.01 100


This should take roughly 3600+3600+3600+3600+3600=18000 seconds. The results are saved in the folder results/. 

Then run the following commands to generate the plot from the results.

>>cd ..
>>cd MATLAB_code
>>matlab 
>>run plot_news20binary_fl






======================
To reproduce Figure 5(a) in our paper, do:

>> cd PUT_PATH_TO_ROOT_WHERE_THIS_README_FILE_IS
>> cd IPALM
>> g++ -o main main.cpp -lgsl -lgslcblas
>> ./main 4 a w4a 1 1 100000 200 1 1000 100 ; ./main 4 b w4a 1 1 100000 200 10 100 20 ; ./main 4 c w4a 1 1 100000 200 100 100 ; ./main 4 d w4a 1 1 100000 200 1 ; ./main 4 e w4a 1 1 1000000 200 0.1 100


This should take roughly 200+200+200+200+200=1000 seconds. The results are saved in the folder results/. 

Then run the following commands to generate the plot from the results.

>>cd ..
>>cd MATLAB_code
>>matlab 
>>run plot_w4a_svm





======================
To reproduce Figure 5(b) in our paper, do:

>> cd PUT_PATH_TO_ROOT_WHERE_THIS_README_FILE_IS
>> cd IPALM
>> g++ -o main main.cpp -lgsl -lgslcblas
>> ./main 4 a a7a 1 1 100000 400 1 1000 100 ; ./main 4 b a7a 1 1 100000 400 1 100 20 ; ./main 4 c a7a 1 1 100000 1200 100 100 ;./main 4 d a7a 1 1 100000 1200 10 ; ./main 4 e a7a 1 1 1000000 1200 0.1 100


This should take roughly 400+400+1200+1200+1200=4400 seconds. The results are saved in the folder results/. 

Then run the following commands to generate the plot from the results.

>>cd ..
>>cd MATLAB_code
>>matlab 
>>run plot_a7a_svm






======================
To reproduce Figure 5(c) in our paper, do:

>> cd PUT_PATH_TO_ROOT_WHERE_THIS_README_FILE_IS
>> cd IPALM
>> g++ -o main main.cpp -lgsl -lgslcblas
>> ./main 4 a a8a 1 1 100000 400 1 1000 100 ; ./main 4 b a8a 1 1 100000 400 1 100 20 ; ./main 4 c a8a 1 1 100000 1200 100 100 ; ./main 4 d a8a 1 1 100000 1200 10 ; ./main 4 e a8a 1 1 1000000 1200 0.1 100


This should take roughly 400+400+1200+1200=1200=4400 seconds. The results are saved in the folder results/. 

Then run the following commands to generate the plot from the results.

>>cd ..
>>cd MATLAB_code
>>matlab 
>>run plot_a8a_svm





======================
To reproduce Figure 5(d) in our paper, do:

>> cd PUT_PATH_TO_ROOT_WHERE_THIS_README_FILE_IS
>> cd IPALM
>> g++ -o main main.cpp -lgsl -lgslcblas
>> ./main 4 a a9a 1 1 100000 400 1 1000 100 ; ./main 4 b a9a 1 1 100000 400 1 100 20 ; ./main 4 c a9a 1 1 100000 1200 100 100 ; ./main 4 d a9a 1 1 100000 1200 10 ; ./main 4 e a9a 1 1 1000000 1200 0.01 100


This should take roughly 400+400+1200+1200+1200=4400 seconds. The results are saved in the folder results/. 

Then run the following commands to generate the plot from the results.

>>cd ..
>>cd MATLAB_code
>>matlab 
>>run plot_a9a_svm






======================
To reproduce Figure 5(e) in our paper, do:

>> cd PUT_PATH_TO_ROOT_WHERE_THIS_README_FILE_IS
>> cd IPALM
>> g++ -o main main.cpp -lgsl -lgslcblas
>> ./main 4 a w6a 1 1 100000 400 1 1000 100 ; ./main 4 b w6a 1 1 100000 400 1 100 20 ; ./main 4 c w6a 1 1 100000 1200 100 100 ; ./main 4 d w6a 1 1 100000 1200 10 ; ./main 4 e w6a 1 1 1000000 1200 0.1 100


This should take roughly 400+400+1200+1200+1200=4400 seconds. The results are saved in the folder results/. 

Then run the following commands to generate the plot from the results.

>>cd ..
>>cd MATLAB_code
>>matlab 
>>run plot_w6a_svm





======================
To reproduce Figure 5(f) in our paper, do:

>> cd PUT_PATH_TO_ROOT_WHERE_THIS_README_FILE_IS
>> cd IPALM
>> g++ -o main main.cpp -lgsl -lgslcblas
>> ./main 4 a w7a 1 1 100000 400 1 1000 100 ; ./main 4 b w7a 1 1 100000 400 1 100 20 ; ./main 4 c w7a 1 1 100000 1200 100 100 ; ./main 4 d w7a 1 1 100000 1200 10 ; ./main 4 e w7a 1 1 1000000 1200 0.1 100


This should take roughly 400+400+1200+1200+1200=4400 seconds. The results are saved in the folder results/. 

Then run the following commands to generate the plot from the results.

>>cd ..
>>cd MATLAB_code
>>matlab 
>>run plot_w7a_svm






======================
To reproduce Figure 5(g) in our paper, do:

>> cd PUT_PATH_TO_ROOT_WHERE_THIS_README_FILE_IS
>> cd IPALM
>> g++ -o main main.cpp -lgsl -lgslcblas
>> ./main 4 a ijcnn1 1 1 100000 400 1 1000 100 ; ./main 4 b ijcnn1 1 1 100000 400 1 100 20 ; ./main 4 c ijcnn1 1 1 100000 1200 100 100 ; ./main 4 d ijcnn1 1 1 100000 1200 10 ; ./main 4 e ijcnn1 1 1 1000000 1200 1 100


This should take roughly 400+400+1200+1200+1200=4400 seconds. The results are saved in the folder results/. 

Then run the following commands to generate the plot from the results.

>>cd ..
>>cd MATLAB_code
>>matlab 
>>run plot_ijcnn_svm





======================
To reproduce Figure 5(h) in our paper, do:

>> cd PUT_PATH_TO_ROOT_WHERE_THIS_README_FILE_IS
>> cd IPALM
>> g++ -o main main.cpp -lgsl -lgslcblas
>> ./main 4 a realsim 1 1 100000 2000 1 100 100 ; ./main 4 b realsim 1 1 100000 2000 1 100 20 ; ./main 4 c realsim 1 1 100000 3600 1 100 ; ./main 4 d realsim 1 1 100000 3600 1 ; ./main 4 e realsim 1 1 1000000 3600 0.01 100


This should take roughly 2000+2000+3600+3600+3600=14800 seconds. The results are saved in the folder results/. 

Then run the following commands to generate the plot from the results.

>>cd ..
>>cd MATLAB_code
>>matlab 
>>run plot_realsim_svm





======================
To reproduce Figure 5(i) in our paper, do:

>> cd PUT_PATH_TO_ROOT_WHERE_THIS_README_FILE_IS
>> cd IPALM
>> g++ -o main main.cpp -lgsl -lgslcblas
>> ./main 4 a w8a 1 1 100000 400 1 1000 100 ; ./main 4 b w8a 1 1 100000 400 1 100 20 ; ./main 4 c w8a 1 1 100000 1200 100 100 ; ./main 4 d w8a 1 1 100000 1200 10 ; ./main 4 e w8a 1 1 1000000 1200 0.1 100


This should take roughly 400+400+1200+1200+1200=4400 seconds. The results are saved in the folder results/. 

Then run the following commands to generate the plot from the results.

>>cd ..
>>cd MATLAB_code
>>matlab
>>run plot_w8a_svm




======================
To reproduce Figure 5(j) in our paper, do:

>> cd PUT_PATH_TO_ROOT_WHERE_THIS_README_FILE_IS
>> cd IPALM
>> g++ -o main main.cpp -lgsl -lgslcblas
>> ./main 4 a covtype 1 1 100000 1200 1 1000 100 ; ./main 4 b covtype 1 1 100000 1200 1 100 20 ; ./main 4 c covtype 1 1 100000 1200 1 100 ; ./main 4 d covtype 1 1 100000 1200 10 ; ./main 4 e covtype 1 1 1000000 1200 0.01 100


This should take roughly 1200+1200+1200+1200+1200=6000 seconds. The results are saved in the folder results/. 

Then run the following commands to generate the plot from the results.

>>cd ..
>>cd MATLAB_code
>>matlab 
>>run plot_covtype_svm


All the plots are saved in .eps format in the folder Matlab_code/myplots/

======================
To reproduce Figure 6 in our paper, do:

>> cd PUT_PATH_TO_ROOT_WHERE_THIS_README_FILE_IS
>> cd IPALM
>> g++ -o main main.cpp -lgsl -lgslcblas
>> ./main 5 b qcqp1 1 1 200 10 100 1 50 ; ./main 5 b qcqp2 1 1 1000 20 100 1 50 ; ./main 5 b qcqp3 1 1 200 30 10 1 50 ; ./main 5 b qcqp4 1 1 190 500 1 1000 50

This should take roughly 10+20+30+500=560 seconds. The results are saved in the folder results/. In the output file, the second column is the running time; the third column is the objective value; and the forth column is the infeasibility. 

Then run the following commands to generate the plot from the realists.

>>cd ..
>>cd MATLAB_code
>>Matlab
>>run plot_qcqp1
>>run plot_qcqp2
>>run plot_qcqp3
>>run plot_qcqp4
