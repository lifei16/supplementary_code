Instruction for LIBSVM

This folder contains the code and data used to solve support vector machine problems by LIBSVM. Readers need download and install LIBSVM package (https://www.csie.ntu.edu.tw/~cjlin/libsvm/#download).

To generate Table 7 in our paper, 
First, do
>> cd ..
>> cd IPALM
>> ./main 4 b a7a 1 1 100000 10 1 100 20; ./main 4 b a8a 1 1 100000 19 1 100 20; ./main 4 b a9a 1 1 100000 41 1 100 20; ./main 4 b w6a 1 1 100000 2 1 100 20; ./main 4 b w7a 1 1 100000 3 1 100 20; ./main 4 b ijcnn1 1 1 100000 31 1 100 20; ./main 4 b w8a 1 1 100000 13 1 100 20
copy the output x from the folder "../IPALM/results/" to current folder:
>> cp ../IPALM/results/L_Katyusha_x_SVM_outer_a7atau_126 L_Katyusha_x_SVM_outer_a7atau_126
>> cp ../IPALM/results/L_Katyusha_x_SVM_outer_a8atau_150 L_Katyusha_x_SVM_outer_a8atau_150
>> cp ../IPALM/results/L_Katyusha_x_SVM_outer_a9atau_180 L_Katyusha_x_SVM_outer_a9atau_180
>> cp ../IPALM/results/L_Katyusha_x_SVM_outer_ijcnn1tau_223 L_Katyusha_x_SVM_outer_ijcnn1tau_223
>> cp ../IPALM/results/L_Katyusha_x_SVM_outer_w6atau_131 L_Katyusha_x_SVM_outer_w6atau_131
>> cp ../IPALM/results/L_Katyusha_x_SVM_outer_w7atau_157 L_Katyusha_x_SVM_outer_w7atau_157
>> cp ../IPALM/results/L_Katyusha_x_SVM_outer_w8atau_223 L_Katyusha_x_SVM_outer_w8atau_223

Second, 
>> matlab
>> addpath PUT_PATH_TO_ROOT_WHERE_LIBSVM_PACKAGE_FILE_IS
>> run svm_compare.m
For different datasets, just change the dataset name to define [b, A], [b2, A2] in the 1st and 2nd line, and x in the 18th line.

In the output variables, "acc" is the accuracy of LIBSVM and "accuracy" is the accuracy of IPALM-Katyuhsa.
