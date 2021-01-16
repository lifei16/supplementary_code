n= 10;
% the number of the constraints, for qcqp2, m= 500; for qcqp3, m= 1000; for qcqp4, m= 10000.
m= 100;
A= sparse(n*(m+ 1),n);
b= zeros(n*(m+ 1),1);
for i= 1:m+1
    R = sprandsym(n,0.5,0.1);
    R = R'*R;
    c= randn(n,1);
    A((i- 1)*n+1:i*n,:)= R;
    b((i- 1)*n+1:i*n,:)= c;
end

% For n= 10, we can directly save.
save qcqp1_A.mat A;
save qcqp1_b.mat b;


% To store in libsvm format for later use of IPALM_Katyusha, we need use libsvmwrite.m in LIBSVM package
%libsvmwrite('matrix_qcqp4',b,A);
