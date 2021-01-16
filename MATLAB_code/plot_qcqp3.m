f_opt1= -1.968035633575473;
f_opt2=-1.092297236140837;
f_opt3= -0.995184307600784;
f_opt4=-1.250143629451023;

t1=2.54;
t2=11.90;

t3=22.14;
t4=580.84;
qcqp1= importdata('../IPALM/results/ALM_QCQP_outer_qcqp1_tau_10',' ');
qcqp2= importdata('../IPALM/results/ALM_QCQP_outer_qcqp2_tau_22',' ');
qcqp3= importdata('../IPALM/results/ALM_QCQP_outer_qcqp3_tau_31',' ');
qcqp4= importdata('../IPALM/results/ALM_QCQP_outer_qcqp4_tau_100',' ');
s1= size(qcqp1,1);
s2= size(qcqp2,1);
s3= size(qcqp3,1);
s4= size(qcqp4,1);
 plot(qcqp3(1:3:s3,2),log10(max( qcqp3(1:3:s3,4), (abs(qcqp3(1:3:s3,3)- f_opt3)/f_opt3))),'-->','LineWidth',2);

hold on;
plot(t3,-8,'r-pentagram','markersize',20);
ylim([-8 2]);
xlim([0 t3]);
xlabel('time');
ylabel('log(E(x))');
title('QCQP n=10 m=1000');
set(gcf,'Position',[10 10 400 400]);
saveas(gcf,'myplots/qcqp3.eps','epsc');





