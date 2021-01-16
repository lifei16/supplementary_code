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
 plot(qcqp4(1:3:s4,2),log10(max( qcqp4(1:3:s4,4), (abs(qcqp4(1:3:s4,3)- f_opt4)/f_opt4))),'-->','LineWidth',2);

hold on;
plot(t4,-8,'r-pentagram','markersize',20);
ylim([-8 2]);
xlim([0 t4]);
xlabel('time');
ylabel('log(E(x))');
title('QCQP n=10 m=10000');
set(gcf,'Position',[10 10 400 400]);
saveas(gcf,'myplots/qcqp4.eps','epsc');





