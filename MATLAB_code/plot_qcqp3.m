f_opt3= -0.995184307600784;

t3=22.14;

qcqp3= importdata('../IPALM/results/ALM_QCQP_outer_qcqp3_tau_31',' ');

s3= size(qcqp3,1);
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





