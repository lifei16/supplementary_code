f_opt2=-1.092297236140837;

t2=11.90;

qcqp2= importdata('../IPALM/results/ALM_QCQP_outer_qcqp2_tau_22',' ');

s2= size(qcqp2,1);
plot(qcqp2(1:4:s2,2),log10(max( qcqp2(1:4:s2,4), (abs(qcqp2(1:4:s2,3)- f_opt2)/f_opt2))),'-->','LineWidth',2);

hold on;
plot(t2,-8,'r-pentagram','markersize',20);
ylim([-8 2]);
xlim([0 t2]);
xlabel('time');
ylabel('log(E(x))');
title('QCQP n=10 m=500');
set(gcf,'Position',[10 10 400 400]);
saveas(gcf,'myplots/qcqp2.eps','epsc');





