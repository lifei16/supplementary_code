f_opt4=-1.250143629451023;

t4=580.84;

qcqp4= importdata('../IPALM/results/ALM_QCQP_outer_qcqp4_tau_100',' ');

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





