function out = sensor_report(samp,out_file)


samp.mag = samp.mag(10:end,:);
h = histogram(samp.mag(:,1));
num_bins = h.NumBins;


sdev = std(samp.mag);
mn = mean(samp.mag);
s = sprintf('std: %f, %f ,%f',sdev(1),sdev(2),sdev(3));
m = sprintf('mean: %f, %f, %f',mn(1),mn(2),mn(3));


figure(1);
subplot(3,1,1);
h = histogram(samp.mag(:,1),num_bins);
title('mag');
xlabel(s);
subplot(3,1,2);
histogram(samp.mag(:,2),num_bins);
xlabel(m);
subplot(3,1,3);
histogram(samp.mag(:,3),num_bins);
xlabel('gaus');

saveas(figure(1), 'fig1.pdf');



sdev = std(samp.acc);
mn = mean(samp.acc);
s = sprintf('std: %f, %f ,%f',sdev(1),sdev(2),sdev(3));
m = sprintf('mean: %f, %f, %f',mn(1),mn(2),mn(3));

figure(2);
subplot(3,1,1);
histogram(samp.acc(:,1));
title('acc');
xlabel(s);
subplot(3,1,2);
histogram(samp.acc(:,2));
xlabel(m);
subplot(3,1,3);
h = histogram(samp.acc(:,3));
xlabel('g');

saveas(figure(2), 'fig2.pdf');


sdev = std(samp.ang);
mn = mean(samp.ang);
s = sprintf('std: %f, %f ,%f',sdev(1),sdev(2),sdev(3));
m = sprintf('mean: %f, %f, %f',mn(1),mn(2),mn(3));


figure(3);
subplot(3,1,1);
histogram(samp.ang(:,1));
title('ang');
xlabel(s);
subplot(3,1,2);
histogram(samp.ang(:,2));
xlabel(m);
subplot(3,1,3);
histogram(samp.ang(:,3));
xlabel('rad/s');

saveas(figure(3), 'fig3.pdf');
close all;

system(['pdftk fig*.pdf cat output ' out_file]);
system('rm fig1.pdf fig2.pdf fig3.pdf');
