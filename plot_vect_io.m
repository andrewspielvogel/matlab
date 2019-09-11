function out = plot_vect_io(out,phins)

err_adap = wrapToPi(unwrap(out.heading_adap) - unwrap(phins.att(:,3)*pi/180));
err_uncal = (unwrap(out.heading_uncal) - unwrap(phins.att(:,3)*pi/180));

T=[1.1,0.2;0.2,0.95];
Tinv = inv(T);
Tinv2 = Tinv*Tinv;
b = [0.06;-0.07];
bTinv2b = b'*Tinv2*b;
m_norm = norm([0.205796;-0.040654]);
phi = m_norm^2 -bTinv2b;

tp = [Tinv2(1,1);Tinv2(2,2);Tinv2(1,2)]/phi;
 phinst = phins.t;
thresh = .25;
err_indx = find(abs(err_adap- mean(err_adap)) > thresh);
phins.t(err_indx) = [];
err_adap(err_indx) = [];

err_indx = find(abs(err_uncal- mean(err_uncal)) > thresh);
phinst(err_indx) = [];
err_uncal(err_indx) = [];



rmse_adap = rms(err_adap - mean(err_adap))*180/pi;
rmse_uncal = rms(err_uncal - mean(err_adap))*180/pi;

tf = 400;
figure;
subplot(2,1,1);
plot(taxis(phins.t(1:tf)-phins.t(1)),unwrap(phins.att(1:tf,3)*pi/180)*180/pi);xlabel(tlabel(phins.t(1:tf)-phins.t(1)));ylabel('Degrees');title('True Instrument Heading');grid on;xlim(taxis([0,phins.t(tf)-phins.t(1)]));ylim([-330;330]);
subplot(2,1,2);
plot(taxis(phinst(1:tf)-phins.t(1)),err_uncal(1:tf)*180/pi-mean(err_adap)*180/pi,'.',taxis(phins.t(1:tf)-phins.t(1)),err_adap(1:tf)*180/pi-mean(err_adap)*180/pi,'.');grid on;xlabel(tlabel(phins.t(1:tf)-phins.t(1)));ylabel('Degrees');title('Heading Error');legend(sprintf('uncal: RMSE= %.2f Deg',rmse_uncal),sprintf('cal: RMSE= %.2f Deg',rmse_adap));xlim(taxis([0,phins.t(tf)-phins.t(1)]));

out.t = out.t(1:10:end);
out.tp = out.tp(:,1:10:end);
out.b = out.b(:,1:10:end);
out.theta = out.theta(:,1:10:end);

out.samp.mag = out.samp.mag(1:20:end,:);
out.mag_adap = out.mag_adap(:,1:20:end);

figure;plot(taxis(out.t-out.t(1)),out.theta);grid on;xlabel(tlabel(out.t-out.t(1)));ylabel('$\hat{\theta}$','Interpreter','latex');title('$\hat{\theta}$','Interpreter','latex');
figure;plot(taxis(out.t-out.t(1)),out.tp-tp);grid on;xlabel(tlabel(out.t-out.t(1)));title('$\hat{t}_p(t)-t_p$','Interpreter','latex');ylabel('$\hat{t}_p-t_p$','Interpreter','latex');legend('a','b','c');
figure;plot(taxis(out.t-out.t(1)),out.b-b);grid on;xlabel(tlabel(out.t-out.t(1)));title('$\hat{b}(t) - b$','Interpreter','latex');ylabel('$\hat{b}-b$ (G)','Interpreter','latex');legend('b_x','b_x');
figure;plot(taxis(out.t-out.t(1)),out.tp);grid on;xlabel(tlabel(out.t-out.t(1)));title('$\hat{t}_p(t)$','Interpreter','latex');ylabel('$\hat{t}_p$','Interpreter','latex');legend('a','b','c');xlim(taxis([0,out.t(end)-out.t(1)]));
figure;plot(taxis(out.t-out.t(1)),out.b);grid on;xlabel(tlabel(out.t-out.t(1)));title('$\hat{b}(t)$','Interpreter','latex');ylabel('$\hat{b}$ (G)','Interpreter','latex');legend('b_x','b_x');xlim(taxis([0,out.t(end)-out.t(1)]));



%figure;plot(out.mag_adap(1,:),out.mag_adap(2,:));grid on; grid minor; axis equal;title('mag_adap');
%figure;plot(out.mag_ls(1,:),out.mag_ls(2,:));grid on; grid minor; axis equal;title('mag_ls');
%figure;plot(out.samp.mag(:,1),out.samp.mag(:,2));grid on;axis equal;grid minor;title('mag_uncal');
figure;plot(out.samp.mag(:,1),out.samp.mag(:,2),'.',out.mag_adap(1,:),out.mag_adap(2,:),'.');grid on;axis equal;legend('uncal mag','cal mag');xlabel('m_x (G)');ylabel('m_y (G)');title('Calibrated Versus Uncalibrated Magnetometer');

