

%mag_bias = calc_mag_bias2(mst,[1;1;0;0;0;0]);
%mag_bias3 = calc_mag_bias(mst,[1;1;1;0;0;0;0;0;0;0]);
%ang_acc_bias = calc_acc_ang_bias(mst,[mst.acc(1,:)';zeros(6,1)]);
%mag_bias_ls = calc_mag_ls(mst.mag);

mag_cal = mag_bias.T*(mst.mag'-mag_bias.b);
mag_cal3 = mag_bias3.T*(mst.mag'-mag_bias3.b);
mag_cal_ls = mag_bias_ls.T2*(mst.mag'-mag_bias_ls.b2);

rph_start = (skew(mst.mag(1,:))*[1;0;0]*180/pi).*[0;0;1];

%att_uncal = calc_att(K,mst.acc',mst.mag',mst.ang',mst.t,rph2R(rph_start*pi/180));
%att_cal = calc_att(K,mst.acc',mag_cal,mst.ang',mst.t,rph2R(rph_start*pi/180));
%att_cal_adap = calc_att(K,mst.acc',mag_bias3.mag,mst.ang',mst.t,rph2R(rph_start*pi/180));
%att_cal3 = calc_att(K,mst.acc',mag_cal3,mst.ang',mst.t,rph2R(rph_start*pi/180));
%att_cal_ls = calc_att(K,mst.acc',mag_cal_ls,mst.ang',mst.t,rph2R(rph_start*pi/180));

att_resamp_uncal = resample2(mst.t,att_uncal.att',phins.t);
att_resamp_cal = resample2(mst.t,att_cal.att',phins.t);
att_resamp_cal_adap = resample2(mst.t,att_cal_adap.att',phins.t);
att_resamp_cal3 = resample2(mst.t,att_cal3.att',phins.t);
att_resamp_cal_ls = resample2(mst.t,att_cal_ls.att',phins.t);

uncal_error = wrapTo180(unwrap(att_resamp_uncal)*180/pi-unwrap(phins.att*pi/180)*180/pi);
cal_error = wrapTo180(unwrap(att_resamp_cal)*180/pi-unwrap(phins.att*pi/180)*180/pi);
cal_error_adap = wrapTo180(unwrap(att_resamp_cal_adap)*180/pi-unwrap(phins.att*pi/180)*180/pi);
%cal_error3 = wrapTo180(unwrap(att_resamp_cal3)*180/pi-unwrap(phins.att*pi/180)*180/pi);
ls_error = wrapTo180(unwrap(att_resamp_cal_ls)*180/pi-unwrap(phins.att*pi/180)*180/pi);

uncal_rms = sqrt(uncal_error(:,3).*uncal_error(:,3));
cal_rms = sqrt(cal_error(:,3).*cal_error(:,3));
ls_rms = sqrt(ls_error(:,3).*ls_error(:,3));

t = taxis(phins.t);
%figure;plot(t,uncal_error(:,3),t,cal_error(:,3),t,ls_error(:,3),t,cal_error3(:,3));grid on;legend('uncal heading','cal heading','ls heading','cal heading3');xlabel(tlabel(phins.t));ylabel('degrees');ylim([-200,200]);xlim([t(1),t(end)]);
figure;plot(t,uncal_error(:,3),t,cal_error(:,3),t,ls_error(:,3),t,cal_error3(:,3),t,cal_error_adap(:,3));grid on;legend('uncal heading','cal heading','ls heading','cal heading new gains','adap mag new gains');xlabel(tlabel(phins.t));ylabel('degrees');ylim([-200,200]);xlim([t(1),t(end)]);title('Heading Error');
figure;
subplot(3,1,1);
plot(taxis(mst.t),mst.mag(:,1));grid on;xlim([t(1),t(end)]);ylabel('gauss');title('X-Component: Magnetometer');
subplot(3,1,2);
plot(taxis(mst.t),mst.mag(:,2));grid on;xlim([t(1),t(end)]);ylabel('gauss');title('Y-Component: Magnetometer');
subplot(3,1,3);
plot(taxis(mst.t),mst.mag(:,3));grid on;xlim([t(1),t(end)]);ylabel('gauss');title('Z-Component: Magnetometer');xlabel(tlabel(phins.t));

figure;
subplot(3,1,1);
plot(taxis(mst.t),mst.ang(:,1));grid on;xlim([t(1),t(end)]);ylabel('rad/s');title('X-Component: Angular-Rate Gyro');
subplot(3,1,2);
plot(taxis(mst.t),mst.ang(:,2));grid on;xlim([t(1),t(end)]);ylabel('rad/s');title('Y-Component: Angular-Rate Gyro');
subplot(3,1,3);
plot(taxis(mst.t),mst.ang(:,3));grid on;xlim([t(1),t(end)]);ylabel('rad/s');title('Z-Component: Angular-Rate Gyro');xlabel(tlabel(phins.t));

figure;
subplot(3,1,1);
plot(taxis(mst.t),mst.acc(:,1));grid on;xlim([t(1),t(end)]);ylabel('m/s^2');title('X-Component: Linear Accelerometer');
subplot(3,1,2);
plot(taxis(mst.t),mst.acc(:,2));grid on;xlim([t(1),t(end)]);ylabel('m/s^2');title('Y-Component: Linear Accelerometer');
subplot(3,1,3);
plot(taxis(mst.t),mst.acc(:,3));grid on;xlim([t(1),t(end)]);ylabel('m/s^2');title('Z-Component: Linear Accelerometer');xlabel(tlabel(phins.t));

phinst = taxis(phins.t);
figure;
subplot(3,1,1);
plot(phinst,phins.att(:,1));grid on;xlim([phinst(1),phinst(end)]);ylabel('deg');title('X-Component: Attitude');
subplot(3,1,2);
plot(phinst,phins.att(:,2));grid on;xlim([phinst(1),phinst(end)]);ylabel('deg');title('Y-Component: Attitude');
subplot(3,1,3);
plot(phinst,phins.att(:,3));grid on;xlim([phinst(1),phinst(end)]);ylabel('deg');title('Z-Component: Attitude');xlabel(tlabel(phinst));
