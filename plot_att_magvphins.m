function out = plot_att_magvphins(att,phins)

phins.att = phins.att*pi/180;
att.att = att.att';

phins.t = phins.t-0.04;

%temp = phins;
%phins = att;
%att = temp;

phins_att = resample2(phins.t,unwrap(phins.att),att.t,'linear');

figure;plot(att.t,att.att*180/pi,phins.t,phins.att*180/pi);grid on;xlabel('Seconds');ylabel('Degrees')

figure;plot(att.t,unwrap(att.att)*180/pi,att.t,unwrap(phins_att)*180/pi);grid on;xlabel('Seconds');ylabel('Degrees')

figure;plot(att.t-att.t(1),unwrap(unwrap(att.att)-unwrap(phins_att))*180/pi);grid on;xlabel('Seconds');ylabel('Error (degrees)');