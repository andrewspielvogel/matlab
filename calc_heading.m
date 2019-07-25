function out = calc_heading(t,samp,phins)

    %samp.mag = rph2R([0;0;pi])*samp.mag;
    R = rph2rot([pi/2;0;-pi/2]);
    R = rph2rot([pi;0;pi]);


    samp.acc = R*samp.acc;
    samp.mag = R*samp.mag;
    
    %samp.acc = my_lowpass(samp.acc',1000,1,10)';
    %samp.mag = my_lowpass(samp.mag',1000,1,10)';

    roll  = atan2(-samp.acc(2,:),-samp.acc(3,:));
    pitch = atan2(samp.acc(1,:),sqrt(samp.acc(2,:).*samp.acc(2,:) + samp.acc(3,:).*samp.acc(3,:)));
    
    mag = zeros(size(samp.mag));
    
    %phins_att = resample2(phins.t,phins.att,t)';
    %roll = phins_att(1,:)*pi/180;
    %pitch = phins_att(2,:)*pi/180;
    
    for i=1:size(samp.t,1)
       
        mag(:,i) =  rph2R([roll(i),pitch(i),0])*samp.mag(:,i);

    end

        
    heading = atan2(-mag(2,:),mag(1,:));

    att = [roll',pitch',heading'];

    att = resample2(t,att,phins.t)';

    out.att = unwrap(att')';
    out.att_error = wrapTo180(att'*180/pi-unwrap(phins.att*pi/180)*180/pi)';
    out.mag_ned = mag;
    figure;plot(out.mag_ned(1,:),out.mag_ned(2,:),'.');grid on;axis equal;
    figure;plot(taxis(phins.t),out.att_error);grid on;xlabel(tlabel(phins.t));
    figure;plot(taxis(phins.t),unwrap(phins.att*pi/180)*180/pi,taxis(phins.t),out.att*180/pi);grid on;legend('px','py','pz','x','y','z');xlabel(tlabel(phins.t));
    figure;plot(taxis(samp.t),samp.b);grid on;xlabel(tlabel(samp.t));legend('x','y','z');
    figure;plot(taxis(samp.t),samp.theta);grid on;legend('a','b','c','d','e','f','a1','a2','a3','btb');xlabel(tlabel(samp.t));

        
        