function out = calc_heading(t,mag_samp,phins)

    num_samps = size(mag_samp,2);
    heading = zeros(num_samps,1);
    
    for i=1:num_samps
       
        mag = (eye(3)-[0,0,0;0,0,0;0,0,1])*mag_samp(:,i);
        mag = mag/norm(mag);
        heading(i) = atan2(mag(1),mag(2));

    end
    
    out.heading = resample2(t,heading,phins.t);
    out.heading_error = wrapTo180(unwrap(out.heading)*180/pi+90-unwrap(phins.att(:,3)*pi/180)*180/pi);
