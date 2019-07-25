function out = calc_mag_ls(samp,phins)

Rmp = rph2rot([0;0;pi]);
Rmp = rph2rot([0;pi/2;-pi/2]);

%samp.acc = (Rmp*samp.acc')';
%samp.mag = (Rmp*samp.mag')';
out.samp = samp.mag;

mag = samp.mag;

num = size(mag,1);

W = zeros(num,10);
W2 = zeros(num,6);

for i=1:num
    
W(i,:) = [mag(i,1)^2,mag(i,2)^2,mag(i,3)^2,2*mag(i,1)*mag(i,2),2*mag(i,1)*mag(i,3),2*mag(i,2)*mag(i,3),-2*mag(i,:),1];
W2(i,:) = [mag(i,1)^2,mag(i,2)^2,2*mag(i,1)*mag(i,2),-2*mag(i,1:2),1];

end

out.W = W;
out.W2 = W2;

[out.V,out.D] = eig(out.W'*out.W);
[out.V2,out.D2] = eig(out.W2'*out.W2);

sign = -1;
%sign =1;

out.V2 = out.V2*sign;
out.V = out.V*sign;
out.T = inv(sqrtm([out.V(1,1),out.V(4,1),out.V(5,1);out.V(4,1),out.V(2,1),out.V(6,1);out.V(5,1),out.V(6,1),out.V(3,1)]));
out.b = [out.V(1,1),out.V(4,1),out.V(5,1);out.V(4,1),out.V(2,1),out.V(6,1);out.V(5,1),out.V(6,1),out.V(3,1)]\out.V(7:9,1);

out.T2 = inv(sqrtm([out.V2(1,1),out.V2(3,1),0;out.V2(3,1),out.V2(2,1),0;0,0,1]));
out.T2 = out.T2/out.T2(1,1);
out.b2 = [out.V2(1,1),out.V2(3,1),0;out.V2(3,1),out.V2(2,1),0;0,0,1]\[out.V2(4:5,1);0];

out.mag2 = real(out.T2\(mag'-out.b2));
out.mag = real(out.T\(mag'-out.b));

out.t = samp.t;

out.acc = samp.acc';

% att = calc_heading(samp.t,out,phins);
% out.att = att.att;
% out.att_error = att.att_error;
% out.mag_ned = att.mag_ned;
% figure;plot(phins.t,out.att_error);grid on;
