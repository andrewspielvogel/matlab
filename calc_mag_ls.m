function out = calc_mag_ls(samp,phins)

Rmp = rph2rot([0;0;pi]);
Rmp = rph2rot([0;pi/2;-pi/2]);


%samp.acc = (Rmp*samp.acc')';
%samp.mag = (Rmp*samp.mag')';

mag = samp.mag;
dmag = samp.dmag;
ang = samp.ang;


num = size(mag,1);

%W = zeros(num*3,9);
W = zeros(num*2,10);
W_ls = zeros(num,9);

for i=1:num
    
Y = [dmag(i,1),0,0,dmag(i,2),dmag(i,3),0;0,dmag(i,2),0,dmag(i,1),0,dmag(i,3);0,0,dmag(i,3),0,dmag(i,1),dmag(i,2)];
X = [0,-mag(i,2)*ang(i,3),mag(i,3)*ang(i,2),-mag(i,1)*ang(i,3),mag(i,1)*ang(i,2),mag(i,2)*ang(i,2)-mag(i,3)*ang(i,3);
    mag(i,1)*ang(i,3),0,-mag(i,3)*ang(i,1),mag(i,2)*ang(i,3),mag(i,3)*ang(i,3)-mag(i,1)*ang(i,1),-mag(i,2)*ang(i,1);
    -mag(i,1)*ang(i,2),mag(i,2)*ang(i,1),0,mag(i,1)*ang(i,1)-mag(i,2)*ang(i,2),-mag(i,3)*ang(i,2),mag(i,3)*ang(i,1)];

%W(i*3-2:i*3,:) = [X+Y,skew(ang(i,:))];
W(i*2-1:i*2,:) = [mag(i,1)^2,mag(i,2)^2,mag(i,3)^2,2*mag(i,1)*mag(i,2),2*mag(i,1)*mag(i,3),2*mag(i,2)*mag(i,3),-2*mag(i,:),1;
                  0,0,1,zeros(1,7)];
W_ls(i,:) = [mag(i,1)^2,mag(i,2)^2,mag(i,3)^2,2*mag(i,1)*mag(i,2),2*mag(i,1)*mag(i,3),2*mag(i,2)*mag(i,3),-2*mag(i,:)];

end
out.W = W;
out.W_ls = W_ls;
[V,D] = eig(out.W'*out.W);
    

K = ones(num,1);

out.theta = (out.W_ls'*out.W_ls)\out.W_ls'*K;

%K =  ones(num*2,1);
%K(1:2:end) = 0;
%out.theta = (out.W'*out.W)\out.W'*K;

Tinv = [out.theta(1),out.theta(4),out.theta(5);out.theta(4),out.theta(2),out.theta(6);out.theta(5),out.theta(6),out.theta(3)];
out.T = inv(sqrtm(Tinv));
out.Tinv=Tinv;
out.b = Tinv\out.theta(7:9);
i = sum(out.theta(1:3));
j = out.theta(1)*out.theta(2)+out.theta(1)*out.theta(3)+out.theta(2)*out.theta(3)-out.theta(4)*out.theta(4)-out.theta(5)*out.theta(5)-out.theta(6)*out.theta(6);
out.i = i;
out.j = j;

out.T = out.T;
[X,Y,Z] = sphere(50);
Cir = [reshape(X,1,51*51);reshape(Y,1,51*51);reshape(Z,1,51*51)];
Ell = out.T*Cir + out.b;

figure;plot3(Ell(1,:),Ell(2,:),Ell(3,:),'b.',mag(:,1),mag(:,2),mag(:,3),'r.');grid on;
out.Cir = Cir;
out.Ell = Ell;

out.acc = samp.acc';
out.mag = out.T\(samp.mag'-out.b);
out.t = samp.t;

out.att = calc_heading(out.t,out,phins);