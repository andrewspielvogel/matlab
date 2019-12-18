function out = vector_io_adap2(samp,Tb,k,phins)

K = diag(k);

Rk=eye(3);
Rm=eye(3);

samp.mag = samp.mag(1:20:end,:);
%samp.mag = (rph2R([0;0;pi])*samp.mag(1:20:end,:)')'; %kvh
%samp.mag = (rph2R([0;0;pi])*samp.mag')'; %kvh
Rk = rph2rot([-90;0;90]*pi/180); %kvh
samp.t = samp.t(1:20:end);
%Rm = rph2rot([0;0;-180]*pi/180); %mst

%Rm = rph2rot([-pi;0;0]);
R = Rm*Rk;
samp.mag = (R*samp.mag')';

T0 = inv([Tb(1),Tb(3);Tb(3),Tb(2)]);
Tinv2 = T0*T0;
bTinv2b = Tb(4:5)'*Tinv2*Tb(4:5);

m_norm = norm([0.205796;-0.040654]);

phi = m_norm^2 -bTinv2b;

theta0 = [Tinv2(1,1);Tinv2(2,2);Tinv2(1,2);Tinv2*Tb(4:5)]/phi;
%theta0 = Tb;


tic
[t,theta] = ode45(@(t,theta) calc_thetadot(t,theta,samp,K),samp.t,theta0);
toc

out.t = t;
out.theta = theta';



mag=samp.mag;
num = size(mag,1);
W = zeros(num,5);


out.tp = zeros(3,num);
out.b = zeros(2,num);


for i=1:num
    
W(i,:) = [mag(i,1)^2,mag(i,2)^2,2*mag(i,1)*mag(i,2),-2*mag(i,1:2)];

T_inv2 = [out.theta(1,i),out.theta(3,i);out.theta(3,i),out.theta(2,i)];

out.tp(:,i) = [out.theta(1,i);out.theta(2,i);out.theta(3,i)];
out.b(:,i) = T_inv2\out.theta(4:5,i);

end

T_inv2 = [out.theta(1,end),out.theta(3,end);out.theta(3,end),out.theta(2,end)];
T = inv(sqrtm(T_inv2));
out.T = T/T(1,1);


out.W  = W;
out.theta_ls = (W'*W)\W'*ones(num,1);
T_inv2 = [out.theta_ls(1,end),out.theta_ls(3,end);out.theta_ls(3,end),out.theta_ls(2,end)];

out.T_ls = inv(sqrtm(T_inv2));
out.b_ls = T_inv2\out.theta_ls(4:5,:);
out.T_ls = out.T_ls/out.T_ls(1,1);


mag = samp.mag;

out.mag_adap = out.T\(mag(:,1:2)' - out.b(:,end));
out.mag_ls = out.T_ls\(mag(:,1:2)' - out.b_ls(:,end));
out.mag_hard = (mag(:,1:2)' - out.b_ls(:,end));


out.heading_adap = atan2(-out.mag_adap(2,:),out.mag_adap(1,:));
out.heading_ls = atan2(-out.mag_ls(2,:),out.mag_ls(1,:));
out.heading_uncal = atan2(-samp.mag(:,2),samp.mag(:,1));

out.heading_adap = resample2(samp.t,out.heading_adap',phins.t);
out.heading_ls = resample2(samp.t,out.heading_ls',phins.t);
out.heading_uncal = resample2(samp.t,out.heading_uncal,phins.t);


out.samp = samp;

plot_vect_io(out,phins);


function theta_dot = calc_thetadot(t,theta,samp,K)

    m = interp1(samp.t,samp.mag,t)';

    u = [m(1)*m(1);m(2)*m(2);2*m(1)*m(2);-2*m(1:2)];

    theta_dot = K*(-u*u'*theta + u);

