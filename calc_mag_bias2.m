function out = calc_mag_bias2(samp,theta0,phins,k)

K.K = zeros(6,6);
% % gains for adaptive observer
K.K(1:3,1:3) = eye(3)*k(1);
K.K(4:5,4:5) = eye(2)*k(2);    
K.K(6,6)=k(3);

K.K

%     R = rph2rot([-pi/2;0;pi/2]);
%     Rmg = rph2rot([pi;0;0]);
% 
% 
%     samp.acc = (R*samp.acc')';
%     samp.mag = (Rmg*R*samp.mag')';
    
    
    R = rph2rot([0;0;pi]);
    samp.acc = (R*samp.acc')';
    samp.mag = (R*samp.mag')';

%K.K=K.K;

% calc initial condition
%T_inv2 = inv(T)*inv(T);

%theta0 = [T_inv2(1,1);T_inv2(2,2);T_inv2(1,2);T_inv2*mb;mb'*T_inv2*mb-norm(m_n)*norm(m_n)];

% run ode45
tic
[t,theta] = ode45(@(t,theta) calc_thetadot(t,theta,samp,K),samp.t,theta0);
toc

% set outputs
out.t = t;
out.theta = theta';

% calc m_b, alpha, and T
T_inv2 = [out.theta(1,end),out.theta(3,end);out.theta(3,end),out.theta(2,end)];

out.T = [inv(sqrtm(T_inv2)),zeros(2,1);0,0,1];

out.b = [T_inv2\out.theta(4:5,end);0];

%out.mag = inv(out.T)*(samp.mag'-out.b);
out.mag = calc_mag_cor(samp,out);
out.acc = samp.acc';
out.ang = samp.ang';
rms = calc_mag_rms(out,phins);

rms_error = rms.rms_error;

out.rmse = rms_error;
out.att = rms.heading.att;
out.att_error = rms.heading.att_error;
out.mag_ned = rms.heading.mag_ned;

%fprintf('rms_error: %f\n',rms_error);


function mag = calc_mag_cor(samp,out)

    mag = zeros(size(samp.mag'));
    
    for i=1:size(samp.t,1)
       
        % calc m_b, alpha, and T
        T_inv2 = [out.theta(1,i),out.theta(3,i);out.theta(3,i),out.theta(2,i)];
        
        T = [inv(sqrtm(T_inv2)),zeros(2,1);0,0,1];
        b = [T_inv2\out.theta(4:5,i);0];
        mag(:,i) = T\(samp.mag(i,:)'-b);

    end


function theta_dot = calc_thetadot(t,theta,samp,K)

    mag = interp1(samp.t,samp.mag,t);

    W = [mag(1)^2,mag(2)^2,2*mag(1)*mag(2),-2*mag(1:2),1];
         
    theta_dot = -K.K*(W'*W)*theta;