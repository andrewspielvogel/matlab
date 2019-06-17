function out = calc_mag_bias(samp,theta0)


% local true magnetic field
m_n = [0.205796;-0.040654;0.468785];

samp.m_n = 1;%m_n;
 
% gains for adaptive observer
K.K = eye(10);
K.K(7:9,7:9) = eye(3)/10;
K.K(10,10) = 1/100;


% gain for scaling theta estimate to correct mag
K.k_alpha = 0.1;


% calc initial condition
%T_inv2 = inv(T)*inv(T);

%theta0 = [T_inv2(1,1);T_inv2(2,2);T_inv2(3,3);T_inv2(1,2);T_inv2(1,3);T_inv2(2,3);T_inv2*mb;mb'*T_inv2*mb-norm(m_n)*norm(m_n)];
%theta0 = [T_inv2(1,1);T_inv2(2,2);T_inv2(3,3);T_inv2(1,2);T_inv2(1,3);T_inv2(2,3);T_inv2*mb;0];


% run ode45
tic
[t,theta] = ode45(@(t,theta) calc_thetadot(t,theta,samp,K),samp.t,theta0);
toc

% set outputs
out.t = t;
out.theta = theta';

b_true = [1;-3;-1.5]/10;
T_true=[.95,0,0;0,1.1,0;0,0,1.05];
%T_true=[.7,0,0;0,1.3,0;0,0,1];


% calc m_b, alpha, and T
T_inv2 = [out.theta(1,end),out.theta(4,end),out.theta(5,end);out.theta(4,end),out.theta(2,end),out.theta(6,end);out.theta(5,end),out.theta(6,end),out.theta(3,end)];

out.T = inv(sqrtm(T_inv2));
out.b = inv(T_inv2)*out.theta(7:9,:);

% plot outputs
plot_calc_mag_bias(out);
%plot_calc_mag_bias_V(out,T_true,b_true,K);

function theta_dot = calc_thetadot(t,theta,samp,K)

    mag = interp1(samp.t,samp.mag,t);

    W = [mag(1)^2,mag(2)^2,mag(3)^2,2*mag(1)*mag(2),2*mag(1)*mag(3),2*mag(2)*mag(3),-2*mag,1];
    
    T_inv2 = [theta(1),theta(4),theta(5);theta(4),theta(2),theta(6);theta(5),theta(6),theta(3)];
    
    b = inv(T_inv2)*theta(7:9);
        
    theta_dot = -K.K*W'*W*theta;% - K.k_alpha*(b'*theta(7:9)-theta(10)-norm(samp.m_n)^2)*theta;


