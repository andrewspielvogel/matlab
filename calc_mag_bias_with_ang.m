function out = calc_mag_bias_with_ang(samp,T,b)


% local true magnetic field
m_n = [0.205796;-0.040654;0.468785];

samp.m_n = 1;%m_n;
 
% gains for adaptive observer
K.m = 1;
K.phi = 1*100;
K.gamma = 1*100;

phi0 = kron(reshape(inv(T),9,1)',T);
gamma0 = kron((inv(T)*b)',T);

theta0 = [samp.mag(1,:)';stack(phi0);stack(gamma0)];

% run ode45
tic
[t,theta] = ode45(@(t,theta) calc_thetadot(t,theta,samp,K),samp.t,theta0);
toc

T=[.95,0,0;0,1.1,0;0,0,1.05];
b = [0.1;-0.3;-0.15];

phi_true   = stack(kron(stack(inv(T))',T));
gamma_true = stack(kron((inv(T)*b)',T));


out.theta = theta';
out.t = t;

figure;plot(out.t,out.theta(4:end,:)-[phi_true;gamma_true]);grid on;

function theta_dot = calc_thetadot(t,theta,samp,K)

    mag = interp1(samp.t,samp.mag,t)';
    ang = interp1(samp.t,samp.ang,t)';

    m_hat     = theta(1:3);
    phi_hat   = reshape(theta(4:84),3,27);
    gamma_hat = reshape(theta(85:end),3,9);
    
    dm = m_hat - mag;
        
    m_dot     = -phi_hat*stack(kron(mag',skew(ang))) + gamma_hat*stack(skew(ang)) - K.m*dm;
    phi_dot   = kron(stack(kron(mag',skew(ang))),dm);
    gamma_dot = -kron(stack(skew(ang)),dm);

    theta_dot = [m_dot;stack(phi_dot);stack(gamma_dot)];
    

