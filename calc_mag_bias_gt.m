function out = calc_mag_bias_gt(samp,theta0,phins)
 
%% gains for adaptive observer
K.a = eye(3);
K.b = eye(3);

%% run ode45
tic
% [t,theta] = ode45(@(t,theta) calc_thetadot(t,theta,samp,K),samp.t,theta0);
% theta = theta';

num_samp = size(samp.t,1);

theta = zeros(6,num_samp);

theta(:,1) = theta0;

for i=2:num_samp
    
    dt = samp.t(i) - samp.t(i-1);

    theta(:,i) = theta(:,i-1) + dt*calc_thetadot(samp.t(i),theta(:,i-1),samp,K);

end

toc

%% set outputs
out.t = samp.t;
out.theta = theta;
out.mag = samp.mag'-out.theta(4:6,end);

out.heading = calc_heading(samp.t,out.mag,phins);
figure;plot(phins.t,out.heading.heading_error);grid on;


function theta_dot = calc_thetadot(t,theta,samp,K)

    mag = interp1(samp.t,samp.mag,t);
    ang = interp1(samp.t,samp.ang,t);

    dx = theta(1:3) - mag';
    
    x_dot = -skew(ang)*(theta(1:3) - theta(4:6)) - K.a*dx;
    b_dot = K.b*skew(ang)*dx;
        
    theta_dot = [x_dot;b_dot];

