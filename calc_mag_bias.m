function out = calc_mag_bias(samp,theta0,phins,k)

K.K = diag(k);

samp.t = samp.t(4:end);
samp.acc = samp.acc(4:end,:);
samp.ang = samp.ang(4:end,:);
samp.mag = samp.mag(4:end,:);

    %R = rph2rot([-pi/2;0;pi/2]);
    %Rmg = rph2rot([pi;0;0]);


    %samp.acc = (R*samp.acc')';
    %samp.mag = (Rmg*R*samp.mag')';
    
    
%     R = rph2rot([0;0;pi]);
%     samp.acc = (R*samp.acc')';
%     samp.mag = (R*samp.mag')';


% run ode45
tic
[t,theta] = ode45(@(t,theta) calc_thetadot(t,theta,samp,K),samp.t,theta0);
toc

% set outputs
out.t = t;
out.theta = theta';
%figure;plot(out.t,out.theta);grid on;


% calc m_b, alpha, and T
T_inv2 = [out.theta(1,end),out.theta(4,end),out.theta(5,end);out.theta(4,end),out.theta(2,end),out.theta(6,end);out.theta(5,end),out.theta(6,end),out.theta(3,end)];

out.T = inv(sqrtm(T_inv2));
out.b = T_inv2\out.theta(7:9,:);

out.mag = calc_mag_cor(samp,out);
out.acc = samp.acc';
out.ang = samp.ang';
% rms = calc_mag_rms(out,phins);
% 
% rms_error = rms.rms_error;
% 
% out.rmse = rms_error;
% out.att = rms.heading.att;
% out.att_error = rms.heading.att_error;
% out.mag_ned = rms.heading.mag_ned;

%fprintf('rms_error: %f\n',rms_error);


function mag = calc_mag_cor(samp,out)

    mag = zeros(size(samp.mag'));
    
    for i=1:size(samp.t,1)
       
        % calc m_b, alpha, and T
        T_inv2 = [out.theta(1,i),out.theta(4,i),out.theta(5,i);out.theta(4,i),out.theta(2,i),out.theta(6,i);out.theta(5,i),out.theta(6,i),out.theta(3,i)];

        T = inv(sqrtm(T_inv2));
        b = T_inv2\out.theta(7:9,i);

        mag(:,i) = T\(samp.mag(i,:)'-b);

    end

function theta_dot = calc_thetadot(t,theta,samp,K)

    mag = interp1(samp.t,samp.mag,t);

    W = [mag(1)^2,mag(2)^2,mag(3)^2,2*mag(1)*mag(2),2*mag(1)*mag(3),2*mag(2)*mag(3),-2*mag,1];
        
    theta_dot = -K.K*(W'*W)*theta;


