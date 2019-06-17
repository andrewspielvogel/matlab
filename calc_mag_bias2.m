function out = calc_mag_bias2(samp,theta0,phins)


% gains for adaptive observer
K.K = eye(6)*0;
K.K(4:5,4:5) = eye(2)/5;    
%K.K(6,6)=1/100;



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

out.b = [inv(T_inv2)*out.theta(4:5,end);0];

figure;plot(out.t,out.theta);grid on;

%out.mag = inv(out.T)*(samp.mag'-out.b);
out.mag = calc_mag_cor(samp,out);
out.heading = calc_heading(samp.t,out.mag,phins);
figure;plot(out.mag(1,:),out.mag(2,:));grid on;axis equal;
figure;plot(phins.t,out.heading.heading_error);grid on;


function mag = calc_mag_cor(samp,out)

    mag = zeros(size(samp.mag'));
    
    for i=1:size(samp.t,1)
       
        % calc m_b, alpha, and T
        T_inv2 = [out.theta(1,i),out.theta(3,i);out.theta(3,i),out.theta(2,i)];

        T = [inv(sqrtm(T_inv2)),zeros(2,1);0,0,1];
        b = [inv(T_inv2)*out.theta(4:5,i);0];

        %mag(:,i) = inv(T)*(samp.mag(i,:)'-b);
        mag(:,i) = (samp.mag(i,:)'-b);

    end


function theta_dot = calc_thetadot(t,theta,samp,K)

    mag = interp1(samp.t,samp.mag,t);

    W = [mag(1)^2,mag(2)^2,2*mag(1)*mag(2),-2*mag(1:2),1];
         
    theta_dot = -K.K*W'*W*theta;