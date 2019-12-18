function out = calc_mag_bias(samp,theta0,k)


K.K = diag(k);


    %R = rph2rot([-pi/2;0;pi/2]);
    %Rmg = rph2rot([pi;0;0]);


    %samp.acc = (R*samp.acc')';
    %samp.mag = (Rmg*R*samp.mag')';
    
    
%     R = rph2rot([0;0;pi]);
%     samp.acc = (R*samp.acc')';
%     samp.mag = (R*samp.mag')';

%theta0 = [reshape(eye(3),9,1);zeros(4,1)];

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
T_inv2 = reshape(out.theta(1:9,end),3,3);
out.T = inv(sqrtm(T_inv2));
out.b = T_inv2\out.theta(7:9,:);
out.b = T_inv2\out.theta(10:12,:);

out.mag = real(out.T\(samp.mag' - out.b(:,end)));
%out.mag = real(calc_mag_cor(samp,out));
out.acc = samp.acc';
out.ang = samp.ang';

%rms = calc_heading(out.t,out,phins);
%out = rms;



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

    %W = [mag(1)^2,mag(2)^2,mag(3)^2,2*mag(1)*mag(2),2*mag(1)*mag(3),2*mag(2)*mag(3),-2*mag,1];
    
    W = [kron(mag,mag), -2*mag, 1];
        
    theta_dot = -K.K*(W'*W)*theta;


