function out = calc_mag_bias2(samp,theta0,k)
    
K.K = diag(k);
    
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
out.T = out.T/out.T(1,1);

out.b = [T_inv2\out.theta(4:5,:);zeros(1,size(out.theta,2))];

%out.mag = inv(out.T)*(samp.mag'-out.b);
out.mag = calc_mag_cor(samp,out);
out.acc = samp.acc';
out.ang = samp.ang';



function mag = calc_mag_cor(samp,out)

    mag = zeros(size(samp.mag'));
    
    for i=1:size(samp.t,1)
       
        % calc m_b, alpha, and T
        T_inv2 = [out.theta(1,i),out.theta(3,i);out.theta(3,i),out.theta(2,i)];
        
        T = [inv(sqrtm(T_inv2)),zeros(2,1);0,0,1];
        T=T/T(1,1);
        b = [T_inv2\out.theta(4:5,i);0];
        mag(:,i) = T\(samp.mag(i,:)'-b);

    end


function theta_dot = calc_thetadot(t,theta,samp,K)

    mag = interp1(samp.t,samp.mag,t);

    W = [mag(1)^2,mag(2)^2,2*mag(1)*mag(2),-2*mag(1:2),1];
         
    theta_dot = -K.K*(W'*W)*theta;