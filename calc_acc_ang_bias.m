function out = calc_acc_ang_bias(samp,wb,ab)


% gains for adaptive observer
K.ka = 1;
K.K_bw = eye(3)/1000;
K.K_ba = eye(3);

%K.K_bw(3,3) = 0;
%K.K_ba(3,3) = 0;

theta0 = [samp.acc(1,:)';wb;ab];
% run ode45
tic
%[t,theta] = ode45(@(t,theta) calc_thetadot(t,theta,samp,K),samp.t,theta0);


theta = zeros(9,size(samp.t,1));
theta(:,1) = theta0;


for i=2:size(samp.t,1)
    
    ang = samp.ang(i,:);
    acc = samp.acc(i,:);
    
    da = theta(1:3,i-1) - acc';

    theta_dot = [-skew(ang)*(theta(1:3,i-1)-theta(7:9,i-1)) + skew(theta(4:6,i-1))*theta(1:3,i-1) - K.ka*da;
                 -K.K_bw*skew(acc)*da;
                  K.K_ba*skew(ang)*da];
   
    theta(:,i) = theta(:,i-1) + (samp.t(i)-samp.t(i-1))*theta_dot;
    
end
toc


% set outputs
out.t = samp.t;%t;
out.theta = theta;%theta';

figure;plot(out.t,out.theta(1:3,:));grid on;
figure;plot(out.t,out.theta(4:6,:));grid on;xlabel('Seconds');ylabel('w_b');
figure;plot(out.t,out.theta(7:9,:));grid on;xlabel('Seconds');ylabel('a_b');


function theta_dot = calc_thetadot(t,theta,samp,K)

    acc = interp1(samp.t,samp.acc,t);
    ang = interp1(samp.t,samp.ang,t);
    
    da = theta(1:3) - acc';
    
    theta_dot = [-skew(ang)*(theta(1:3)-theta(7:9)) + skew(theta(4:6))*theta(1:3) - K.ka*da;
                 -K.K_bw*skew(acc)*da;
                  K.K_ba*skew(ang)*da];
    