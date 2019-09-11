function out = vector_io_adap(samp,theta0,k)

K = diag(k);

tic
[t,theta] = ode45(@(t,theta) calc_thetadot(t,theta,samp,K),samp.t,theta0);
toc

out.t = t;
out.theta = theta';
T_inv2 = [out.theta(1,end),out.theta(4,end),out.theta(5,end);out.theta(4,end),out.theta(2,end),out.theta(6,end);out.theta(5,end),out.theta(6,end),out.theta(3,end)];

out.T = inv(sqrtm(T_inv2));
out.b = T_inv2\out.theta(7:9,:);

figure;plot(out.t,out.theta);grid on;
figure;plot(out.t,out.b);grid on;



function theta_dot = calc_thetadot(t,theta,samp,K)

    m = interp1(samp.t,samp.mag,t)';
    

    u = [m(1)*m(1);m(2)*m(2);m(3)*m(3);2*m(1)*m(2);2*m(1)*m(3);2*m(2)*m(3);-2*m];

    theta_dot = K*(-u*u'*theta + u);

