function f = plot_calc_mag_bias_V(out,T_true,b_true,K)

T_inv2 = inv(T_true)*inv(T_true);
m_n = [0.205796;-0.040654;0.468785];

theta_true = [T_inv2(1,1);T_inv2(2,2);T_inv2(3,3);T_inv2(1,2);T_inv2(1,3);T_inv2(2,3);T_inv2*b_true;b_true'*T_inv2*b_true-norm(m_n)*norm(m_n)]';


dtheta = out.theta - theta_true;

K_inv = inv(K.K);

V = sum(dtheta.*repmat(diag(K_inv)',size(dtheta,1),1).*dtheta/2,2);

figure;plot(out.t,V);grid on;xlabel('Time (Seconds)');ylabel('V');
figure;plot(out.t,dtheta);grid on;xlabel('Time (Seconds)');ylabel('d\theta');
