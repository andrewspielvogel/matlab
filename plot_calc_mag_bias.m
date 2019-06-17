function f = plot_calc_mag_bias(out)

T_inv2 = [out.theta(1,end),out.theta(4,end),out.theta(5,end);out.theta(4,end),out.theta(2,end),out.theta(6,end);out.theta(5,end),out.theta(6,end),out.theta(3,end)];

figure;plot(out.t,out.theta(7:9,:));grid on;xlabel('seconds');ylabel('T^{-2}m_b');
figure;plot(out.t,out.theta(1:6,:));grid on;xlabel('Seconds');ylabel('T^{-2}');
figure;plot(out.t,out.theta(10,:));grid on;xlabel('Seconds');ylabel('b^TT^{-2}b-||m_t||^2');
figure;plot(out.t,inv((T_inv2))*out.theta(7:9,:));grid on;xlabel('seconds');ylabel('m_b');


