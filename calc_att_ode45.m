function out = calc_att_ode45(K,x,y,w,sampt,R0)


tic
[out.t,R] = ode45(@(t,R) calc_dR(t,R,x,y,w,sampt,K),sampt,R0);
toc
    
out.R = R;

function dR = calc_dR(t,R,acc,mag,ang,sampt,K)
% generate a_n

    x = interp1(sampt,acc',t);
    y = interp1(sampt,mag',t);
    w = interp1(sampt,ang',t);
    
x = x/norm(x);
y = y/norm(y(1:2));

P = reshape(R,3,3)'*[0;0;1]*[0,0,1]*reshape(R,3,3);

x_tilde = K.x*skew(x)*reshape(R,3,3)'*[0;0;-1];
y_tilde = K.y*P*(skew(y)*reshape(R,3,3)'*[1;0;0]);

dR = skew(x_tilde + y_tilde + w);
dR = reshape(dR,9,1);