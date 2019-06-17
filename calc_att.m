function out = calc_att(K,x,y,w,t,R0)

R_prev = R0;

out.att = zeros(3,size(t,1));

out.att(:,1) = rot2rph(R0)';

tic
for i=2:size(t,1)
    
    dR = calc_dR(K,R_prev,x(:,i),y(:,i),w(:,i));
    
    dt = t(i) - t(i-1);
    
    R = R_prev*expm(dR*dt);
    
    out.att(:,i) = rot2rph(R*rph2R([0;0;pi])')';
    %out.att(:,i) = rot2rph(R*rph2R([0;0;0])')';

    %out.att(:,i) = rot2rph(R{i})';
    
    R_prev = R;

end
toc

out.t = t;

function dR = calc_dR(K,R,x,y,w)
% generate a_n

lat = 39.32*pi/180;
r = 6371*1000;
Ren = [-sin(lat),0,-cos(lat);0,1,0;cos(lat),0,-sin(lat)];
a_e = [cos(lat);0;sin(lat)] - (15.04*pi/180/3600)^2*cos(lat)*[r;0;0]/9.81;
a_n = Ren'*a_e;


w_E_e = [0;0;1]*15.04*pi/180/3600;
w_E_n = Ren'*w_E_e;

P = R'*[0;0;1]*[0,0,1]*R;
P = R'*(a_n*a_n')*R/(a_n'*a_n);

x = x/norm(x);
y = (eye(3)-P)*y;
y = y/norm(y);



x_tilde = K.x*skew(x)*R'*[0;0;-1];
y_tilde = K.y*(skew(y)*R'*[1;0;0]);

dR = skew(x_tilde + y_tilde + w - R'*w_E_n);