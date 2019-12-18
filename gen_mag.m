function mag = gen_mag(t_end,dt,D,b)

m_n = [0.205796;-0.040654;0.468785];

t_begin = 0;
%t_end   = 600;
%dt      = 0.1;
t = (t_begin:dt:t_end)';

eul = [ 180*(pi/180)*sin( (pi/60)*t'); 45*(pi/180)*sin((pi/13)*t');  15*(pi/19)*sin((pi/9.6)*t')]';
I = eye(3);

for k = 1:size(t)

    e(1:3,k) = 0.002 * randn(3,1);
        
end

for k = 1:size(t)

  % compute maggue measurement
  A = eul2rotm(eul(k,:), 'zyx'); 
  mag(1:3,k) = (I + D) \ ( A * m_n  +  b + e(1:3,k));
end