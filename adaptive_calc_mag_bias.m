function out = adaptive_calc_mag_bias(samp,T,wb,mb)

Km = [1,0,0;0,1,0;0,0,1];
Ka = .1;
Kg = .1;

m_n = [0.205796;-0.040654;0.468785];

alpha = reshape(kron(reshape(inv(T),9,1)',T),81,1);
gamma = reshape(T*skew(inv(T)*mb),9,1);

out.theta(:,1) = [alpha;gamma];


out.mag(:,1) = samp.mag(1,:)';
out.mag_cor(:,1) = inv(T)*(samp.mag(1,:)' - mb);


for i=2:size(samp.t,1)
    
    dm = out.mag(:,i-1) - samp.mag(i,:)';
    dt = samp.t(i) - samp.t(i-1);
    
    alpha_coef = kron(reshape(kron(samp.mag(i,:),skew(samp.ang(i,:)')),27,1)',eye(3));

    gamma_coef = kron(samp.ang(i,:),eye(3));
    

    dmag = -[alpha_coef,gamma_coef]*out.theta(:,i-1) - Km*dm;
    
    out.mag(:,i) = out.mag(:,i-1) + dt*dmag;
    
    out.theta(:,i) = out.theta(:,i-1) + dt*[Ka*alpha_coef';Kg*gamma_coef']*dm;
    
    %out.mag_cor(:,i) = inv(reshape(out.theta(1:9,i),3,3))*(samp.mag(i,:)' - out)
   
end
