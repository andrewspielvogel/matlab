function acc = calc_g(samp)

    k = .1;
    
    [t,acc] = ode45(@(t,acc) calc_accdot(t,acc,samp,k),samp.t,samp.acc(:,1));
    
    

    
function theta_dot = calc_accdot(t,acc,samp,k)

    acc_samp = interp1(samp.t,samp.acc',t)';
    ang_samp = interp1(samp.t,samp.ang',t)';
    
    theta_dot = -skew(ang_samp)*acc - k*(acc-acc_samp);
