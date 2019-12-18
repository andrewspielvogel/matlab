function [D_est,b_est,n,Cov_est]=TWOSTEP_estimate(B,H,Sigma_noise)
%TWOSTEP ESTIMATE estimation of magnetometer errors.
% [D est,b est,n,Cov est]=TWOSTEP estimate(B,H,Sigma noise)
% B − 3xN matrix, where N is the number of signal samples
% H − Earth's magnetic field in Earth frame coordinates
% Sigma noise − 3xN noise covariance array, where each column is a
% diagonal of the noise covariance matrix.
%
% Implemented based on the work presented in:
% "Complete Linear Attitude−Independent Magnetometer Calibration"
% by R.Alonso and M. Shuster
% in "The Journal of the Astronautical Sciences" Vol. 50 No.4
% October−December 2002.
%
% Sections of code and optimisation adapted from the TWOSTEP code written
% by J.F.Vasconcelos for the validation of geometric calibration presented
% in "A geometric approach to strapdown magnetometer calibratio in sensor
% frame" by J.F. Vasconcelos et al. in Aerospace and Electronic Systems,
% IEEE Transactions 47(2) 2011. Many thanks go to Dr. J o s Vasconcelos for
% his assistance.
%
% Slight variation in the calculation of the centered J function, was used
% for clarity in derivations.
%
% Justin Dinale, DSTO Department of Defence, Australia
% $Revision: 1.0.0 $ $Date: 2012/11/05$
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Define stop condition and perform housekeeping
stop_tol=1e-24; % Stop Condition from Alonso paper
max_n=200; % Max number of iterations
k_samples=size(B,2); % Number of signal samples


%% TWOSTEP Centered estimate
% Set initial guess for b and D.
b0=[0 0 0]';
D0=zeros(3);


% Calculate L
L_k=[2*B;-(B.^2);-2*B(1,:).*B(2,:);-2*B(1,:).*B(3,:);-2*B(2,:).*B(3,:)];

% Calculate scalar measurement z
z_k=total_field(B).^2-total_field(H).^2;

% Calculate mean and variance of noise
mu_k=-sum(Sigma_noise,1); % Mean
sigma_sq_k=4*sum(((eye(3)+D0)*B-b0*ones(1,k_samples)).*(Sigma_noise.*((eye(3)+D0)*B-b0*ones(1,k_samples))),1)+5*(sum(Sigma_noise.^2,1))-(sum(Sigma_noise,1).^2);

% Calculate centered sigma squared
sigma_sq_bar=1/sum(1./sigma_sq_k);

% Center the data
[mu_bar,mu_k_tilde]=center_data(mu_k,sigma_sq_k,sigma_sq_bar);
[z_bar,z_k_tilde]=center_data(z_k,sigma_sq_k,sigma_sq_bar);
[L_bar,L_k_tilde]=center_data(L_k,sigma_sq_k,sigma_sq_bar);

% Calculate fisher information matrix
[I_fisher_tilde,I_fishinv_tilde]=TS_fisher_centered(sigma_sq_k,L_k_tilde);

% Calculate centered estimate
theta_0_tilde=I_fishinv_tilde*(L_k_tilde*((z_k_tilde-mu_k_tilde)'./sigma_sq_k'));


%% TWOSTEP Center correction
theta_n=theta_0_tilde; % Initiare theta for first iteration
n=0; % Initialise iteration counter
TS_err=Inf; %Initial condition for error.

% ABC is used to remove intensive calculations out of for loop.
ABC=-(((z_k_tilde-mu_k_tilde)./sigma_sq_k)*L_k_tilde')';

while (TS_err>stop_tol && n<max_n)
    
    if (n>0) % If we are not in the first iteration
        theta_n=theta_np1;
    end
    
    % Extract c and E components
    [c,E]=theta_to_c_E(theta_n);
    
    % Calculate second derivative of bˆ2 wrt theta
    tmp=((eye(3)+E)\c)*((eye(3)+E)\c)';
    dbsqdtheta_p=[2*((eye(3)+E)\c); -diag(tmp);-2*tmp(1,2); -2*tmp(1,3); -2*tmp(2,3)]';

    %Calculate gradient of J
    dJdThetap_tilde=ABC + I_fisher_tilde*theta_n;
    dJdThetap_bar=(-(1/sigma_sq_bar)*(L_bar' - dbsqdtheta_p)*(z_bar -L_bar'*theta_n +c'*((eye(3)+E)\c) - mu_bar))';
    dJdTheta=dJdThetap_tilde+dJdThetap_bar;

    % Calculate Fisher matrix
    [I_fisher_bar] =TS_fisher_center(sigma_sq_bar,L_bar,dbsqdtheta_p);
    
    % Update theta
    theta_np1=theta_n-(I_fisher_tilde+I_fisher_bar)\dJdTheta;

    % Calculate Error
    TS_err=(theta_np1-theta_n)'*(I_fisher_tilde+I_fisher_bar)*(theta_np1-theta_n);
    n=n+1; %Increase counter
end

[b_est,D_est]=theta_to_b_D(theta_np1);

% Extract covariance matrix
b=b_est;
D=D_est;
M_cD=[b(1) 0 0 b(2) b(3) 0; 0 b(2) 0 b(1) 0 b(3); 0 0 b(3) 0 b(1) b(2)];

M_ED=[2*D(1), 0, 0, 2*D(4), 2*D(5), 0; 0, 2*D(2), 0, 2*D(4), 0, 2*D(6); 0, 0, 2*D(3), 0, 2*D(5), 2*D(6); ...
D(4), D(4), 0, D(1)+D(2), D(6), D(5); D(5), 0, D(5), D(6), D(1)+D(3), D(4); 0, D(6), D(6), D(5), D(4), D(6)+D(5)]; 

dbD_dcE=([(eye(3)+D), M_cD; zeros(6,3), M_ED])\eye(9);
Cov_est=dbD_dcE*(I_fisher_tilde+I_fisher_bar)\dbD_dcE';
%END TWOSTEP

function [B_total]=total_field(B_in)
% Calculates total field of B in (3xn matrix).
B_total=sum(B_in.^2,1).^0.5;

function [X_bar,X_tilde]=center_data(X,sigma_sq_k,sigma_sq_bar)
% Calculates the centered and center components of X
% Center component
X_bar = sigma_sq_bar*(X*(1./sigma_sq_k)');
% Centered component
X_tilde = X - X_bar*ones(1,size(X,2));

function [I_fisher_tilde,I_fishinv_tilde] = TS_fisher_centered(sigma_sq,L_tilde)
% Calculates the fisher information matrix for the centered estimate, when
% given variance sigma sq and centered vectors of L tilde
% Calculate fisher information matrix
I_fisher_tilde=((L_tilde.*(ones(size(L_tilde,1),1)*(1./sigma_sq)))*L_tilde');
% Calculate inverse
I_fishinv_tilde=(I_fisher_tilde)\eye(9);

function [I_fisher_bar] = TS_fisher_center(sigma_sq_bar,L_bar,dbsqdtheta_p)
% Calculates center informaiton matris. Used for readability
I_fisher_bar=((L_bar'-dbsqdtheta_p)'*(L_bar'-dbsqdtheta_p))/sigma_sq_bar;

function [b,D]= theta_to_b_D(theta)
% Converts a value of theta to usable physical values
[c,E]=theta_to_c_E(theta);
[U,S]=eig(E);
W=-eye(3)+(eye(3)+S).^(0.5);
D=U*W*U';
% Calculate b using the inversr of (I+D)
b=(eye(3)+D)\c;

function [c,E]=theta_to_c_E(theta)
% Extracts c and E elements from theta as per Alonso paper.
c=theta(1:3);
E=[theta(4), theta(7), theta(8); theta(7), theta(5), theta(9); theta(8), theta(9), theta(6)];