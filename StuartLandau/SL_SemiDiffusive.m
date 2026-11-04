%% MSF landscape for semi-diffusive Stuart-Landau oscillators
% This code calculates the maximum Lyapunov exponent corresponding to
% complex eigenvalues of the Laplacian matrix that couples the system and
% plots a colormap where blue represents the negative (stable) Lyapunov 
% exponents and red the positive ones (unstable)

% clear all; close all; clc;

%% Parameters
lambda = 0.1;
omega = 0.25; 
gamma = -4.42; 

tau = input('tau: '); 
d_star = input('d_star: ');  


%%%%%--- Setting up the paralelizing ---%%%%%
poolobj = gcp('nocreate');
delete(poolobj);
ncores = input('Number of cores: ');
parpool(ncores);

%% Parameters vector
param.lambda = lambda;
param.omega = omega;
param.gamma = gamma;
param.tau = tau;
param.d_star = d_star;


% finding the sync solution
f_omega = @(x) f_omega_aux(x,param);

%% Fixed point
Omega0 = fzero(f_omega, omega);  % omega is the initial guess

% === Compute phi and r2 ===
phi = -Omega0 * tau;
r2 = lambda + d_star * (cos(phi) - 1);

r0 = sqrt(r2);

%% Matrices 
P = [lambda-3*r2 0;
    -2*gamma*r0 0];
R0 = [-1 r0*sin(phi);
      -(1/r0)*sin(phi) -cos(phi)];
Rtau = [cos(phi) -r0*sin(phi);
        (1/r0)*sin(phi) cos(phi)];

J1 = P + d_star*R0;

re_range = 25;
im_range = 25;
m = 200;
alpha = linspace(-re_range,re_range,m); % Real part of \nu
beta = linspace(-im_range,im_range,m); % Imaginary part of \nu
mtle = zeros(m,m);

tic
for ii = 1:m
    % ii
    parfor jj = 1:m
        % if (ii==91 && jj==99)
        %     disp('stop')
        % end
        J2 = (alpha(ii) + 1i*beta(jj))*Rtau;
        mtle(jj,ii) = dde_rightmost_eig(J1,J2,tau); 
    end
end

fraction_negatives = sum(mtle(:) < 0) / numel(mtle);

disp(fraction_negatives);
filename = "filename.mat";
save(filename)


%%
figure();

imagesc(mtle(end:-1:1, 1:end)) % Flip vertically and exclude extra column
xlabel('Re(\nu)')
ylabel('Im(\nu)')
% Compose the title with variables
title(sprintf('Lyapunov Exponent (\\indegree = %.3f, fracNeg = %.3f)', d_star, fraction_negatives));


% Apply colormap and colorbar
colormap(bluewhitered);


%% Finding the sync solution
function y = f_omega_aux(x,param)
    omega = param.omega;
    gamma = param.gamma;
    lambda = param.lambda;
    d_star = param.d_star;
    tau = param.tau;

    y = omega - gamma*(lambda + d_star*(cos(-x*tau) - 1))...
	    + d_star*(sin(-x*tau)) - x;
end