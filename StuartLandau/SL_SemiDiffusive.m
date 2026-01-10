%% MSF landscape for semi-diffusive Stuart-Landau oscillators
% This code calculates the maximum Lyapunov exponent corresponding to
% complex eigenvalues of the Laplacian matrix that couples the system and
% plots a colormap where blue represents the negative (stable) Lyapunov 
% exponents and red the posive ones (unstable)

% clear all; close all; clc;

%% Parameters
lambda = 0.1;
omega = 0.25; 
gamma = -4.42; 
K = input('kappa: ');
tau = input('tau: '); 
mu = input('mu: ');  
kmu = K*mu;

%%%%%--- Setting up the paralelizing ---%%%%%
poolobj = gcp('nocreate');
delete(poolobj);
ncores = 6; % input('Number of cores: ')
parpool(ncores);

%% Parameters vector
param.lambda = lambda;
param.omega = omega;
param.gamma = gamma;
param.tau = tau;
param.mu = mu;
param.kmu = kmu;

% finding the sync solution
f_omega = @(x) f_omega_aux(x,param);

%% Fixed point
Omega0 = fzero(f_omega, omega);  % omega is the initial guess

% === Compute phi and r2 ===
phi = -Omega0 * tau;
r2 = lambda + kmu * (cos(phi) - 1);

r0 = sqrt(r2);

%% Matrices 
P = [lambda-3*r2 0;
    -2*gamma*r0 0];
R0 = [-1 r0*sin(phi);
      -(1/r0)*sin(phi) -cos(phi)];
Rtau = [cos(phi) -r0*sin(phi);
        (1/r0)*sin(phi) cos(phi)];

J1 = P + K*mu*R0;

re_range = 10;
im_range = 10;
m = 200;
alpha = linspace(-re_range,re_range,m);
beta = linspace(-im_range,im_range,m);
mtle = zeros(m,m);

tic
for ii = 1:m
    % ii
    parfor jj = 1:m
        % if (ii==91 && jj==99)
        %     disp('stop')
        % end
        J2 = K*(alpha(ii) + 1i*beta(jj))*Rtau;
        mtle(jj,ii) = dde_rightmost_eig(J1,J2,tau,1); 
    end
end

fraction_negatives = sum(mtle(:) < 0) / numel(mtle);

disp(fraction_negatives);
filename = "filename.mat";
save(filename)
runtime = toc


%%
figure();

imagesc(mtle(end:-1:1, 1:end)) % Flip vertically and exclude extra column
xlabel('\alpha')
ylabel('\beta')
% Compose the title with variables
title(sprintf('Lyapunov Exponent (\\indegree = %.3f, fracNeg = %.3f)', mu, fraction_negatives));


% Apply colormap and colorbar
colormap(bluewhitered);


%% Finding the sync solution
function y = f_omega_aux(x,param)
    omega = param.omega;
    gamma = param.gamma;
    lambda = param.lambda;
    kmu = param.kmu;
    tau = param.tau;

    y = omega - gamma*(lambda + kmu*(cos(-x*tau) - 1))...
	    + kmu*(sin(-x*tau)) - x;
end