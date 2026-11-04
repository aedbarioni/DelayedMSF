% SL_DelayedDiffusive
%% MSF landscape for delayed diffusive Stuart-Landau oscillators
% This code calculates the maximum Lyapunov exponent corresponding to
% complex eigenvalues of the Laplacian matrix that couples the system and
% plots a colormap where blue represents the negative (stable) Lyapunov 
% exponents and red the positive ones (unstable)

% clear all; close all; clc;

%% Parameters
lambda = 0.1;
omega = 0.25; % 0.12; % -0.28; % 0.12; %
gamma = -4.42; 
tau =  input('tau: '); %  %0.5; % 1.8*pi; %%% Modified

    
    
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


%% Fixed point
Omega0 = omega - lambda*gamma;  

% === Compute phi and r2 ===
phi = -Omega0 * tau;
r2 = lambda;

r0 = sqrt(r2);

%% Matrices 
P = [-2*lambda 0;
    -2*gamma*r0 0];

Rtau = [cos(phi) -r0*sin(phi);
        (1/r0)*sin(phi) cos(phi)];

J1 = P;

re_range = 25;
im_range = 25;
m = 200;
alpha = linspace(-re_range,re_range,m);
beta = linspace(-im_range,im_range,m);
mtle = zeros(m,m);


for ii = 1:m
    ii
    parfor jj = 1:m
        J2 = -(alpha(ii) + 1i*beta(jj))*Rtau;
        mtle(jj,ii) = dde_rightmost_eig(J1,J2,tau); 
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
xlabel('Re(\nu)')
ylabel('Im(\nu)')
% Compose the title with variables
title(sprintf('Lyapunov Exponent (fracNeg = %.3f)', fraction_negatives));


% Apply colormap and colorbar
colormap(bluewhitered);

