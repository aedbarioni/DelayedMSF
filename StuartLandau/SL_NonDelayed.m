%% MSF landscape for non-delayed Stuart-Landau oscillators
% This code calculates the maximum Lyapunov exponent corresponding to
% complex eigenvalues of the Laplacian matrix that couples the system and
% plots a colormap where blue represents the negative (stable) Lyapunov 
% exponents and red the posive ones (unstable)

% clear all; close all; clc;

%% Parameters
lambda = 0.1;
omega = -0.28;
gamma = -4.42;
K = 0.3;

%% Fixed point
r0 = sqrt(lambda);
Omega0 = omega - gamma*lambda;

%% Matrices 
P = [-2*lambda 0;
    -2*gamma*lambda 0];
DH = eye(2);

re_range = 10;
im_range = 10;
m = 200;
alpha = linspace(-re_range,re_range,m);
beta = linspace(-im_range,im_range,m);

for ii = 1:m
    ii
    for jj = 1:m
        J = P + (alpha(ii) + 1i*beta(jj))*DH;
        evals = eig(J);
        % Remove eigenvalues close to zero
        tol = 1e-6;
        evals(abs(evals) < tol) = [];
        mtle(jj,ii) = max(real(evals));
    end
end

fraction_negatives = sum(mtle(:) < 0) / numel(mtle);
disp(fraction_negatives);
filename = "filename.mat";
save(filename)

%%
figure();

imagesc(mtle(end:-1:1, 1:end)) % Flip vertically and exclude extra column
xlabel('\alpha')
ylabel('\beta')
% Compose the title with variables
title(sprintf('Lyapunov Exponent (fracNeg = %.3f)', fraction_negatives));


% Apply colormap and colorbar
colormap(bluewhitered);
