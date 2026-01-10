%% MSF dependence of coupling strength and degree
% This code calculates the MSF and determine the deepness and size of the
% stable region in the complex plane for varying values of coupling
% strength and degree

% clear all; close all; clc;

%% Parameters
lambda = 0.1;
omega =  0.25;
gamma = -4.42; 
K = 3.5;
tau = 0.1; 
M = 3;
a_max = 6;

%%%%%--- Defining topology of coupling ---%%%%%
adj = ones(M,M) - eye(M);
sum_deg = 0;
for jjj = 1:M
    sum_deg = sum_deg + adj(1,jjj);
end
deg = sum_deg;

%% Parameters vector
mu = deg;
kmu = K*mu;
param.lambda = lambda;
param.omega = omega;
param.gamma = gamma;
param.K = K;
param.tau = tau;
param.mu = mu;
param.deg = mu;
param.kmu = kmu;
param.M = M;

m = 5000;
re_fin = 50;
im_fin = 25;
deg_v = linspace(0,30,m);
kappa_v = linspace(0,5,m);
mtle = zeros(m,m);
sizeXi_m = zeros(m,m);
deepXi_m = zeros(m,m);
Omega_m = zeros(m,m);
x00 = [5.3413, 0];

tic 
Omega_ini = 0;
for ii = 1:1 
    for jj = 1:m
        jj
        deg = mu; 
        coup_str = kappa_v(jj);

        % finding the sync solution
        f_omega = @(x) f_omega_aux(x,param,coup_str,deg);
        
        %% Fixed point
        Omega0 = fminsearch(@(x)f_omega_aux(x,param,coup_str,deg),Omega_ini); % frequency of the sync solution
        Omega_ini = Omega0;
        phi = -Omega0*tau;
        r2 = lambda + kmu*(cos(phi) - 1); % square amplitude of the sync solution
        r0 = sqrt(r2);
        
        %% Matrices 
        P = [lambda-3*r2 0;
            -2*gamma*r0 0];
        R0 = [-1 r0*sin(phi);
              -(1/r0)*sin(phi) -cos(phi)];
        Rtau = [cos(phi) -r0*sin(phi);
                (1/r0)*sin(phi) cos(phi)];
        
        [sizeXi,deepXi] = XiChar(mu,coup_str,tau,M,P,R0,Rtau);
        sizeXi_m(jj,ii) = sizeXi;
        deepXi_m(jj,ii) = deepXi;
        Omega_m(jj,ii) = Omega0;  
    end    
end

runtime = toc
filename = "filename.mat";
save(filename)


%% Plotting figure for deepness and size of stability region 
figure();
sizeXi_m(1,1) = 1;
deepXi_m(1,1) = 0;
plot(kappa_v,10*sizeXi_m(:,1),'LineWidth',2.0,'color',[0.4660 0.6740 0.1880]);
hold on 
plot(kappa_v,deepXi_m(:,1),'LineWidth',2.0,'color',[0.8500 0.3250 0.0980]);
yline(0, 'k', 'LineWidth', 0.5); 
xlabel('\kappa'); ylabel('\Lambda')

axis square

%% Calculating the deepness and size of the stable region of the MSF in the complex plane
function [sizeXi,deepXi] = XiChar(mu,K,tau,M,P,R0,Rtau)
    J1 = P + K*mu*R0;

    m = 201;
    re_fin = 25;
    im_fin = 10;
    Re_alpha = linspace(-re_fin,re_fin,m);
    Im_beta = linspace(-im_fin,im_fin,m);
    mtle = zeros(m,m);
    
    for ii = 1:m    
        for jj = 1:m
            J2 = (Re_alpha(ii) + 1i*Im_beta(jj))*K*Rtau;
            mtle(jj,ii) = dde_rightmost_eig(J1,J2,tau,M);
        end
    end
    deepXi = min(min(mtle));
    sizeXi = sum(mtle(:) < 0)/(m^2);

end

%% Finding the synchronous solution 
function y = f_omega_aux(x,param,k,mu)
    omega = param.omega;
    gamma = param.gamma;
    lambda = param.lambda;
    kmu = k*mu; 
    tau = param.tau;

    F = omega - gamma*(lambda + kmu*(cos(-x*tau) - 1))...
	    + kmu*(sin(-x*tau)) - x;
    y = log10(norm(F));
end
