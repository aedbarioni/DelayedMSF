%% MSF landscape for non-diffusive Lang-Kobayashi oscillators
% This code calculates the maximum Lyapunov exponent corresponding to
% complex eigenvalues of the Laplacian matrix that couples the system and
% plots a colormap where blue represents the negative (stable) Lyapunov 
% exponents and red the posive ones (unstable)

% clear all; close all; clc;

%%%%%--- Defining the parameters ---%%%%%
N_0 = 1.5;        %1.5e08;  (rescaled N_0_tilde = N_0*1e-8)
g = 1.5e03;       %1.5e-05; (rescaled g_tilde = g*1e8)
s = 1e-03;        %1e-07;   (rescaled s_tilde = s*1e4)
gamma = 500;      
alpha0 = 5;       
gamma_n = 0.5;    
sigma0 = 3;       
omega0 = 0;       
epsilon = 0.1;    
tau = input('tau: '); 
tau_tilde = 0;    
J_0 = 4*gamma_n*(N_0-(gamma/g)); 
M = 10; 
dx = 0.95;
neig = M;
a = J_0/(gamma_n*(N_0+(gamma/g)));

%%%%%--- Defining topology of coupling ---%%%%%
adj = ones(M,M);
x00 = [5.3413, 0]; 

coup_str = input('kappa: '); 

%%%%%--- Setting up the paralelizing ---%%%%%
poolobj = gcp('nocreate');
delete(poolobj);
ncores = input('Number of cores: ')
parpool(ncores);

deg = input('deg: '); 
neig = 10;

%%%%--- Parameters vector ---%%%%%
param.N_0 = N_0;       
param.g = g;       
param.s = s;        
param.gamma = gamma;      
param.alpha0 = alpha0;       
param.gamma_n = gamma_n;     
param.sigma0 = sigma0;       
param.omega0 = omega0;      
param.epsilon = epsilon;   
param.a = a;            
param.tau = tau;  
param.tau_tilde = tau_tilde;
param.J_0 = J_0;
param.M = M; 
param.dx = dx;
param.adj = adj;
param.deg = deg; 
param.neig = neig;

% finding the sync solution
x= fminsearch(@(x)sync(x,param),x00);

r = x(1);
Omega = x(2);
N =((((g*r^2)/(1+s*r^2))*1e-4 + gamma_n)^(-1))*(J_0 + ((g*r^2)/(1+s*r^2))*N_0*1e-4);

r0 = repmat(r,M,1); N00 = repmat(N,M,1); Omega0 = Omega; delta = zeros(M,1);

%%%%%--- Calculating Matrices ---%%%%%
M1 = zeros(3,3);
M2 = zeros(3,3);
M3 = zeros(3,3);

% Defining DF matrix
M1(1,1) = (1/2)*((g*(N-N_0)/((1+s*r^2)^2))*(1-s*r^2) - gamma);
M1(1,3) = (g/2)*(r/(1+s*r^2));
M1(2,1) = -(alpha0*g*s*r)*((N-N_0)/((1+s*r^2)^2));
M1(2,3) = (alpha0*g/2)*(1/(1+s*r^2));
M1(3,1) = -2*g*r*((N-N_0)/((1+s*r^2)^2))*1e-4;
M1(3,3) = -(gamma_n + 1e-4*g*((r^2)/(1+s*r^2)));

% Defining DH
M2(1,2) = -r*sin(Omega*tau);
M2(2,1) = (1/r)*sin(Omega*tau);
M2(2,2) = -cos(Omega*tau);

% Defining D(tau)H
M3(1,1) = cos(Omega*tau);
M3(1,2) = r*sin(Omega*tau);
M3(2,1) = -(1/r)*sin(Omega*tau);
M3(2,2) = cos(Omega*tau);

J1 = M1 + (coup_str/neig)*deg*M2;

m = 200;
re_fin = 25;
im_fin = 10;
Re_alpha = linspace(-re_fin,re_fin,m);
Im_beta = linspace(-im_fin,im_fin,m);
mtle = zeros(m,m);

for ii = 1:m
    ii
    parfor jj = 1:m
        J2 = (Re_alpha(ii) + 1i*Im_beta(jj))*(coup_str/neig)*M3;
        mtle(jj,ii) = dde_rightmost_eig(J1,J2,tau,M);
    end
end

filename = "filename.mat";
save(filename)


%%
figure();

imagesc(mtle(end:-1:1, 1:end)) % Flip vertically and exclude extra column
xlabel('\alpha')
ylabel('\beta')
% Compose the title with variables
title(sprintf('Lyapunov Exponent (\\indegree = %.3f, fracNeg = %.3f)', deg, fraction_negatives));


% Apply colormap and colorbar
colormap(bluewhitered);


%% Finding the sync solution
function y = sync(x,param)
N_0 = param.N_0;        
g = param.g;       
s = param.s;        
gamma = param.gamma;      
alpha0 = param.alpha0;       
gamma_n = param.gamma_n;    
kappa = param.kappa;      
sigma0 = param.sigma0;       
omega0 = param.omega0;       
a = param.a;            
tau = param.tau;              
deg = param.deg;
neig = param.neig;

F = zeros(2,1);
    
    F(2) = (1/2)*(((1+((g/gamma_n)*1e-4 +s)*(x(1))^2)^(-1))*((a-1)*g*N_0+ a*gamma)-gamma) + (kappa/neig)*deg*cos(x(2)*tau);
    F(1) = -x(2) + (alpha0/2)*(((1+((g/gamma_n)*1e-4 +s)*(x(1))^2)^(-1))*((a-1)*g*N_0+ a*gamma)-gamma) + sigma0*omega0 - (kappa/neig)*deg*sin(x(2)*tau);   

    y = log10(norm(F));
end

