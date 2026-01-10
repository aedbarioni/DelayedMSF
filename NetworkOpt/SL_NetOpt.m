% SLOptAdj
% clear all; close all; clc;
%% Parameters
lambda = 0.1;
omega =  0.25; 
gamma = -4.42; 
K = 0.1;
tau = input('tau: ');
M = 5;
a_max = 6;

%%%%--- Defining the all-2-all network ---%%%%%
adj = ones(M,M);
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

%%%%%--- Setting up the paralelizing ---%%%%%
poolobj = gcp('nocreate');
delete(poolobj);
ncores = input('Number of cores: ')
parpool(ncores);

% finding the sync solution
f_omega = @(x) f_omega_aux(x,param);

%% Fixed point
Omega0 = fzero(f_omega,omega); % frequency of the sync solution
phi = -Omega0*tau;
r2 = lambda + kmu*(cos(phi) - 1); % square amplitude of the sync solution
r0 = sqrt(r2);

r0_ini = zeros(2*M,1);
r0_ini(1:M) = r0*ones(M,1); r0_ini(M+1) = Omega0; r0_ini(M+2:2*M) = zeros(M-1,1); 
r = r0_ini(1:M); Omega = Omega0; delta = zeros(M,1); 

% Initial mtle 
[J1,J2] = Jac_SL_hetnet(param,adj,r,Omega,delta,tau);
mtle_0 = dde_rightmost_eig(J1,J2,tau,M);
disp(mtle_0)

%%%%%--- Defining the perturbation variables ---%%%%%
trial = 100; % number of realizations of random oscillator heterogeneity
m = 500; % maximum number of steps
sstep = 0.5; % size of step 
mtle = Inf(trial,m);
numSteps = zeros(trial,1);

adj_best_v = cell(trial, m);
r0_best_v = cell(trial, m);
mtle_best_v = Inf(trial,m);

mtle_best_v(:,1) = mtle_0*ones(trial,1);
for jj = 1:trial
    adj_best_v{jj,1} = adj;
end

%%%%%--- Main loop ---%%%%%
tic
parfor jj = 1:trial
    disp(['jj= ',num2str(jj)])
    adj0 = adj;
    mtle_best = mtle_0;
    r0 = r0_ini;
    for ii = 2:m
        % Define the maximum number of retries and initial sstep value
        max_retries = 10; 
        retry_count = 0;
        success = false; % Flag to check if optimization is successful

        % a random realization of heterogeneity
        szhet = M^2;
        sample = normrnd(0,1,szhet,1); 
        sample = (sample - mean(sample))/std(sample);
        sstep_current = sstep; % Start with the initial sstep value

        % Set optimization options
        optionss = optimoptions('fmincon', 'Algorithm', 'sqp', 'Display', 'iter','MaxIterations', 5);
    
        % Define the objective function
        obj_fun = @(x) ObjAdj(x,adj0,r0,param);

        % Loop to attempt optimization with adjusted sstep
        while ~success && retry_count < max_retries
            x0 = (sstep_current / (5 * norm(sample))) * sample; % Initial guess

            % Define the nonlinear constraints function
            nonlcon = @(x) constraints(x, adj0, M, sstep_current, a_max);

            try
                % Run the optimization
                [x_opt, fval, exitflag] = fmincon(obj_fun, x0, [], [], [], [], -ones(M*(M-1), 1), [], nonlcon, optionss);
                success = true;

            catch ME
                % Check for specific error message
                if contains(ME.message, 'Objective function is undefined at initial point')
                    disp(['1 - Optimization failed at iteration ii = ', num2str(ii), ' with sstep = ', num2str(sstep_current)]);
                    disp(ME.message);
                    retry_count = retry_count + 1; % Increment retry count
                    success = false;
                    % Adjust sstep value for the next attempt
                    sstep_current = sstep_current * 0.8; % Reduce sstep by 20%
                elseif contains(ME.message, 'Finite difference derivatives at initial point contain Inf, NaN, or complex values')
                    disp(['2 - Optimization failed at iteration ii = ', num2str(ii), ' with sstep = ', num2str(sstep_current)]);
                    % Add your custom logic here to handle this case if needed
                    retry_count = retry_count + 1; % Increment retry count
                    success = false;
                    % Optionally adjust sstep or take another action
                    sstep_current = sstep_current * 0.8; % Reduce sstep by 20%
                else
                    success = true;
                    % Handle other errors
                    rethrow(ME);                    
                end
            end
        end

        if ~success
            disp('Optimization failed after maximum retries. Moving to next iteration.');
            continue; % Skip to the next iteration of the loop if unsuccessful
        end 
          
        % Run the optimization
        [mtle,r0_new,adj_h] = ObjAdj_full(x_opt,adj0,r0,param);

        % Sanity check
        if abs(mtle - fval)>1e-1
            disp("problem")
        end
        
        if abs(mtle - mtle_best)<1e-8
            disp("stopping because the improvement is too little")
            adj0 = adj_h;
            r0 = r0_new;
            mtle_best = mtle;
    
            mtle_best_v(jj,ii) = mtle_best;
            adj_best_v{jj,ii} = adj0;
            r0_best_v{jj,ii} = r0;
            numSteps(jj) = ii;
            break
        end

        adj0 = adj_h;
        r0 = r0_new;
        mtle_best = mtle;

        mtle_best_v(jj,ii) = mtle_best;
        adj_best_v{jj,ii} = adj0;
        r0_best_v{jj,ii} = r0;
        numSteps(jj) = ii;
        
    end
    
    

end

[mtle_best_tot,indmin] = min(mtle_best_v(:));
[jm, im] = ind2sub(size(mtle_best_v), indmin);
adj_best_tot = adj_best_v{jm,im};

% 
filename = "filename.mat";
save(filename)
runtime = toc


function [c, ceq] = constraints(x, adj0, M, sstep, a_max)
    % Tolerance for floating-point comparisons
    tol = 1e-8;

    % Ensure the norm of x is less than or equal to sstep
    c_norm = norm(x) - sstep;

    % Compute updated adjacency matrix
    adj_h = assemble_dif(x, adj0, M);

    % Inequality: adj_h >= 0  -->  -adj_h <= 0
    c_nonneg = -adj_h(:) + tol;

    % Inequality: adj_h <= a_max  -->  adj_h - a_max <= 0
    c_max = adj_h(:) - a_max;

    % Combine all inequality constraints
    c = [c_norm; c_nonneg; c_max];

    % No equality constraints
    ceq = [];
end

function [mtle,r0_new,adj_h] = ObjAdj_full(dif,adj0,r0,param)
% Inputs:
%        dif ---- vector with nondiagonal perturbation
%        adj0   ---- initial (unperturbed) matrix
%        r0     ---- sync sol for the initial matrix
% Outputs:
%        mtle   ---- Lyapunov exponent of solution of perturbed matrix
%        r0_new ---- sync sol for the perturbed matrix
%        adj_h  ---- perturbed matrix

lambda = param.lambda;      
gamma = param.gamma;        
K = param.K;          
mu = param.mu;    
deg = param.deg;        
kmu = param.kmu;   
M = param.M;
omega = param.omega;
tau = param.tau;


adj_h = assemble_dif(dif,adj0,M);

r = r0(1:M); Omega = r0(M+1); Delta = r0(M+2:2*M);

% limit of double-precision numbers ~1e-15
tol = 1e-15;
found = 0;
options = optimoptions('fsolve','Display','off','MaxIterations',1000,'MaxFunctionEvaluations',100000);
for count = 1:10
    		[sol,~,~,~] = fsolve(@(x)limit_cycle_sol(x,omega,gamma,M,K,tau,lambda,kmu,adj_h),[r;Omega;Delta],options);

            if norm(sol(1:M)-r)<sqrt(M)*1e-1 && abs(sol(M+1)- Omega)<0.2
                found = 1;
                r = sol(1:M); Omega = sol(M+1); delta = [0;sol(M+2:2*M)]; Delta = sol(M+2:2*M);    
                break
            else
                found = 2;
            end
end

if found == 1 % sync state identified successfully
    [J1,J2] = Jac_SLtau_hetnet(param,adj_h,r,Omega,delta,tau);
    mtle = dde_rightmost_eig(J1,J2,tau,M);
else % no sync state found
    mtle = Inf;
end

r0_new = zeros(2*M,1);
r0_new(1:M) = r; r0_new(M+1) = Omega; r0_new(M+2:2*M) = Delta;

end

function [adj_h] = assemble_dif(dif,adj0,M)

adjdif = zeros(M,M);
idx = 1;
for i = 1:M
    for j = 1:M
        adjdif(i,j) = dif(idx);
        idx = idx+1;
    end
end

adj_h = adj0 + adjdif;
end

function [mtle] = ObjAdj(dif,adj0,r0,param)
% Inputs:
%        dif ---- vector with nondiagonal perturbation
%        adj0   ---- initial (unperturbed) matrix
%        r0     ---- sync sol for the initial matrix
% Outputs:
%        mtle   ---- Lyapunov exponent of solution of perturbed matrix
%        r0     ---- sync sol for the perturbed matrix

lambda = param.lambda;      
gamma = param.gamma;  
omega = param.omega;
K = param.K;          
mu = param.mu;    
deg = param.deg;        
kmu = param.kmu;   
M = param.M; 
tau = param.tau;

adj_h = assemble_dif(dif,adj0,M);

r = r0(1:M); Omega = r0(M+1); Delta = r0(M+2:2*M);

% limit of double-precision numbers ~1e-15
tol = 1e-15;
found = 0;
options = optimoptions('fsolve','Display','off','MaxIterations',1000,'MaxFunctionEvaluations',100000);
for count = 1:10
    		[sol,~,~,~] = fsolve(@(x)limit_cycle_sol(x,omega,gamma,M,K,tau,lambda,kmu,adj_h),[r;Omega;Delta],options);

            if norm(sol(1:M)-r)<sqrt(M)*1e-2 && abs(sol(M+1)- Omega)<0.2
                found = 1;
                r = sol(1:M); Omega = sol(M+1); delta = [0;sol(M+2:2*M)]; Delta = sol(M+2:2*M);    
                break
            else
                found = 2;
            end
%             found = 1;
%             r = sol(1:M); Omega = sol(M+1); delta = [0;sol(M+2:2*M)]; 
end

if found == 1 % sync state identified successfully
    [J1,J2] = Jac_SLtau_hetnet(param,adj_h,r,Omega,delta,tau);
    mtle = dde_rightmost_eig(J1,J2,tau,M);
else % no sync state found
    mtle = Inf;
end

% r0(1:M) = r; r0(M+1) = Omega; r0(M+2:2*M) = Delta;

end

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