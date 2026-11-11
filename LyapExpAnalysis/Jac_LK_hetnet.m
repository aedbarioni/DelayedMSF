function [J1,J2] = Jac_LK_hetnet(param,adj,r,Omega,delta,N,tau)

%%%%--- Parameters ---%%%%%
N_0 = param.N_0;       
g = param.g;       
s = param.s;        
gamma = param.gamma;      
alpha0 = param.alpha0;       
gamma_n = param.gamma_n;    
kappa = param.kappa;      
M = param.M; 
phase_lag = delta;
neig = param.neig;

J1 = zeros(3*M);
J2 = zeros(3*M);

for kk = 1:M
    J1(kk,kk) = (1/2)*(g*((N(kk)-N_0)/(1+s*(r(kk))^2)^2)*(1-s*(r(kk))^2) - gamma);          % (\partial F1)/(\partial r)
    J1(kk,kk+M) = (kappa/neig)*sum_sin(M,adj,r,phase_lag,Omega,tau,kk,1);                % (\partial F1)/(\partial phi)
    J1(kk,kk+2*M) = (1/2)*((g*r(kk))/(1+s*(r(kk))^2));                                      % (\partial F1)/(\partial N)
    
    J1(kk+M,kk) = -alpha0*g*s*((N(kk)-N_0)/((1+s*(r(kk))^2)^2))*r(kk) - (kappa/neig)*sum_sin(M,adj,r,phase_lag,Omega,tau,kk,2);        % (\partial F2)/(\partial r)     
    J1(kk+M,kk+M) = -(kappa/neig)*sum_cos(M,adj,r,phase_lag,Omega,tau,kk);                                                                  % (\partial F2)/(\partial phi)
    J1(kk+M,kk+2*M) = (alpha0/2)*(g/(1+s*(r(kk))^2));                                                                                     % (\partial F2)/(\partial N)
    
    J1(kk+2*M,kk) = -2*g*r(kk)*((N(kk)-N_0)/((1+s*r(kk)^2)^2))*(1e-4);                 % (\partial F3)/(\partial r)
    J1(kk+2*M,kk+2*M) = -gamma_n - g*1e-4*((r(kk))^2/(1+s*r(kk)^2));                   % (\partial F3)/(\partial N)

    for kkk=1:M
        J2(kk,kkk) = (kappa/neig)*adj(kk,kkk)*cos(phase_lag(kkk)-phase_lag(kk)-Omega*tau);                       % (\partial F1)/(\partial r)
        J2(kk,kkk+M) = -(kappa/neig)*adj(kk,kkk)*r(kkk)*sin(phase_lag(kkk)-phase_lag(kk)-Omega*tau);             % (\partial F1)/(\partial phi)

        J2(kk+M,kkk) = (kappa/(neig*r(kk)))*adj(kk,kkk)*sin(phase_lag(kkk)-phase_lag(kk)-Omega*tau);                       % (\partial F2)/(\partial r)
        J2(kk+M,kkk+M) = (kappa/neig)*adj(kk,kkk)*(r(kkk)/r(kk))*cos(phase_lag(kkk)-phase_lag(kk)-Omega*tau);              % (\partial F2)/(\partial phi)
    end            
end

end


function sum = sum_sin(M,adj,r,phase_lag,Omega,tau,i,type)
    sum = 0;
    if type == 1
        for jkj=1:M
            sum = sum + adj(i,jkj)*r(jkj)*sin(phase_lag(jkj)-phase_lag(i)-Omega*tau);
        end
    elseif type == 2
        for jkj=1:M
            sum = sum + adj(i,jkj)*(r(jkj)/(r(i))^2)*sin(phase_lag(jkj)-phase_lag(i)-Omega*tau);
        end
    end
end

function sum = sum_cos(M,adj,r,phase_lag,Omega,tau,i)
    sum = 0;
    for kj=1:M
    sum = sum + adj(i,kj)*(r(kj)/r(i))*cos(phase_lag(kj)-phase_lag(i)-Omega*tau);
    end
end