function [J1,J2] = Jac_SLtau_hetnet(param,adj,r,Omega,delta,tau)

%%%%--- Parameters ---%%%%%
lambda = param.lambda;      
gamma = param.gamma;        
K = param.K;          
mu = param.mu;    
deg = param.deg;        
kmu = param.kmu;   
M = param.M; 

% neig = M*ones(M,1);
% for iii = 1:M
%     for jjj = 1:M
%         if adj(iii,jjj)==0
%             neig(iii) = neig(iii)-1;
%         end
%     end
% end


phase_lag = delta;
phi = -Omega*tau;


J1 = zeros(2*M);
J2 = zeros(2*M);

for kk = 1:M
    J1(kk,kk) = lambda - 3*r(kk)^2 - kmu ;          % (\partial F1)/(\partial r)
    J1(kk,kk+M) = K*sum_sin(M,adj,r,phase_lag,phi,kk,1);                % (\partial F1)/(\partial phi)
    
    J1(kk+M,kk) = -2*gamma*r(kk) - K*sum_sin(M,adj,r,phase_lag,phi,kk,2);        % (\partial F2)/(\partial r)     
    J1(kk+M,kk+M) = -K*sum_cos(M,adj,r,phase_lag,phi,kk);                                                                  % (\partial F2)/(\partial phi)
        
    
    for kkk=1:M
        J2(kk,kkk) = K*adj(kk,kkk)*cos(phase_lag(kkk)-phase_lag(kk) + phi);                       % (\partial F1)/(\partial r)
        J2(kk,kkk+M) = -K*adj(kk,kkk)*r(kkk)*sin(phase_lag(kkk)-phase_lag(kk) + phi);             % (\partial F1)/(\partial phi)

        J2(kk+M,kkk) = K*adj(kk,kkk)*(1/r(kk))*sin(phase_lag(kkk)-phase_lag(kk) + phi);                       % (\partial F2)/(\partial r)
        J2(kk+M,kkk+M) = K*adj(kk,kkk)*(r(kkk)/r(kk))*cos(phase_lag(kkk)-phase_lag(kk) + phi);              % (\partial F2)/(\partial phi)
    end            
end

end

function sum = sum_sin(M,adj,r,phase_lag,phi,i,type)
    sum = 0;
    if type == 1
        for jkj=1:M
            sum = sum + adj(i,jkj)*r(jkj)*sin(phase_lag(jkj)-phase_lag(i) + phi);
        end
    elseif type == 2
        for jkj=1:M
            sum = sum + adj(i,jkj)*(r(jkj)/(r(i))^2)*sin(phase_lag(jkj)-phase_lag(i) + phi);
        end
    end
    
end

function sum = sum_cos(M,adj,r,phase_lag,phi,i)
    sum = 0;
    
    for jkj=1:M
        sum = sum + adj(i,jkj)*(r(jkj)/r(i))*cos(phase_lag(jkj)-phase_lag(i)-phi);
    end
   

    
end