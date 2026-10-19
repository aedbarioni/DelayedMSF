function y = log_norm_det(x,J,k,P,n,AR,tau)
	%% characteristic equation for Floquet exponents
	y = log10(norm(det(J - k*P - (x(1)+x(2)*i)*eye(2*n) + k*AR*exp(-(x(1)+x(2)*i)*tau))));
end
