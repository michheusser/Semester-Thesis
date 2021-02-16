function R_yu = cross_correl(u,y,tau_vec)
R_yu = zeros(length(tau_vec),1);
N = length(y);

for k = 1 : length(tau_vec)
    tau = tau_vec(k);
    
    u_shifted = circshift(u,tau);
    R_yu(k) = (1/N)*sum(u_shifted.*y);
    
end

