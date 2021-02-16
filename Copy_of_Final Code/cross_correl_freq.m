function R_yu = cross_correl_freq(U,Y,tau_vec)
R_yu = zeros(length(tau_vec),1);
N = length(U);
freq = (2*pi/N)*(0:1:N-1)';

for k = 1 : length(tau_vec)
    tau = tau_vec(k);
    
    exp_term = exp(1i*tau*freq);
    R_yu(k) = (1/N)*sum(Y.*conj(U).*exp_term);
    
end

