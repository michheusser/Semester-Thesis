function R_u = autocorr(u,tau)

R_u = 0;
N = length(u);

for k = 1 : N
    
    if((k-tau)>=1 & (k-tau)<=N)
    R_u = R_u + (1/N)*u(k)*u(k-tau);
    end
    
end

end