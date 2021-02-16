function R_yu = correl(y,u,tau)

R_yu = 0;
N = length(u);

for k = 1 : N
    
    if((k-tau)>=1 & (k-tau)<=N)
    R_yu = R_yu + (1/N)*y(k)*u(k-tau);
    end
    
end

end