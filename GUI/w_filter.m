function w = w_filter(tau_vec,gamma,type)

w = ones(length(tau_vec),1);
switch type
    case 'Bartlett'
        w = (1 - abs(tau_vec)/gamma)';
    case 'Parzen'
        for k = 1 : length(tau_vec)
            tau = tau_vec(k);
            if(abs(tau)<= gamma/2)
                w(k) = 1- 6 *(tau^2)/(gamma^2) * (1 - abs(tau)/gamma);
            elseif(abs(tau)<= gamma)
                w(k) = 2 * (1 - abs(tau)/gamma)^3;
            end
        end
        
    case 'Hamming'
        w = (0.5*(1+cos(pi*tau_vec/gamma)))';
    case 'None'
        w = 1;
        
end