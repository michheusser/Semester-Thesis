function W = w_filter_freq(freq_vec,gamma,type)
N = length(freq_vec);
W = ones(N,1);

k_zero = find(freq_vec == 0);

switch type
    case 'Bartlett'
        
        W = (1/gamma)*(sin(gamma*freq_vec/2)./sin(freq_vec/2))';
        W(k_zero) = gamma;
        
    case 'Parzen'
        
        W = ((4/gamma^3)*(2 + cos(freq_vec)).*(sin(gamma*freq_vec/4)./sin(freq_vec/2)).^4)';
        W(k_zero) = gamma*3/4;
        
    case 'Hamming'
        
         %W = 0.5*sin((gamma+0.5)*freq_vec)./sin(freq_vec/2) + 0.25*sin((gamma+0.5)*(freq_vec-pi/5))./sin((freq_vec-pi/5)/2) + 0.25*sin((gamma+0.5)*(freq_vec+pi/5))./sin((freq_vec+pi/5)/2)
%         
%         D_1 = sin((gamma + 0.5)*freq_vec)./sin(freq_vec/2);
%         D_2 = sin((gamma + 0.5)*(freq_vec - pi/gamma))./sin((freq_vec - pi/gamma)/2);
%         D_3 = sin((gamma + 0.5)*(freq_vec + pi/gamma))./sin((freq_vec + pi/gamma)/2);
%         W = (0.5*D_1 + 0.25*D_2 + 0.25*D_3)';
%         
%         %W(k_zero) = gamma + 0.5;
        
    case 'None'
        W = 2*pi;
end

W = W/(2*pi);