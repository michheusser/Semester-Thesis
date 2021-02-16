function [psy_yu, psy_u] = spect_filtered_freq(U_data,Y_data,freq_vec,gamma,delta,filter)

disp('Calculating filtered spectra')

N = length(freq_vec);

psy_yu = zeros(N,1);
psy_u = zeros(N,1);

tau = -delta : 1 : delta;
w = w_filter(tau,gamma,filter);
disp('Calculating cross-correlations...')
R_uy = cross_correl_freq(U_data,Y_data,tau);
R_u = cross_correl_freq(U_data,U_data,tau);


disp('Calculating filtered spectra...')
w_Ruy = w.*R_uy;
w_Ru = w.*R_u;

progress = 0;
for k = 1 : N
    freq = freq_vec(k);
    
    exp_term = exp(-1i*tau*freq)';
      
    psy_yu(k) = sum(w_Ruy.*exp_term);
    psy_u(k) = sum(w_Ru.*exp_term);
    
    new_progress = round(k*100/N/10)*10;
    if(new_progress > progress)
        progress = new_progress;
        disp([num2str(new_progress) '%']);
    end
    
end

disp('Calculations finished.')



