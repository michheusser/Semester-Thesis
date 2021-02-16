G = G_model_smooth_av;
N = length(G);

freq = (2*pi/N)*(0:1:N-1)';
freq_herz = (Fs/N)*(0:1:N-1)';
Fs = 44.1e3;
Ts = 1/Fs;
n_order = 20;

idfrd_Data = idfrd(G,freq,'Ts',Ts,'FrequencyUnit','rad/s', 'TimeUnit','seconds')
sys = n4sid(idfrd_Data,n_order,'Ts',Ts)

% frd_Data = frd(G,freq_herz,'Ts',-1,'FrequencyUnit','Hz')
% sys = n4sid(frd_Data,n_order)
% 
[num , den] = tfdata(sys);

G_est = zeros(length(freq),1);
for k = 1 : length(freq)
    G_est(k) = polyval(num{1,1},exp(i*freq(k)))/polyval(den{1,1},exp(i*freq(k)));
end

% sys = d2c(sys);
% [num , den] = tfdata(sys);
% 
% G_est = zeros(length(freq),1);
% for k = 1 : length(freq)
%     G_est(k) = polyval(num{1,1},i*freq(k))/polyval(den{1,1},i*freq(k));
% end



f1 = figure(1);
    
 s1 = subplot(2,1,1)
    loglog(freq_herz,abs(G_model_smooth_av),freq_herz,abs(G_est))
    %axis tight
    grid on
 xlim(s1,[f_min f_max])
    s2 =  subplot(2,1,2)
    semilogx(freq_herz,unwrap(angle(G_model_smooth_av))*180/pi,freq_herz, unwrap(angle(G_est))*180/pi)
    %axis tight
    grid on
  xlim(s2,[f_min f_max])
 
    
f3 = figure('Name','Bode Plots')
bode(sys)
tf(sys)