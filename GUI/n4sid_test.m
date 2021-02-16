G = G_model_smooth_av;
N = length(G_model_smooth_av);

freq = (2*pi/N)*(0:1:N-1);
Fs = 44.1e3;
Ts = 1/Fs;
n_order = 1 : 10;

ResponseData(1,1,:) =  G;
idfrd_Data = idfrd(ResponseData,freq,Ts,'FrequencyUnit','rad/s')
sys = n4sid(idfrd_Data,n_order,'Ts',Ts)
%%
[num , den] = tfdata(sys);

G_est = zeros(length(freq),1);
for k = 1 : length(freq)
    G_est(k) = polyval(num{1,1},i*freq(k))/polyval(den{1,1},i*freq(k));
end

 f1 = figure(1);
    
 s1 = subplot(2,1,1)
    loglog(freq,abs(G_model_smooth_av),freq,abs(G_est))
    axis tight
    grid on
  s2 =  subplot(2,1,2)
    semilogx(freq,angle(G_model_smooth_av),freq, angle(G_est))
    axis tight
    grid on
  
f2 = figure(2)
bode(sys)