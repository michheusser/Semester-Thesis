close all
clear all
clc
format long




%% Dataset loading
load('Data_Struct.mat');
experiment = 'FF'; %To chose: 'FF' or 'FB'
user = 'ricardo'; %To chose: 'ricardo', 'patricio' or 'demian'

u_data = Data.(experiment).(user).('Noise').Data(:,2);
y_data = Data.(experiment).(user).('Noise').Data(:,1);
Fs = Data.(experiment).(user).('Noise').Fs;

t = (0 : length(u_data)-1)/Fs;

%% Analysis Variables
Bode_Mode = 'Amplitude'
thinning_ratio = 1;
I_beg = 226138;
I_end = 1543430;

f_min = 20; % Hz
f_max = 5e3; % Hz
N_sample = Fs/f_min;
n_sub = floor((I_end - I_beg + 1)/N_sample);
method = 'Average - Smooth'; %To chose: 'Average - Smooth', 'Smooth - Average'
gamma = 400;
delta = gamma;
filter = 'Parzen'; %To chose: 'Hamming', 'Bartlett', 'Parzen', 'None'

%% Cropping
u_model = u_data(I_beg : I_end);
y_model = y_data(I_beg : I_end);

%% Thinning:
u_model = matrix_thinner(u_model,thinning_ratio);
y_model = matrix_thinner(y_model,thinning_ratio);
Fs = Fs/thinning_ratio;

%% Validation Model
val = 0;
N_model = floor((1-val)*length(u_model));
u_val = u_model(N_model + 1 : length(u_model));
y_val = y_model(N_model + 1: length(y_model));

u_model = u_model(1: N_model);
y_model = y_model(1: N_model);

%% Batch Generation
N_val = length(u_val);
n_sub_val = floor(N_val/N_sample);
u_model_batch = reshape(u_model(1:N_sample*n_sub),N_sample,n_sub);
y_model_batch = reshape(y_model(1:N_sample*n_sub),N_sample,n_sub);
u_val_batch = reshape(u_val(1 : N_sample*n_sub_val),N_sample,n_sub_val);
y_val_batch = reshape(y_val(1 : N_sample*n_sub_val),N_sample,n_sub_val);

%% Fourier Transform
U_model_batch = zeros(size(u_model_batch));
Y_model_batch = zeros(size(y_model_batch));
for k = 1 : n_sub
    U_model_batch(:,k) = (1/N_sample)*fft(u_model_batch(:,k));
    Y_model_batch(:,k) = (1/N_sample)*fft(y_model_batch(:,k));
end

U_model_av = zeros(N_sample,1);
Y_model_av = zeros(N_sample,1);
for k = 1 : n_sub
    U_model_av = U_model_av + U_model_batch(:,k);
    Y_model_av = Y_model_av + Y_model_batch(:,k);
end
U_model_av = (1/n_sub)*U_model_av;
Y_model_av = (1/n_sub)*Y_model_av;

U_val_batch = zeros(size(u_val_batch));
Y_val_batch = zeros(size(y_val_batch));
for k = 1 : n_sub_val
    U_val_batch(:,k) = (1/N_sample)*fft(u_val_batch(:,k));
    Y_val_batch(:,k) = (1/N_sample)*fft(y_val_batch(:,k));
end

U_val_av = zeros(N_sample,1);
Y_val_av = zeros(N_sample,1);

for k = 1 : n_sub_val
    U_val_av = U_val_av + U_val_batch(:,k);
    Y_val_av = Y_val_av + Y_val_batch(:,k);
end
U_val_av = (1/n_sub_val)*U_val_av;
Y_val_av = (1/n_sub_val)*Y_val_av;


%% TF Estimate
freq = (2*pi/N_sample)*(0:1:N_sample-1)';

G_model_batch = zeros(size(U_model_batch));
G_model_smooth_batch = zeros(size(u_model_batch));

if(strcmp(method,'Smooth - Average'))
    % Generate G_batch
    for k = 1 : n_sub
        G_model_batch(:,k) = Y_model_batch(:,k)./U_model_batch(:,k);
    end
    
    % Average G_batch
    G_model_av = zeros(N_sample,1);
    for k = 1 : n_sub
        G_model_av = G_model_av + G_model_batch(:,k);
    end
    G_model_av = (1/n_sub)*G_model_av;
    
    
    % Smoothing of every G
    for k = 1 : n_sub
        [Psy_yu, Psy_u] = spect_filtered(u_model_batch(:,k),y_model_batch(:,k),freq,delta,gamma,filter);
        G_model_smooth_batch(:,k) = Psy_yu./Psy_u;
        disp([num2str(k) '/' num2str(n_sub) ' Finished.'])
    end
    
    %Averaging of smoothed G
    G_model_smooth_av = zeros(N_sample,1);
    for k = 1 : n_sub
        G_model_smooth_av = G_model_smooth_av + G_model_smooth_batch(:,k);
    end
    G_model_smooth_av = (1/n_sub)*G_model_smooth_av;
                    
else
    %Calculate G_model out of averaged Y_model and U_model
    G_model_av = Y_model_av./U_model_av;
    
    %Smooth averaged U and Y
    [Psy_yu, Psy_u] = spect_filtered_freq(U_model_av,Y_model_av,freq,gamma,delta,filter);
    G_model_smooth_av = (Psy_yu./Psy_u);
end


%% Subspace Identification
freq_herz = freq*(Fs/(2*pi));
freq_herz_relevant = freq_herz(freq_herz >=0 & freq_herz <= f_max);
freq_relevant = freq_herz_relevant*2*pi/Fs;
%G_model_smooth_av_relevant = G_model_smooth_av(freq >=0 & freq <= f_max);

G_model_smooth_av_relevant = G_model_smooth_av(freq_herz >=0 & freq_herz <= f_max); 

%h = ifft(G_model_smooth_av_relevant); 
%H = hankel(h);
q = 50;
r = 50;
%H = H(1:q,1:r);
%SV = svd(H);
%log_SV = log(SV);
%semilogy(SV,'o')
%plot(SV,'o')
%xlim([0 30])
%ylim([0 max(SV)])
%axis('tight')

% 
%G_1 = ones(q,1)*G_model_smooth_av_relevant';
%G_2 = ones(q,1)*exp(j*freq_relevant');
%G_3 = (0:q-1)'*ones(1,length(G_model_smooth_av_relevant));

%G = 1/sqrt(length(G_model_smooth_av_relevant))*G_1.*(G_2.^G_3);
%W = 1/sqrt(length(G_model_smooth_av_relevant))*G_2.^G_3;

G = zeros(q,length(G_model_smooth_av_relevant));
W = zeros(q,length(G_model_smooth_av_relevant));

for k = 1 : q
    for kk = 1 : length(G_model_smooth_av_relevant)
        G(k,kk) = 1/sqrt(length(G_model_smooth_av_relevant))*G_model_smooth_av_relevant(kk)*exp(j*freq_relevant(kk)*(k-1));
    end
end

for k = 1 : q
    for kk = 1 : length(G_model_smooth_av_relevant)
        W(k,kk) = 1/sqrt(length(G_model_smooth_av_relevant))*exp(j*freq_relevant(kk)*(k-1));
    end
end


[Q,R_T] = qr([real(W) imag(W) ; real(G) imag(G)]);

R = R_T';
R_22 = R(length(G_model_smooth_av_relevant)+1 : 2*length(G_model_smooth_av_relevant),  q + 1 : 2*q  )';
%
%KKT_eig = eig(real(W*cov*conj(W')))
%K = chol(vpa(real(W*conj(W)'),40),'lower');

s = svd(R_22)
plot(s)

%% Validation

%Sweep 1
validation_signal = 'Sweep1';

u_val_1 = Data.(experiment).(user).(validation_signal).Data(:,2);
y_val_1 = Data.(experiment).(user).(validation_signal).Data(:,1);

U_val_1 = (1/(length(u_val_1)))*fft(u_val_1);
Y_val_1 = (1/(length(y_val_1)))*fft(y_val_1);

[Psy_yu, Psy_u] = spect_filtered_freq(U_model_av,Y_model_av,(2*pi/length(u_val_1))*(0:1:length(u_val_1)-1),gamma,delta,filter);
Y_val_est_1 = (Psy_yu./Psy_u).*U_val_1;

y_val_est_1 = ifft(length(u_val_1)*Y_val_est_1);


%Sweep 2
validation_signal = 'Sweep2';

u_val_2 = Data.(experiment).(user).(validation_signal).Data(:,2);
y_val_2 = Data.(experiment).(user).(validation_signal).Data(:,1);

U_val_2 = (1/(length(u_val_2)))*fft(u_val_2);
Y_val_2 = (1/(length(y_val_2)))*fft(y_val_2);

[Psy_yu, Psy_u] = spect_filtered_freq(U_model_av,Y_model_av,(2*pi/length(u_val_2))*(0:1:length(u_val_2)-1),gamma,delta,filter);
Y_val_est_2 = (Psy_yu./Psy_u).*U_val_2;

y_val_est_2 = ifft(length(u_val_2)*Y_val_est_2);


%% RMS

% I_freq_herz_crop = find(freq_herz >= f_min & freq_herz <= f_max);
% freq_herz_crop = freq_herz(I_freq_herz_crop);
RMS = sqrt(  1/(length(u_val_1) + length(u_val_2)) * (sum(abs(y_val_1 - y_val_est_1).^2) + sum(abs(y_val_2 - y_val_est_2).^2)));   


%% Plots

%Definition
f(1) = figure('Name','Bode plot: Empirical and smoothed estimates'); %Bode plot: G_ETFE and G_smoothed
f(2) = figure('Name','Real and estimated outputs (Sweep1.wav)'); %y_val_1 and y_val_est_1
f(3) = figure('Name','Real and estimated outputs (Sweep2.wav)'); %y_val_2 and y_val_est_2
f(4) = figure('Name','Measurement (Time Domain)','Position',[0 0 1300 400]);
f(5) = figure('Name','Measurement (Frequency Domain)','Position',[0 0 1300 400]);
f(6) = figure('Name','Parzen Window');

s(1,1) = subplot(2,1,1,'Parent',f(1)); %Amplitude
s(1,2) = subplot(2,1,2,'Parent',f(1)); %Phase

s(2,1) = subplot(2,1,1,'Parent',f(2)); %y_val_1
s(2,2) = subplot(2,1,2,'Parent',f(2)); %y_val_est_1

s(3,1) = subplot(2,1,1,'Parent',f(3)); %y_val_2
s(3,2) = subplot(2,1,2,'Parent',f(3)); %y_val_est_2

s(4,1) = subplot(2,1,1,'Parent',f(4)); %Input
s(4,2) = subplot(2,1,2,'Parent',f(4)); %Output

s(5,1) = subplot(2,1,1,'Parent',f(5)); %Input
s(5,2) = subplot(2,1,2,'Parent',f(5)); %Output

s(6,1) = axes('Parent',f(6));

%Plotting
loglog(s(1,1),freq_herz,abs(G_model_av),freq_herz,abs(G_model_smooth_av))
semilogx(s(1,2),freq_herz,angle(G_model_av)*180/pi,freq_herz,angle(G_model_smooth_av)*180/pi)

plot(s(2,1),(0 : length(u_val_1)-1)/Fs,y_val_1)
plot(s(2,2),(0 : length(u_val_1)-1)/Fs,y_val_est_1)

plot(s(3,1),(0 : length(u_val_2)-1)/Fs,y_val_2)
plot(s(3,2),(0 : length(u_val_2)-1)/Fs,y_val_est_2)

u_to_show = Data.('FF').('ricardo').('Sweep2').Data(:,2);
y_to_show = Data.('FF').('ricardo').('Sweep2').Data(:,1);
plot(s(4,1), (0 : length(u_to_show)-1)/Fs,u_to_show)
plot(s(4,2), (0 : length(u_to_show)-1)/Fs,y_to_show)

U_to_show = 1/length(u_to_show)*fft(u_to_show);
Y_to_show = 1/length(y_to_show)*fft(y_to_show);

loglog(s(5,1),Fs/length(U_to_show)*(0:length(U_to_show)-1),abs(U_to_show),[4e3 4e3],[min(abs(U_to_show)) max(abs(U_to_show))],'r')
loglog(s(5,2),Fs/length(Y_to_show)*(0:length(Y_to_show)-1),abs(Y_to_show),[4e3 4e3],[min(abs(Y_to_show)) max(abs(Y_to_show))],'r')

freq_parzen = (-1000 : 1000)*pi/1000;
plot(s(6,1),freq_parzen,w_filter_freq(freq_parzen,5,'Parzen'))
set(s(6,1),'XTick',[-pi -pi/2 0 pi/2 pi])
set(s(6,1),'XLim',[-pi pi])
set(s(6,1),'XTickLabel',{})
xlabel(s(6,1),'\omega [rad/s]')
ylabel(s(6,1),'W_\gamma(\omega)')
title(s(6,1),'Parzen Frequency Window')

%Properties
for k = 1 : 5
    
    for kk = 1 : 2
        set(s(k,kk),'XGrid','on')
        set(s(k,kk),'YGrid','on')
        set(s(k,kk),'XMinorGrid','on')
        set(s(k,kk),'YMinorGrid','on')
        axis(s(k,kk),'tight')
    end
end

title(s(1,1),'|G(e^{j\omega})|')
xlabel(s(1,1),'Frequency [Hz]')
ylabel(s(1,1),'Amplitude [dB]')
set(s(1,1),'XTick',[20:10:90 100:100:900 1000:1000:f_max])
set(s(1,1),'XTickLabel',{'20';'';'';'';'';'';'';'';'100';'';'';'';'';'';'';'';'';'1000';'';'';'';'5000'})
xlim(s(1,1),[f_min f_max])
legend(s(1,1),'Empirical','Smoothed','Location','Southwest')

title(s(1,2),'\angleG(e^{j\omega})')
xlabel(s(1,2),'Frequency [Hz]')
ylabel(s(1,2),'Phase [°]')
set(s(1,2),'XTick',[20:10:90 100:100:900 1000:1000:f_max])
set(s(1,2),'XTickLabel',{'20';'';'';'';'';'';'';'';'100';'';'';'';'';'';'';'';'';'1000';'';'';'';'5000'})
xlim(s(1,2),[f_min f_max])
legend(s(1,2),'Empirical','Smoothed','Location','Southwest')

title(s(2,1),'Real')
xlabel(s(2,1),'Time [s]')
ylabel(s(2,1),'Voltage')

title(s(2,2),'Estimated')
xlabel(s(2,2),'Time [s]')
ylabel(s(2,2),'Voltage')

title(s(3,1),'Real')
xlabel(s(3,1),'Time [s]')
ylabel(s(3,1),'Voltage')

title(s(3,2),'Estimated')
xlabel(s(3,2),'Time [s]')
ylabel(s(3,2),'Voltage')

title(s(4,1),'Chin Canal Signal')
xlabel(s(4,1),'Time [s]')
ylabel(s(4,1),'Voltage')

title(s(4,2),{'','Ear Canal Signal'})
xlabel(s(4,2),'Time [s]')
ylabel(s(4,2),'Voltage')

title(s(5,1),'|U(e^{j\omega})|')
xlabel(s(5,1),'Frequency [Hz]')
ylabel(s(5,1),'Amplitude [dB]')
%set(s(5,1),'XTick',[20:10:90 100:100:900 1000:1000:f_max])
%set(s(5,1),'XTickLabel',{'20';'';'';'';'';'';'';'';'100';'';'';'';'';'';'';'';'';'1000';'';'';'';'5000'})
xlim(s(5,1),[0 Fs/2])


title(s(5,2),'|Y(e^{j\omega})|')
xlabel(s(5,2),'Frequency [Hz]')
ylabel(s(5,2),'Amplitude [dB]')
%set(s(5,1),'XTick',[20:10:90 100:100:900 1000:1000:f_max])
%set(s(5,1),'XTickLabel',{'20';'';'';'';'';'';'';'';'100';'';'';'';'';'';'';'';'';'1000';'';'';'';'5000'})
xlim(s(5,2),[0 Fs/2])


% 
% xlabel(s)
%     title_G = '|G(\omega)|';
%     
%       
%      title_G = '\angle G(\omega)';
%    
% 
% 
% xlabel(axes1,'Time [s]')
% title(axes1,'u(t)')
% 
% set(axes2,'XGrid','on')
% set(axes2,'YGrid','on')
% set(axes2,'XMinorGrid','on')
% set(axes2,'YMinorGrid','on')
% axis(axes2,'tight')
% xlabel(axes2,'Time [s]')
% title(axes2,'y(t)')
% 
% set(axes3,'XGrid','on')
% set(axes3,'YGrid','on')
% set(axes3,'XMinorGrid','on')
% set(axes3,'YMinorGrid','on')
% axis(axes3,'tight')
% xlabel(axes3,'Frequency [Hz]')
% title(axes3,title_U)
% xlim(axes3,[x_min x_max])
% legend(axes3,'System','Validation','Location','Best')
% 
% set(axes4,'XGrid','on')
% set(axes4,'YGrid','on')
% set(axes4,'XMinorGrid','on')
% set(axes4,'YMinorGrid','on')
% axis(axes4,'tight')
% xlabel(axes4,'Frequency [Hz]')
% title(axes4,title_Y)
% xlim(axes4,[x_min x_max])
% legend(axes4,'System','Validation','Location','Best')
% 
% set(axes5,'XGrid','on')
% set(axes5,'YGrid','on')
% set(axes5,'XMinorGrid','on')
% set(axes5,'YMinorGrid','on')
% axis(axes5,'tight')
% xlabel(axes5,'Frequency [Hz]')
% title(axes5,title_G)
% xlim(axes5,[x_min x_max])
% legend(axes5,'Non-smoothed','Smoothed','Location','Best')