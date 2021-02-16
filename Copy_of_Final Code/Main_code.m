close all
clear all
clc
format long


load('Data_Struct.mat');

%% Dataset loading
gamma_vec = [45 45]';
n_order = 3;

for kkk = 1 : length(gamma_vec)
    
    experiment = 'FB'; %To chose: 'FF' or 'FB'
    user = 'ricardo'; %To chose: 'ricardo', 'patricio' or 'demian'
    
    u_data = Data.(experiment).(user).('Noise').Data(:,2);
    y_data = Data.(experiment).(user).('Noise').Data(:,1);
    Fs = Data.(experiment).(user).('Noise').Fs;
    
    t = (0 : length(u_data)-1)'/Fs;
    
    % Analysis Variables
    Bode_Mode = 'Amplitude';
    thinning_ratio = 1;
    
    if(strcmp(experiment,'FF'))
    I_beg = 349318; %For ricardo FF: 349318, For ricardo FB: 353474 (ca. 3s transient)
    I_end = 1541386; % For ricardo FF: 1541386, For ricardo FB: 1542509 (ca. 3s transient)
    else
    I_beg = 353474; %For ricardo FF: 349318, For ricardo FB: 353474 (ca. 3s transient)
    I_end = 1542509; % For ricardo FF: 1541386, For ricardo FB: 1542509 (ca. 3s transient)
    end
    
    
    f_min = 20; % Hz
    f_max = 5e3; % Hz
    N_sample = Fs/f_min;
    n_sub = floor((I_end - I_beg + 1)/N_sample);
    
    gamma = gamma_vec(kkk);
    delta = gamma;
    filter = 'Parzen'; %To chose: 'Hamming', 'Bartlett', 'Parzen', 'None'
    
    % Cropping
    u_model = u_data(I_beg : I_end);
    y_model = y_data(I_beg : I_end);
    
    % Thinning:
    u_model = matrix_thinner(u_model,thinning_ratio);
    y_model = matrix_thinner(y_model,thinning_ratio);
    Fs = Fs/thinning_ratio;
    
    % Validation Model
    val = 0;
    N_model = floor((1-val)*length(u_model));
    u_val = u_model(N_model + 1 : length(u_model));
    y_val = y_model(N_model + 1: length(y_model));
    
    u_model = u_model(1: N_model);
    y_model = y_model(1: N_model);
    
    % Batch Generation
    N_val = length(u_val);
    n_sub_val = floor(N_val/N_sample);
    u_model_batch = reshape(u_model(1:N_sample*n_sub),N_sample,n_sub);
    y_model_batch = reshape(y_model(1:N_sample*n_sub),N_sample,n_sub);
    u_val_batch = reshape(u_val(1 : N_sample*n_sub_val),N_sample,n_sub_val);
    y_val_batch = reshape(y_val(1 : N_sample*n_sub_val),N_sample,n_sub_val);
    
    % Fourier Transform
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
    
    
    % TF Estimate
    freq = (2*pi/N_sample)*(0:1:N_sample-1)';
    freq_herz = freq*(Fs/(2*pi));
    
    G_model_batch = zeros(size(U_model_batch));
    G_model_smooth_batch = zeros(size(u_model_batch));
    
    
    %Calculate G_model out of averaged Y_model and U_model
    G_model_av = Y_model_av./U_model_av;
    
    %Smooth averaged U and Y
    [Psy_yu, Psy_u] = spect_filtered_freq(U_model_av,Y_model_av,freq',gamma,delta,filter);
    G_model_smooth_av = (Psy_yu./Psy_u);
    
    
    % Subspace Identification
    
    if(kkk == 1)
        G_analize = G_model_av;
        G_estimates(:,kkk) = G_model_av;
    else
        G_analize = G_model_smooth_av;
        G_estimates(:,kkk) = G_model_smooth_av;
    end
    
    freq_herz = freq*Fs/(2*pi); % Converts frequencies in Hz
    freq_relevant = freq_herz(freq_herz >= f_min & freq_herz <= f_max); %Saves only relevant frequencies (20 Hz to 5 Hz)
    G_analize_relevant = G_analize(freq_herz >= f_min & freq_herz <= f_max); %Saves only relevant datapoints of the transfer function
    
    
    q = 100; %Controllability/Reachability parameters
    r = 100;
    
    W = zeros(q,length(G_analize_relevant)); %Pre-allocation of W and G matrices
    G = zeros(q,length(G_analize_relevant));
    for k = 1 : q
        for kk = 1 : length(G_analize_relevant)
            W(k,kk) = 1/sqrt(length(G_analize_relevant))*exp(1i*freq_relevant(kk)*(k-1));
            G(k,kk) = W(k,kk)*G_analize_relevant(kk);
        end
    end
    
    %I try to calculate the QR factorization of the G_special (the one in the
    %weird font) and W_special_orthogonal (The one with the weird font and
    %orthogonal superindex) shown between formulas 61 and 62 in McKelvey
    
    W_sp = [real(W) imag(W)];
    G_sp = [real(G) imag(G)];
    W_orth = eye(2*length(G_analize_relevant)) - transpose(W_sp)*((W_sp*transpose(W_sp))\W_sp);
    [Q,R] = qr(transpose(G_sp*W_orth));
    R_22 = transpose(R);
    
    %Formula 60 in McKelvey
    K = chol(real(W*W'),'lower');
    sigma(:,kkk) = svd(K\R_22);
    % figure(1)
    % plot(sigma)
    % xlim([0 200])
    % title(['\gamma = ' num2str(gamma)])
    
%N4Sid
  
    Ts = 1/Fs;
    
    freq = (2*pi/length(G_analize))*(0:length(G_analize)-1)'; 
    freq_herz = Fs/(2*pi)*freq;
    freq_relevant = freq(freq_herz >= f_min & freq_herz <= f_max); 
    freq_herz_relevant = Fs/(2*pi)*freq_relevant;
    G_analize_relevant = G_analize(freq_herz >= f_min & freq_herz <= f_max);
    
    idfrd_Data = idfrd(G_analize_relevant,freq_herz_relevant,'Ts',Ts,'FrequencyUnit','Hz', 'TimeUnit','seconds');
 
    %G_n4sid_sys = n4sid(idfrd_Data,n_order,'Ts',0);
    opt = n4sidOptions('InitialState','zero');
    G_n4sid_sys = n4sid(idfrd_Data,n_order,'Ts',0,'DisturbanceModel','none',opt);
    %G_est = freqresp(G_n4sid_sys,freq,'rad/s');
    G_est = freqresp(G_n4sid_sys,freq_herz,'Hz');
    
%     G_est = zeros(length(freq),1);
%     for k = 1 : length(freq_relevant)
%         G_est(k) = evalfr(G_n4sid_sys,exp(1j*freq_relevant(k)));
%     end
    G_estimates(:,kkk + length(gamma_vec)) = G_est;
        
    
end

impulse(G_n4sid_sys)
% Plots


%1
%f(1) = figure('Name','Bode plot: Empirical and smoothed estimates','Position',[0 0 2600 800]); %Bode plot: G_ETFE and G_smoothed
f(1) = figure(1);
%set(f(1),'Position',[0 0 1600 800]); %Bode plot: G_ETFE and G_smoothed

tol_unwrap = 1.6*pi;
s(1,1) = subplot(2,1,1,'Parent',f(1)); %Amplitude

loglog(s(1,1),freq_herz_relevant,abs(G_estimates(freq_herz >= f_min & freq_herz <= f_max,1)),'LineWidth',1.5,'Color', hsv2rgb(1*(1/length(gamma_vec)),0,0.2));
hold(s(1,1),'on')
%loglog(s(1,1),freq_herz_relevant,abs(G_estimates(freq_herz >= f_min & freq_herz <= f_max,length(gamma_vec) + 1)),'LineWidth',1.5,'Color', hsv2rgb(1*(1/length(gamma_vec)),0,0.2),'LineStyle','--');
for kkk = 2 : length(gamma_vec)
    loglog(s(1,1),freq_herz_relevant,abs(G_estimates(freq_herz >= f_min & freq_herz <= f_max,kkk)),'LineWidth',1.5,'Color', hsv2rgb(kkk*(1/length(gamma_vec)),1,0.7));
    loglog(s(1,1),freq_herz_relevant,abs(G_estimates(freq_herz >= f_min & freq_herz <= f_max,length(gamma_vec) + kkk)),'LineWidth',1.5,'Color', hsv2rgb(kkk*(1/length(gamma_vec)),1,0.7),'LineStyle','--');
end
hold(s(1,1),'off')


title(s(1,1),'|G(e^{j\omega})|')
xlabel(s(1,1),'Frequency [Hz]')
ylabel(s(1,1),'Gain')
set(s(1,1),'XTick',[20:10:90 100:100:900 1000:1000:f_max])
set(s(1,1),'XTickLabel',{'20';'';'';'';'';'';'';'';'100';'';'';'';'';'';'';'';'';'1000';'';'';'';'5000'})
%axis(s(1,1),'tight')
xlim(s(1,1),[f_min f_max])
%legend(s(1,1), 'Non-Smoothed',  ['N4SID: Non-Smoothed, n = ' num2str(n_order)], 'Smoothed: \gamma = 500', ['N4SID: \gamma = 500, n = ' num2str(n_order)], 'Smoothed: \gamma = 200',  ['N4SID: \gamma = 200, n = ' num2str(n_order)], 'Smoothed: \gamma = 500', ['N4SID: \gamma = 500, n = ' num2str(n_order)], 'Smoothed: \gamma = 1000', ['N4SID: \gamma = 1000, n = ' num2str(n_order)], 'Location','Southwest')
legend(s(1,1), 'Non-Smoothed',  ['Smoothed: \gamma =' num2str(gamma_vec(2))], ['N4SID: \gamma =' num2str(gamma_vec(2)) ', n = ' num2str(n_order)])

s(1,2) = subplot(2,1,2,'Parent',f(1)); %Phase

semilogx(s(1,2),freq_herz_relevant,angle(G_estimates(freq_herz >= f_min & freq_herz <= f_max,1))*180/pi,'LineWidth',1.5,'Color', hsv2rgb(1*(1/length(gamma_vec)),0,0.2));
hold(s(1,2),'on')
%semilogx(s(1,2),freq_herz_relevant,angle(G_estimates(freq_herz >= f_min & freq_herz <= f_max,length(gamma_vec) + 1))*180/pi,'LineWidth',1.5,'Color', hsv2rgb(1*(1/length(gamma_vec)),0,0.2),'LineStyle','--');
for kkk = 2 : length(gamma_vec)
    semilogx(s(1,2),freq_herz_relevant,angle(G_estimates(freq_herz >= f_min & freq_herz <= f_max,kkk))*180/pi,'LineWidth',1.5,'Color', hsv2rgb(kkk*(1/length(gamma_vec)),1,0.7));
    semilogx(s(1,2),freq_herz_relevant,angle(G_estimates(freq_herz >= f_min & freq_herz <= f_max,length(gamma_vec) + kkk))*180/pi,'LineWidth',1.5,'Color', hsv2rgb(kkk*(1/length(gamma_vec)),1,0.7),'LineStyle','--');
end
hold(s(1,2),'off')

title(s(1,2),'\angleG(e^{j\omega})')
xlabel(s(1,2),'Frequency [Hz]')
ylabel(s(1,2),'Phase [°]')
set(s(1,2),'XTick',[20:10:90 100:100:900 1000:1000:f_max])
set(s(1,2),'XTickLabel',{'20';'';'';'';'';'';'';'';'100';'';'';'';'';'';'';'';'';'1000';'';'';'';'5000'})

xlim(s(1,2),[f_min f_max])
%legend(s(1,2), 'Non-Smoothed',  ['N4SID: Non-Smoothed, n = ' num2str(n_order)], 'Smoothed: \gamma = 500', ['N4SID: \gamma = 500, n = ' num2str(n_order)], 'Smoothed: \gamma = 200',  ['N4SID: \gamma = 200, n = ' num2str(n_order)], 'Smoothed: \gamma = 500', ['N4SID: \gamma = 500, n = ' num2str(n_order)], 'Smoothed: \gamma = 1000', ['N4SID: \gamma = 1000, n = ' num2str(n_order)], 'Location','Southwest')
legend(s(1,2), 'Non-Smoothed',  ['Smoothed: \gamma =' num2str(gamma_vec(2))], ['N4SID: \gamma =' num2str(gamma_vec(2)) ', n = ' num2str(n_order)])

for kk = 1 : 2
    set(s(1,kk),'XGrid','on')
    set(s(1,kk),'YGrid','on')
    set(s(1,kk),'XMinorGrid','on')
    set(s(1,kk),'YMinorGrid','on')
end


%% Validation

%Sweep 1
validation_signal = 'Sweep1';

u_val_1 = Data.(experiment).(user).(validation_signal).Data(:,2);
y_val_1 = Data.(experiment).(user).(validation_signal).Data(:,1);

U_val_1 = (1/(length(u_val_1)))*fft(u_val_1);
Y_val_1 = (1/(length(y_val_1)))*fft(y_val_1);

freq_herz = (Fs/length(u_val_1))*(0:1:length(u_val_1)-1)';
G_est = freqresp(G_n4sid_sys,freq_herz,'Hz');
%[Psy_yu, Psy_u] = spect_filtered_freq(U_model_av,Y_model_av,(2*pi/length(u_val_1))*(0:1:length(u_val_1)-1),gamma,delta,filter);
Y_val_est_1 = G_est(:).*U_val_1;

y_val_est_1 = ifft(length(u_val_1)*Y_val_est_1);


%% Sweep 2
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



%%

%1
%f(1) = figure('Name','Bode plot: Empirical and smoothed estimates','Position',[0 0 2600 800]); %Bode plot: G_ETFE and G_smoothed
f(1) = figure(1);
set(f(1),'Position',[0 0 1400 600]); %Bode plot: G_ETFE and G_smoothed

tol_unwrap = 1.6*pi;
s(1,1) = subplot(2,1,1,'Parent',f(1)); %Amplitude

p(1,1) = loglog(s(1,1),freq_herz_relevant,abs(G_estimates(freq_herz >= f_min & freq_herz <= f_max,1)),'LineWidth',1.4,'Color', hsv2rgb(1*(1/length(gamma_vec)),0,0.4));
hold(s(1,1),'on')
for kkk = 2 : length(gamma_vec)
    p(1,kkk) = loglog(s(1,1),freq_herz_relevant,abs(G_estimates(freq_herz >= f_min & freq_herz <= f_max,kkk)),'LineWidth',1.7,'Color', hsv2rgb(kkk*(1/length(gamma_vec)),1,0.8));
end
hold(s(1,1),'off')

title(s(1,1),{'Feedback Model','|G(e^{j\omega})|'})
xlabel(s(1,1),'Frequency [Hz]')
ylabel(s(1,1),'Gain')
set(s(1,1),'XTick',[20:10:90 100:100:900 1000:1000:f_max])
set(s(1,1),'XTickLabel',{'20';'';'';'';'';'';'';'';'100';'';'';'';'';'';'';'';'';'1000';'';'';'';'5000'})
%axis(s(1,1),'tight')
xlim(s(1,1),[f_min f_max])
legend(s(1,1),'Non-Smoothed','Smoothed: \gamma = 100', 'Smoothed: \gamma = 200', 'Smoothed: \gamma = 500',  'Smoothed: \gamma = 1000',   'Location','Southwest')

s(1,2) = subplot(2,1,2,'Parent',f(1)); %Phase

p(2,1) = semilogx(s(1,2),freq_herz_relevant,angle(G_estimates(freq_herz >= f_min & freq_herz <= f_max,1))*180/pi,'LineWidth',1.4,'Color', hsv2rgb(1*(1/length(gamma_vec)),0,0.4));
hold(s(1,2),'on')
for kkk = 2 : length(gamma_vec)
   p(2,kkk) = semilogx(s(1,2),freq_herz_relevant,angle(G_estimates(freq_herz >= f_min & freq_herz <= f_max,kkk))*180/pi,'LineWidth',1.7,'Color', hsv2rgb(kkk*(1/length(gamma_vec)),1,0.8));
end
hold(s(1,2),'off')

title(s(1,2),'\angleG(e^{j\omega})')
xlabel(s(1,2),'Frequency [Hz]')
ylabel(s(1,2),'Phase [°]')
set(s(1,2),'XTick',[20:10:90 100:100:900 1000:1000:f_max])
set(s(1,2),'XTickLabel',{'20';'';'';'';'';'';'';'';'100';'';'';'';'';'';'';'';'';'1000';'';'';'';'5000'})

xlim(s(1,2),[f_min f_max])
legend(s(1,2),'Non-Smoothed','Smoothed: \gamma = 100', 'Smoothed: \gamma = 200', 'Smoothed: \gamma = 500',  'Smoothed: \gamma = 1000',  'Location','Southwest')

for kk = 1 : 2
    set(s(1,kk),'XGrid','on')
    set(s(1,kk),'YGrid','on')
    set(s(1,kk),'XMinorGrid','on')
    set(s(1,kk),'YMinorGrid','on')
end

set(p(1,1),'Marker','o')
set(p(2,1),'Marker','o')


%% Plots
close all

%1
%f(1) = figure('Name','Bode plot: Empirical and smoothed estimates','Position',[0 0 2600 800]); %Bode plot: G_ETFE and G_smoothed
f(1) = figure(1);
set(f(1),'Position',[0 0 1400 600]); %Bode plot: G_ETFE and G_smoothed

tol_unwrap = 1.5*pi;
s(1,1) = subplot(2,1,1,'Parent',f(1)); %Amplitude

loglog(s(1,1),freq_herz_relevant,abs(G_estimates(freq_herz >= f_min & freq_herz <= f_max,1)),'LineWidth',1.5,'Color', hsv2rgb(1*(1/length(gamma_vec)),0,0.2));
hold(s(1,1),'on')
%loglog(s(1,1),freq_herz_relevant,abs(G_estimates(freq_herz >= f_min & freq_herz <= f_max,length(gamma_vec) + 1)),'LineWidth',1.5,'Color', hsv2rgb(1*(1/length(gamma_vec)),0,0.2),'LineStyle','--');
for kkk = 2 : length(gamma_vec)
    loglog(s(1,1),freq_herz_relevant,abs(G_estimates(freq_herz >= f_min & freq_herz <= f_max,kkk)),'LineWidth',1.5,'Color', hsv2rgb(kkk*(1/length(gamma_vec)),1,0.7));
    loglog(s(1,1),freq_herz_relevant,abs(G_estimates(freq_herz >= f_min & freq_herz <= f_max,length(gamma_vec) + kkk)),'LineWidth',1.5,'Color', hsv2rgb(kkk*(1/length(gamma_vec)),1,0.7),'LineStyle','--');
end
hold(s(1,1),'off')


title(s(1,1),{'Feedforward model','|G(e^{j\omega})|'})
xlabel(s(1,1),'Frequency [Hz]')
ylabel(s(1,1),'Gain')
set(s(1,1),'XTick',[20:10:90 100:100:900 1000:1000:f_max])
set(s(1,1),'XTickLabel',{'20';'';'';'';'';'';'';'';'100';'';'';'';'';'';'';'';'';'1000';'';'';'';'5000'})
%axis(s(1,1),'tight')
xlim(s(1,1),[f_min f_max])
legend(s(1,1), 'Non-Smoothed', 'Smoothed: \gamma = 45', ['N4SID: \gamma = 45, n = ' num2str(n_order)], 'Smoothed: \gamma = 200',  ['N4SID: \gamma = 200, n = ' num2str(n_order)], 'Smoothed: \gamma = 500', ['N4SID: \gamma = 500, n = ' num2str(n_order)], 'Smoothed: \gamma = 1000', ['N4SID: \gamma = 1000, n = ' num2str(n_order)], 'Location','Southwest')

s(1,2) = subplot(2,1,2,'Parent',f(1)); %Phase

semilogx(s(1,2),freq_herz_relevant,unwrap(angle(G_estimates(freq_herz >= f_min & freq_herz <= f_max,1)),1.05*tol_unwrap)*180/pi,'LineWidth',1.5,'Color', hsv2rgb(1*(1/length(gamma_vec)),0,0.2));
hold(s(1,2),'on')
%semilogx(s(1,2),freq_herz_relevant,angle(G_estimates(freq_herz >= f_min & freq_herz <= f_max,length(gamma_vec) + 1))*180/pi,'LineWidth',1.5,'Color', hsv2rgb(1*(1/length(gamma_vec)),0,0.2),'LineStyle','--');
for kkk = 2 : length(gamma_vec)
    semilogx(s(1,2),freq_herz_relevant,unwrap(angle(G_estimates(freq_herz >= f_min & freq_herz <= f_max,kkk)),0.6*tol_unwrap)*180/pi,'LineWidth',1.5,'Color', hsv2rgb(kkk*(1/length(gamma_vec)),1,0.7));
    semilogx(s(1,2),freq_herz_relevant,unwrap(angle(G_estimates(freq_herz >= f_min & freq_herz <= f_max,length(gamma_vec) + kkk)),tol_unwrap)*180/pi,'LineWidth',1.5,'Color', hsv2rgb(kkk*(1/length(gamma_vec)),1,0.7),'LineStyle','--');
end
hold(s(1,2),'off')

title(s(1,2),'\angleG(e^{j\omega})')
xlabel(s(1,2),'Frequency [Hz]')
ylabel(s(1,2),'Phase [°]')
set(s(1,2),'XTick',[20:10:90 100:100:900 1000:1000:f_max])
set(s(1,2),'XTickLabel',{'20';'';'';'';'';'';'';'';'100';'';'';'';'';'';'';'';'';'1000';'';'';'';'5000'})

xlim(s(1,2),[f_min f_max])
legend(s(1,2), 'Non-Smoothed', 'Smoothed: \gamma = 45', ['N4SID: \gamma = 45, n = ' num2str(n_order)], 'Smoothed: \gamma = 200',  ['N4SID: \gamma = 200, n = ' num2str(n_order)], 'Smoothed: \gamma = 500', ['N4SID: \gamma = 500, n = ' num2str(n_order)], 'Smoothed: \gamma = 1000', ['N4SID: \gamma = 1000, n = ' num2str(n_order)], 'Location','Southwest')

for kk = 1 : 2
    set(s(1,kk),'XGrid','on')
    set(s(1,kk),'YGrid','on')
    set(s(1,kk),'XMinorGrid','on')
    set(s(1,kk),'YMinorGrid','on')
end


%% 2
f(2) = figure('Name','Real and estimated outputs (Sweep1.wav)','Position',[0 0 1000 400]); %y_val_1 and y_val_est_1

s(2,1) = subplot(2,1,1,'Parent',f(2)); %y_val_1
plot(s(2,1),(0 : length(u_val_1)-1)/Fs,y_val_1)
title(s(2,1),{'Real vs Estimated Outputs (Feedforward)','','Real'})
xlabel(s(2,1),'Time [s]')
ylabel(s(2,1),'Voltage')
axis(s(2,1),'tight')

s(2,2) = subplot(2,1,2,'Parent',f(2)); %y_val_est_1
plot(s(2,2),(0 : length(u_val_1)-1)/Fs,y_val_est_1,'Color',[0 0.5 0])
title(s(2,2),'Estimated')
xlabel(s(2,2),'Time [s]')
ylabel(s(2,2),'Voltage')
axis(s(2,2),'tight')

%% 3
f(3) = figure('Name','Real and estimated outputs (Sweep2.wav)'); %y_val_2 and y_val_est_2


s(3,1) = subplot(2,1,1,'Parent',f(3)); %y_val_2
plot(s(3,1),(0 : length(u_val_2)-1)/Fs,y_val_2)
title(s(3,1),'Real')
xlabel(s(3,1),'Time [s]')
ylabel(s(3,1),'Voltage')
axis(s(3,1),'tight')

s(3,2) = subplot(2,1,2,'Parent',f(3)); %y_val_est_2
plot(s(3,2),(0 : length(u_val_2)-1)/Fs,y_val_est_2)
title(s(3,2),'Estimated')
xlabel(s(3,2),'Time [s]')
ylabel(s(3,2),'Voltage')
axis(s(3,2),'tight')


%% 4
f(4) = figure('Name','Measurement (Time Domain)','Position',[0 0 1300 400]);
% u_to_show = Data.('FF').('ricardo').('Noise').Data(:,2);
% y_to_show = Data.('FF').('ricardo').('Noise').Data(:,1);

 u_to_show = Data.(experiment).('ricardo').('Sweep1').Data(:,2);
 y_to_show = y_val_est_1;

s(4,1) = subplot(2,1,1,'Parent',f(4)); %Input
plot(s(4,1), (0 : length(u_to_show)-1)/Fs,u_to_show)
title(s(4,1),'Headphone Canal Signal')
xlabel(s(4,1),'Time [s]')
ylabel(s(4,1),'Voltage')
axis(s(4,1),'tight')

s(4,2) = subplot(2,1,2,'Parent',f(4)); %Output
plot(s(4,2), (0 : length(u_to_show)-1)/Fs,y_to_show,'Color',[0 0.5 0])
title(s(4,2),{'','Ear Canal Signal'})
xlabel(s(4,2),'Time [s]')
ylabel(s(4,2),'Voltage')
axis(s(4,2),'tight')

U_to_show = 1/length(u_to_show)*fft(u_to_show);
Y_to_show = 1/length(y_to_show)*fft(y_to_show);

%% 5
f(5) = figure('Name','Measurement (Frequency Domain)','Position',[0 0 1300 400]);


s(5,1) = subplot(2,1,1,'Parent',f(5)); %Input
loglog(s(5,1),Fs/length(U_to_show)*(0:length(U_to_show)-1),abs(U_to_show))
title(s(5,1),{'Real vs Estimated Output (Frequency Domain)','Real |Y(e^{j\omega})|'})
xlabel(s(5,1),'Frequency [Hz]')
ylabel(s(5,1),'Amplitude [dB]')
%set(s(5,1),'XTick',[20:10:90 100:100:900 1000:1000:f_max])
%set(s(5,1),'XTickLabel',{'20';'';'';'';'';'';'';'';'100';'';'';'';'';'';'';'';'';'1000';'';'';'';'5000'})
axis(s(5,1),'tight')
xlim(s(5,1),[f_min f_max])

s(5,2) = subplot(2,1,2,'Parent',f(5)); %Output
loglog(s(5,2),Fs/length(Y_to_show)*(0:length(Y_to_show)-1),abs(Y_to_show),'Color',[0 0.5 0])
title(s(5,2),'Estimated |Y(e^{j\omega})|')
xlabel(s(5,2),'Frequency [Hz]')
ylabel(s(5,2),'Amplitude [dB]')
%set(s(5,1),'XTick',[20:10:90 100:100:900 1000:1000:f_max])
%set(s(5,1),'XTickLabel',{'20';'';'';'';'';'';'';'';'100';'';'';'';'';'';'';'';'';'1000';'';'';'';'5000'})
axis(s(5,2),'tight')
xlim(s(5,2),[f_min f_max])


%% 6
f(6) = figure('Name', 'Nyquist Plots of estimated systems');
s(6,1) = axes('Parent',f(6));
polar(s(6,1),(angle(G_model_av(freq_herz>=f_min & freq_herz <=f_max))),abs(G_model_av(freq_herz>=f_min & freq_herz <=f_max)))

% hold(s(1,1),'on')
% loglog(s(1,1),freq_herz,abs(G_model_smooth_av_100),'LineWidth',1.2,'Color',[0 0 0.7])
% loglog(s(1,1),freq_herz,abs(G_model_smooth_av_300),'LineWidth',1.2,'Color',[0 0.7 0])
% loglog(s(1,1),freq_herz,abs(G_model_smooth_av_500),'LineWidth',1.2,'Color',[0.7 0 0])
% %loglog(s(1,1),[f_max f_max],[min(abs([G_model_av ; G_model_smooth_av])) max(abs([G_model_av ; G_model_smooth_av]))],'Color','r', 'LineWidth',2)
% hold(s(1,1),'off')
% title(s(1,1),'|G(e^{j\omega})|')
% xlabel(s(1,1),'Frequency [Hz]')
% ylabel(s(1,1),'Amplitude [dB]')
% set(s(1,1),'XTick',[20:10:90 100:100:900 1000:1000:f_max])
% set(s(1,1),'XTickLabel',{'20';'';'';'';'';'';'';'';'100';'';'';'';'';'';'';'';'';'1000';'';'';'';'5000'})
% axis(s(1,1),'tight')
% xlim(s(1,1),[f_min f_max])
% legend(s(1,1),'Non-Smoothed','Smoothed: \gamma = 100','Smoothed: \gamma = 300','Smoothed: \gamma = 500','Location','Southwest')


%% 7
f(7) = figure('Name','Parzen Window');
freq_parzen = (-1000 : 1000)*pi/1000;

s(7,1) = axes('Parent',f(7));
plot(s(7,1),freq_parzen,w_filter_freq(freq_parzen,5,'Parzen'))
axis(s(7,1),'tight')
set(s(7,1),'XTick',[-pi -pi/2 0 pi/2 pi])
set(s(7,1),'XLim',[-pi pi])
set(s(7,1),'XTickLabel',{})
xlabel(s(7,1),'\omega [rad/s]')
ylabel(s(7,1),'W_\gamma(\omega)')
title(s(7,1),'Parzen Frequency Window')

%% 8

f(8) = figure('Name','Order reduction');
index = (1 : q )';

s(8,1) = axes('Parent',f(8));
curr_axes = s(8,1);

hold(curr_axes,'on')
for kkk = 1 : length(gamma_vec)
plot(curr_axes,index, sigma(:,kkk),'Color',hsv2rgb(kkk*(1/length(gamma_vec)),1,0.7));
end
hold(curr_axes,'off')

axis(curr_axes,'tight')
set(curr_axes,'XLim',[0 100])
xlabel(curr_axes,'Index of singular value')
ylabel(curr_axes,'Magnitude of singular value')
title(curr_axes,{'Singular Values (ordered in decreasing value)',['User: ' user ' (' experiment ')']})
legend(curr_axes,'\gamma = 50','\gamma = 100', '\gamma = 200', '\gamma = 500', '\gamma = 1000', 'Non-smoothed')

set(curr_axes,'XGrid','on')
set(curr_axes,'YGrid','on')
set(curr_axes,'XMinorGrid','on')
set(curr_axes,'YMinorGrid','on')



%% Properties
for k = 2 : 6
    
    for kk = 1 : 2
        set(s(k,kk),'XGrid','on')
        set(s(k,kk),'YGrid','on')
        set(s(k,kk),'XMinorGrid','on')
        set(s(k,kk),'YMinorGrid','on')
    end
end


%%

close all

%1
%f(1) = figure('Name','Bode plot: Empirical and smoothed estimates','Position',[0 0 2600 800]); %Bode plot: G_ETFE and G_smoothed
f(1) = figure(1);
set(f(1),'Position',[0 0 2600 800]); %Bode plot: G_ETFE and G_smoothed

tol_unwrap = 1.6*pi;
s(1,1) = subplot(2,1,1,'Parent',f(1)); %Amplitude

loglog(s(1,1),freq_herz_relevant,abs(G_estimates(freq_herz >= f_min & freq_herz <= f_max,1)),'LineWidth',1.5,'Color', hsv2rgb(1*(1/length(gamma_vec)),0,0.2));
hold(s(1,1),'on')
loglog(s(1,1),freq_herz_relevant,abs(G_estimates(freq_herz >= f_min & freq_herz <= f_max,length(gamma_vec) + 1)),'LineWidth',1.5,'Color', hsv2rgb(1*(1/length(gamma_vec)),0,0.2),'LineStyle','--');
for kkk = 2 : length(gamma_vec)
    loglog(s(1,1),freq_herz_relevant,abs(G_estimates(freq_herz >= f_min & freq_herz <= f_max,kkk)),'LineWidth',1.5,'Color', hsv2rgb(kkk*(1/length(gamma_vec)),1,0.7));
    loglog(s(1,1),freq_herz_relevant,abs(G_estimates(freq_herz >= f_min & freq_herz <= f_max,length(gamma_vec) + kkk)),'LineWidth',1.5,'Color', hsv2rgb(kkk*(1/length(gamma_vec)),1,0.7),'LineStyle','--');
end
hold(s(1,1),'off')


title(s(1,1),'|G(e^{j\omega})|')
xlabel(s(1,1),'Frequency [Hz]')
ylabel(s(1,1),'Gain')
set(s(1,1),'XTick',[20:10:90 100:100:900 1000:1000:f_max])
set(s(1,1),'XTickLabel',{'20';'';'';'';'';'';'';'';'100';'';'';'';'';'';'';'';'';'1000';'';'';'';'5000'})
%axis(s(1,1),'tight')
xlim(s(1,1),[f_min f_max])
legend(s(1,1), 'Non-Smoothed',  ['N4SID: Non-Smoothed, n = ' num2str(n_order)], 'Smoothed: \gamma = 500', ['N4SID: \gamma = 500, n = ' num2str(n_order)], 'Smoothed: \gamma = 200',  ['N4SID: \gamma = 200, n = ' num2str(n_order)], 'Smoothed: \gamma = 500', ['N4SID: \gamma = 500, n = ' num2str(n_order)], 'Smoothed: \gamma = 1000', ['N4SID: \gamma = 1000, n = ' num2str(n_order)], 'Location','Southwest')

s(1,2) = subplot(2,1,2,'Parent',f(1)); %Phase

semilogx(s(1,2),freq_herz_relevant,angle(G_estimates(freq_herz >= f_min & freq_herz <= f_max,1))*180/pi,'LineWidth',1.5,'Color', hsv2rgb(1*(1/length(gamma_vec)),0,0.2));
hold(s(1,2),'on')
semilogx(s(1,2),freq_herz_relevant,angle(G_estimates(freq_herz >= f_min & freq_herz <= f_max,length(gamma_vec) + 1))*180/pi,'LineWidth',1.5,'Color', hsv2rgb(1*(1/length(gamma_vec)),0,0.2),'LineStyle','--');
for kkk = 2 : length(gamma_vec)
    semilogx(s(1,2),freq_herz_relevant,angle(G_estimates(freq_herz >= f_min & freq_herz <= f_max,kkk))*180/pi,'LineWidth',1.5,'Color', hsv2rgb(kkk*(1/length(gamma_vec)),1,0.7));
    semilogx(s(1,2),freq_herz_relevant,angle(G_estimates(freq_herz >= f_min & freq_herz <= f_max,length(gamma_vec) + kkk))*180/pi,'LineWidth',1.5,'Color', hsv2rgb(kkk*(1/length(gamma_vec)),1,0.7),'LineStyle','--');
end
hold(s(1,2),'off')

title(s(1,2),'\angleG(e^{j\omega})')
xlabel(s(1,2),'Frequency [Hz]')
ylabel(s(1,2),'Phase [°]')
set(s(1,2),'XTick',[20:10:90 100:100:900 1000:1000:f_max])
set(s(1,2),'XTickLabel',{'20';'';'';'';'';'';'';'';'100';'';'';'';'';'';'';'';'';'1000';'';'';'';'5000'})

xlim(s(1,2),[f_min f_max])
legend(s(1,1), 'Non-Smoothed',  ['N4SID: Non-Smoothed, n = ' num2str(n_order)], 'Smoothed: \gamma = 500', ['N4SID: \gamma = 500, n = ' num2str(n_order)], 'Smoothed: \gamma = 200',  ['N4SID: \gamma = 200, n = ' num2str(n_order)], 'Smoothed: \gamma = 500', ['N4SID: \gamma = 500, n = ' num2str(n_order)], 'Smoothed: \gamma = 1000', ['N4SID: \gamma = 1000, n = ' num2str(n_order)], 'Location','Southwest')

for kk = 1 : 2
    set(s(1,kk),'XGrid','on')
    set(s(1,kk),'YGrid','on')
    set(s(1,kk),'XMinorGrid','on')
    set(s(1,kk),'YMinorGrid','on')
end



