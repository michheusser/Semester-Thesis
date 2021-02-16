close all
clear all
clc
format long


load('Data_Struct.mat');

%% Dataset loading
%gamma_vec = [10:10:800]';
gamma_vec = [10:10:1000]';
n_order_vec = [1:20]';


    experiment = 'FF'; %To chose: 'FF' or 'FB'
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
    
    
    %Validation preparations
    validation_signal = 'Sweep1';
    
    
    u_val_1 = Data.(experiment).(user).(validation_signal).Data(:,2);
    y_val_1 = Data.(experiment).(user).(validation_signal).Data(:,1);
    
    U_val_1 = (1/(length(u_val_1)))*fft(u_val_1);
    Y_val_1 = (1/(length(y_val_1)))*fft(y_val_1);
    
    freq_1 = 2*pi/length(u_val_1)*(0:length(u_val_1)-1);
    freq_herz_1 = freq_1*(Fs/(2*pi));
    
    freq_relevant_1 = freq_1(freq_herz_1>=f_min & freq_herz_1 < f_max);
    freq_herz_relevant_1 = freq_herz_1(freq_herz_1>=f_min & freq_herz_1 < f_max);
    
    
    validation_signal = 'Sweep2';
    
    u_val_2 = Data.(experiment).(user).(validation_signal).Data(:,2);
    y_val_2 = Data.(experiment).(user).(validation_signal).Data(:,1);
    
    U_val_2 = (1/(length(u_val_2)))*fft(u_val_2);
    Y_val_2 = (1/(length(y_val_2)))*fft(y_val_2);
    
    freq_2 = 2*pi/length(u_val_2)*(0:length(u_val_2)-1);
    freq_herz_2 = freq_2*(Fs/(2*pi));
    
    freq_relevant_2 = freq_2(freq_herz_2>=f_min & freq_herz_2 < f_max);
    freq_herz_relevant_2 = freq_herz_2(freq_herz_2>=f_min & freq_herz_2 < f_max);
    
    STABILITY_MATRIX = zeros(length(gamma_vec),length(n_order_vec));
    RMS_MATRIX = NaN(length(gamma_vec),length(n_order_vec));
    
    for kk = 1 : length(gamma_vec)
        gamma = gamma_vec(kk);
        delta = gamma;
        [Psy_yu, Psy_u] = spect_filtered_freq(U_model_av,Y_model_av,freq',gamma,delta,filter);
        G_model_smooth_av = (Psy_yu./Psy_u);
        
        G_analize = G_model_smooth_av;
        % Subspace Identification
        
        Ts = 1/Fs;
        
        freq = (2*pi/length(G_analize))*(0:length(G_analize)-1)';
        freq_herz = Fs/(2*pi)*freq;
        freq_relevant = freq(freq_herz >= f_min & freq_herz <= f_max);
        freq_herz_relevant = Fs/(2*pi)*freq_relevant;
        G_analize_relevant = G_analize(freq_herz >= f_min & freq_herz <= f_max);
        
        idfrd_Data = idfrd(G_analize_relevant,freq_herz_relevant,'Ts',Ts,'FrequencyUnit','Hz', 'TimeUnit','seconds');
        
        opt = n4sidOptions('InitialState','zero');
        
        
        
        for kkk = 1 : length(n_order_vec)
            gamma
            n_order = n_order_vec(kkk)
            G_n4sid_sys = n4sid(idfrd_Data,n_order,'Ts',0,'DisturbanceModel','none',opt);
            
            if(sum(~(real(pole(G_n4sid_sys))<0)) == 0)
                STABILITY_MATRIX(kk,kkk) = 1;
                
                
                %Validation
                
                %Sweep 1
                              
                %[Psy_yu, Psy_u] = spect_filtered_freq(U_model_av,Y_model_av,freq_relevant_1',gamma,delta,filter);
                %G_est_1 = (Psy_yu./Psy_u);
                G_est_1 = freqresp(G_n4sid_sys,freq_herz_relevant_1,'Hz');
                                
                Y_val_est_1 = G_est_1(:).*U_val_1(freq_herz_1>=f_min & freq_herz_1 < f_max);
                
                Y_val_relevant_1 = Y_val_1(freq_herz_1>=f_min & freq_herz_1 < f_max);
                
                
                %Sweep 2
                                
                %[Psy_yu, Psy_u] = spect_filtered_freq(U_model_av,Y_model_av,freq_relevant_2',gamma,delta,filter);
                %G_est_2 = (Psy_yu./Psy_u);
                
                G_est_2 = freqresp(G_n4sid_sys,freq_herz_relevant_2,'Hz');
                Y_val_est_2 = G_est_2(:).*U_val_2(freq_herz_2>=f_min & freq_herz_2 < f_max);
                
                Y_val_relevant_2 = Y_val_2(freq_herz_2>=f_min & freq_herz_2 < f_max);
                
                %RMS
                RMS = sqrt(1/(length(Y_val_relevant_1) + length(Y_val_relevant_2)) * (sum(abs(Y_val_relevant_1 - Y_val_est_1).^2) + sum(abs(Y_val_relevant_2 - Y_val_est_2).^2)));
                
                
                
                RMS_MATRIX(kk,kkk) = RMS;

            end
            
        end
    end


%% Plots
%save('GAMMA_n_RMS_FF.mat','gamma_vec','n_order_vec','STABILITY_MATRIX','RMS_MATRIX')

%%
 [RMS_MIN,ind] = min(RMS_MATRIX(:));
 [m,n] = ind2sub(size(RMS_MATRIX),ind)
n_order_opt = n_order_vec(n)
gamma_opt = gamma_vec(m)

surf(n_order_vec,gamma_vec,RMS_MATRIX)
