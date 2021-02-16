close all
clear all
clc
format long


load('Data_Struct.mat');

%% Dataset loading


experiment_vec = {'FF', 'FB'};
I_beg_vec = [349318 353474];
I_end_vec = [1541386 1542509];

%for kkkk = 1 : 2
    kkkk = 2;
    
    experiment = experiment_vec{kkkk}; %To chose: 'FF' or 'FB'
    user = 'ricardo'; %To chose: 'ricardo', 'patricio' or 'demian'
    
    %Sweep 1
    validation_signal = 'Sweep1';
    
    u_val_1 = Data.(experiment).(user).(validation_signal).Data(:,2);
    y_val_1 = Data.(experiment).(user).(validation_signal).Data(:,1);
    
    U_val_1 = (1/(length(u_val_1)))*fft(u_val_1);
    Y_val_1 = (1/(length(y_val_1)))*fft(y_val_1);
    
    validation_signal = 'Sweep2';
    
    u_val_2 = Data.(experiment).(user).(validation_signal).Data(:,2);
    y_val_2 = Data.(experiment).(user).(validation_signal).Data(:,1);
    
    U_val_2 = (1/(length(u_val_2)))*fft(u_val_2);
    Y_val_2 = (1/(length(y_val_2)))*fft(y_val_2);
    
    
    
    
    u_data = Data.(experiment).(user).('Noise').Data(:,2);
    y_data = Data.(experiment).(user).('Noise').Data(:,1);
    Fs = Data.(experiment).(user).('Noise').Fs;
    
    t = (0 : length(u_data)-1)'/Fs;
    
    % Analysis Variables
    Bode_Mode = 'Amplitude';
    thinning_ratio = 1;
    I_beg = I_beg_vec(kkkk); %For ricardo FF: 349318, For ricardo FB: 353474 (ca. 3s transient)
    I_end = I_end_vec(kkkk); % For ricardo FF: 1541386, For ricardo FB: 1542509 (ca. 3s transient)
    
    f_min = 20; % Hz
    f_max = 5e3; % Hz
    N_sample = Fs/f_min;
    n_sub = floor((I_end - I_beg + 1)/N_sample);
    method = 'Average - Smooth'; %To chose: 'Average - Smooth', 'Smooth - Average'
    
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
    
    
    gamma_vec = [450 : 50 : 1000];
    
    
    %%
    %for kkk = 1 : length(gamma_vec)
        
        kkk = 1
        %gamma = gamma_vec(kkk)
        gamma = 45;
        delta = gamma;
        
        % TF Estimate
        freq = (2*pi/N_sample)*(0:1:N_sample-1)';
        freq_herz = freq*(Fs/(2*pi));
        
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
                [Psy_yu, Psy_u] = spect_filtered(u_model_batch(:,k),y_model_batch(:,k),freq',delta,gamma,filter);
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
            [Psy_yu, Psy_u] = spect_filtered_freq(U_model_av,Y_model_av,freq',gamma,delta,filter);
            G_model_smooth_av = (Psy_yu./Psy_u);
        end
        
        G_analize = G_model_smooth_av;
        
        
        Ts = 1/Fs;
        %n_order_vec = [1:20 25:5:50 60:10:100];
        n_order_vec = [1:1:10];
        
        freq = (2*pi/length(G_analize))*(0:length(G_analize)-1)';
        freq_herz = Fs/(2*pi)*freq;
        freq_relevant = freq(freq_herz >= f_min & freq_herz <= f_max);
        freq_herz_relevant = Fs/(2*pi)*freq_relevant;
        G_analize_relevant = G_analize(freq_herz >= f_min & freq_herz <= f_max);
        
        for kk = 1 : length(n_order_vec)
            %kk = 1;
            n_order = n_order_vec(kk)
            step1 = 1
            idfrd_Data = idfrd(G_analize_relevant,freq_herz_relevant,'Ts',Ts,'FrequencyUnit','Hz', 'TimeUnit','seconds');
            step2 = 1
            G_n4sid_sys = n4sid(idfrd_Data,n_order,'Ts',0,'DisturbanceModel','none','Display','on');
            step3 = 1
            
            G_est_1 = freqresp(G_n4sid_sys,(Fs/length(U_val_1))*(0:1:length(U_val_1)-1),'Hz');
            G_est_1 = G_est_1(:);
            step4 = 1
            Y_val_est_1 = G_est_1.*U_val_1;
            step5 = 1
            
%             G_est_2 = freqresp(G_n4sid_sys,(Fs/length(U_val_2))*(0:1:length(U_val_2)-1),'Hz');
%             step6 = 1
%             G_est_2 = G_est_2(:);
%             Y_val_est_2 = G_est_2.*U_val_2;
%             step7 = 1
            %RMS = sqrt(  1/(length(U_val_1) + length(U_val_2)) * (sum(abs(Y_val_1 - Y_val_est_1).^2) + sum(abs(Y_val_2 - Y_val_est_2).^2)));
            
            RMS = sqrt(  1/length(U_val_2) * (sum(abs(Y_val_1 - Y_val_est_1).^2)));
            
            RMS_vec(kkk,kk) = RMS;
            step8 = 1
        end
        
    %end
    
    
    save(['GAMMA_N_RMS_' experiment_vec{kkkk} '.mat'],'n_order_vec','RMS_vec', 'gamma_vec')
    
%end
