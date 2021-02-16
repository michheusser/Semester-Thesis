
gammavec = [1000 900 800 700 600 500 400 300 200 100 10];
RMS_abs_vec = zeros(length(gammavec),1);
RMS_angle_vec = zeros(length(gammavec),1);
u_model_batch = Data.Model.u_model_batch;
y_model_batch = Data.Model.y_model_batch;
U_model_av = Data.Model.U_model_av;
U_val_av = Data.Validation.U_val_av;
Y_val_av = Data.Validation.Y_val_av;
n_sub = 250;

for kk = 1 : length(gammavec)
    disp([num2str(kk) '/' num2str(length(gammavec))])
gamma = gammavec(kk);
delta = gamma;
filter = 'Hamming';

handles = [];

N_sample = length(U_model_av);
freq = (2*pi/N_sample)*(0:1:N_sample-1);

G_model_smooth_batch = zeros(size(u_model_batch));
    for k = 1 : n_sub
        [Psy_yu, Psy_u] = spect_filtered(u_model_batch(:,k),y_model_batch(:,k),freq,delta,gamma,filter,handles);
        G_model_smooth_batch(:,k) = Psy_yu./Psy_u;
        disp([num2str(k) '/' num2str(n_sub) ' Finished.'])
    end
    
    %Averaging of smoothed G
    G_model_smooth_av = zeros(N_sample,1);
    for k = 1 : n_sub
        G_model_smooth_av = G_model_smooth_av + G_model_smooth_batch(:,k);
    end
    G_model_smooth_av = (1/n_sub)*G_model_smooth_av;
    
    Y_est = G_model_smooth_av.*U_val_av;
    RMS_abs = sqrt((1/N_sample)*sum((abs(Y_est)-abs(Y_val_av)).^2));
    RMS_angle = sqrt((1/N_sample)*sum((angle(Y_est)-angle(Y_val_av)).^2));
    
    RMS_abs_vec(kk) = RMS_abs;
    RMS_angle_vec(kk) = RMS_angle;
    
    
end


Data = struct('RMS_abs_vec',RMS_abs_vec,'RMS_angle_vec',RMS_angle_vec,'gammavec',gammavec);

save('RMS_Gamma.mat','-struct','Data')

% subplot(2,1,1)
% plot(gammavec,RMS_abs_vec)
% subplot(2,1,2)
% plot(gammavec,RMS_angle_vec)