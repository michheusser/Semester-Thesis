close all
clear all
clc

load('G.mat'); %Loads smoothed transfer function and frequencies in rad

Fs = 44.1e3; %Hz
f_min = 20; %Hz
f_max = 5e3; %Hz

freq_herz = freq*Fs/(2*pi); % Converts frequencies in Hz
freq_relevant = freq_herz(freq_herz >= f_min & freq_herz <= f_max); %Saves only relevant frequencies (20 Hz to 5 Hz)
G_model_smooth_av_relevant = G_model_smooth_av(freq_herz >= f_min & freq_herz <= f_max); %Saves only relevant datapoints of the transfer function


q = 200; %Controllability/Reachability parameters
r = 200;

W = zeros(q,length(G_model_smooth_av_relevant)); %Pre-allocation of W and G matrices
G = zeros(q,length(G_model_smooth_av_relevant));
for k = 1 : q
    k
    for kk = 1 : length(G_model_smooth_av_relevant)
        W(k,kk) = 1/sqrt(length(G_model_smooth_av_relevant))*exp(1i*freq_relevant(kk)*(k-1));
        G(k,kk) = W(k,kk)*G_model_smooth_av_relevant(kk);
    end
end

%I try to calculate the QR factorization of the G_special (the one in the
%weird font) and W_special_orthogonal (The one with the weird font and
%orthogonal superindex) shown between formulas 61 and 62 in McKelvey

W_sp = [real(W) imag(W)]; 
G_sp = [real(G) imag(G)];
W_orth = eye(2*length(G_model_smooth_av_relevant)) - transpose(W_sp)*((W_sp*transpose(W_sp))\W_sp);
[Q,R] = qr(transpose(G_sp*W_orth));
R_22 = transpose(R);

%Formula 60 in McKelvey
K = chol(real(W*W'),'lower');
s= svd(K\R_22)
figure(1)
plot(s)
xlim([0 200])