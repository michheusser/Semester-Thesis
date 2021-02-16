close all
clear all
clc

disp('Loading data...')
load('Data_Struct.mat');
disp('Data_Struct.mat loaded.');
load('Data_Matrix.mat');
disp('Data_Matrix.mat loaded');
load('Data_Matrix_Length.mat');
disp('Data_Matrix_Length.mat loaded');
disp('Finished loading files.')

Number_of_datasets = 18;


%% Visualisation - Fourier Transform


DTA = Data.FB.demian.Noise; %Data to analyse
N = length(DTA(:,1));
freq = (2*pi/N)*(0:1:ceil(N/2));  %rad/s
Fs = 44.1e3; %Hz
t = (0:1:N-1)/Fs;

U = (1/N)*abs(fft(DTA(:,1)));
Y = (1/N)*abs(fft(DTA(:,2)));
G = Y./U;

f1 = figure(1);
subplot(5,1,1)
plot(t,DTA(:,1))
axis tight
title('Input (Time-Domain)')
xlabel('t [s]')
ylabel('u(t)')
subplot(5,1,2)
plot(t,DTA(:,2))
axis tight
title('Output (Time-Domain)')
xlabel('Time [s]')
ylabel('y(t)')
subplot(5,1,3)
plot(freq,U(1:length(freq)))
axis tight
title('Input (Frequency-Domain)')
xlabel('\omega [rad/s]')
ylabel('|U(f)|')
subplot(5,1,4)
plot(freq,Y(1:length(freq)))
axis tight
title('Output (Frequency-Domain)')
xlabel('\omega [rad/s]')
ylabel('|Y(f)|')
subplot(5,1,5)
plot(freq,G(1:length(freq)))
axis tight
title('Empirical Transfer-function (|G(f)|) Estimate')
xlabel('\omega [rad/s]')
ylabel('|G(f)|')


%Measurements
%     1   -   RFFS1
%     2   -   RFFS2
%     3   -   RFFN
%     4   -   RFBS1
%     5   -   RFBS2
%     6   -   RFBN
%     7   -   PFFS1
%     8   -   PFFS2
%     9   -   PFFN
%     10  -   PFBS1
%     11  -   PFBS2
%     12  -   PFBN
%     13  -   DFFS1
%     14  -   DFFS2
%     15  -   DFFN
%     16  -   DFBS1
%     17  -   DFBS2
%     18  -   DFBN


%% Estimates of G

f2 = figure(2);
Data_Struct_Cell = cell(1,Number_of_datasets);
                    
for k = 1:Number_of_datasets;
    
    DTA = Data_Matrix(1:Data_Matrix_Length(k),2*k-1:2*k);
    N = length(DTA(:,1));
    freq = (2*pi/N)*(0:1:ceil(N/2)); %rad/sec
    Fs = 44.1e3; %Hz
    t = (0:1:N-1)/Fs;
    
    disp([num2str(k) '/' num2str(Number_of_datasets) ' FFT Calculation...'])
    U = (1/N)*abs(fft(DTA(1:N,1)));
    Y = (1/N)*abs(fft(DTA(1:N,2)));
    G = Y./U;
    Data_Struct_Cell{1,k} = struct('U',U,'Y',Y,'G',G);
    
    disp(['Calculation finished.'])
        
    subplot(6,3,k)
    plot(freq,G(1:length(freq)))
    axis tight
    title([num2str(k) '. (|G(f)|) Estimate'])
    xlabel('\omega [rad/s]')
    ylabel('|G(f)|')
    
end
   save('Data_Struct_Cell.mat','Data_Struct_Cell');

    %UYG = struct('M1',Data_Struct_Cell{1,1})

    

%% Averaging G

% N = max(Data_Matrix_Length);
% G_average = zeros(N,1);
% 
% p=0;
% for k = 1 : N
% 
%     f = (k-1)*2*pi/N;
%     
% sum_U_Sq = 0;
% for l = 1 : Number_of_datasets
% 
%     U = Data_Struct_Cell{1,l}.U;
%     sum_U_Sq = sum_U_Sq + abs(nearestval(U,f))^2;    
%     
% end
% 
% 
% for l = 1 : Number_of_datasets
% 
%     U = Data_Struct_Cell{1,l}.U;
%     G = Data_Struct_Cell{1,l}.G;
%     
%     alpha = abs(nearestval(U,f))^2/sum_U_Sq;
%     
%  G_average(k) = G_average(k) + nearestval(G,f)*alpha;
% 
% end
% 
% 
% disp([num2str(p) '% Completed'])
% 
% if(p ~= floor(k/N*100*10)/10)
% p = floor(k/N*100*10)/10;
% disp([num2str(p) '% Completed'])
% end
% 
% 
% end

%% Smoothing Estimates
disp('Loading data...');
load('Data_Struct_Cell.mat');
load('Data_Struct.mat')
disp('Data loaded.')
%%
u = Data.FF.ricardo.Sweep1(:,1);
y = Data.FF.ricardo.Sweep1(:,2);

k = 1;
U = Data_Struct_Cell{1,k}.U;
Y = Data_Struct_Cell{1,k}.Y;
G = Data_Struct_Cell{1,k}.G;

N = length(U);
omega = (2*pi/N)*(0:1:N-1);

%Filters (Fourier Transformed)
w_Bartlett = @(t,g) 1 - abs(t)/g;
w_Hamming = @(t,g) 0.5*(1 + cos(pi*t/g));

%Correlations
gamma = 10;
delta = 5; 

%Psi Estimates
Psi_u = zeros(N,1);
disp('Calculating Psi_u...')
for kk = 1 : N
for k = 1: delta*2+1
    
    tau = k-1-delta;
    R_u = autocorr(u,tau);
    Psi_u(kk) = Psi_u(kk) + w_Bartlett(tau,gamma)*...
        R_u*exp(-j*tau*omega(kk));
    
end
end

disp('Calculating Psi_yu...')
Psi_yu = zeros(N,1);
for kk = 1 : N
for k = 1: delta*2+1
    
    tau = k-1-delta;
    R_yu = correl(y,u,tau);
    Psi_yu(kk) = Psi_yu(kk) + w_Bartlett(tau,gamma)*...
        R_yu*exp(-j*tau*omega(kk));
    
end
end

G_smooth = Psi_yu./Psi_u;
freq = (2*pi/N)*(0:1:ceil(N/2));


subplot(2,1,1)
plot(omega,abs(G(1:length(omega))))
axis tight
title('(|G(f)|) Estimate')
xlabel('\omega [rad/s]')
ylabel('|G(f)|')

subplot(2,1,2)
plot(omega,abs(G_smooth(1:length(omega))))
axis tight
title('(|G(f)|) Smoothed Estimate')
xlabel('\omega [rad/s]')
ylabel('|G(f)| Smoothed')

















