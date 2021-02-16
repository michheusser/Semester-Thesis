
load('FB_ricardo_Noise_44100_Hamming_100_100_SA_10.mat')
G1 = G_sm;
load('FB_ricardo_Noise_44100_Hamming_200_200_SA_10.mat')
G2 = G_sm;
load('FB_ricardo_Noise_44100_Hamming_500_500_SA_10.mat')
G3 = G_sm;
%%
freq    =  (Fs/length(G))*(0:1:ceil(length(G)/2));
freq1 = (Fs/length(G1))*(0:1:ceil(length(G1)/2));
freq2 = (Fs/length(G2))*(0:1:ceil(length(G2)/2));
freq3 = (Fs/length(G3))*(0:1:ceil(length(G3)/2));

title_cell =     {  ['Comparison between \gamma'],...
                    ['Fs: ' num2str(Fs) ' Hz    Averaging: ' num2str(Averaging)], ...
                    ['Filter: ' Filter '   \gamma = ' num2str(Gamma) '   \delta = ' num2str(Delta)],...
                    ['Method: ' Method]}';



f1 = figure(1);
set(f1,'Position', [0 0 1000 500])
set(f1,'PaperPositionMode','auto')
set(f1,'PaperType','A4')
set(f1,'PaperOrientation','landscape');
                
%|G(w)|
subplot(2,1,1)
loglog(freq,abs(G(1:length(freq))),'Color',[0.6 0.6 0.6])
hold on
h = loglog(freq1,abs(G1(1:length(freq1))),freq2,abs(G2(1:length(freq2))),freq3,abs(G2(1:length(freq3))));
hold off
set(h(2),'LineWidth',1.1);
xlabel('Frequency [Hz]')
ylabel('|G(\omega)|')
grid on
grid minor
axis tight
legend('ETFE','\gamma = 100','\gamma = 200','\gamma = 500')
title(title_cell)

subplot(2,1,2)
semilogx(freq,angle(G(1:length(freq))),'Color',[0.6 0.6 0.6])
hold on
h = semilogx(freq1,angle(G1(1:length(freq1))),freq2,angle(G2(1:length(freq2))),freq3,angle(G2(1:length(freq3))));
set(h(2),'LineWidth',1.1);
xlabel('Frequency [Hz]')
ylabel('\angle G(\omega)')
grid on
grid minor
axis tight
legend('ETFE','\gamma = 100','\gamma = 200','\gamma = 500')

print(f1,'Comparison gamma','-dpdf')
disp('File Printed')

% %arg(G(w))
% subplot(8,1,8)
% h = semilogx(freq,angle(G_data(1:length(freq))),freq,angle(G_smooth(1:length(freq))),'r');
% set(h(2),'LineWidth',1.1);
% xlabel('Frequency [Hz]')
% ylabel('\angle G(\omega)')
% grid on
% grid minor
% axis tight
% 
% 
% 
% 
% title_cell =     {  [experiment '\_' user '\_' signal],...
%                     ['Fs: ' num2str(Fs) ' Hz    Averaging: ' num2str(n_sub)], ...
%                     ['Filter: ' filter '   \gamma = ' num2str(gamma) '   \delta = ' num2str(delta)],...
%                     ['Method: ' method]}';
% title(title_cell)
% 
% disp(title_cell)

% %y(t)
% subplot(8,1,2)
% plot(t,y_data(1:length(t)))
% xlabel('Time [s]')
% ylabel('y(t)')
% grid on
% grid minor
% axis tight
% 
% %|U(w)|
% subplot(8,1,3)
% loglog(freq,abs(U_data(1:length(freq))))
% xlabel('Frequency [Hz]')
% ylabel('|U(\omega)|')
% grid on
% grid minor
% axis tight
% 
% %arg(U(w))
% subplot(8,1,4)
% semilogx(freq,angle(U_data(1:length(freq))))
% xlabel('Frequency [Hz]')
% ylabel('\angle U(\omega)')
% grid on
% grid minor
% axis tight
% 
% %|Y(w)|
% subplot(8,1,5)
% loglog(freq,abs(Y_data(1:length(freq))))
% xlabel('Frequency [Hz]')
% ylabel(' |Y(\omega)|')
% grid on
% grid minor
% axis tight
% 
% %arg(Y(w))
% subplot(8,1,6)
% semilogx(freq,angle(Y_data(1:length(freq))))
% xlabel('Frequency [Hz]')
% ylabel('\angle Y(\omega)')
% grid on
% grid minor
% axis tight
% 
% %|G(w)|
% subplot(8,1,7)
% h = loglog(freq,abs(G_data(1:length(freq))),freq,abs(G_smooth(1:length(freq))),'r');
% set(h(2),'LineWidth',1.1);
% xlabel('Frequency [Hz]')
% ylabel('|G(\omega)|')
% grid on
% grid minor
% axis tight
% 
% %arg(G(w))
% subplot(8,1,8)
% h = semilogx(freq,angle(G_data(1:length(freq))),freq,angle(G_smooth(1:length(freq))),'r');
% set(h(2),'LineWidth',1.1);
% xlabel('Frequency [Hz]')
% ylabel('\angle G(\omega)')
% grid on
% grid minor
% axis tight