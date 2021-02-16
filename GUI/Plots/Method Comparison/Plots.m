
load('FF_ricardo_Noise_44100_Hamming_500_500_SA_50.mat')
G1sm = G_sm;
G1 = G;
load('FF_ricardo_Noise_44100_Hamming_500_500_AS_50.mat')
G2sm = G_sm;
G2 = G;
%%
freq1sm = (Fs/length(G1sm))*(0:1:ceil(length(G1sm)/2));
freq1 = (Fs/length(G1))*(0:1:ceil(length(G1)/2));
freq2sm = (Fs/length(G2sm))*(0:1:ceil(length(G2sm)/2));
freq2 = (Fs/length(G2))*(0:1:ceil(length(G2)/2));



title_cell =     {  ['Comparison between methods'],...
                    ['Fs: ' num2str(Fs) ' Hz    Averaging: ' num2str(Averaging)], ...
                    ['Filter: ' Filter '   \gamma = ' num2str(Gamma) '   \delta = ' num2str(Delta)]}';



f1 = figure(1);
set(f1,'Position', [0 0 1000 500])
set(f1,'PaperPositionMode','auto')
set(f1,'PaperType','A4')
set(f1,'PaperOrientation','landscape');
                
%|G(w)|
subplot(2,1,1)
h = loglog(freq1,abs(G1(1:length(freq1))),'Color',[0 0 0.3]);
hold on
loglog(freq2,abs(G2(1:length(freq2))),'Color',[0 0.3 0]);
hold off
set(h(1),'LineWidth',1.1);
xlabel('Frequency [Hz]')
ylabel('|G(\omega)|')
grid on
grid minor
axis tight
h = legend('Smooth - Average: ETFE','Average - Smooth: ETFE')
title(title_cell)

subplot(2,1,2)
h = semilogx(freq1,angle(G1(1:length(freq1))),'Color',[0 0 0.3]);
hold on
semilogx(freq2,angle(G2(1:length(freq2))),'Color',[0 0.3 0]);
hold off
set(h(1),'LineWidth',1.1);
xlabel('Frequency [Hz]')
ylabel('|G(\omega)|')
grid on
grid minor
axis tight
h = legend('Smooth - Average: ETFE','Average - Smooth: ETFE')

print(f1,'Comparison method ETFE','-dpdf')
disp('File Printed')

f1 = figure(1);
set(f1,'Position', [0 0 1000 500])
set(f1,'PaperPositionMode','auto')
set(f1,'PaperType','A4')
set(f1,'PaperOrientation','landscape');


%%

%|G(w)|
subplot(2,1,1)
h = loglog(freq1sm,abs(G1sm(1:length(freq1sm))),'Color',[0 0 0.3]);
hold on
loglog(freq2sm,abs(G2sm(1:length(freq2sm))),'Color',[0 0.3 0]);
hold off
set(h(1),'LineWidth',1.1);
xlabel('Frequency [Hz]')
ylabel('|G(\omega)|')
grid on
grid minor
axis tight
h = legend('Smooth - Average: ETFE','Average - Smooth: ETFE')
title(title_cell)

subplot(2,1,2)
h = semilogx(freq1sm,angle(G1sm(1:length(freq1sm))),'Color',[0 0 0.3]);
hold on
semilogx(freq2sm,angle(G2sm(1:length(freq2sm))),'Color',[0 0.3 0]);
hold off
set(h(1),'LineWidth',1.1);
xlabel('Frequency [Hz]')
ylabel('|G(\omega)|')
grid on
grid minor
axis tight
h = legend('Smooth - Average: Smoothed','Average - Smooth: Smoothed')

print(f1,'Comparison method Smoothed G','-dpdf')
disp('File Printed')


