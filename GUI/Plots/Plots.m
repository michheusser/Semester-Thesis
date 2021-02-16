

freq1 = (Fs/length(G1))*(0:1:ceil(length(G1)/2));
freq2 = (Fs/length(G2))*(0:1:ceil(length(G2)/2));
freq3 = (Fs/length(G3))*(0:1:ceil(length(G3)/2));

%|G(w)|
subplot(2,1,1)
h = loglog(freq1,abs(G1(1:length(freq1))),freq2,abs(G2(1:length(freq2))),freq3,abs(G2(1:length(freq3))));
set(h(2),'LineWidth',1.1);
xlabel('Frequency [Hz]')
ylabel('|G(\omega)|')
grid on
grid minor
axis tight

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