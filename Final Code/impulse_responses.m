
impulse(G_n4sid_sys_FB,3)
[y,t] = impulse(G_n4sid_sys_FB,3);
[y_max, I_max] = max(y)
grid on
title('Impulse Response of FB model (3 Seconds)')
text(t(I_max),y(I_max),['Amplitude: ' num2str(y_max)])
text(t(end),y(end),['Amplitude: ' num2str(y(end))])
hold on
plot(t(I_max),y(I_max),'LineStyle','None','Marker','o','MarkerEdgeColor','r','MarkerSize',7)
plot(t(end),y(end),'LineStyle','None','Marker','o','MarkerEdgeColor','r','MarkerSize',7)
hold off