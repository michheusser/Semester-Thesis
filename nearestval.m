function val = nearestval(G, search_freq);

N = length(G);
%available_frequencies = (0:1:N-1)*(2*pi/N);
k = mod(round(search_freq/(2*pi/N))+1,N);
val = G(k); 


end