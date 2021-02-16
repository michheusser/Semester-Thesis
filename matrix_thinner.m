function output_matrix = matrix_thinner(input_matrix,thinning_ratio)

[M , N] = size(input_matrix);
output_matrix = zeros(floor(M/thinning_ratio), N);

for k = 1 : floor(M/thinning_ratio)
    
   output_matrix(k,:) = input_matrix(thinning_ratio*k,:);

end

