%function struct_thinner(struct_name,thinning_ratio)

close all
clear all
clc

struct_name = 'Data_Matrix.mat';
struct_name_2 = 'Data_Matrix_Length.mat';
thinning_ratio = 4;

load(struct_name);
load(struct_name_2);

[M N] = size(Data_Matrix);

Data_Matrix_thin = zeros(floor(M/thinning_ratio), N);

for k = 1 : length(Data_Matrix_thin)
    
    Data_Matrix_thin(k,:) = Data_Matrix(thinning_ratio*k,:);
    
end

Data_Matrix_Length_thin = floor(Data_Matrix_Length/thinning_ratio);

Data_Matrix = Data_Matrix_thin;
Data_Matrix_Length = Data_Matrix_Length_thin;

save('Data_Matrix_thin.mat','Data_Matrix');
save('Data_Matrix_Length_thin.mat','Data_Matrix_Length');
