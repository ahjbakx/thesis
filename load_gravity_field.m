clear all; close all; clc;

file = '/Users/aaron/thesis/Data/moon_gravity/gggrx_1200b_sha.tab.txt';

data = importdata(file);

Rm = data(1,1);
mu = data(1,2);
dmu = data(1,3);
degree = data(1,4);
order = data(1,5);

coefficients = data(2:end,1:6);

coefficients_new = zeros(size(coefficients));

index = 1;
for l = 0:1200
    l_rows = find(coefficients(:,2) == l);
    for i = 1:length(l_rows)
        coefficients_new(index,:) = coefficients(l_rows(i),:);
        index = index + 1;
    end
end

%save('Data/moon_gravity/SH_coefficients_GSHformat.mat', 'coefficients_new');

clear data