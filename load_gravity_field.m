clear all; close all; clc;

file = '/Users/aaron/thesis/Data/moon_gravity/sha.grgm1200b_rm1_1e1_sigma.txt';

data = importdata(file);

mu = 4.9028001224452998E+12;
Rm = 1.7380000000000000E+06;
degree = 1200;
order = 1200;

coefficients = data(2:end,:);

coefficients_new = zeros(size(coefficients));

index = 1;
for l = 0:1200
    l_rows = find(coefficients(:,2) == l);
    for i = 1:length(l_rows)
        coefficients_new(index,:) = coefficients(l_rows(i),:);
        index = index + 1;
    end
end

%save('/Users/aaron/thesis/Data/moon_gravity/SH_coefficients_GSHformat_1200B-RM1.mat', 'coefficients_new');

clear data