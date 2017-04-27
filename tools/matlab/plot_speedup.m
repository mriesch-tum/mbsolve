function plot_speedup()

[filename, folder] = uigetfile('*.csv', 'Select performance data');

data = csvread(fullfile(folder, filename));

[unique_x, ~, subscripts] = unique(data(:, 1));
ymean = accumarray(subscripts, data(:, 2), [], @mean);


figure;
plot(data(:, 1), data(:, 2), 'o');
hold on;
plot(unique_x, ymean);
plot(unique_x, ymean(1) .* unique_x);

end
