clear all;
close all;

[filename, folder] = uigetfile('../../*.mat', 'Select result data');
load(fullfile(folder, filename));

t = 76;

x = 0:GridPointSize:XDim;

[filename, folder] = uigetfile('../../*.csv', 'Select compare trace');
comp_data = csvread(fullfile(folder, filename));

% explicitely non-complex quantities to real numbers
d11 = real(d11);
d22 = real(d22);
e = real(e);

figure;
plot(x, d11(:, t) + d22(:, t));
xlim([0, XDim]);
ylabel('Trace as sanity check');

figure;
plot(x, d11(:, t) - d22(:, t));
xlim([0, XDim]);
hold on;
plot(comp_data(:, 1) * 1e-6, comp_data(:, 2));
ylabel('Population  inversion');




figure;
plot(x, d11(:, t));
ylabel('dm11');
xlim([0, XDim]);

figure;
plot(x, d22(:, t));
ylabel('dm22');
xlim([0, XDim]);

figure;
plot(x, e(:, t));
ylabel('E-Field');
xlim([0, XDim]);
