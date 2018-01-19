clear all;
close all;

[filename, folder] = uigetfile('../../*.mat', 'Select result data');
load(fullfile(folder, filename));

t = 76;

x = 0:GridPointSize:XDim;

[filename, folder] = uigetfile('../../*.csv', 'Select compare trace');
use_reference = 0;
if (filename ~= 0)
    comp_data = csvread(fullfile(folder, filename));
    use_reference = 1;
end

% explicitely non-complex quantities to real numbers
inv12 = real(inv12);
e = real(e);

figure;
plot(x, inv12(:, t));
xlim([0, XDim]);
if (use_reference == 1)
    hold on;
    plot(comp_data(:, 1) * 1e-6, comp_data(:, 2));
end
ylabel('Population  inversion');

figure;
plot(x, e(:, t));
ylabel('E-Field');
xlim([0, XDim]);
