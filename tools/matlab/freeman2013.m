clear all;
close all;

[filename, folder] = uigetfile('../../*.mat', 'Select result data');
load(fullfile(folder, filename));

t = 116;

x = 0:GridPointSize:XDim;

% explicitely non-complex quantities to real numbers
inv12 = real(inv12);
e = real(e);

figure;
plot(x, inv12(:, t));
ylabel('Population  inversion');
xlim([0, XDim]);

figure;
plot(x, e(:, t));
ylabel('E-Field');
xlim([0, XDim]);

x_idx = size(e, 1);
figure;
plot(linspace(0, SimEndTime, size(e, 2)), e(x_idx, :));
ylabel('E-Field');
xlabel('Time');
