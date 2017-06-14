clear all;
close all;

[filename, folder] = uigetfile('../../*.mat', 'Select result data');
load(fullfile(folder, filename));

t = 116;

x = 0:GridPointSize:XDim;

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
ylabel('Population  inversion');
xlim([0, XDim]);

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

x_idx = size(e, 1);
figure;
plot(linspace(0, SimEndTime, size(e, 2)), e(x_idx, :));
ylabel('E-Field');
xlabel('Time');
