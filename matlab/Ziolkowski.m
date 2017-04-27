clear all;
close all;

%[filename, folder] = uigetfile('~/Desktop/*.mat', 'Select performance data');
%[filename, folder] = uigetfile('~/CPH/Articles/riesch2017a/Data/*.mat', 'Select performance data');
[filename, folder] = uigetfile('~/CPH/Work/mbsolve/build-openmp/*.mat', 'Select performance data');
load(fullfile(folder, filename));

t = 1;

x = 0:GridPointSize:XDim;

figure;
plot(x, dm11(:, t) + dm22(:, t));
xlim([0, XDim]);
ylabel('Trace as sanity check');

figure;
plot(x, dm11(:, t) - dm22(:, t));
ylabel('Population  inversion');
xlim([0, XDim]);
figure;
plot(x, dm11(:, t));
ylabel('dm11');
xlim([0, XDim]);

figure;
plot(x, dm22(:, t));
ylabel('dm22');
xlim([0, XDim]);

figure;
plot(x, e(:, t));
ylabel('E-Field');
xlim([0, XDim]);
