clear all;
close all;

%  choose hdf5 file
[filename, folder] = uigetfile('../../*.hdf', 'Select result data');
f = fullfile(folder, filename);

[filename, folder] = uigetfile('../../*.csv', 'Select compare trace');
use_reference = 0;
if (filename ~= 0)
    comp_data = csvread(fullfile(folder, filename));
    use_reference = 1;
end

% read global attributes
d_x = h5readatt(f, '/', 'gridpoint_size');
d_t = h5readatt(f, '/', 'timestep_size');
t_e = h5readatt(f, '/', 'sim_endtime');
L_x = h5readatt(f, '/', 'dev_length');

% grid
x = 0:d_x:L_x;

% data
e = h5read(f, '/e/real');
inv12 = h5read(f, '/inv12/real');

% select time index
i = 76;

figure;
plot(x, inv12(:, i));
xlim([0, L_x]);
if (use_reference == 1)
    hold on;
    plot(comp_data(:, 1) * 1e-6, comp_data(:, 2));
end
ylabel('Population  inversion');

figure;
plot(x, e(:, i));
ylabel('E-Field');
xlim([0, L_x]);
