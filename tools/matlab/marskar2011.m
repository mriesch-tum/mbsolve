% marskar2011  Displays result of marskar2011 setup.
%              Reference data available in tools/reference-data.

% mbsolve: Framework for solving the Maxwell-Bloch/-Lioville equations
%
% Copyright (c) 2016, Computational Photonics Group, Technical University 
% of Munich.
%
% This program is free software; you can redistribute it and/or modify
% it under the terms of the GNU General Public License as published by
% the Free Software Foundation; either version 3 of the License, or
% (at your option) any later version.
%
% This program is distributed in the hope that it will be useful,
% but WITHOUT ANY WARRANTY; without even the implied warranty of
% MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
% GNU General Public License for more details.
%
% You should have received a copy of the GNU General Public License
% along with this program; if not, write to the Free Software Foundation,
% Inc., 51 Franklin Street, Fifth Floor, Boston, MA 02110-1301  USA

clear all;
close all;

%  choose hdf5 file
[filename, folder] = uigetfile('../../*.hdf', 'Select result data');
f = fullfile(folder, filename);
if (filename == 0)
   return;
end

comp_data = {};
[filenames, folder] = uigetfile('../reference-data/*.csv', ...
    'Select compare trace', 'MultiSelect', 'on');
if iscell(filenames)
    for i = 1:length(filenames)
        comp_data{i} = csvread(fullfile(folder, filenames{i}));
    end
elseif (filenames == 0)
    % noop
else
    comp_data{1} = csvread(fullfile(folder, filenames));
end

% read global attributes
d_x = h5readatt(f, '/', 'gridpoint_size');
d_t = h5readatt(f, '/', 'timestep_size');
t_e = h5readatt(f, '/', 'sim_endtime');
L_x = h5readatt(f, '/', 'dev_length');

% time grid
t = 0:d_t:t_e;

% data
e = h5read(f, '/e/real');
rho11 = h5read(f, '/d11/real');
rho22 = h5read(f, '/d22/real');
rho33 = h5read(f, '/d33/real');
rho44 = h5read(f, '/d44/real');
rho55 = h5read(f, '/d55/real');
rho66 = h5read(f, '/d66/real');

figure;
plot(t, e);
xlabel('Time/s');
ylabel('E-Field/Vm^{-1}');
xlim([0, t_e]);

figure;
plot(t, rho11, 'DisplayName', '\rho_{11}');
hold on;
plot(t, rho22, 'DisplayName', '\rho_{22}');
plot(t, rho33, 'DisplayName', '\rho_{33}');
plot(t, rho44, 'DisplayName', '\rho_{44}');
plot(t, rho55, 'DisplayName', '\rho_{55}');
plot(t, rho66, 'DisplayName', '\rho_{66}');
for i = 1:length(comp_data)
    plot(comp_data{i}(:, 1) * 1e-12, comp_data{i}(:, 2), ...
        'DisplayName', 'Ref');
end
xlabel('Time/s');
ylabel('Populations');
xlim([0, t_e]);
ylim([0, 1]);
legend show;
