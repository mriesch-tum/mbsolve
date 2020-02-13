% ziolkowski1995  Displays result of ziolkowski1995 setup.
%                 Reference data available in tools/reference-data.

% mbsolve: An open-source solver tool for the Maxwell-Bloch equations.
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

[filename, folder] = uigetfile('../reference-data/*.csv', ...
    'Select compare trace');
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
xlabel('Propagation direction/m');
ylabel('Population  inversion');

figure;
plot(x, e(:, i));
xlabel('Propagation direction/m');
ylabel('E-Field/Vm^{-1}');
xlim([0, L_x]);
