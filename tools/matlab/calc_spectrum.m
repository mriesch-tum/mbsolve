function [ spectrum, freqs ] = calc_spectrum(field, d_t)
% calc_spectrum  Calculate single-sided spectrum of parameter using FFT. 

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

% use power of two as size for efficiency
size_fft = 2^nextpow2(length(field));

% window function
win = hanning(length(field));
win = reshape(win, size(field));

% amplitude correction factor for hanning window
corr = 2;

%figure;
%plot(field);

%figure;
%plot(field .* win);

% transform field
Y = fft(field .* win, size_fft);

% single-sided spectrum
P2 = abs(Y/size_fft * corr);
P1 = P2(1:size_fft/2 + 1);
P1(2:end-1) = 2 * P1(2:end-1);
spectrum = P1;

% frequencies
freqs = 1/d_t * (0:(size_fft/2))/size_fft;



% TODO return phase as well?

end

