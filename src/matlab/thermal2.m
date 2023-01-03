function map = thermal2(varargin)
%THERMAL  Purple-red-orange-yellow colormap
%
% Examples:
%   map = thermal;
%   map = thermal(len);
%   B = thermal(A);
%   B = thermal(A, lims);
%
% A colormap designed to replicate the tones of thermal images.
% IN:
%   len - Scalar length of the output colormap. If len == Inf the concise
%         table is returned. Default: len = size(get(gcf, 'Colormap'), 1);
%   A - Non-scalar numeric array of real values to be converted into
%       truecolor.
%   lims - 1x2 array of saturation limits to be used on A. Default:
%          [min(A(:)) max(A(:))].
%
% OUT:
%   map - (len)x3 colormap table.
%   B - size(A)x3 truecolor array.

map =   [0.1500        0    0.6000;
        0.4000         0    0.5286;
        0.6714    0.0643    0.3714;
        0.9143    0.1857    0.1857;
        0.9857    0.3643    0.1143;
        1.0000    0.5714    0.0714;
        1.0000    0.7857    0.0357;
        1.0000    1.0000         0];
% map =  [0.15  0.0  0.6; % 0.75
%         0.5   0.0  0.5;  % 1.0
%         0.9   0.15 0.2; % 1.25
%         1.0   0.4  0.1;  % 1.5
%         1.0   0.7  0.05; % 1.75
%         1.0   1.0  0.0 ];% 2.0
map = colormap_helper(map, varargin{:});
