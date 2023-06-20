function varargout = vec2grid(varargin)
%VEC2GRID Grid ND data without interpolation
%
% [xg, yg] = vec2grid(x1, x2, ..., y);
% [xg, yg] = vec2grid(xy);
% [xg1, xg2, ... yg] = vec2grid(...);
%
% This function reformats a list of n-dimensional datpoints to a grid,
% similar to the results of the griddata function but without any
% interpolation. It reshapes the data using all unique values in the 1st,
% 2nd, 3rd, etc. as rows, columns, pages, etc. of the final data grid.     
%
% Input variables:
%
%   x#:     coordinates of data points, vectors.  Each input corresponds to
%           one dimension of the data 
%
%   y:      y coordinates of data points, vector the same size as x(s)
%
%   xy:     npoint x ndim+1 array, where each row represents one point.
%           The first ndim columns are the coordinates of the points, and
%           the final column holds the data value at that coordinate
%
% Output variables:
%
%   xg#:    npoint x 1 cell array of unique values of the #th dimension of
%           yg 
%
%   yg:     ndim-dimensional array, y values rearranged into a grid
%
%   xg:     1 x ndim cell array holding all xg# vectors
%
% Example:data
%
% [x,y] = meshgrid(1:5,0:2:6);
% x = x(:); y = y(:);
% holes = rand(size(x)) < .2;
% x(holes) = []; y(holes) = [];
% v = 1:length(x);
% [xg,yg,vg] = vec2grid(x,y,v)
% xg =
%      1
%      2
%      3
%      4
%      5
% yg =
%      0
%      2
%      4
%      6
% vg =
%      1   NaN   NaN   NaN    13
%      2     5     8    11    14
%      3     6     9   NaN    15
%      4     7    10    12    16
% Copyright 2008 Kelly Kearney
%-------------------------
% Check input
%-------------------------
if nargin == 1
    x = num2cell(varargin{1}(:,1:end-1), 1);
    y = varargin{1}(:,end);
else
    if ~all(cellfun(@isvector, varargin)) || any(diff(cellfun(@length, varargin)))
        error('Inputs must be vectors of the same length');
    end
    x = varargin(1:end-1);
    y = varargin{end};
end
   
% error(nargchk(3, 4, nargin));
%-------------------------
% Create grid
%-------------------------
ndim = length(x);
[xunique, blah, xidx] = cellfun(@unique, x, 'uni', 0);
nx = cellfun(@max, xidx);
ynew = nan(nx);
ynew(sub2ind(nx, xidx{:})) = y;
if nargout == (ndim+1)
    varargout = cell(ndim+1,1);
    tmp = [xunique ynew];
    [varargout{:}] = deal(tmp{:});
elseif nargout == 2
    varargout{1} = xunique;
    varargout{2} = ynew;
end
% if nargin == 3
%     
%     x = varargin{1};
%     y = varargin{2};
%     z = varargin{3};
% 
%     [xunique, blah, xidx] = unique(x);
%     [yunique, blah, yidx] = unique(y);
% 
%     znew = nan(max(yidx),max(xidx));
%     znew(sub2ind(size(znew), yidx, xidx)) = z;
%     
%     varargout = {xunique, yunique, znew};
%     
% elseif nargin == 4
%     
%     x = varargin{1};
%     y = varargin{2};
%     z = varargin{3};
%     v = varargin{4};
%     
%     [xunique, blah, xidx] = unique(x);
%     [yunique, blah, yidx] = unique(y);
%     [zunique, blah, zidx] = unique(z);
%     
%     vnew = nan(max(zidx), max(yidx), max(xidx));
%     vnew(sub2ind(size(vnew), zidx, yidx, xidx)) = v;
%     
%     varargout = {xunique, yunique, zunique, vnew};
%     
% end
%     