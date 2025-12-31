function XYZ = iBT_sphere_roi(r, voxsize, center)
% 
%   Returns an array of co-ordinates that samples a sphere 
%   of radius r centered on center. 
%
% FORMAT XYZ = iBT_sphere_roi(r, voxsize, center)
%
%   The returned array is sampled at, and returned in,
%   the same units as the center co-ordinate if one of 
%   the following calling conventions is followed:
%
%   r       = radius (mm)
%   voxsize = voxel sizes (mm)
%   center  = center coord (voxels)
%   Returns XYZ as an array of co-ordinates inside the 
%   sphere in units of voxels and sampled at 1-voxel spacing.
% or
%   r       = radius (voxels)
%   voxsize = [1 1 1]
%   center  = center coord (voxels)
%   Returns XYZ as an array of co-ordinates inside the 
%   sphere in units of voxels and sampled at 1-voxel spacing.
% or
%   r       = radius (mm)
%   voxsize = [1 1 1]
%   center  = center coord (mm)
%   Returns XYZ as an array of co-ordinates inside the 
%   sphere in units of mm and sampled at 1mm spacing.
% or
%   r       = radius (voxels)
%   voxsize = inverse voxel sizes (/mm)
%   center  = center coord (mm)
%   Returns XYZ as an array of co-ordinates inside the 
%   sphere in units of mm and sampled at 1mm spacing.
%
%___________________________________________________________________________
% Copyright 2009-2011 The Florey Institute of Neuroscience and Mental Health
%
% This file is part of iBT (the Integrated Brain Analysis Toolbox for SPM).
% See iBT_SPM.m for more information.

% iBT is free software: you can 
% redistribute it and/or modify it under the terms of the GNU General 
% Public License as published by the Free Software Foundation, either
% version 3 of the License, or (at your option) any later version.
% 
% iBT is distributed in the hope 
% that it will be useful, but WITHOUT ANY WARRANTY; without even the  
% implied warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR 
% PURPOSE.  See the GNU General Public License for more details.
%
% You should have received a copy of the GNU General Public License
% along with this program. If not, see <http://www.gnu.org/licenses/>.
%___________________________________________________________________________
%
% Recent Version History:
% 2011-08-26: (dfa) Updated for public release
%

% 2004-2011: iBT analysis scripts were developed by David Abbott based
%            upon the who & what structure concept of Tony Waites, with 
%            substantial contributions from Matt Harvey, Kaushik Bhaganagarapu,
%            and Richard Masterton.
%            The script that ultimately became iBT_sphere_roi
%            was originally written by Richard Masterton


    try
        center;
    catch
        center = [0 0 0];
    end;

    try
        voxsize = abs(voxsize);
    catch
        voxsize = [1 1 1];
    end;

    try
        r;
    catch
        r = 5;
    end

    extent = ceil(r./voxsize);
    [x y z] = meshgrid([-extent(1):extent(1)], ...
                        [-extent(2):extent(2)], ...
                        [-extent(3):extent(3)]);
    s = r - sqrt((x.*voxsize(1)).^2 + (y.*voxsize(2)).^2 + (z.*voxsize(3)).^2);

    X = x(s > 0);
    Y = y(s > 0);
    Z = z(s > 0);
    XYZ = [X+center(1) Y+center(2) Z+center(3)]';

end
