function result = ibt_nanmean(x,dim)
% Works like Matlab's mean() function but excludes NaNs.
% Note: In Matlab R2015a the MathWorks finally (belatedly for our purposes)
%	included the ability to exclude NaN's from the standard descriptive 
%	statistics functions such as mean, std, etc. Previously MathWorks 
%	only offered functions in the Matlab Statistics toolbox to do this.
%	Our own implementation is designed to work in both recent and older 
%	versions of Matlab without requiring the Statistics toolbox. It can
%	therefore, for example, still be used on older MATLAB versions that
%       are required when running on Red Hat Enterprise Linux 5 / CenOS 5, 
%	and Max OS X Mountain Lion (OS X 10.8) or earlier, which cannot run
%	Matlab R2015a or later.
%	We also avoid buggy nanmean replacements in older versions of SPM, 
%	and divide by zero warnings in current versions of SPM.
%___________________________________________________________________________
% Copyright 2015-2016 The Florey Institute of Neuroscience and Mental Health
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
% 2016-04-27: (dfa) Updated documentation.
% 2015-11-11: (dfa) If dim is specified but is invalid, cause an error.
% 2015-11-29: (dfa) Original version by David Abbott

NaN_indicies = isnan(x); % Find location of the NaNs
x(NaN_indicies) = 0; % Set them all to zero

dim_defined = 0; % Initialise
try, dim;
	dim_defined = 1;
	N = sum(~NaN_indicies,dim); %How many real numbers are there (sum across dim)?
	N(N==0) = NaN; % Avoid dividing by zero if there are no real numbers
	result = sum(x,dim) ./ N;
catch
	if dim_defined == 1,
		% We must have arrived here due to a crash above AFTER verifying that the dim argument was supplied.
		error(['Error(ibt_nanmean): Invalid dimension (dim) argument given for supplied array. Correct syntax is: ibt_nanmean(x,dim)']);
	else
		N = sum(~NaN_indicies); %How many real numbers are there (sum across first non-singleton)?
		N(N==0) = NaN; % Avoid dividing by zero if there no real numbers - result will be NaN.
		result = sum(x) ./ N;
	end % if dim_defined == 1
end
