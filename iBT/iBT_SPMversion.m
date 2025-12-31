function [SpmVersion, SpmRelease, SpmVersionReleaseString, MatlabVersion, MatlabDate, MatlabArch ] = iBT_SPMversion(display_flag)
%
% Get version of SPM in a consistent format.
%
% FORMAT [SpmVersion,SpmRelease, SpmVersionReleaseString] = iBT_SPMversion(display_flag)
%
% Where display_flag = 	2 or missing to print a message with the SpmVersionReleaseString and MatlabVersion
%			1 to print a message with the SpmVersionReleaseString,
% 			0 otherwise
%
%___________________________________________________________________________
% Copyright 2011-2016 The Florey Institute of Neuroscience and Mental Health
% David Abbott
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
% 2016-05-04: (dfa) Added Matlab version information
% 2012-03-12: (dfa) Permit call without arguments
% 2011-08-26: (dfa) Updated for public release
%

try tmp = display_flag; catch display_flag = 2; end 

MatlabArch = computer('arch');
[MatlabVersion, MatlabDate] = version;

try
    [SpmVersion,SpmRelease]   = spm('Ver','',1); SpmRelease = str2double(SpmRelease);
catch
    error('SPM cannot be found in MATLAB path.');
end


if  isnan(SpmRelease)
	SpmVersionReleaseString = SpmVersion;
else
	SpmVersionReleaseString = sprintf('%s (release %i)',SpmVersion,SpmRelease);
end

if ( display_flag == 1 ), disp(sprintf('SPM version %s',SpmVersionReleaseString)); end
if ( display_flag == 2 ), disp(sprintf('SPM version %s; MATLAB version %s (%s)',SpmVersionReleaseString,MatlabVersion,MatlabArch)); end
