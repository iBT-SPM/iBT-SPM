function [iBTversion,iBTrelease,iBTversionReleaseString,iBTversionReleaseSentence,iBTmenuTitle] = iBT_version(display_flag)
%
% Get version of the Integrated Brain Analysis Toolbox for SPM (iBT).
%
% FORMAT [iBTversion,iBTrelease, iBTversionReleaseString,iBTversionReleaseSentence,iBTmenuTitle] = iBT_version(display_flag)
%
% Where display_flag = 1 to print the iBTversionReleaseSentence, 0 otherwise
%
%___________________________________________________________________________
% Copyright 2011-2025 The Florey Institute 
%                     of Neuroscience and Mental Health
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
% 2012-03-12: (dfa) Moved version info from iBT_start.

iBTversion = '4.0.1';
iBTrelease = '2026.01.03';

iBTversionReleaseString = sprintf('%s (release %s)',iBTversion,iBTrelease);

iBTmenuTitle        = sprintf('iBT version %s',iBTversion);
iBTversionReleaseSentence = sprintf('Integrated Brain Analysis Toolbox (iBT) version %s for SPM',iBTversionReleaseString);

try tmp = display_flag; catch display_flag = 1; end 

if ( display_flag == 1 ), disp(iBTversionReleaseSentence); end
