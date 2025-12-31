function  iBT_custom_stats(iBT)
% Perform some user-customised analysis for the session specified in iBT.
% This script is designed to be called by iBT_analysis for each
% subject and session for which stats output would normally be generated.
% It is called after the iBT.what.stats.do and iBT.what.cons.do sections
% of the processing pipeline have been completed. It is enabled by setting
% the iBT.what.custom_stats.do =1 flag in your runinfo.
% Users can modify the present file to perform customised  post-processing
% on the stats results.
%
% FORMAT iBT_custom_stats(iBT)
%  where 'iBT.what' is a structure indicating what is to be done,
%    and 'iBT.who' is a structure indicating subject specific details.
% For any extra variables you wish to pass to this script from a runinfo,
% please use a "custom_stats" field in these structures. e.g. 
% 	iBT.who.sub{sub}.custom_stats.yourvariable1
% and	iBT.what.custom_stats.yourvariable2
%___________________________________________________________________________
% Copyright 2013 The Florey Institute of Neuroscience and Mental Health
%
% This file is part of iBT (the Integrated Brain Analysis Toolbox for SPM).
% See iBT_SPM.m for more information.
%
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
% Version History:
% 2013-08-13: (dfa) First version by David Abbott.


disp(sprintf('Entered custom analysis script...'));

[iBTversion,iBTrelease,iBTversionReleaseString] = iBT_version(1); % Display version so it gets logged.

pwd_orig=pwd;

try
    [iBT.what.SpmVersion,iBT.what.SpmRelease]   = spm('Ver','',1); iBT.what.SpmRelease = str2double(iBT.what.SpmRelease);
catch
    error('SPM cannot be found in MATLAB path.');
end

if  isnan(iBT.what.SpmRelease)
	iBT.what.SpmVersionReleaseString = iBT.what.SpmVersion;
else
	iBT.what.SpmVersionReleaseString = sprintf('%s (release %i)',iBT.what.SpmVersion,iBT.what.SpmRelease);
end
disp(sprintf('SPM version %s',iBT.what.SpmVersionReleaseString))

if isempty(who), error('Error: "who" structure is not defined.'), end

spm('defaults','fmri'); % Ensure we are starting with standard defaults

format compact

first_ses = iBT.what.first_ses;  % start with this session (counting from 1) (set in iBT_start.m)
last_ses = iBT.what.last_ses;  % end with this session (set in iBT_start.m)
nsess = 1 + last_ses - first_ses;  % number of sessions to analyse in this instance of iBT_custom_stats

sub = iBT.what.do_sub;

for ses=first_ses:last_ses 

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
%  Add your custom analysis code below this line
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%












%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
%  End your custom analysis code before this line
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

end %for ses

% Go back to initial directory in case it was changed above
cd(pwd_orig);

return
% End of iBT_custom_stats
