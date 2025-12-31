function ibt(varargin)
%
%  Command-line script to start the Integrated Brain Analysis Toolbox for SPM.
%  Place this script somewhere in your Matlab path.
%
%  The script will add the iBT folder to the path if necessary, 
%	 (assuming it is in the current working directory in a folder called
%  iBT or bin, or in the iBT folder in the toolbox 
%  folder in your SPM distribution).
%
%  It will then either prompt for a runinfo to run (if called with no 
%  arguments), or just print the toolbox version if called as ibt('version').
%
% FORMAT ibt
%	Print the version of the toolbox that will be used and then
% 	promt the user for runinfo file(s) to run
%       Multiple runinfos can be specified, in which case they are
%       run as separate analyses. There is presently no interactive
%       mechanism to chain runinfo pairs or tupples in a single analysis.
%
% FORMAT ibt version | ver | v
%	Just print the version of the toolbox that will be used.
%
% FORMAT ibt('version' | 'ver' | 'v')
%	Just print the version of the toolbox that will be used.
%
% FORMAT ibt(filespec)
%	Print the version of the toolbox that will be used and then
%	run with runinfo filenames matching filespec. Each matching
%       runinfo file is run as a separate analysis.
%
% FORMAT ibt(filespec_part1, filespec_part2, ...)
%	Multi-file runinfo parts are chained together: Filenames matching each
%       argument are resolved and then the first analysis proceeds by chaining
%       together the first match from filespec_part1 and the first match
%       from filespec_part2 etc. That is, it acts as if the contents of
%       these first matching runinfos are concatenated. Once that analysis 
%       is complete, the next analysis proceeds by chaining the second match 
%       in each filespec (if any) and so on. 
%       If any of the filespecs match only one file, then that
%       file is used in every analysis. In this way, one can
%       utilise a single master runinfo with settings in common to
%       several analyses, specified for example as filespec_part1, and
%       then chain to it a much simpler runinfo containing only the relevant 
%       overrides to the master, specified as filespec_part2. Each filename
%       matching filespec_part2 would then execute as a separate analysis
%       chained to the single master runinfo.
%       Another potential use case would be to provde a common override to
%       multiple runinfos. For example one could indicate
%       pre-processing of all runinfos has already been completed, without
%       having to edit every single runinfo.
%
%  See README_iBT.txt 
%  folder for details.
%______________________________________________________________________________
% Copyright 2011,2019-21 The Florey Institute of Neuroscience and Mental Health
% David Abbott
%
% This file is part of iBT (the Integrated Brain Analysis Toolbox for SPM).
%
% This toolbox was presented at the 2011 OHBM meeting:
%   Abbott DF, Masterton RAJ, Waites AB, Bhaganagarapu K, Pell GS, 
%   Harvey MR, Sharma GJ, Jackson GD. "The iBrain(TM) Analysis Toolbox 
%   for SPM". Proceedings of the 17th Annual Meeting of the Organisation
%   for Human Brain Mapping, Quebec City, Canada. 364 WTh. (2011).
%
% Please refer to the following for a description of the laterality
% methods implemented in this toolbox. If you publish laterality
% results using these methods or variations of them, please cite 
% this paper: 
%     Abbott DF, Waites AB, Lillywhite LM, Jackson GD.
%     fMRI assessment of language lateralization: An objective approach. 
%     Neuroimage 50(4):1446-1455 (2010).
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
% Recent Version History:
% 2021-08-06: (dfa) Updated usage documentation
% 2021-07-20: (dfa) Support multiple runinfo filespec arguments
% 2020-04-30: (dfa) Support runinfo filespec argument
% 2019-04-03: (dfa) If the toolbox is not already in the path, add
%                   the appropriate folder from the current working 
%                   directory (if one exists there), otherwise add
%                   from the toolbox folder of the SPM installation.
% 2011-09-19: (dfa) Don't add toolbox to path if a version already there.
% 2011-08-26: (dfa) Updated for public release
%


% Determine best option we have for the iBT version to use
  newpath = fullfile(pwd,'iBT');
  if exist(fullfile(newpath,'iBT_start.m'), 'file') == 2,
  	ibt_folder = 1;
  else
    newpath = fullfile(pwd,'bin');
    if exist(fullfile(newpath,'iBT_start.m'), 'file') == 2,
  	ibt_folder = 2;
    else
      newpath = fullfile(spm('dir'),'toolbox','iBT');
      if exist(fullfile(newpath,'iBT_start.m'), 'file') == 2,
  	ibt_folder = 3;
      else
  	ibt_folder = -1;
      end
    end
  end

if ibt_folder > 0

  % To avoid comparison mismatches if symbolic links are involved,
  % first ensure we have the real physical path for newpath
  cwd = pwd;
  cd(newpath);
  newpathReal = pwd;
  cd(cwd)

end

start = which('iBT_start');
currentStartFile = fullfile(start);
currentStartPath=fileparts(currentStartFile);

if strcmp(start,'') % Not yet in path - so put our preferred option in

  if ibt_folder < 0
  	disp(['Unable to find iBT. Please check your Matlab path.']);
  else
        addpath(newpathReal);
 	disp(['Added ' newpathReal ' to Matlab path.']);
  end

else

  if ibt_folder > 0	
    newStartFile = fullfile(newpathReal,'iBT_start.m');

    if ~strcmp(currentStartFile,newStartFile) % then switch to the preferred path
        rmpath(currentStartPath);
        addpath(newpathReal);
	if strcmp(newpath,newpathReal)
  		disp(['Switched to using iBT in ' newpathReal]);
	else
  		disp(['Switched to using iBT in ' newpath ' (a link to "' newpathReal '").']);
	end
    else % prefered option matches current, so continue with existing path
  	disp(['Using iBT in ' currentStartPath]); 
    end
  else % ~ ibt_folder > 0, so no prefered option, so continue with existing path
  	disp(['Using iBT in ' currentStartPath]); 
  end
end

if nargin>0
 switch lower(varargin{1}), 
 case 'version'
    iBT_version(1);
 case 'ver'
    iBT_version(1);
 case 'v'
    iBT_version(1);
 otherwise
    	%iBT_version(1);
	% disp(['Warning(iBT): Unknown argument provided to ibt ']);
	iBT_start(varargin{:}); % This will also print the toolbox version used
		% The {:} breaks out the individual arguments so that the varargin of 
		% the called routine will maintain the same dimensions as here in the caller.
 end %switch
else
	iBT_start; % This will also print the toolbox version used
end %if
return  
