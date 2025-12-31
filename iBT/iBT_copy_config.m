function  iBT_copy_config(source,destination)
% Create runinfo and study list files from pre-defined templates
% FORMAT iBT_copy_config(source,destination)
%  where source indicates which runinfo template file to copy 
% from spm('dir')/toolbox/iBT/ to destination.
%
%___________________________________________________________________________
% Copyright 2006-2011 The Florey Institute of Neuroscience and Mental Health
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
% 2011-09-01: (dfa) Updated for public release
%

  full_source = fullfile(spm('dir'),'toolbox','iBT',source);
  target = destination;

  [SUCCESS,MESSAGE,MESSAGEID] = fileattrib(target);

  if  SUCCESS == 0 % Good, target does not exist
	[SUCCESS,MESSAGE,MESSAGEID] = copyfile(full_source,target);
	if  SUCCESS == 0 % Bad, copy failed
		disp(MESSAGE);
	else
		disp(['Successfully created: ' target ' in ' pwd]);
		disp(['Please edit this file to suit your requirements.']);
	end
  else
	disp(['Action aborted, file already exists: ' target])
  end


return
 
