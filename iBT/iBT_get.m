function files=iBT_get(dir,wild)
% Similar to spm_select('FPList',dir,spm_wildconvert(wild)) for SPM5 or later
% or spm_get('Files',dir,wild) for SPM2, but also allows wild to have
% a relative path component.
% 
% FORMAT iBT_preprocess(dir,wild)
%  where dir is a directory and wild is a human-readable filespec that can
%  contain wildcard characters. This function internally converts the
%  wildcards to regex format for SPM5 and later.
%___________________________________________________________________________
% Copyright 2004-2025 The Florey Institute of Neuroscience and Mental Health
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
% 2025-12-26: (dfa) More robust path element concatenation
% 2025-11-11: (dfa) Allow wild to contain relative path elements
% 2025-11-04: (dfa) Adapted code duplicated many times elsewhere in iBT
%   so can eventually tidy up / consolidate.
%

  try SPMver = iBT.what.SpmVersion; catch; SPMver=''; end % Will default to conversion to regex
  
  [pth,nam,ext] = fileparts(wild);
  path = fullfile(dir, pth); % Put the path elements together
  filespec = [nam ext]; % Reassemble the filespec, without any path element

  switch (SPMver),
    case 'SPM2'
      files=spm_get('Files',path, filespec);
    otherwise
      files=spm_select('FPList',path, spm_wildconvert(filespec));
    end


function new_str = spm_wildconvert(orig_str)
%
% Sub function to convert normal wildcards to a regexp format for spm
%
% Modification History
% ---------------------
%
% 2008-08-18 - matt: Function creation

new_str =  strcat('^', orig_str, '$');
new_str = strrep(new_str,'.','\.');
new_str = strrep(new_str,'?','.');
new_str = strrep(new_str,'*','(.*)');
