function TR = iBT_extract_TR(input_dir, wild)
% Extract the repetition time (TR) from JSON sidecar
% FORMAT iBT_extract_TR(input_dir, wild)
%  where input_dir is a string name of the folder containing the JSON sidecar
%    and wild is a string to match the filename 
%    (e.g. for fMRI prep it might be ["*task-*" "pseusoword" "*_bold.json"] )
% Returns the TR in seconds, or 0 if error (including if the TR is inconsistent)
%
% Note: The TR in a BIDS JSON sidecar should be in seconds according to the BIDS
%   standard. However SPM12's JSON sidecar is not a BIDS sidecar, and its field
%   acqpar.RepetitionTime has been observed to be in milliseconds - possibly
%   vendor units which might vary? Therefore, as a somewhat blunt heuristic, 
%   if the TR read from an SPM sidecar is 100 or more it is treated as in ms
%   and the value is divided by 1000 before being returned.
%
% Requires MATLAB R2016b or later, for jsondecode()
%
%___________________________________________________________________________
% Copyright 2022-2025 The Florey Institute of Neuroscience and Mental Health
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
% 2025-12-26: (dfa) Change on 2025-11-06 requires processing iBT_get()
%                   filenames with deblank()
% 2025-11-06: (dfa) Changed task argument to full wildcard; changed strings
%    to legacy to maintain compatibility with MATLAB r2016b.
%    Support SPM JSON files.
%    Return to caller with TR set to 0 if any error.
% 2024-10-31: (dfa) Improved diagnostic output formatting
% 2024-06-25: (dfa) Document and incorpotate into public iBT release
% 2023-10-10: (ac) extract_TR() by Aaron Capon

  sub_dir = fullfile(input_dir);

  if not(exist(sub_dir, 'dir') == 7)
    disp('Could not locate input dir to determine TR');
    disp(sprintf('Searched for: %s\n', sub_dir));
    TR = 0; return;
  end

  %in_files = dir(fullfile(sub_dir, wild));
  in_files=iBT_get(sub_dir,wild);
  if length(in_files) == 0
    disp('Could not locate json files to determine TR');
    TR = 0; return;
  end

  TR = [];

  for i = 1 : size(in_files,1)
    json_file = fileread( deblank(in_files(i,:)) );
    [~, name, ext] = fileparts( deblank(in_files(i,:)) );
    fprintf('Reading repetition time from: %s ', [name ext]);
    json_data = jsondecode(json_file);
    if isfield(json_data,'RepetitionTime') % Hopefully a BIDS JSON sidecar
      TR(i) = json_data.RepetitionTime;
    elseif isfield(json_data,'acqpar') && isfield(json_data.acqpar,'RepetitionTime') % Could be an SPM JSON sidecar
       TR(i) = json_data.acqpar.RepetitionTime;
       if TR(i) >= 100 TR(i) = TR(i) / 1000.0; end  % An arbitrary heuristic for when to assume ms and convert to s
    else
      disp('WARNING: Could not locate RepetitionTime in this json side car.');
    end
    disp(sprintf(': %g', TR(i)));
  end

  if isempty(TR)
    disp('Could not locate json side cars to read in repetition time.');
    TR = 0; return;
  end

  % Consistency check
  if all(TR == TR(1))
    TR = TR(1);
    if i > 1
    	disp('Repetition time is consistent between available json files');
    end
  else
    disp('!!! ERROR !!!');
    disp('Repetition time varies in json files: ');
    disp(sprintf('%s ', string(TR), ' '));
    TR = 0; return;
  end

  

end
