function SPM=ibt_update_spm_paths(SPM, old_path, new_path)
% Function to update paths in an SPM.mat file or structure. 
% Inputs:
%   SPM - An SPM structure, or a full path to an SPM.mat file
%   old_path - Old base path to replace
%   new_path - New base path to use
% Outputs:
%   If successful, this function returns the modified SPM structure. 
%   If SPM is a filename, a backup file is first created (if it
%     does not already exist) that includes the original's last  
%     modification timestamp within the name. The SPM structure is  
%     then loaded, modified, and saved to replace the original file.
% Operation:
%   A search and replace operation is performed, replacing old_path 
%     with new_path in relevant path strings, including paths to 
%     pre-processed and stats files (which are likely in different,
%     folders, so please specify only the common base path to change).
% Examples:
%   The following would update an SPM.mat file after creating a backup:
%        ibt_update_spm_paths('./SPM.mat', '/data1/', '/data2/');
%   The following does the same,  and also retains the new SPM structure
%       in memory via the function return value:
%        SPM=ibt_update_spm_paths('./SPM.mat', '/data1/', '/data2/');
%   The following would update an already loaded SPM structure,
%   without reading or writing any files:
%        SPM=ibt_update_spm_paths(SPM, '/data1/', '/data2/');
%
% This function is part of iBT (integrated Brain Analysis Toolbox for SPM),
% and can be run stand-alone.
%___________________________________________________________________________
% Copyright 2025 The Florey Institute of Neuroscience and Mental Health

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

% Version History:
% 2025-03-26 (dfa): Original version by David Abbott
%

    if isstruct(SPM)
        file=0; % We will modify an already loaded SPM structure
    else
        file=1; % We will load and modify an SPM.mat file
	spm_input = SPM;
        % Verify SPM.mat file existence:
        if ~exist(spm_input, 'file')
            error(['The specified SPM.mat file is not found: ' spm_input]);
        end

        % Get the last modification timestamp of the SPM.mat file
        file_info = dir(spm_input);
        timestamp = datestr(file_info.datenum, 'yyyy-mm-dd_HH-MM-SS');
        backup_file = [spm_input, '_backup_', timestamp];
    
        % Check if a backup already exists with the same timestamped name
        if ~exist(backup_file, 'file')
            try
                copyfile(spm_input, backup_file);
                fprintf('Backup created successfully: %s\n', backup_file);
            catch ME
                error('Failed to create backup. Aborting. Error: %s', ME.message);
            end
        else
            fprintf('Backup already exists: %s\n', backup_file);
        end

        % Load the SPM.mat file
        load(spm_input, 'SPM');
    end %if isstruct(SPM) else

    % Update paths in SPM.xY.P 
    if isfield(SPM, 'xY') && isfield(SPM.xY, 'P')
      if iscell(SPM.xY.P), 
        SPM.xY.P = replace_paths(SPM.xY.P, old_path, new_path);
      else % will be a char array
	paths = cellstr(SPM.xY.P); % convert from char array to cellstr as easier to deal with
        paths = replace_paths(paths, old_path, new_path);
     	SPM.xY.P = char(paths); % convert back from cellstr to char array
      end% if iscell    
    end% if isfield

    % Check and update SPM.xY.VY().fname
    if isfield(SPM, 'xY') && isfield(SPM.xY, 'VY')
        for i = 1:numel(SPM.xY.VY)
	    if isfield(SPM.xY.VY(i), 'fname')
                SPM.xY.VY(i).fname = replace_paths(SPM.xY.VY(i).fname, old_path, new_path);
	    end
	end
    end

    % Update paths in SPM.swd
    if isfield(SPM, 'swd')
        SPM.swd = strrep(SPM.swd, old_path, new_path);
    end

    % Check and update SPM.xM.VM.fname
    if isfield(SPM, 'xM') && isfield(SPM.xM, 'VM') && isfield(SPM.xM.VM, 'fname')
        SPM.xM.VM.fname = replace_paths(SPM.xM.VM.fname, old_path, new_path);
    end

    %% The rest below might not contain a path (they didn't in the SPM.mat that I was
    %% using to guide me), but no harm done if not (they'll stay the same):
    
    % Check and update SPM.xVol.VRpv.fname
    if isfield(SPM, 'xVol') && isfield(SPM.xVol, 'VRpv') && isfield(SPM.xVol.VRpv, 'fname')
        SPM.xVol.VRpv.fname = replace_paths(SPM.xVol.VRpv.fname, old_path, new_path);
    end

    % Check and update SPM.Vbeta().fname
    if isfield(SPM, 'Vbeta') 
        for i = 1:numel(SPM.Vbeta)
	    if isfield(SPM.Vbeta(i), 'fname')
                SPM.Vbeta(i).fname = replace_paths(SPM.Vbeta(i).fname, old_path, new_path);
	    end
	end
    end

    % Check and update SPM.VM.fname
    if isfield(SPM, 'VM') && isfield(SPM.VM, 'fname')
        SPM.VM.fname = replace_paths(SPM.VM.fname, old_path, new_path);
    end

    % Check and update SPM.VResMS.fname
    if isfield(SPM, 'VResMS') && isfield(SPM.VResMS, 'fname')
        SPM.VResMS.fname = replace_paths(SPM.VResMS.fname, old_path, new_path);
    end

    % Check and update SPM.xCon().fname
    if isfield(SPM, 'xCon')
        for i = 1:numel(SPM.xCon)
            if isfield(SPM.xCon(i), 'Vspm') && isfield(SPM.xCon(i).Vspm, 'fname')
                SPM.xCon(i).Vspm.fname = replace_paths(SPM.xCon(i).Vspm.fname, old_path, new_path);
            end
            if isfield(SPM.xCon(i), 'Vcon') && isfield(SPM.xCon(i).Vcon, 'fname')
                SPM.xCon(i).Vcon.fname = replace_paths(SPM.xCon(i).Vcon.fname, old_path, new_path);
            end
        end
    end
    
    if file,
        % Save the updated SPM.mat file
        save(spm_input, 'SPM');
        disp(['Paths updated successfully in ' spm_input]);
    end % if file

    return

end

function updated_paths = replace_paths(paths, old_path, new_path)
    % Helper function to replace old paths with new paths
    if iscell(paths)
        updated_paths = cellfun(@(p) strrep(p, old_path, new_path), paths, 'UniformOutput', false);
    else
        updated_paths = strrep(paths, old_path, new_path);
    end
end

