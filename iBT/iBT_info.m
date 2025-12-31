function [present, value] = iBT_info(infoFile, token, varargin)
% Read a value from a previously saved header info file or the like
%
% FORMAT [present, value] = iBT_TR(infoFile, token, [type])
% where infoFile is the filename of the information file,
% 	token 	is the string preceding a colon for which to return the 
% 		value after the colon.
%	type	is an integer code indicating the return type as follows:
%		0 = string (default)
%		1 = number (via call to str2num() )
%
%	present is set to 1 if the token is found in the file, 0 otherwise
%
% Example: If a text file contains, amongst other things, 
%	   a line as follows:
%
% TR(s):3.2
%
% then [present, value] =  iBT_info(infoFile, 'TR(s)') 
% will return [1, '3.2']
% and [present, value] = iBT_info(infoFile, 'TR(s)', 1)
% will return [1, 3.2000]
%
%___________________________________________________________________________
% Copyright 2011-2012 The Florey Institute of Neuroscience and Mental Health
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
% 2012-02-14: (dfa) Issue an informative error if can't open file.
% 2011-08-26: (dfa) Updated for public release
%
  optargin = size(varargin,2);
  if optargin < 1, type = 0; else 
 	type = varargin{1};
   end

  fsettings_is_OK = 0; % Initialise
  fsettings = fopen(infoFile); % open info file and read the data
  if fsettings == -1
  		throw(MException('iBT_info', ['Unable to find file matching ' infoFile]));
  end
  
  % Most of the file consists of label:value pairs. Keep reading until we find desired label
   buffer = fgets(fsettings);
   while buffer ~= -1 % Continue unless we are at the end of the file
	[tok,value] = strtok(buffer,':');
	if strcmp(tok,token)
		fsettings_is_OK = 1;
        	break; %we have found what we want
	end
   	buffer = fgets(fsettings);
   end % While
 if fsettings_is_OK == 1,
	value = strrep(value,':',' ');
	if type > 0 value = str2num(value);	end
	present = 1;
 else
 	present = 0; % Unable to find value in file.
 end
 fclose(fsettings);
  
% end of function 
