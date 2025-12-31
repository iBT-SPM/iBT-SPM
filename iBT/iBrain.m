function status = iBrain(job, dir, varargin)
% Function to call an external iBrain command.
% Note: iBrain is a separate, self-contained program, not included in this toolbox, 
%   It is free for non-commercial use. See florey.edu.au/imaging-software
%
% FORMAT iBrain(job, dir [,use_server [, halt_on_error [, timeout [, del ]]]])
% where:
%   job is a string consisting of arguments compatible with 
%	the iBrain command-line syntax.
%
%   dir should be set to the fully qualified directory from which 
%	it is desired the iBrain command be executed. This must be
%	writable by the present user.
%
%   use_server 	= 1 to use an iBrain Server; jobs are placed in 
%			files in a queue folder whose location
%			depends on the computer operating system: 
%			On Windows the queue is %LOCALAPPDATA%\iBrain_queue\
% 			or if %LOCALAPPDATA% is not set, then it is 
%				%APPDATA%\iBrain_queue\
%			On other platforms the queue is $HOME/.iBrain_queue/
%		= n > 1 to use an iBrain Server and attempt to 
%			start a new iBrain Server process if the
%			job is not accepted within the specified
%			time n in seconds (default 10, minimum 5).
%		= 0 to call iBrain from the system command line
%	 (only use 0 if you have a script called iBrainr
%	  that will start iBrain with a runtime or full 
%	  IDL licence to avoid the IDL splash screen).
%
%  halt_on_error = 1 to throw a matlab error on error (default), 
%		otherwise 0 to return error status value
% The following optional arguments only apply if use_server == 1:
%  timeout = time in seconds to wait for iBrain to accept the job 
%		before giving up; if not specified then wait a year!
%  del = 1 to delete job from iBrain queue if timeout occurs, or
%	 0 to leave job in queue after timeout (default = 1)
%
%
%
%  Returns status: 	 0 = command completed successfully
%			-1 = failed to create iBrain queue folder (a file
%			     of that name exists and could not be deleted).
%			-2 = failed to create iBrain queue folder
%			-3 = failed to create iBrain command script
%			-4 = timeout (iBrain server did not accept command)
%			-5 = failed to run iBrainr system command
%			-10 = an error occured when iBrain server attempted 
%				to complete the command.
%___________________________________________________________________________
% Copyright 2009-2011,2017 The Florey Institute of Neuroscience and Mental Health
% David Abbott and Gagan Sharma
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

% Version History:
% 2017-10-30: (dfa) Updated website reference to florey.edu.au/imaging-software
% 2011-09-30: (dfa) Fixed syntax error when looking for iBrain in non-standard
%        	locations: I'd used IDL's ; for a comment instead of Matlab's %
% 2011-09-05: (dfa) Tries harder to find and start an iBrain server if job not accepted.
% 2011-08-31: (dfa) Added halt_on_error, timeout & status. Changed queue
%		folder name to more standard location for various platforms
%		(this now requires use of iBrain v5.4 or later).
%		Switched to iBT structure and updated comments for public release
% 2011-01-13: (dfa) Removed obsolete 'versn' from call to fileparts
% 2010-08-23: (dfa) Fixed when queue folder does not yet exist or is a file
% 2010-07-20: (dfa) Updated comments
% 2010-01-30: (dfa) Added additional message if waiting a while for completion
% 2009-02-27: Original version by David Abbott and Gagan Sharma

optargin = size(varargin,2);
if optargin < 1, use_server = 10; else 
 	if varargin{1} > 0, use_server = max( [varargin{1}, 5] ); else use_server = 0; end
end
if optargin < 2, halt_on_error = 1; else halt_on_error = varargin{2} ; end

start=pwd;
cd(dir);
 
if use_server == 0 
 
 command = ['iBrainr ' job];
 disp(command); % system(['echo ' command ]);
 [status] = system(command);
 if status ~= 0
 	disp(['Error starting external program: ' command]);
	status = -5;
	if halt_on_error, error('Halting due to error.'); exit; else cd(start); return; end
 end
  
else % use_server 

 if optargin < 3, timeout = 3153600000; else timeout = varargin{3}; end % Default 100 years!
 if optargin < 4, del_job = 1; else del_job = varargin{4}; end % Default deletes job on timeout
 
 command = ['iBrain ' job];
 disp(command);
 timeout = ceil(timeout); % Integers only
 if timeout < 1,timeout = 3153600000; end % Zero or negative effectively no timeout (100 years).
 
 polling_interval = min([5,timeout]); % seconds between check of iBrainServer status
 num_polls = fix ( timeout / polling_interval);
 timeout_residual = timeout - ( num_polls * polling_interval );
 if use_server > 1,  num_prestart_polls = ceil ( use_server / polling_interval);,  else num_prestart_polls = -1; end
 
   working_directory = pwd;

   % Check queue directory exists (make it if not)  

   if ispc(),
   	appdata = getenv('LOCALAPPDATA');
	if strcmp(appdata, ''), appdata = getenv('APPDATA'); end
   	queue_dir_name = fullfile(appdata,'iBrain_queue'); 
   else   
   	queue_dir_name = fullfile(getenv('HOME'),'.iBrain_queue');
   end 
   
   queue_dir_exist = exist(queue_dir_name,'dir');
   
   if queue_dir_exist ~= 7
   	% If the queue folder does not exsit already, make it.
	if exist(queue_dir_name,'file')
   		disp(['Removing existing file ' queue_dir_name]);
		delete(queue_dir_name);
		if exist(queue_dir_name,'file')
   			disp(['Error: Unable to removing existing file ' queue_dir_name]);
			status = -1;
			if halt_on_error, error('Halting due to error.'); exit; else cd(start); return; end
		end
	end
	[SUCCESS,MESSAGE,MESSAGEID] = mkdir(queue_dir_name);
	if SUCCESS == 0
   			disp(['Error: Unable to make job queue directory: ' queue_dir_name]);
			status = -2;
			if halt_on_error, error('Halting due to error.'); exit; else cd(start); return; end
	else
   			disp(['Made job queue directory: ' queue_dir_name]);
	end
   else
    			disp(['Using job queue directory: ' queue_dir_name]);
   end

   % Construct an iBrain pre-processing script, then move it to the queue directory
   % so that iBrain server will run it   

   [pathstr, name, ext] = fileparts(tempname);
   script_filename_ibc = [ name '_iBrainCommand.ibc'];
   script_filename_ibs = [ name '_iBrainCommand.ibs']; 

   [script_file_ibc, message] = fopen(script_filename_ibc,'w');
   if script_file_ibc < 0 
   		disp(['Error: Unable to create temporary command script ' script_filename_ibc ' in ' pwd]);
		status = -3;
		if halt_on_error, error('Halting due to error.'); exit; else cd(start); return; end
   end
   fprintf(script_file_ibc,'iBrain_DIR=%s\r\n',dir);
   
   fprintf(script_file_ibc,'%s\r\n',command);
   fclose(script_file_ibc);
   
   ibc_exist = exist(script_filename_ibc);
   
   if ibc_exist == 2 
   	copyfile(script_filename_ibc,queue_dir_name);
   elseif ibc_exist ~=2 
   	disp('Missing *.ibc file. Please check the script.');
	status = -3;
	if halt_on_error, error('Halting due to error.'); exit; else cd(start); return; end
   end
  
   delete( script_filename_ibc )
      
   cd(queue_dir_name);
   ibc_exist = exist( script_filename_ibc,'file' );
   i = 0;
   server_just_started = 0; % Initialise

   while ibc_exist == 2
   	if num_prestart_polls == 0 
		% iBrain Server has not accepted a job in time, so attempt to start a new iBrain server process
	  	disp(['An iBrain Server has not yet accepted the job: ' script_filename_ibc ] );
	  	disp(['Now attempting to start a new iBrain server process... ' ] );
		% Return to our start directory (in case there is an iBrain startup script there)
		cd(start);

		if ispc()
			iBrain_prog = 'iBrain.bat';
		else
			iBrain_prog = 'iBrain';
		end
		
		disp('  looking for iBrain in default path...');
				
		if ispc()
			% No "which" command on Windows, but can work around this as follows.
			% First check the current directory:
			if exist(iBrain_prog, 'file') == 2
				status = 0; % Found
			else % check the path
			    [status,result] = system('for %i in (iBrain.bat) do @echo.%~$PATH:i');
			    if status == 1
				disp('An unexpected error occured trying to check for the presence of iBrain.bat');
			    elseif fix(result(1)) == 10, % A line-feed return indicates nothing found. 
				status = 1; % Not found
			    else 
				status = 0; % Found	
			    end
			end % if exist('iBrain.bat' ...
		else
			[status,result] = system('which iBrain');
		end

		if status == 0 
			disp(['found iBrain command in default path.']);
			iBrain_prog_exists = 2;
		else % That didn't work, so see if we can find an iBrain installation...
		    if status == 1
			disp('  iBrain not in default path; checking other standard locations...');
			iBrain_prog_exists = 0; % initialise
			if ispc()
				progs=getenv('ProgramFiles');
				iBrain_prog = fullfile(progs, 'iBrain', 'bin', 'Windows', 'iBrain.bat');
				if exist(iBrain_prog,'file') 
				  iBrain_prog_exists = 1;
				else
				  % Non-standard
				  progs=getenv('SystemDrive');
				  iBrain_prog = fullfile(progs, 'programs', 'iBrain', 'bin', 'Windows', 'iBrain.bat');
				  if exist(iBrain_prog,'file') 
				    iBrain_prog_exists = 1;
				  else
				    % home directory (perhaps a non-administrative user has installed it here)
				    progs=getenv('UserProfile');
				    iBrain_prog = fullfile(progs, 'iBrain', 'bin', 'Windows', 'iBrain.bat');
				    if exist(iBrain_prog,'file') iBrain_prog_exists = 1; end  
				  end
				end
			elseif ismac()
				progs='/Applications';
				iBrain_prog = fullfile(progs, 'iBrain', 'bin', 'MacOSX', 'iBrain');
				if exist(iBrain_prog,'file') 
				  iBrain_prog_exists = 1;
				else
				  % Non-standard
				  progs='/usr/local';
				  iBrain_prog = fullfile(progs, 'iBrain', 'bin', 'MacOSX', 'iBrain');
				  if exist(iBrain_prog,'file') 
				    iBrain_prog_exists = 1;
				  else
				    % home directory (perhaps a non-administrative user has installed it here)
				    progs=getenv('HOME');
				    iBrain_prog = fullfile(progs, 'iBrain', 'bin', 'MacOSX', 'iBrain');
				    if exist(iBrain_prog,'file') iBrain_prog_exists = 1; end  
				  end
				end
			else % UNIX or GNU+Linux
				progs='/usr/local';
				iBrain_prog = fullfile(progs, 'iBrain', 'bin', 'linux', 'iBrain');
				if exist(iBrain_prog,'file') 
				  iBrain_prog_exists = 1;
				else
				  progs='/opt';
				  iBrain_prog = fullfile(progs, 'iBrain', 'bin', 'linux', 'iBrain');
				  if exist(iBrain_prog,'file') 
				    iBrain_prog_exists = 1;
				  else
				    % home directory (perhaps a non-administrative user has installed it here)
				    progs=getenv('HOME');
				    iBrain_prog = fullfile(progs, 'iBrain', 'bin', 'linux', 'iBrain');
				    if exist(iBrain_prog,'file') iBrain_prog_exists = 1; end  
				  end
				end
			end %if ispc() else...
		    else % if status == 1 else...
		    	disp(['  check returned unrecognised error code: ' num2str(status)]);
			iBrain_prog_exists = 0; % If you want to continue anyway, set this to 1
		    end % if status == 1 else...
		end %if status == 0 

		if iBrain_prog_exists == 1
%		  % Upon reflection, the following is somewhat excessive and 
% 		  % breaks backward compatibility so we won't use it.
%		  % Instead we'll just assume that what we have found works.
%		    command = ['"' iBrain_prog '" -exists']; % Check to see if the iBrain we found actually runs
% 		    [status] = system(command);
		    status = 0; % If ever re-enable the two lines above, then remove this line.
		    if status == 0, 
		  	iBrain_prog_exists = 2; 
			disp(['Found ' iBrain_prog ]);
		    end		  
		end% if iBrain_prog_exists == 1
		
		if iBrain_prog_exists == 1
	  		disp(['An iBrain startup file was found in ' iBrain_prog ] );
	  		disp(['however it failed to execute the command as expected.' ] );
	  		disp(['Please ensure that it is installed correctly ' ] );
	  		disp(['(the startup script changed in iBrain v5.4c so' ] );
	  		disp([' please be sure you are using a recent startup' ] );
	  		disp([' script and you have permission to run it).' ] );
		elseif iBrain_prog_exists == 0
	  		disp(['An iBrain startup file was not found.' ] );
		end
		
		if iBrain_prog_exists < 2
	  		disp(['The iBrain software package needs to be downloaded and ' ] );
	  		disp(['installed separately from the iBrain Toolkit for SPM. ' ] );
	  		disp(['See florey.edu.au/imaging-software. If you already have ' ] );
			disp(['iBrain, perhaps it is installed in a non-standard location. '] );
			disp(['In that case, please ensure the iBrain folder is in '] );
			disp(['your default path, OR make a symbolic link (not a copy) '] );
			disp(['of the iBrain startup script in a folder that is already in your path.'] );
			disp(['Alternatively, if you do not want to install iBrain, you may '] );
			disp(['restrict the selection of options in the runinfo to those that '] );
			disp(['do not require iBrain . ' ] );
			if halt_on_error, error('Halting due to error.'); exit; else cd(start); return; end
		end
		
		command = ['"' iBrain_prog '" -server &']; % Must be run in background so this script can continue.
 		disp(command); 
 		[status] = system(command);
		% It is pointless checking the status here because the system shell is run in 
		% the background so we don't get the return. 

		server_just_started = 1;

		num_prestart_polls = -1; % Don't try this again
		num_polls = num_polls + 1; % Add an extra polling interval to allow new 
					% server time to start and user to click spash screen
		cd(queue_dir_name); % Return to the queue directory
	end   
	
	if i == 0
	  disp(['Waiting for the iBrain Server to accept the job: ' script_filename_ibc ] );
	  i=10;
	end
	if num_polls < 1, 
		if timeout_residual > 0, 
			polling_interval = timeout_residual; 
			timeout_residual = 0;
		else
	  		disp(['Warning: iBrain Server failed to accept the job: ' script_filename_ibc ] );
	  		disp(['         in iBrain queue folder : ' pwd ] );
			if del_job == 1, delete( script_filename_ibc ); end
			cd(start);
			status = -4;
			if halt_on_error, error('Halting due to error.'); exit; else cd(start); return; end
			return
		end
	end; % if num_polls < 1	
	
  	 pause(polling_interval);
	 if num_prestart_polls > 0, num_prestart_polls = num_prestart_polls - 1;
	 else num_polls = num_polls - 1; end
	 ibc_exist = exist(script_filename_ibc,'file');
	 if (ibc_exist == 2 ) && (server_just_started == 1) 
	 	disp(['Please click on the IDL splash screen... ' ] );
		server_just_started = 0;
	 end
	 i = i-1;
   end %end while 

   disp('iBrain command accepted for processing');
   ibs_exist = exist(script_filename_ibs,'file');

   i = 0;
   j = 16;
   while ibs_exist ~= 2 % Keep polling until status file appears
  	 pause(polling_interval);        
	 ibs_exist = exist(script_filename_ibs,'file');
	 
    	 if i >= 10
   		disp('still waiting for iBrain command to complete...'); 
		i = 0;
		j = j + 1;
    	 	if j >= 20
   			disp(['   ( polling for status file: ' script_filename_ibs ' )' ] ); 
			% If iBrain crashes and you want this script to continue, you can look at the last command,
			% and its log, execute it manually if it failed, and if you do that successfully 
			% then echo 0 to the status file name above to tell this script to continue 
			% (and restart an iBrain server to accept new jobs)
			j = 0;
	 	end	
	 end	
	 i = i+1;
   end %end while 

   copyfile(script_filename_ibs,working_directory);
   delete(script_filename_ibs) % from queue directory
   
   cd(working_directory);
   fid = fopen(script_filename_ibs,'r');
   line = fgetl(fid);
   fclose(fid);
   if ~strcmp(line,'0') 
   	disp(['iBrain command failed with status code ' line]);
			status = -10;
			if halt_on_error, error('Halting due to error.'); exit; else cd(start); return; end
   else
    	disp('iBrain command completed successfully.');
   end 
   delete(script_filename_ibs)

end % else use_server

cd(start);
