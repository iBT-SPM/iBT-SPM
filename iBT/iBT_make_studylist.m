function iBT_make_studylist(make_new_path)
%
% Create a studylist from subject details contained in a runinfo file
%
% FORMAT iBT_make_studylist([make_new_path])
%
% This script can be called on its own to create a studylist.txt
% file from a given runinfo file. It doesn't matter if the runinfo
% has non-contiguous (or missing) subjects and/or sessions - the
% studylist will be created for all valid subjects and sessions. 
% This can be useful if you wish to migrate from
% an existing study stream specified within a runinfo.m file.
%
% make_new_path = 0 (default) or 1 % If 1, this script will regenerate 
%	the path to the raw location file by reading the header of the
%	first image in the specified rawloc and determining its
% 	new path location in our image database (using the same 
%	program that generates that filename). This currently works 
%	for GE images only and requires an external program 
%	(ximg_filename) that is not presently included in the public
%	release of this toolbox.
%___________________________________________________________________________
% Copyright 2009-2011 The Florey Institute of Neuroscience and Mental Health
%
% This file is included with the Integrated Brain Analysis Toolbox for SPM.
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
% 2011-10-14: (dfa) Improved compatibility with older runinfo files
% 2011-08-26: (dfa) Updated for compatibility with the public release of
%	the Integrated Brain Analysis Toolbox for SPM
%
% 2004-2009: iBT analysis scripts were developed by David Abbott
%     based upon the who & what structure concept of Tony Waites, with 
%     substantial contributions from Matt Harvey, Kaushik Bhaganagarapu,
%     and Richard Masterton.
%     In 2009 make_studylist.m was created by David Abbott.
%

dbstop if error % enter debug mode if an error occurs
format compact

if ~exist('make_new_path','var'), make_new_path = 0; end

% store current path
pwd_orig = pwd;

try
    [SpmVersion,SpmRelease]   = spm('Ver','',1); SpmRelease = str2double(SpmRelease);
catch
    error('SPM cannot be found in MATLAB path.');
end

spm('defaults','fmri'); % Start with standard defaults

clear iBT who what;

  switch (SpmVersion),
    case 'SPM2'
      P = spm_get(1,'runinfo*.m',{'Select runinfo file specifying WHO and WHAT'});
    otherwise
      [P, status] = spm_select(1,'^runinfo.*\.m$',{'Select runinfo file specifying WHO and WHAT'});
  end

if isempty(P), 
    disp('Make studylist job cancelled: No runinfo file specified.')
    return
end

switch (SpmVersion),
  case 'SPM2'
    runinfo_filename = P{1};
  otherwise
    runinfo_filename = P;
end
[pathstr filen] = fileparts(runinfo_filename);  %create new analysis dir where results will be stored

cd(pathstr)
%filen=[filen]
eval(filen);

if ~exist('iBT','var'), % Perhaps we have an old runinfo that doesn't use iBT...
	if ( ~exist('who','var') ) , % Nope, problem is more serious
		error('Unable to determine instructions from the specified runinfo file.')
	else % looks like an old-format runinfo - add required details to iBT structure:
		iBT.who.sub = who.sub;
		% note that we don't need the WHAT when creating a studylist
		clear who what;
	end
end

outfilename = [runinfo_filename '_studylist.txt']

sz_who = size(iBT.who.sub);
nsub = sz_who(2)
subs_done = 0; % Initialise
sub = 1; % Initialise
FID = fopen(outfilename,'w');

if FID == -1,

    disp(['Unable to create ' outfilename]);
    disp(['Please ensure you have permission to create this file, ']);
    disp(['or first copy ' runinfo_filename ]);
    disp(['into a folder where you do have write permission.']);
    
else 
    while subs_done < nsub
	try 
		ID = iBT.who.sub{sub};
		sub_valid = 1;
	catch
		sub_valid = 0;
		sub = sub + 1; % sub not valid so try incrementing
	end
	if sub_valid == 1
		sz_who_sub = size(iBT.who.sub{sub}.sess);
		nses = sz_who_sub(2);
		sessions_done = 0; % Initialise
		ses = 1; % Initialise
		while sessions_done < nses
			try 
				loc = iBT.who.sub{sub}.sess{ses}.loc;
				ses_valid = 1;
			catch
				ses_valid = 0;
				ses = ses + 1; % ses not valid so try incrementing
			end
			if ses_valid == 1
				try 
					MasterSession_flag = (iBT.who.sub{sub}.MasterSession == ses);
				catch 
					MasterSession_flag = 0;
                end
                rawloc =  iBT.who.sub{sub}.sess{ses}.rawloc;
				if make_new_path == 1 
					[status result] = unix(['./ximg_filename ' rawloc '/I.001']);
					if status == 0, rawloc = ['/data/mri/GE/' strtrim(result)]; end;
                end
				
				a = sprintf('%s %i %s %s %s %i %i\n', iBT.who.sub{sub}.ID, MasterSession_flag, rawloc, iBT.who.sub{sub}.sess{ses}.loc,...
					 iBT.who.sub{sub}.sess{ses}.task, iBT.who.sub{sub}.sess{ses}.pick(1), iBT.who.sub{sub}.sess{ses}.pick(end) );
				disp(a);
				fprintf(FID,'%s', a);
				sessions_done = sessions_done + 1;
				ses = ses + 1;
			end; % if ses_valid == 1
		end; %while sesions_done < nses
 		subs_done = subs_done + 1;
		sub = sub + 1;
	end; % if sub_valid == 1
    end; %while subs_done < nsub

    fclose(FID);
    disp(['Created ' outfilename]);
    disp(['iBT_make_studylist completed.']);
end; %if FID == -1, else...
