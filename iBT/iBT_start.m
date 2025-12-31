function iBT_start(varargin)
% This is the command that the GUI calls to run a processing job. 
%
% Whilst iBT_start can be called directly from the Matlab
% command line, it is probably better to use the ibt command instead,
% as ibt ensures the path is correct, and is fewer characters to type.
%
% FORMAT iBT_start
%	Prompt the user for runinfo file(s) to run.
%       Multiple runinfos can be specified, in which case they are
%       run as separate analyses. There is presently no interactive
%       mechanism to chain runinfo pairs or tupples in a single analysis.
%
% FORMAT iBT_start(filespec)
%	Run with runinfo filenames matching filespec. Each matching
%       runinfo file is run as a separate analysis.
%
% FORMAT iBT_start(filespec_part1, filespec_part2, ...)
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
%
%  See README_iBT.txt 
%  folder for more details.
%
%___________________________________________________________________________
% Copyright 2007-2016, 2020-21 The Florey Institute 
%                           of Neuroscience and Mental Health
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
% 2021-07-20: (dfa) Support chaining of pairs or tupples of runinfos
%                   Added diary of iBT_start
% 2020-11-06: (dfa) Support iBT.what.pre.norm.do = 3 (rigid)
% 2020-04-30: (dfa) Support passing runinfo filename to this function
% 2020-04-19: (dfa) Support iBT.what.li.own.do & iBT.what.li.norm.do
%                   Write to fully-qualified log filename so when toggle
%                     diary off/on to flush, continues to write to the same
%                     file even if we've changed directories during processing.
% 2016-05-04: (dfa) Log Matlab Version
% 2015-11-11: (dfa) Bugfix: when selecting multiple runinfos, settings were
%		    cleared only prior to the first job, so contamination
%		    of subsequent jobs could occur for variables not
%		    explicitly set in those jobs.
% 2014-09-01: (dfa) Log the studylist name also to the analysis log file. 
% 2014-05-23: (dfa) Log the studylist name to pre-processing log file.
% 2013-12-09: (dfa) Flush and turn off diary logfile when finished.
% 2013-10-18: (dfa) Support iBT.what.disp.OTHER
% 2013-10-04: (dfa) Support sequential execution of multiple runinfo files.
% 2013-10-03: (dfa) Forgot to log subject number for multisession analyses.
% 2013-08-12: (dfa) Support iBT.what.custom_stats.do
% 2013-08-07: (dfa) Tidied up a bit
% 2013-03-25: (dfa) Added iBT.iBTpath
% 2012-03-12: (dfa) Move version info to iBT_version();
% 2012-02-11: (dfa) Added iBT.what.debug.justReadRuninfo for debugging
%               & renamed non-user variable .what.nses to .what.last_ses
% 2011-10-19: (dfa) Only use addpath if not running pre-compiled Matlab.
% 2011-10-18: (dfa) Improved compatibility with older runinfo files
% 2011-08-26: (dfa) Updated for public release
%


[iBTversion,iBTrelease,iBTversionReleaseString] = iBT_version(1);
start = which('iBT_start');
[iBTpath,NAME,EXT] = fileparts(start);
clear start NAME EXT

dbstop if error % enter debug mode if an error occurs

format compact

diary( fullfile(pwd,['matlab_iBT_start' sprintf('-%.2i',fix(clock)) '.log'] ) )

% store current path
pwd_orig = pwd;
if (~isdeployed) % If we are not running in pre-compiled (stand-alone) mode...
  addpath(pwd_orig); % Add to path so that any custom analysis functions are defined in Matlab path.
end

[SpmVersion, SpmRelease, SpmVersionReleaseString, MatlabVersion, MatlabDate, MatlabArch] = iBT_SPMversion(2);

spm('defaults','fmri'); % Start with standard defaults

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
% section 1: read WHO and WHAT from runinfo*.m
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

if nargin==0,

 switch (SpmVersion),
   case 'SPM2'
	P = spm_get([0,Inf],'runinfo*.m',{'Select runinfo file(s) specifying WHO and WHAT'});
	num_runs = numel(P);
   otherwise
	[P, status] = spm_select(Inf,'^runinfo.*\.m$',{'Select runinfo file(s) specifying WHO and WHAT'});
	num_runs = size(P,1);
 end; %switch

 if isempty(P), 
    disp('Integrated Brain Analysis Toolbox job cancelled: No runinfo file specified.')
    return
 end % if isempty(P),

 n_filespecs = 1;
 
else

 n_filespecs = size(varargin,2);

% disp('varargin:')
% varargin
 
 for fsn = 1:n_filespecs
   filespec{fsn} = char(varargin{fsn});
   D{fsn}=dir(filespec{fsn}); % Find all files matching input filespec
   if isempty(D{fsn}),
      disp(['Integrated Brain Analysis Toolbox job cancelled: Unable to find a file matching ' filespec{fsn}])
      return
   end; %if isempty(D{fsn})
   num_runs_in(fsn) = numel(D{fsn});
 end % for fsn = 1:n_filespecs
 
 max_num_runs = max(num_runs_in);
 min_num_runsGT = min(num_runs_in(find(num_runs_in - 1)));
 if max_num_runs > ( min_num_runsGT + 1)
      disp(['Integrated Brain Analysis Toolbox job cancelled: Inconsistent number of file matches across filespecs; '])
      disp(sprintf('  each should be either 1 or the maximum number found = %i, whereas found:', max_num_runs))
      disp(num_runs_in)
      return
 end; %if max_num_runs
 
 num_runs = max_num_runs;
 
end % if nargin==0, ...else...

for irun = 1:num_runs

 clear iBT who what; % Clear iBT and in case a legacy runinfo has been used, also clear who and what.

 if nargin==0,
   switch (SpmVersion),
    case 'SPM2'
	runinfo_filename = P{irun};
    otherwise
	runinfo_filename = P(irun,:);
   end %switch
 else
	if num_runs_in(1) == 1,
 		runinfo_filename = fullfile(D{1}(1).folder, D{1}(1).name);
	else
 		runinfo_filename = fullfile(D{1}(irun).folder, D{1}(irun).name);
	end
 end; % nargin==0,

 if num_runs > 1,
 	disp(sprintf('Commencing iBT job %i of %i (%s)',irun,num_runs,runinfo_filename));
 else
 	disp(sprintf('Commencing iBT job %s',runinfo_filename));
 end

 [pathstr filen] = fileparts(runinfo_filename); 
 iBT.what.runinfo_filename = runinfo_filename;

 cd(pathstr)
 eval(filen);
 
 % If we are chaining runinfos, evaluate the additional ones here
 for fsn = 2:n_filespecs
	if num_runs_in(fsn) == 1,
 		runinfo_filename = fullfile(D{fsn}(1).folder, D{fsn}(1).name);
	else
 		runinfo_filename = fullfile(D{fsn}(irun).folder, D{fsn}(irun).name);
	end
 	disp(sprintf(' including additional instructions from  %s',runinfo_filename));

        [pathstr2 filen ext] = fileparts(runinfo_filename); 
	if strcmp(pathstr, pathstr2),
		iBT.what.runinfo_filename = [iBT.what.runinfo_filename ' & ' filen ext]
	else
		iBT.what.runinfo_filename = [iBT.what.runinfo_filename ' & ' runinfo_filename]
		pathstr = pathstr2;
	end
	
        cd(pathstr)
        eval(filen);

 end % for fsn = 2:n_filespecs

 if ~exist('iBT','var'), % Perhaps we have an old runinfo that doesn't use iBT...
	if ( ~exist('who','var') ) || ( ~exist('what','var') ), % Nope, problem is more serious
		error('Unable to determine instructions from the specified runinfo file.')
	else % Looks like an old-format runinfo - add details to iBT structure:
	     % We need to be careful here - when running from within a function some Matlab
	     % versions will call the what and who functions for the following assignments:
	     % iBT.what = what; iBT.who = who;
	     % even though an identical assignment command typed at the Matlab prompt will
	     % use the what and who variables when they exist.
	     % To work around this, we set the subfields so there is no ambiguity
		iBT.what.SpmVersion = SpmVersion;
		iBT.what.MatlabVersion = MatlabVersion;
		try iBT.what.prefer = what.prefer; catch, end
		try iBT.what.where = what.where; catch, end
		try iBT.what.pre = what.pre; catch, end
		try iBT.what.stats = what.stats; catch, end
		try iBT.what.cons = what.cons; catch, end
		try iBT.what.MELODIC = what.MELODIC; catch, end
		try iBT.what.li = what.li; catch, end
		try iBT.what.disp = what.disp; catch, end
		try iBT.what.fc = what.fc; catch, end
		try iBT.what.debug = what.debug; catch, end
		iBT.who.sub = who.sub; % This must exist, otherwise we can't do anything.
		clear who what;
	end
 end

 iBT.iBTversion = iBTversion;
 iBT.iBTrelease = iBTrelease;
 iBT.iBTversionReleaseString = iBTversionReleaseString;
 iBT.iBTpath = iBTpath;
 iBT.SpmVersion = SpmVersion;
 iBT.SpmRelease = SpmRelease;
 iBT.SpmVersionReleaseString = SpmVersionReleaseString;

 try
	temp = iBT.what.disp.analysis;
   	 disp(sprintf('WARNING: iBT.what.disp.analysis is set by user - this will be ignored.'))
 catch; end
 try
	temp = iBT.what.disp.type;
   	 disp(sprintf('WARNING: iBT.what.disp.type is set by user - this will be ignored.'))
 catch; end
 try
	temp = iBT.what.cons.type;
   	 disp(sprintf('WARNING: iBT.what.cons.type is set by user - this will be ignored.'))
 catch; end
 try
	temp = iBT.what.disp.loc;
   	 disp(sprintf('WARNING: iBT.what.disp.loc is set by user - this will be ignored.'))
 catch; end

 try
	temp = iBT.what.prefer.SPM;
 catch; 
	iBT.what.prefer.SPM = ''; % Default if not set is no preference.
 end

 if strcmp(iBT.what.prefer.SPM,'') == 0
   if strcmp(iBT.what.prefer.SPM,iBT.what.SpmVersion) == 0
   	  disp(sprintf('WARNING: iBT.what.prefer.SPM is set to %s but you are actually running %s.',iBT.what.prefer.SPM, iBT.what.SpmVersion))
   	  disp(sprintf('         Please either re-run with the preferred SPM version, or edit'))
   	  disp(sprintf('         the iBT.what.prefer.SPM setting in:'))
   	  disp(sprintf('         %s',iBT.what.runinfo_filename))
	  return
   end
 end

 % Backwards compatibility - set variables that may not exist (or have changed) in old or stipped-down runinfo files:
 try, iBT.what.pre.do; catch iBT.what.pre.do = 0; end
 try, iBT.what.stats.do; catch iBT.what.stats.do = 0; end
 try, iBT.what.cons.do; catch iBT.what.cons.do = 0; end
 try, iBT.what.custom_stats.do; catch iBT.what.custom_stats.do = 0; end
 
 try, global_li_flag=iBT.what.li.do; catch global_li_flag = 0; end % We previoulsy could not distinguish between own and norm space when requesting LI analysis.
 try, iBT.what.li.own.do; catch iBT.what.li.own.do = global_li_flag; end % Maintain backward compatibility
 try, iBT.what.li.norm.do; catch iBT.what.li.norm.do = global_li_flag; end % Maintain backward compatibility
 iBT.what.li.do = iBT.what.li.norm.do | iBT.what.li.own.do; 
 
 try, iBT.what.disp.do; catch iBT.what.disp.do = 0; end
 try, iBT.what.pre.norm.do; catch try iBT.what.pre.norm.do = iBT.what.pre.SpatialNorm; catch iBT.what.pre.norm.do = 0; end; end
 try, iBT.what.pre.norm.intra.do; catch try iBT.what.pre.norm.intra.do = iBT.what.pre.SpatialNorm_multiscan; catch iBT.what.pre.norm.intra.do = 0; end; end
 try, iBT.what.pre.norm.intra.niterations; catch try iBT.what.pre.norm.intra.niterations = iBT.what.pre.SpatialNorm_niterations; catch iBT.what.pre.norm.intra.niterations = 0; end; end
 try, iBT.what.pre.norm.intra.strategy; catch iBT.what.pre.norm.intra.strategy = 2; end
 try, iBT.what.pre.norm.BiasCorrect; catch try iBT.what.pre.norm.BiasCorrect = iBT.what.pre.BiasCorrect; catch iBT.what.pre.norm.BiasCorrect = 0; end; end
 try, iBT.what.pre.norm.WriteNorm; catch try iBT.what.pre.norm.WriteNorm = iBT.what.pre.WriteNorm; catch iBT.what.pre.norm.WriteNorm = 0; end; end
 try, iBT.what.stats.sessions.do; catch iBT.what.stats.sessions.do = 1; end % Default is to analyse each subject independently
 try, iBT.what.stats.subjects.do; catch iBT.what.stats.subjects.do = 0; end 
 try, iBT.what.stats.group.do; catch iBT.what.stats.group.do = 0; end 
 try, iBT.what.stats.own_space.do; catch iBT.what.stats.own_space.do = 0; end 
 try, iBT.what.stats.norm_space.do; catch iBT.what.stats.norm_space.do = 0; end 
 try, iBT.what.stats.sessions.which; catch iBT.what.stats.sessions.which = 0; end 
 try, iBT.what.stats.tasks_text_max_length; catch iBT.what.stats.tasks_text_max_length = 30; end 
 try, iBT.what.stats.number_dirs; catch iBT.what.stats.number_dirs = 1; end 
 try, iBT.what.stats.data_link; catch iBT.what.stats.data_link = 0; end 
 try, iBT.what.cons.do; catch iBT.what.cons.do = 0; end;

 if abs(iBT.what.pre.norm.do) == 2,
 	iBT.what.stats.normName = 'affine';
 elseif abs(iBT.what.pre.norm.do) == 3,
	iBT.what.stats.normName = 'rigid';
 else
	iBT.what.stats.normName = 'norm';
 end

 % Count subjects and sessions
 iBT.who.nsub = length(iBT.who.sub);			% number of subjects to process
 missing_sessions = 0;
 for sub=1:iBT.who.nsub
   iBT.who.sub{sub}.nses=length(iBT.who.sub{sub}.sess);	% number of sessions for this subject
   for ses=1:iBT.who.sub{sub}.nses % Set names for info files for each session for each subject
  	try
  	  iBT.who.sub{sub}.sess{ses}.infoFile = fullfile(iBT.who.sub{sub}.sess{ses}.loc, '_header_info.txt');
	catch
	  missing_sessions = missing_sessions + 1;
	end
   end; % for ses
 end; % for sub


 %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
 %
 % section n: call routines
 %
 % does the work
 %
 %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

 try, iBT.what.li.li_only; catch, 
	% Work out if not specified
	iBT.what.li.li_only = (iBT.what.li.do~=0) && ( (iBT.what.pre.do ==0) && (iBT.what.stats.do ==0) && (iBT.what.cons.do ==0) && (iBT.what.custom_stats.do ==0)  && (iBT.what.disp.do ==0) )
 end 

 try tmp = iBT.what.debug.justReadRuninfo; catch tmp = 0; end
 if tmp == 1
	disp('The variable iBT.what.debug.justReadRuninfo is set...' );
	throw(MException('iBT_start:debug', 'Execution stopped so that you can examine the iBT structure.' ));
 end

 if iBT.what.li.li_only == 0
  spm fmri; % Starts the user interface unless we are only running the laterality stuff.
  	  % We could also run the other stuff without it, but having it open 
          % lets us automatically generate and save postscript pages
	  % of output such as realignment parameter plot and normalisation summary.
 end

 % call preprocess routine
 if (iBT.what.pre.do == 1)  %
    diary( fullfile(pwd,['matlab_iBT_preprocess' sprintf('-%.2i',fix(clock)) '.log'] ) )
    disp(['Pre-processing data with parameters specified in ' iBT.what.runinfo_filename ]);
    if length(studylistfile) > 0
    	disp(['and studylist file: ' studylistfile]);
    else
    	disp(['(without a studylist file).']);   
    end
        iBT_preprocess(iBT);
 end   

 % call analysis
 if (iBT.what.stats.do == 1) || (iBT.what.cons.do == 1) || (iBT.what.custom_stats.do == 1) || (iBT.what.li.do == 1) || (iBT.what.disp.do == 1) %
  diary( fullfile(pwd,['matlab_iBT_analysis' sprintf('-%.2i',fix(clock)) '.log'] ) )
  disp(['Analysing data with parameters specified in ' iBT.what.runinfo_filename ]);
  if length(studylistfile) > 0
    	disp(['and studylist file: ' studylistfile]);
  else
    	disp(['(without a studylist file).']);   
  end
 
  if (iBT.what.stats.sessions.do ~= 0 ) %Analyse each session of each subject independently (single-subject, single-session design)
    iBT.what.disp.analysis = 1; % Used to inform iBT_display.m of the type of analysis to display

    if iBT.what.stats.sessions.do < 0  % if negative, -iBT.what.stats.sessions.do specifies a particular subject to analyse
   	start = - iBT.what.stats.sessions.do;
	finish = start;
    elseif iBT.what.stats.sessions.do == 1 
   	start = 1;
	finish = iBT.who.nsub;
    else
   	disp(sprintf('Error: iBT.what.stats.sessions.do = %d is greather than 1. Use negative numbers to specify a particular subject.',iBT.what.stats.sessions.do))
	return
    end
    disp(sprintf('Commencing single-subject single-session analyses...'))
    for sub=start:finish
      nses=length(iBT.who.sub{sub}.sess);		% number of sessions per subject
      disp(sprintf('Analysing subject %d (%s)...',sub,iBT.who.sub{sub}.ID))
      iBT.what.do_sub = sub; % just do this subject
      if iBT.what.stats.sessions.which == 0 
     	StartSes = 1; EndSes = nses; 
      else 
     	StartSes = iBT.what.stats.sessions.which; EndSes = StartSes; 
      end
      for ses=StartSes:EndSes
	iBT.what.first_ses = ses; % In loops, start with this session....
   	iBT.what.last_ses = ses; % .. and do no more.
	% First see if we're using direct paths to the analysis directory (e.g. as used in stripped-down runinfo_li script)
	try
		temp = iBT.who.sub{sub}.sess{ses}.spmloc; % If we are using runinfo_li
		empty_session = 0;
		unknown_do = 1;
	catch
		iBT.who.sub{sub}.sess{ses}.spmloc = ''; % Not specified - we determine later.
		empty_session = 1;
		unknown_do = 0;
	end %try
	
	if empty_session == 1
	  try
		temp = iBT.who.sub{sub}.sess{ses}.rawloc; % If this fails the session is skipped
		empty_session = 0;
	  catch
		empty_session = 1;
	  end %try
	end % if empty_session

	try, if iBT.what.disp.OTHER.do > 0,
			empty_session = 0; % i.e. Don't skip
			unknown_do = 1;
	     end
	catch;
	end

	if empty_session == 1 disp(sprintf('Ignoring empty session %d',ses)); end

	if (empty_session ~= 1) && ((iBT.what.stats.own_space.do == 1) || (iBT.what.stats.norm_space.do == 1) || (unknown_do == 1 ))

	  try
		iBT.what.stats.tasks = iBT.who.sub{sub}.sess{ses}.task{1}; % if not a cell this fails.
		ntasks = length(iBT.who.sub{sub}.sess{ses}.task);
    		for t=2:ntasks % concatenate all the task names into a string separated by underscores
			iBT.what.stats.tasks = [iBT.what.stats.tasks '_' iBT.who.sub{sub}.sess{ses}.task{t}];
		end % for		
	  catch
		iBT.what.stats.tasks = iBT.who.sub{sub}.sess{ses}.task; % failed above so must not be a cell.
	  end % try

	  if length(iBT.what.stats.tasks) > iBT.what.stats.tasks_text_max_length; % Shorten 
		  iBT.what.stats.tasks = [iBT.what.stats.tasks(1:(iBT.what.stats.tasks_text_max_length-1)) '_'];
	  end		  

	  % Own space and norm_space require separate calls to iBT_analysis:
	  if iBT.what.stats.own_space.do == 1			
	  	if iBT.what.stats.number_dirs == 1
   		  iBT.what.stats.dir = fullfile(iBT.what.where.awd,iBT.what.stats.subdir, 'own', iBT.who.sub{sub}.ID, sprintf('%02d_%s',ses,iBT.what.stats.tasks), '');
		else
   		  iBT.what.stats.dir = fullfile(iBT.what.where.awd,iBT.what.stats.subdir, 'own', iBT.who.sub{sub}.ID, iBT.what.stats.tasks, '');
		end
		clear iBT.what.MELODIC.loc
		try
		  iBT.what.MELODIC.loc = iBT.who.sub{sub}.sess{ses}.MELODIC.loc;
		end % try
		iBT.what.stats.norm = 0;
   		disp(sprintf('Analysing session %d of %d (own space)',ses,nses));
	  	iBT.what.disp.ses=ses; iBT_analysis(iBT);
	  end % if iBT.what.stats.own_space.do
	
	  if iBT.what.stats.norm_space.do == 1
	  	if iBT.what.stats.number_dirs == 1
  		  iBT.what.stats.dir = fullfile(iBT.what.where.awd,iBT.what.stats.subdir, iBT.what.stats.normName, iBT.who.sub{sub}.ID, sprintf('%02d_%s',ses,iBT.what.stats.tasks), '');
		else
  		  iBT.what.stats.dir = fullfile(iBT.what.where.awd,iBT.what.stats.subdir, iBT.what.stats.normName, iBT.who.sub{sub}.ID, iBT.what.stats.tasks, '');
		end
		try
		  iBT.what.MELODIC.loc = iBT.who.sub{sub}.sess{ses}.MELODIC.loc;
		end % try
		iBT.what.stats.norm = 1;
   		disp(sprintf('Analysing session %d of %d (spatialy normalised)',ses,nses))
	  	iBT.what.disp.ses=ses; iBT_analysis(iBT);
	  end %if iBT.what.stats.norm_space.do
	  
	  if unknown_do == 1 
  		iBT.what.stats.dir = iBT.who.sub{sub}.sess{ses}.spmloc; 
		try
		  iBT.what.MELODIC.loc = iBT.who.sub{sub}.sess{ses}.MELODIC.loc;
		end % try
   		disp(sprintf('Analysing session %d of %d',ses,nses))
	  	iBT.what.disp.ses=ses; iBT_analysis(iBT);
	  end %if ~strcmp(iBT.who.sub{sub}.sess{ses}.spmloc, '')
	  
	end % if empty_session != 1		

      end % for ses
    end % for sub
    disp(sprintf('Completed single-subject single-session analyses.'))

  end % if (iBT.what.stats.sessions.do ~= 0 )


  if (iBT.what.stats.subjects.do ~= 0) %Analyse each subject independently (single-subject, multi-session design)
    iBT.what.disp.analysis = 2; % Used to inform iBT_display.m of the type of analysis to display
    if iBT.what.stats.subjects.do < 0  
   	start = - iBT.what.stats.subjects.do;
	finish = start;
    elseif iBT.what.stats.subjects.do == 1 
   	start = 1;
	finish = iBT.who.nsub;
    else
   	disp(sprintf('Error: iBT.what.stats.subjects.do = %d is greather than 1. Use negative numbers to specify a particular subject.',iBT.what.stats.subjects.do))
	return
    end
    disp(sprintf('Commencing single-subject multi-session analyses...'))
    iBT.what.stats.tasks = 'multisession'; % failed above so must not be a cell.
    for sub=start:finish     
        disp(sprintf('Analysing subject %d (%s)...',sub,iBT.who.sub{sub}.ID))
   	iBT.what.last_ses=length(iBT.who.sub{sub}.sess);	% number of sessions for this subject
	iBT.what.first_ses = 1; % In loops, start with first session and include them all
	iBT.what.do_sub = sub;
	% Own space and norm_space require separate calls to iBT_analysis:
	if iBT.what.stats.own_space.do == 1 
    		iBT.what.stats.dir = fullfile(iBT.what.where.awd,iBT.what.stats.subdir, 'own', iBT.who.sub{sub}.ID, 'multi_session', '');
		iBT.what.stats.norm = 0;
		iBT.what.disp.ses=iBT.what.first_ses; iBT_analysis(iBT);
	end
	if iBT.what.stats.norm_space.do == 1 
    		iBT.what.stats.dir = fullfile(iBT.what.where.awd,iBT.what.stats.subdir, iBT.what.stats.normName, iBT.who.sub{sub}.ID, 'multi_session', '');
		iBT.what.stats.norm = 1;
		iBT.what.disp.ses=iBT.what.first_ses; iBT_analysis(iBT);
	end
    end % for sub
    disp(sprintf('Completed single-subject multi-session analyses.'))

  end % if (iBT.what.stats.subjects.do ~= 0)

  if (iBT.what.stats.group.do > 0) % Analyse subjects together (multi-subject, single-session fixed-effects design)
    iBT.what.disp.analysis = 4; % Used to inform iBT_display.m of the type of analysis to display.

    % For each session, we fake a multi-session design where each "fake session" is actually a different subject 
    fake_iBT = iBT;
    fake_iBT.what.do_sub = 1; % We pretend our fake sessions are for subject 1
    fake_iBT.what.first_ses = 1; % In loops, start with this subject
    fake_iBT.what.last_ses = iBT.who.nsub;    % and include them all
    nses=length(iBT.who.sub{1}.sess);	% number of sessions - assume same for all subjects
   		% Perhaps here we should find max number of sessions and elsewhere elegantly ignore subject sessions that don't exist? 
    if iBT.what.stats.group.do < 0  
   	start = - iBT.what.stats.group.do;
	finish = start;
    elseif iBT.what.stats.group.do == 1 
   	start = 1;
	finish = nses;
    else
   	disp(sprintf('Error: iBT.what.stats.subjects.do = %d is greather than 1. Use negative numbers to specify a particular subject.',iBT.what.stats.sessions.do));
	return
    end 
   
    if start > nses 
   	disp(sprintf('Error: Specified session (%d) does not exist for subject 1.',start));
	return
    end
   
    disp(sprintf('Commencing multi-subject single-session analyses...'))
    for ses=start:finish
   	clear fake_iBT.who.sub;
	for sub=1:iBT.who.nsub
		fake_iBT.who.sub{1}.sess{sub} = iBT.who.sub{sub}.sess{ses};
	end
	% Own space and norm_space require separate calls to iBT_analysis:
	if iBT.what.stats.own_space.do == 1 
     		fake_iBT.what.stats.dir = fullfile(iBT.what.where.awd,iBT.what.stats.subdir, 'own', 'group', '');
		fake_iBT.what.stats.norm = 0;
		iBT.what.disp.ses=fake_iBT.what.first_ses; iBT_analysis(fake_iBT);
	end
	if iBT.what.stats.norm_space.do == 1 
     		fake_iBT.what.stats.dir = fullfile(iBT.what.where.awd,iBT.what.stats.subdir, iBT.what.stats.normName, 'group', '');
		fake_iBT.what.stats.norm = 1;
		iBT.what.disp.ses=fake_iBT.what.first_ses; iBT_analysis(fake_iBT);
	end
    end % for ses
    disp(sprintf('Completed multi-subject single-session analyses.'))

  end % if (iBT.what.stats.group.do > 0)
 
 end % (iBT.what.stats.do == 1) 

 cd(pwd_orig)

 % Reset the following flag so if now run a different runinfo iBT doesn't get confused.
 iBT.what.li.li_only = 0;

 if num_runs > 1,
 	disp(sprintf('iBT job %i of %i (%s) completed!',irun,num_runs,iBT.what.runinfo_filename));
 else
 	disp(sprintf('iBT job completed!'));
 end
 
 diary off;
 
end % for irun = 1:num_runs

return  

