function iBT_preprocess(iBT)
% Perform preprocessing.
% This script is designed to be called by iBT_start to do preprocessing
%  (e.g. image conversion, header information extraction, slice timing,
%   realignment, bias correction, spatial normalisation and creation of
%   within-brain image masks.)
% FORMAT iBT_preprocess(iBT)
%  where 'iBT.what' is a structure indicating what is to be done,
%    and 'iBT.who' is a structure indicating subject specific details.%
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
% 2026-01-03: (dfa)  Restore accidentally deleted line finding .n_inputAnat
% 2025-12-28: (dfa)  Support locating and copying anatomy file from fMRIPrep
%                    via new iBT.who.sub{sub}.sess{ses}.wild.[r,n]_inputAnat
% 2025-11-04: (dfa)  Support option iBT.what.pre.DICOMconvert = 2, to have
%                      SPM DICOM conversion also write JSON sidecar files.
% 2024-10-17: (dfa)  Support option iBT.what.pre.header_info = 4, to read
%                      TR from NIfTI file.
% 2024-09-06: (dfa)  Support option iBT.what.pre.header_info = 3, to read
%                      TR from JSON sidecar
%                    Reformatted with spaces rather than mixed tabs & spaces
% 2021-10-14: (dfa)  Fix backwards compatibility for legacy
%                      iBT.what.pre.wild.[slice realign]
% 2021-09-21: (dfa)  Support   iBT.who.sub{sub}.sess{ses}.wild.raw
%                      and     iBT.who.sub{sub}.sess{ses}.wild.input
%                      and     iBT.who.sub{sub}.sess{ses}.wild.slice
%                      and     iBT.who.sub{sub}.sess{ses}.wild.realign
%                    If these are not set, then fall back to the legacy
%                    global iBT.what.pre.wild.[raw input slice realign]
%      Also now use iBT.who.sub{sub}.sess{ses}.wild.confounds_fmriprep
%        instead of iBT.who.sub{sub}.sess{ses}.fmriprep.wild
% 2021-06-17: (dfa)  Support iBT.who.sub{sub}.sess{ses}.fmriprep.wild
%                    Bug fix for basename provided to fslsplit (the bug was
%                      mitigated when using fsl versions earlier than 6.0.4)
% 2020-11-06: (dfa)  Support iBT.what.pre.norm.do = 3 (rigid-body).
%                    Support iBT.what.pre.norm.tempLabel
%                    Don't write lines to header info file that pertain to
%                    options for processing steps that are turned off (=0).
% 2020-05-21: (dfa)  Support iBT.what.pre.fmriprep.done. Also copy relevant
%                      .JSON and .TSV files to preproc dir if present in
%                      raw data folder when iBT.what.pre.FSLsplit == 1
% 2020-04-19: (dfa)  Bugfix: If BT.what.pre.mean is nonzero then spatially
%                      normalise to mean even if no realignment undertaken.
%                    Bugfix: Don't look for mean file ending with
%                      iBT.what.pre.wild.ref if realignment was not done.
% 2020-04-09: (dfa)  Switch to using iBT.what.pre.wild.input and
%                    new iBT.what.pre.FSLsplit to support original data
%                    in fslsplit supported formats including 4D NIfTI_GZ.
% 2019-07-08: (dfa)  Flush log file more frequently
%                    Rename slice timing 'correction' to 'interpolation'
% 2019-06-19: (dfa)  Quote first_dicom_file so robust to spaces. This now
%                    works when using iBrain in server mode (the default and
%                    highly recommended); it might not work otherwise (relying
%                    on iBrain in command-line mode) as that is UNIX shell
%                    dependent which often does not handle spaces well.
% 2019-05-24: (dfa)  If we will call iBrain, then call it right aways so the
%                     user can click on the IDL OK button immediately.
% 2018-09-19: (dfa)  Support iBT.what.pre.mean
% 2018-09-19: (dfa)  Support iBT.what.ext.SPM, iBT.what.ext.FSL and
%			iBT.what.ext.iB (image filename extensions)
% 2016-05-04: (dfa)  Log Matlab Version
% 2014-05-06: (dfa)  Better error reporting if no files found from which to
%                    read header
% 2013-11-08: (dfa)  iBT.what.pre.norm.templates_folder to override default.
%		     Support SPM12b (only in legacy "old normalisation" mode
%                    so presently no advantage in using SPM12b over SPM8 for
%		    preprocessing).
% 2013-09-26: (dfa)  Bugfix: multi-session intra-subject-normalisation with
%		     strategy=2 included an intermediate transform back to
%		     the final session rather than master session (which would
%		     not have been much of a problem except that the image to
%		     base inverse on is the master session so truncation
%		     of the field of view occurred if the two sessions
%		     were not originally in close alignment).
% 2013-09-23: (dfa)  Fixed cleanup when '?' wildcards used in filespec.
% 2013-09-09: (dfa)  Fixed crash when cleaning up after own-space analysis
% 2013-04-05: (dfa)  Support in-plane realignment (iBT.what.pre.Realign == 2)
% 		     Improved the smoothing log file name when
%		     an anisotropic smooth is specified.
% 2013-03-26: (dfa)  Improved support for bet on different platforms.
% 2013-03-21: (dfa)  Fix backward compatibility for runinfo's without
%			iBT.what.pre.slice_time_extra_gap set.
% 2013-03-04: (dfa)  Support iBT.what.pre.interleaved = 4 and added
%                    some descending slice-order options.
%		     Support iBT.what.pre.slice_time_extra_gap
%		     Added logging of inter-slice time spacing
%		     Fixed reference slice used when
%                    iBT.what.pre.target_time is set to 'first'
%                    & iBT.what.pre.interleaved is 3 (bug introduced 2013-02-19).
% 2013-02-20: (dfa & ct) Support BETwithinBrainMask flags to call FSL's brain
%  extraction tool to create within-brain mask (Reference: S.M. Smith. "Fast and
%  robust automated brain extraction. Human Brain Mapping, 2002:17(3):143-155)
%  See also http://fsl.fmrib.ox.ac.uk/fsl/fslwiki/BET
% 2013-02-19: (dfa) Support iBT.what.pre.interleaved = 3
% 2012-09-11: (dfa) Maintain backward compatibility when ImageSWIdrawnUpon
%		    is not set.
% 2012-08-16: (dfa) If specified, adjust SWI mask by coregistering image
%		    upon which mask was drawn to the current pre-processed
%		    mean image and applying transform to mask.
% 2012-06-14: (dfa) Work around bug in SPM if back-transformed SWI has
%                   values outside range [0:1].
% 2012-03-12: (dfa) Display version info from iBT_version()
% 2012-02-12: (dfa) if iBT.what.pre.sessions_together == 1 then now
%		  the realignment target used is from the session indicated
% 		  as the master session (i.e. iBT.who.sub{sub}.MasterSession ),
%		  rather than the 1st session; also slice-time logs,
%		  realignment parameters and masks are created in
%		  individual session folders rather than just the first.
%		No longer change some iBT fields here.
%		Recoded some loops to reduce code duplication.
% 2011-12-13: (dfa) Permit iBT.who.sub{}.sess{}.swi = '' or '0' to indicate none.
% 2011-12-06: (dfa) Correct n_iterations in header_info if
%                  re-generating it with iBT.what.pre.norm.do == -2.
% 2011-11-29: (dfa) Improved error reporting
% 2011-09-22: (dfa) Added option to tidy up within-brain mask
%		using milxview software from CSIRO:
%		(iBT.what.pre.Create_n_iBrainWithinBrainMask = 1)
% 2011-09-12: (dfa) Fixed failure of iBT.what.pre.norm.intra.do = 2
% 2011-09-12: (dfa) Fixed typo exists -> exist that caused a crash if
% 		checking for an existing header info file.
% 2011-09-05: (dfa) Updated for public release.
%

% 2004-2011: iBT analysis scripts were developed by David Abbott based
%            upon the who & what structure concept of Tony Waites, with
%            substantial contributions from Matt Harvey, Kaushik Bhaganagarapu,
%            Richard Masterton and some assistance from Mark Griffin.
%            In 2004 the script that ultimately became iBT_preprocess.m began as
%            an adaptation of a much simpler script originally by R. Henson 2003-05-20.
%
%
% Wish list : 	* Implement sub-sampling to enable retention of resolution
%                 when large voxels are used.
%               * Make iBrain conversion recognise the input filetype
%                 automatically (so could then run a mix of GE and DICOM
%                 conversions from the one script).
%

disp(sprintf('Entered preprocessing script'));

if isempty(iBT), error('Error: "iBT" structure is not defined.'), end

[iBTversion,iBTrelease,iBTversionReleaseString] = iBT_version(1); % Display version so it gets logged.
[iBT.what.SpmVersion, iBT.what.SpmRelease, iBT.what.SpmVersionReleaseString, MatlabVersion, MatlabDate, MatlabArch] = iBT_SPMversion(2);

spm('defaults','fmri'); % Ensure we are starting with standard defaults
if (strcmp(iBT.what.SpmVersion, 'SPM2'))
  global defaults;  % Give us access to the SPM defaults variable (requied only to change defaults.printstr later)
end

format compact

% Note that as listed below, most pre-processing items depend on all previous items
% having been completed. Exceptions are:
%   * only one of iBrainConvert or DICOMconvert is required
%   * WriteRealign and SmoothRealign can be omitted if you don't want to do
%     individual-space analysis.
%   * norm.do can be omitted if you don't want to do group analysis.

% Maintain compatibility with old runinfo files that may not have
% some new fields that are used in this function.
try, iBT.what.pre.hideName; catch, iBT.what.pre.hideName = 0; end % Default if not specified
try, iBT.what.pre.hideID; catch, iBT.what.pre.hideID = 0; end % Default if not specified
try, iBT.what.pre.Create_r_BETwithinBrainMask; catch, iBT.what.pre.Create_r_BETwithinBrainMask = 0; end % Default if not specified
try, iBT.what.pre.Create_n_BETwithinBrainMask; catch, iBT.what.pre.Create_n_BETwithinBrainMask = 0; end % Default if not specified
try, iBT.what.pre.slice_time_extra_gap; catch, iBT.what.pre.slice_time_extra_gap = 0; end % Default if not specified
try, templates_folder=iBT.what.pre.norm.templates_folder; catch templates_folder=''; end; % '' indicates use default location
try, iBT.what.ext.SPM; catch iBT.what.ext.SPM = '.img'; end % Used for images generated by SPM
try, iBT.what.ext.FSL; catch iBT.what.ext.FSL = '.img'; end % Used for images generated by MELODIC
try, iBT.what.ext.iB;  catch iBT.what.ext.iB  = '.img'; end % Used for images generated by iBrain
try, iBT.what.pre.mean;  catch iBT.what.pre.mean = 1; end %
try, iBT.what.pre.FSLsplit; catch iBT.what.pre.FSLsplit = 0; end % Input data is in single file NIfTI compressed format
try, iBT.what.pre.fmriprep.done; catch iBT.what.pre.fmriprep.done = 0; end % Input data is from fmriprep

try
  use_iBrainServer = iBT.what.pre.iBrainServer;
catch
  use_iBrainServer = 10; % 0 = No, 1= yes, >1= yes and attempt to start a new server if job not accepted after this many seconds.
end

% If only global wildcards were specified, set the equivalent subject/session specific variables:
for sub = 1:iBT.who.nsub
  for ses = 1:iBT.who.sub{sub}.nses
    
    try, iBT.who.sub{sub}.sess{ses}.wild.raw; catch
      try, iBT.who.sub{sub}.sess{ses}.wild.raw = iBT.what.pre.wild.raw; catch
        iBT.who.sub{sub}.sess{ses}.wild.raw = '*.nii'; %Default if unspecified
      end % inner try wild.raw
    end % outer try wild.raw
    
    
    % We originally used iBT.what.pre.wild.DICOM. Then we switched to a more
    % generic name iBT.what.pre.wild.input as the input images might not always be DICOM.
    % More recently we added the more flexible session-specific BT.who.sub{sub}.sess{ses}.wild.input:
    try, iBT.who.sub{sub}.sess{ses}.wild.input; catch
      try, iBT.who.sub{sub}.sess{ses}.wild.input = iBT.what.pre.wild.input; catch
        try, iBT.who.sub{sub}.sess{ses}.wild.input = iBT.what.pre.wild.DICOM; catch
          iBT.who.sub{sub}.sess{ses}.wild.input = '*.nii'; %Default if unspecified
        end % inner try wild.input.DICOM exists
      end % mid try pre.wild.input exists
    end % outer try iBT.who.sub{sub}.sess{ses}.wild.input exists
    
    try, iBT.who.sub{sub}.sess{ses}.wild.slice; catch
      try, iBT.who.sub{sub}.sess{ses}.wild.slice = iBT.what.pre.wild.slice; catch
        iBT.who.sub{sub}.sess{ses}.wild.slice = cat(2, 'a', iBT.who.sub{sub}.sess{ses}.wild.raw); %Default if unspecified
      end % inner try wild.slice
    end % outer try wild.slice
    
    try, iBT.who.sub{sub}.sess{ses}.wild.realign; catch
      try, iBT.who.sub{sub}.sess{ses}.wild.realign = iBT.what.pre.wild.realign; catch
        iBT.who.sub{sub}.sess{ses}.wild.realign = cat(2, 'ra', iBT.who.sub{sub}.sess{ses}.wild.raw); %Default if unspecified
      end % inner try wild.realign
    end % outer try wild.realign
    
  end % for ses
end % for sub

if ispc()
  iBT_BET = fullfile(iBT.iBTpath,'iBT_BET.bat');
else
  iBT_BET = fullfile(iBT.iBTpath,'iBT_BET');
end

do_iBrainConvert = iBT.what.pre.iBrainConvert;
do_DICOMconvert = iBT.what.pre.DICOMconvert;
do_SliceTiming = iBT.what.pre.SliceTiming;
do_iBrainRealignTarget = iBT.what.pre.iBrainRealignTarget;
do_BiasCorrect = iBT.what.pre.norm.BiasCorrect;
do_norm_intra = iBT.what.pre.norm.intra.do;
intra_niterations = iBT.what.pre.norm.intra.niterations;

if do_BiasCorrect == 1
  switch(iBT.what.SpmVersion),
    case {'SPM2','SPM5'}
      error('Error: This script implements iBT.what.pre.norm.BiasCorrect only with SPM8 or later.')
    otherwise %  SPM8 or later should be OK
      ;
  end % switch(iBT.what.SpmVersion)
  if do_norm_intra == 0
    do_norm_intra = 1;         % We do the BC step in the intra-norm section even if only one session
    intra_niterations = 1; % ...but we don't need to continue to iterate after the first pass.
  end
end

do_WriteRealign = iBT.what.pre.WriteRealign;
do_SmoothRealign = iBT.what.pre.SmoothRealign;
if iBT.what.pre.norm.do > 0, do_SpatialNorm = 1; else do_SpatialNorm = 0; end

do_WriteNorm = iBT.what.pre.norm.WriteNorm;
if do_SpatialNorm <= 0 do_WriteNorm = 0; end % Curently can't write if don't calculate
do_Cleanup = iBT.what.pre.Cleanup;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% subject specific parameters, change for each subject
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
cwd = iBT.what.where.awd;	% analysis working directory

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% analysis specific parameters, change if necessary
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

r_smth = iBT.what.pre.r_smooth_fwhm ;		% smoothing kernel(mm) for sr*
n_smth = iBT.what.pre.n_smooth_fwhm;		% smoothing kernel(mm) for snr*

target_time = iBT.what.pre.target_time;		% Set to 'first' or 'middle' (for slice timing, middle picks the slice acquired nearest to the middle of the TR)

pp_defaults.normalise.write = iBT.what.pre.norm.write; % The settings for normalisation that we will use in this pre-processing script

to_template.normalise.estimate = iBT.what.pre.norm.estimate;
intra_session.normalise.estimate = iBT.what.pre.norm.estimate;

if abs(iBT.what.pre.norm.do) == 2
  to_template.normalise.estimate.nits = 0; % Number of nonlinear iterations set to zero for affine-only normalisation to template
end

if abs(iBT.what.pre.norm.do) == 3
  to_template.normalise.estimate.nits = 0; % Number of nonlinear iterations set to zero for rigid-only normalisation to template
  to_template.normalise.estimate.regtype = 'rigid'; %  Normalisation to template constrained to be almost rigid (see spm_normalise(), spm_affreg() and spm_affine_priors(typ))
end

if abs(do_norm_intra) == 2
  intra_session.normalise.estimate.nits = 0; % Set to affine-only for intra-subject inter-session normalisation
end

if strcmp(templates_folder,'')
  switch(iBT.what.SpmVersion), %% Set default location of normalisation templates
    case {'SPM2','SPM5','SPM8'}
      templates_folder=fullfile(spm('Dir'),'templates');
    otherwise %  SPM12 or later
      templates_folder=fullfile(spm('Dir'),'toolbox','OldNorm');
  end % switch(iBT.what.SpmVersion)
end; %if strcmp(templates_folder,'')


VG = fullfile(templates_folder,iBT.what.pre.norm.temp);  % register to specified template located in spm templates directory

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Hopefully you won't need to change anything below this line
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

nsub = iBT.who.nsub;	% number of subjects to pre-process

% First check that all the specified directories exist
fprintf('Checking raw data directories...');
for sub = 1:iBT.who.nsub
  for ses = 1:iBT.who.sub{sub}.nses
    rwd = cd(iBT.who.sub{sub}.sess{ses}.rawloc);
    cd(rwd);
  end
end
clear rwd
fprintf('OK.\n');


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Will we need to use iBrain?
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% If we will need iBrain server, check that it is running now,
% (this will start it if necessary, so the user can click OK on the IDL splash screen and leave after that)

% First need to know if we will need to create header info:
if iBT.what.pre.header_info < 0 % check for existence
  if exist(iBT.who.sub{sub}.sess{ses}.infoFile,'file'),
    do_header_info = 0; % No need to re-do it
  else
    do_header_info = -iBT.what.pre.header_info; % Re-do it
  end
else
  do_header_info = iBT.what.pre.header_info;
end

% Start iBrain server if necessary, so the user can click OK on the IDL splash screen and leave the job alone after that...
if (use_iBrainServer ~= 0) && ( ((do_iBrainConvert == 1) && ~( do_DICOMconvert < 0)) || (do_header_info == 1) || (do_iBrainRealignTarget == 1) )
  job=['-version'];
  iBrain( job, pwd, use_iBrainServer );
end %if (use_iBrainServer...
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

cd(cwd);
%warning off MATLAB:MKDIR:DirectoryExists
%mkdir(iBT.what.where.pre.subdir);
%cwd = fullfile(cwd,iBT.what.where.pre.subdir,'');

fprintf('Analysis dir: %s\n',cwd)

cd(cwd);

for sub=1:nsub
  disp(sprintf('****** Preprocessing subject %d of %d...',sub,nsub))
  diary off; diary on; %Flush the current diary logfile
  
  nses=iBT.who.sub{sub}.nses;		% number of sessions per subject
  clear P
  clear V
  
  if iBT.what.pre.sessions_together == 1 % How big we'll define P...
    pnses = nses; % deal together with sessions from the same subject
  else
    pnses = 1; % We'll deal with each session independently
  end
  
  for ses=1:nses
    disp(sprintf('*** Converting and organising data for subject %d of %d, session %d of %d...',sub,nsub,ses,nses))
    
    fprintf('Unprocessed data dir: %s\n',iBT.who.sub{sub}.sess{ses}.rawloc)
    fprintf('Processed   data dir: %s\n',iBT.who.sub{sub}.sess{ses}.loc)
    
    % define wildcards
    switch (iBT.what.SpmVersion),
      case 'SPM2',
        input_wild.sub{sub}.sess{ses} = iBT.who.sub{sub}.sess{ses}.wild.input; % wildcard for raw DICOM images
        raw_wild.sub{sub}.sess{ses} = iBT.who.sub{sub}.sess{ses}.wild.raw; % wildcard for raw anz images
        ref_wild = iBT.what.pre.wild.ref; % wildcard for realignment target
        if (iBT.what.pre.iBrainRealignTarget == 0) || (iBT.what.pre.Realign == 0) || ( (iBT.what.pre.Realign < 0) && (iBT.what.pre.fmriprep.done == 1) ),
          realigned_mean_wild = ['mean' iBT.what.ext.SPM];
        else
          realigned_mean_wild = ['mean' iBT.what.pre.wild.ref];
        end
        if do_SliceTiming == 0
          slice_wild.sub{sub}.sess{ses}  = raw_wild.sub{sub}.sess{ses}; % wildcard for images to realign
          wslice_wild.sub{sub}.sess{ses}  = ['w' raw_wild.sub{sub}.sess{ses}]; % wildcard for normalised images
        else
          slice_wild.sub{sub}.sess{ses}  = iBT.who.sub{sub}.sess{ses}.wild.slice; % wildcard for slice timing corrected images (images to realign)
          wslice_wild.sub{sub}.sess{ses}  = ['w' iBT.who.sub{sub}.sess{ses}.wild.slice]; % wildcard for normalised images
        end
        if iBT.what.pre.Realign == 0
          realign_wild.sub{sub}.sess{ses} = slice_wild.sub{sub}.sess{ses}; % wildcard for images to normalise
        else
          realign_wild.sub{sub}.sess{ses} = iBT.who.sub{sub}.sess{ses}.wild.realign; % wildcard for realigned images (images to normalise)
        end
        mean_wild = ['mean*' iBT.what.ext.SPM]; % wildcard for mean image
      otherwise
        input_wild.sub{sub}.sess{ses} = spm_wildconvert(iBT.who.sub{sub}.sess{ses}.wild.input);
        raw_wild.sub{sub}.sess{ses} = spm_wildconvert(iBT.who.sub{sub}.sess{ses}.wild.raw); % wildcard for raw anz images
        ref_wild = spm_wildconvert(iBT.what.pre.wild.ref); % wildcard for realignment target
        if (iBT.what.pre.iBrainRealignTarget == 0) || (iBT.what.pre.Realign == 0) || ( (iBT.what.pre.Realign < 0) && (iBT.what.pre.fmriprep.done == 1) ),
          realigned_mean_wild = spm_wildconvert(['mean' iBT.what.ext.SPM]);
        else
          realigned_mean_wild = spm_wildconvert(['mean' iBT.what.pre.wild.ref]);
        end
        if do_SliceTiming == 0
          slice_wild.sub{sub}.sess{ses}  = raw_wild.sub{sub}.sess{ses} ; % wildcard for images to realign
          wslice_wild.sub{sub}.sess{ses}  = spm_wildconvert(['w' iBT.who.sub{sub}.sess{ses}.wild.raw]); % wildcard for normalised images
        else
          slice_wild.sub{sub}.sess{ses}  = spm_wildconvert(iBT.who.sub{sub}.sess{ses}.wild.slice); % wildcard for slice timing corrected images (images to realign)
          wslice_wild.sub{sub}.sess{ses}  = spm_wildconvert(['w' iBT.who.sub{sub}.sess{ses}.wild.slice]); % wildcard for normalised images
        end
        if iBT.what.pre.Realign == 0
          realign_wild.sub{sub}.sess{ses} = slice_wild.sub{sub}.sess{ses}; % wildcard for images to normalise
        else
          realign_wild.sub{sub}.sess{ses} = spm_wildconvert(iBT.who.sub{sub}.sess{ses}.wild.realign); % wildcard for realigned images (images to normalise)
        end
        mean_wild = spm_wildconvert(['mean*' iBT.what.ext.SPM]); % wildcard for mean image
    end % case 'SPM2'
    
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % Convert images to NIfTI or Analyze format
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    
    if ( do_iBrainConvert == 1 ) && ~( do_DICOMconvert < 0)
      
      if ( do_DICOMconvert == 0 )
        conv_format = 'Signa';
      elseif ( do_DICOMconvert == 1 ) || ( do_DICOMconvert == 2 )
        conv_format = 'DICOM';
      else
        error('Unknown DICOM conversion format code specified: %i',do_DICOMconvert);
      end
      
      [SUCCESS,MESSAGE,MESSAGEID] = mkdir(iBT.who.sub{sub}.sess{ses}.loc);
      if SUCCESS ~=1,
        pause(2) % in case the problem is a slow automounting filesystem
        [SUCCESS,MESSAGE,MESSAGEID] = mkdir(iBT.who.sub{sub}.sess{ses}.loc);
      end
      if SUCCESS ~=1, %OK, really failed.
        error(MESSAGE);
      end
      
      if use_iBrainServer == 0
        job = sprintf('-Convert %s Analyze_TRANSAXIAL "''%s''" "''%s%s''" \n',conv_format, iBT.who.sub{sub}.sess{ses}.rawloc ,iBT.who.sub{sub}.sess{ses}.loc, filesep());
      else
        job=['-Convert ' conv_format ' Analyze_TRANSAXIAL "' iBT.who.sub{sub}.sess{ses}.rawloc '" "' iBT.who.sub{sub}.sess{ses}.loc filesep() '"'];
      end
      iBrain( job, pwd, use_iBrainServer );
      
    end % if do_iBrainConvert
    
    if ( do_iBrainConvert == 0 ) && (( do_DICOMconvert == 1 ) || ( do_DICOMconvert == 2 )) %%%%%%%%%%%%%%
      switch (iBT.what.SpmVersion),
        case 'SPM2'
          dicom.fname=spm_get('Files',iBT.who.sub{sub}.sess{ses}.rawloc, input_wild.sub{sub}.sess{ses});
        otherwise
          dicom.fname=spm_select('FPList',iBT.who.sub{sub}.sess{ses}.rawloc, input_wild.sub{sub}.sess{ses});
      end
      
      dicom.task=iBT.who.sub{sub}.sess{ses}.task;
      d_hdr = spm_dicom_headers(dicom.fname);
      if(size(d_hdr,1) > 0)
        switch (iBT.what.SpmVersion),
          case 'SPM2'
            % In spm2 one can use a hacked spm_dicom_convert, although the
            % iBrain converter is usually still much better.
            spm_dicom_convert_new(d_hdr,dicom,iBT.who.sub{sub}.sess{ses}.loc,'all');
          otherwise
            % Should work for spm5, although iBrain is often better.
            % Recent releases of spm8 (r4010 or later) now seem to do a much
            % better job for Siemens DICOM files than before.
            mkdir(iBT.who.sub{sub}.sess{ses}.loc);
            tmp_dir = pwd;
            cd(iBT.who.sub{sub}.sess{ses}.loc);
            % We assume here that the DICOM files are already arranged into sensible
            % subdirectories (i.e. one image set per folder - at our site we have an
            % automated database system that does this).
            % See spm_dicom_convert.m from SPM8 for details of other possible conversion options.
            do_spm_meta = do_DICOMconvert-1;
            spm_dicom_convert(d_hdr,'all','flat','img',pwd,do_spm_meta);% the last value flag here, if true, exports 
                  %  metadata as JSON during DICOM to NIfTI conversion for SPM in r7201 (public release r7219 on 16-Nov-2017) 
                  % or later.
            cd(tmp_dir);
        end
      end
    end % if do_dicomConvert  %%%%%%%%%%%%%%%
    
    if ( iBT.what.pre.FSLsplit == 1 ) % The original data is already in a format supported by fslsplit (e.g. could be a single compressed NIfTI time series file (4D volume))
      switch (iBT.what.SpmVersion),
        case 'SPM2'
          input_fname=spm_get('Files',iBT.who.sub{sub}.sess{ses}.rawloc, input_wild.sub{sub}.sess{ses});
        otherwise
          input_fname=spm_select('FPList',iBT.who.sub{sub}.sess{ses}.rawloc, input_wild.sub{sub}.sess{ses});
      end %switch
      
      num_input_files = size(input_fname,1);
      if (num_input_files) == 0,
        error(['ERROR: No input image found matching ' input_wild.sub{sub}.sess{ses} ' in ' iBT.who.sub{sub}.sess{ses}.rawloc]);
      end; %if (size(input_fname,1)) == 0 else...
      disp(sprintf('Found %0i input files matching %s in %s', num_input_files, iBT.who.sub{sub}.sess{ses}.wild.input, iBT.who.sub{sub}.sess{ses}.rawloc));
      disp('Checking location of fslsplit command...');
      [status,result] = system('which fslsplit');
      if status ~=0,
        error(['ERROR: ' result]);
      else
        disp(['Found ' result]);
      end
      
      for i = 1:num_input_files
        disp(['Extracting volumes from ' input_fname(i,:)  ' into ' iBT.who.sub{sub}.sess{ses}.loc ' ...']);
        [SUCCESS,MESSAGE,MESSAGEID] = mkdir(iBT.who.sub{sub}.sess{ses}.loc);
        [dir,fname,ext] = fileparts(input_fname(i,:));
        if strcmpi(deblank(ext),'.gz'), [dir,fname,ext] = fileparts(fname); end; % ext was .gz so need to do it again to get real .nii extension off
        command = ['export FSLOUTPUTTYPE=NIFTI; fslsplit ' input_fname(i,:) ' ' fullfile(iBT.who.sub{sub}.sess{ses}.loc,fname)  ];
        disp(command)
        [status,result] = system(command);
        if status ~=0,
          error(['ERROR: ' result]);
        end
        % Also copy related BIDS JSON file if it exists
        jsonFilename = fullfile(dir,[fname,'.json']);
        if exist(jsonFilename, 'file')==2,
          [SUCCESS,MESSAGE,MESSAGEID] = copyfile(jsonFilename,iBT.who.sub{sub}.sess{ses}.loc);
          if SUCCESS,
            disp(['  also copied ' jsonFilename ]);
          else
            disp(['   WARNING:  Failed to copy ' jsonFilename ]);
          end; %if SUCCESS
        end;
      end; %for
      
      try
        iBT.what.stats.mc.wild;
        look_for_mc = 1;
      catch
        look_for_mc = 0;
      end; %catch
      if look_for_mc == 1,
        % Also copy related motion correction file if it already exists in the rawloc  (e.g. if this is an fmriprep output folder)
        if iBT.what.stats.mc.fmriprep == 1,
          try,
            mc_wild=iBT.who.sub{sub}.sess{ses}.wild.confounds_fmriprep;
          catch,
            try,
              mc_wild=iBT.who.sub{sub}.sess{ses}.fmriprep.wild; % We introduced this only in v3.9alpha03 and replaced with above in v3.9alpha06
            catch,
              mc_wild = iBT.what.stats.mc.wild;
            end
          end
        else
          mc_wild = iBT.what.stats.mc.wild; 	% wildcard for mc parameters
        end
        switch (iBT.what.SpmVersion),
          case 'SPM2'
            tsvFilename=spm_get('Files',iBT.who.sub{sub}.sess{ses}.rawloc, mc_wild);
          otherwise
            tsvFilename=spm_select('FPList',iBT.who.sub{sub}.sess{ses}.rawloc, spm_wildconvert(mc_wild));
        end %switch
        num_tsv_files = size(tsvFilename,1);
        if (num_tsv_files) == 1,
          if exist(tsvFilename(1,:), 'file')==2, [SUCCESS,MESSAGE,MESSAGEID] = copyfile(tsvFilename,iBT.who.sub{sub}.sess{ses}.loc); end;
          if SUCCESS,
            disp(['  also copied ' tsvFilename ]);
          else
            disp(['   WARNING:  Failed to copy ' tsvFilename ]);
          end; %if SUCCESS
        elseif (num_tsv_files) > 1,
          disp(['Error: Multiple files matching ' mc_wild ' found in ' iBT.who.sub{sub}.sess{ses}.rawloc ]);
          error(['Please adjust your wildcard specification to pick only the relevant file for this analysis stream']);
        end %if num_tsv_files
      end % if look_for_mc ==1,

      look_for_r_anat = 0; % initialise
      look_for_n_anat = 0; % initialise
      try tmp=iBT.who.sub{sub}.sess{ses}.wild.r_inputAnat;
          if ~strcmp(iBT.who.sub{sub}.sess{ses}.wild.r_inputAnat, ''), look_for_r_anat = 1; end
      catch; end
      try tmp=iBT.who.sub{sub}.sess{ses}.wild.n_inputAnat;
          if ~strcmp(iBT.who.sub{sub}.sess{ses}.wild.n_inputAnat, ''), look_for_n_anat = 1; end
      catch; end

      if look_for_r_anat == 1,
        % Copy (via fslsplit) related anatomical scan if we can find it in the usual place (../anat/)
        anat_rawloc = fullfile(iBT.who.sub{sub}.sess{ses}.rawloc,'..','anat');
        anatFilename = iBT_get( anat_rawloc, [iBT.who.sub{sub}.sess{ses}.wild.r_inputAnat '*']); % extra * wildcard to allow for .gz
        num_anat_files = size(anatFilename,1);
        if (num_anat_files) == 1,
          [dir,fname,ext] = fileparts(anatFilename);
          if strcmpi(deblank(ext),'.gz'), [dir,fname,ext] = fileparts(fname); end; % ext was .gz so need to do it again to get real .nii extension off
          disp(['Extracting anatomical image ' anatFilename  ' into ' iBT.who.sub{sub}.sess{ses}.loc ' ...']);
          command = ['export FSLOUTPUTTYPE=NIFTI; fslsplit ' anatFilename ' ' fullfile(iBT.who.sub{sub}.sess{ses}.loc,fname)  ];
          disp(command)
          [status,result] = system(command);	
        elseif (num_anat_files) > 1,
          disp(['Warning: Multiple files matching ' anat_wild ' found in ' anat_rawloc ]);
          error(['Please adjust your r_inputAnat wildcard specification to pick only the relevant file for this analysis stream']);
	end % if num_anat_files
      end % if look_for_r_anat ==1,

      if look_for_n_anat == 1,
        % Copy (via fslsplit) related anatomical scan if we can find it in the usual place (../anat/)
        anat_rawloc = fullfile(iBT.who.sub{sub}.sess{ses}.rawloc,'..','anat');
        anatFilename = iBT_get( anat_rawloc, [iBT.who.sub{sub}.sess{ses}.wild.n_inputAnat '*']); % extra * wildcard to allow for .gz
        num_anat_files = size(anatFilename,1);
        if (num_anat_files) == 1,
          [dir,fname,ext] = fileparts(anatFilename);
          if strcmpi(deblank(ext),'.gz'), [dir,fname,ext] = fileparts(fname); end; % ext was .gz so need to do it again to get real .nii extension off
          disp(['Extracting anatomical image ' anatFilename  ' into ' iBT.who.sub{sub}.sess{ses}.loc ' ...']);
          command = ['export FSLOUTPUTTYPE=NIFTI; fslsplit ' anatFilename ' ' fullfile(iBT.who.sub{sub}.sess{ses}.loc,fname)  ];
          disp(command)
          [status,result] = system(command);	
        elseif (num_anat_files) > 1,
          disp(['Warning: Multiple files matching ' anat_wild ' found in ' anat_rawloc ]);
          error(['Please adjust your n_inputAnat wildcard specification to pick only the relevant file for this analysis stream']);
	end % if num_anat_files
      end % if look_for_n_anat ==1,

      disp('Done.');
      
      if (do_header_info == 1), error('ERROR: iBrain header file creation requested from NIfTI files. This is not currenty supported.');end
      
    end % if iBT.what.pre.FSLsplit  %%%%%%%%%%%%%%%
    
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % Get information from image headers if we've just done conversion
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    
    if do_header_info > 0 %%%%%%%%%%%%%%
      
      dest_dir = iBT.who.sub{sub}.sess{ses}.loc;
      header_is_OK = 0; % initialise
      
      if do_header_info == 2 %%%%%%%%%%%%%%
        disp(sprintf(' Adding TR specified in iBT.who.sub{sub}.sess{ses}.TR to image information for subject %i session %i...',sub, ses))
        disp(sprintf(' WARNING: Manually specifying this information is subject to human error and so is not recommended.'))
        header_is_OK =1;
      elseif do_header_info == 3 %%%%%%%%%%%%%%
        disp(sprintf(' Adding TR specified in JSON sidecar to image information for subject %i session %i...',sub, ses))
        if iBT.what.pre.fmriprep.done ~= 0
          TR_wildcard = ['*task-*' iBT.who.sub{sub}.sess{ses}.name '*_bold.json']; 
          JSON_loc = iBT.who.sub{sub}.sess{ses}.rawloc; % fMRIprep's output folder
       else
          [filepath, name, ~] = fileparts(iBT.what.pre.wild.raw);
          TR_wildcard = fullfile(filepath, [name, '.json']);
          JSON_loc = iBT.who.sub{sub}.sess{ses}.loc; % SPM's output folder
        end
        TR = iBT_extract_TR(JSON_loc, TR_wildcard);  % Extract TR from JSON file(s)
        if TR == 0
          disp(['ERROR: Unable to determine TR from files matching ' fullfile(JSON_loc, TR_wildcard)]);
          error('Unable to continue without a specified TR.');
        else
          iBT.who.sub{sub}.sess{ses}.TR = TR;
          header_is_OK =1;
        end
      elseif do_header_info == 4 %%%%%%%%%%%%%%
        disp(sprintf(' Adding TR specified in NIfTI header to image information for subject %i session %i...',sub, ses))
        % Find an appropriate image from which to read header information
        source_dir = iBT.who.sub{sub}.sess{ses}.loc; % Look in the pre-processed folder
        nifti.fname=spm_select('FPList', source_dir, raw_wild.sub{sub}.sess{ses}); % Look for SPM-compatible images
        found_file = 0; % Initialise
        if ~isempty(nifti.fname)
          first_nifti_file = nifti.fname(1,:);
          if exist(first_nifti_file, 'file'),
            found_file = 1;
            Vol = spm_vol(first_nifti_file);  % Extract information from the NIfTI file
            if isfield(Vol.private,'timing') && isfield(Vol.private.timing,'tspace')
              TR = Vol.private.timing.tspace;
            else
              TR = 0; 
            end
            if TR == 0,
              disp('ERROR: Unfortunately it looks like your pre-processed images do not contain TR information.');
              disp('       Consider using an improved conversion tool, or choose another option for obtaining');
              disp('       the TR - see the iBT.what.pre.header_info flag.');
              error('Unable to continue without a specified TR.');
            end
            iBT.who.sub{sub}.sess{ses}.TR = TR;
            header_is_OK =1;
          end
        end % if length(nifti.fname) > 0
        if found_file == 0
          error(['Unable to find an image file matching' raw_wild.sub{sub}.sess{ses} ' in ' source_dir ]);
        end
      else % Use iBrain to read header information
        temp_header_file_m = fullfile(iBT.who.sub{sub}.sess{ses}.loc, '__temp_header.txt'); % Single-quoted version for use within matlab
        temp_header_file = ['"' temp_header_file_m '"']; % Double quoted version to pass to Unix commands
        
        disp(sprintf(' Reading image header for subject %i session %i...',sub, ses))
        [SUCCESS,MESSAGE,MESSAGEID] = mkdir(iBT.who.sub{sub}.sess{ses}.loc);
        if SUCCESS ~=1,
          error(MESSAGE);
        end
        header_info = 'raw';
        
        if iBT.what.pre.hideID ~= 0
          hh_switch = '-hhai ';
        elseif iBT.what.pre.hideName ~= 0
          hh_switch = '-hha ';
        else
          hh_switch = '-hh ';
        end
        
        while strcmp(header_info, 'done') ~= 1
          
          if strcmp(header_info, 'raw') == 1
            source_dir = iBT.who.sub{sub}.sess{ses}.rawloc;
          else
            source_dir = iBT.who.sub{sub}.sess{ses}.loc;
          end
          
          if exist(temp_header_file_m, 'file'),  %Make sure there is not a left-over present.
            delete(temp_header_file_m);
            if exist(temp_header_file_m, 'file'),
              error(['Unable to remove existing file: ' temp_header_file_m])
            end
          end % Make sure there is not a left-over present.
          
          % Find an appropriate image from which to read header information
          first_GE_file_m = fullfile(source_dir,'I.001'); % Single-quoted version for use within matlab
          first_GE_file = ['"' first_GE_file_m '"']; % Double quoted version to pass to Unix commands
          if exist( first_GE_file_m, 'file'),
            [SUCCESS,MESSAGE,MESSAGEID] = copyfile( first_GE_file_m, dest_dir,'f');
            if SUCCESS ~=1,
              error(MESSAGE);
            else
              job = [hh_switch first_GE_file ' -o ' temp_header_file ];
              iBrain( job, pwd, use_iBrainServer );
            end
          else %'I.001' does not exist
            first_dicom_file = ' ';
            switch (iBT.what.SpmVersion),
              case 'SPM2',
                dicom.fname=spm_get('Files', iBT.who.sub{sub}.sess{ses}.rawloc, input_wild.sub{sub}.sess{ses});
              otherwise
                dicom.fname=spm_select('FPList', iBT.who.sub{sub}.sess{ses}.rawloc, input_wild.sub{sub}.sess{ses});
            end
            found_file = 0 % Initialise
            if ~isempty(dicom.fname)
              first_dicom_file =  dicom.fname(1,:)
              if exist(first_dicom_file, 'file'),
                found_file = 1;
              end
            end
            
            if found_file == 1,
              [SUCCESS,MESSAGE,MESSAGEID] = copyfile(first_dicom_file, dest_dir);
              if SUCCESS ~=1,
                error(MESSAGE);
              else
                job = [hh_switch '"' first_dicom_file '" -o ' temp_header_file ];
                iBrain( job, pwd, use_iBrainServer );
              end % if SUCCESS ~=1
            else
              if strcmp(header_info, 'raw') == 1
                header_info = 'not_found_in_raw';
              else
                disp(sprintf('Error: Unable to find I.001 or %s in %s or %s.',first_dicom_file,iBT.who.sub{sub}.sess{ses}.rawloc, iBT.who.sub{sub}.sess{ses}.loc));
                error('Now terminating script execution prematurely - please fix problem and try again.');
              end
            end % if exist(first_dicom_file, 'file'), else
          end %if exist('I.001', 'file'), else
          
          if strcmp(header_info, 'not_found_in_raw') == 1
            header_info = 'preprocessed';
          else
            header_info = 'done';
          end
        end; % while strcmp(header_info, 'done') ~= 1
        
        ftemp_header_file = fopen(temp_header_file_m, 'r'); % open for reading the temporary header file just created
        % Most of the file consists of label:value pairs but there may be some junk at the start that we don't need.
        % Keep reading until we find first useful label =  "Image filename"
        buffer = fgets(ftemp_header_file);
        while buffer ~= -1 % Continue unless we are at the end of the file
          [token,TR_String] = strtok(buffer,':');
          if strcmp(token,'Image filename')
            header_is_OK = 1;
            break; %we have what we want, so break and start writing from here
          end
          buffer = fgets(ftemp_header_file);
        end % While
      end % if do_header_info
      
      if header_is_OK == 1,
        disp(sprintf(' Saving header info for subject %i session %i to %s...',sub, ses, iBT.who.sub{sub}.sess{ses}.infoFile))
        % Store all the clean info into information file (filename set in iBT_start.m)
        freal_header_file = fopen(iBT.who.sub{sub}.sess{ses}.infoFile, 'w'); % create final file (trashing any existing)
        % The first two lines are a fixed format identification and with the exception
        % of the minor version number should not be changed unless you change the reader.
        fprintf(freal_header_file, 'Info_version=\r\n'); % Version of the info file follows.
        if strcmp(to_template.normalise.estimate.regtype,'mni'), % This is the SPM default value
          fprintf(freal_header_file, '2.3\r\n');
        else % iBT began supporting a non-default regtype value only in v3.8h, so we increased the header
          %  major version so a warning is issued if this is encountered by an earlier iBT version.
          fprintf(freal_header_file, '3.0\r\n');
        end
        % Note: Increment major version (i.e. integer portion) only if make changes that require changes
        % when reading, including if changes could result in erroneous information display if interpreted
        % by older iBT versions.
        % If just adding extra discretionary fields, can increment the minor version number
        %   Note: From 2020-11-06 (iBT 3.8h) onwards, iBT_display will halt if an unfamiliar
        %         major version is encoundered, or will log a warning and continue if only the minor
        %         version is unfamiliar. Any unfamiliar fields will be ignored.
        %         Prior to that (versions of iBT_display between 2008-07-24 and 2020-11-05)
        %         if the version is >=3.0 then only a warning is logged and execution continues. Again,
        %         any unfamiliar fields are ignored.
        fprintf(freal_header_file, '--------------------------------------------\r\n');
        if do_header_info == 2
          fprintf(freal_header_file, 'Data obtained from iBT runinfo file\r\n');
          fprintf(freal_header_file, '--------------------------------------------\r\n');
          fprintf(freal_header_file, 'TR(s):%g\r\n',iBT.who.sub{sub}.sess{ses}.TR);
        elseif do_header_info == 3
          fprintf(freal_header_file, 'Data obtained from JSON sidecar\r\n');
          fprintf(freal_header_file, '--------------------------------------------\r\n');
          fprintf(freal_header_file, 'TR(s):%g\r\n',iBT.who.sub{sub}.sess{ses}.TR);
        elseif do_header_info == 4
          fprintf(freal_header_file, 'Data obtained from NIfTI header\r\n');
          fprintf(freal_header_file, '--------------------------------------------\r\n');
          fprintf(freal_header_file, 'TR(s):%g\r\n',iBT.who.sub{sub}.sess{ses}.TR);
        else
          fprintf(freal_header_file, 'Data obtained by iBrain from original image\r\n');
          fprintf(freal_header_file, '--------------------------------------------\r\n');
          while buffer ~= -1 % Continue until we reach the end of the file
            fprintf(freal_header_file, buffer);
            buffer = fgets(ftemp_header_file);
          end % While
          fclose(ftemp_header_file);
        end %if do_header_info == 2 else
        % Append some other useful information about the processing we are about to undertake
        fprintf(freal_header_file, '--------------------------------------------\r\n');
        fprintf(freal_header_file, 'Processing by iBT\r\n');
        fprintf(freal_header_file, '--------------------------------------------\r\n');
        fprintf(freal_header_file, 'iBT_version:%s\r\n', iBT.iBTversionReleaseString);
        fprintf(freal_header_file, 'SPM_version:%s\r\n', iBT.SpmVersionReleaseString);
        if iBT.what.pre.hideName
          %Write subject ID as specified in runinfo, instead of name obtained from image header.
          fprintf(freal_header_file, 'Subject name:%s\r\n', iBT.who.sub{sub}.ID);
        end
        if iBT.what.pre.SmoothRealign ~= 0,
          rSmooth = strcat('r_smth:', num2str(iBT.what.pre.r_smooth_fwhm));
          fprintf(freal_header_file, '%s\r\n', rSmooth);
        end
        if iBT.what.pre.norm.WriteNorm ~= 0,
          nSmooth = strcat('n_smth:', num2str(iBT.what.pre.n_smooth_fwhm));
          fprintf(freal_header_file, '%s\r\n', nSmooth);
        end
        if iBT.what.pre.norm.do ~= 0,
          nTemplate = strcat('n_template:', iBT.what.pre.norm.temp);
          fprintf(freal_header_file, '%s\r\n', nTemplate);
          if isfield(iBT.what.pre.norm,'tempLabel'),
            nTemplateLabel = strcat('n_templateLabel:', iBT.what.pre.norm.tempLabel);
            fprintf(freal_header_file, '%s\r\n', nTemplateLabel);
          end
          nRegtype = (['n_regtype:' to_template.normalise.estimate.regtype])
          fprintf(freal_header_file, '%s\r\n', nRegtype);
          nIts = strcat('n_iterations:', num2str(to_template.normalise.estimate.nits));
          fprintf(freal_header_file, '%s\r\n', nIts);
        end
        fclose(freal_header_file);
      else
        disp(sprintf('Unable to find "Image filename:" label in temporary header information file: %s', temp_header_file))
        fclose(ftemp_header_file);
        error('Now terminating script execution prematurely - please fix problem and try again.');
      end %if header_is_OK
      
      if do_header_info == 1, delete(temp_header_file_m); end
      
    end % Get information from image headers
    
    % Now check that the preprocessing data directory exists in case we reached here without doing
    % anything above (e.g. if variables set to indicate above had been done already).
    
    [success,message,messageid] = fileattrib(iBT.who.sub{sub}.sess{ses}.loc);
    if success
      if message.directory == 0
        error(['Error: Not a directory: ' iBT.who.sub{sub}.sess{ses}.loc]);
      end
    else
      error(['Error:  directory does not exist: ' iBT.who.sub{sub}.sess{ses}.loc]);
    end
    
    
    if do_SliceTiming == 1 % we'll need to know the TR, so read it now.
      
      % Now read TR from saved settings file:
      disp(sprintf(' Reading TR from saved header info... '))
      [present, value] = iBT_info(iBT.who.sub{sub}.sess{ses}.infoFile,'TR(s)',1);
      if present,
        TRs{ses} = value;
      else
        error(['Unable to find "TR(s):" in: ' iBT.who.sub{sub}.sess{ses}.infoFile]);
      end
      disp(sprintf(' TR is %g seconds.',TRs{ses}));
      
    end % if do_SliceTiming == 1 % we'll need to know the TR, so read it now.
    
    % If there is a header override file, copy this to the data directory for subesquent use.
    try
      copyfile('_header_override.txt', fullfile(iBT.who.sub{sub}.sess{ses}.loc,''));
    catch
    end
    
    cd(cwd)
    
    disp(sprintf('*** Converting and organising data completed for subject %d of %d, session %d of %d',sub,nsub,ses,nses))
    
  end %for ses %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  
  
  if iBT.what.pre.sessions_together == 1
    loop_ses = 1; % we loop within sections below to collect data together
  else
    loop_ses = nses; % one big loop for separate processing
  end
  
  for ses=1:loop_ses
    
    if iBT.what.pre.sessions_together == 1
      % These are used to facilitate code that loops ses in
      % small sections below when doing sessions_together but
      % just does a single value of ses otherwise.
      ses_start=1;
      ses_end=nses;
      % We are already in a for ses for loop but it is really
      % just a one-value loop. It is therefore OK to enter
      % and leave another ses for loop so log as the value at the
      % end isn't less than the value at the start.
      % This allows us to reduce code dupliction
      % between sessions_together and not
      % sessions_together cases.
    else
      ses_start=ses;
      ses_end=ses;
      % We are already in a for loop that increments ses,
      % but it is OK to enter another ses for loop when we
      % start and stop on the same value of ses. This reduces
      % code duplication between sessions_together and
      % not sessions_together cases.
    end
    
    if (iBT.what.pre.sessions_together == 1) && (ses_end > ses_start)
      disp(sprintf('*** Preprocessing subject #%d (all sessions %d - %d)...',sub,ses_start,ses_end))
    else
      disp(sprintf('*** Preprocessing subject #%d session #%d...',sub,ses))
    end
    
    diary off; diary on; %Flush the current diary logfile
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % Slice-timing interpolation - must do before spatial realignment.
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    
    if do_SliceTiming == 1
      
      cd(iBT.who.sub{sub}.sess{ses}.loc);			% where spm*.ps is output
      
      disp(sprintf('Subject #%d, retrieving filenames for slice timing interpolation...',sub))
      
      diary off; diary on; %Flush the current diary logfile
      
      P = cell(1,1); % Irrespective of sessions_together, we always do slice timing separately
      %                since result of one image does not depend on any other, and by doing separately the
      %                log files for each session remain in that sessions directory.
      
      for ses=ses_start:ses_end % Just does one session unless iBT.what.pre.sessions_together; Using the existing loop variable ses within another loop here is intentional,
      %                            see comments near "if iBT.what.pre.sessions_together" above.
        switch (iBT.what.SpmVersion),
          case 'SPM2',
            P{1}   = spm_get('Files', iBT.who.sub{sub}.sess{ses}.loc, raw_wild.sub{sub}.sess{ses});
          otherwise
            P{1}   = spm_select('FPList', iBT.who.sub{sub}.sess{ses}.loc, raw_wild.sub{sub}.sess{ses});
        end
        if isempty(P{1})
          throw(MException('iBT_preprocess:slice_timing',['No files found matching "' raw_wild.sub{sub}.sess{ses} '" in ' iBT.who.sub{sub}.sess{ses}.loc] ));
        end
        
        V = spm_vol(P);
        
        num_slices = V{1}(1).dim(3);
        
        slice_time = (TRs{ses} - iBT.what.pre.slice_time_extra_gap) / num_slices; % time between slices
        slice_gap  = slice_time + iBT.what.pre.slice_time_extra_gap; % time between last slices and next volume; continuous scanning if iBT.what.pre.slice_time_extra_gap is 0.
        
        if iBT.what.pre.interleaved == 4,
          if ( (num_slices/2) - fix(num_slices/2) ) == 0.5, %then odd number of slices
            interleaved = 1; % Use this scheme when an odd total number of slices acquired
          else
            interleaved = 3; % Use this scheme when an even total number of slices acquired
          end
        elseif iBT.what.pre.interleaved == 14,
          if ( (num_slices/2) - Fix(num_slices/2) ) == 0.5, %then odd number of slices
            interleaved = 11; % Use this scheme when an odd total number of slices acquired
          else
            interleaved = 13; % Use this scheme when an even total number of slices acquired
          end
        else
          interleaved = iBT.what.pre.interleaved;
        end
        
        if interleaved == -1
          slice_ord = iBT.what.pre.slice_order;
        elseif interleaved == 0
          slice_ord  = [1:1:num_slices];	% Ascending not interleaved: 1,2,3,...
        elseif interleaved == 10
          slice_ord  = [num_slices:-1:1];	% descending not interleaved: ...,3,2,1
        elseif interleaved == 1
          slice_ord  = [1:2:num_slices 2:2:num_slices];	% single ascending interleaved. e.g. 5 slices would be: 1,3,5,2,4
        elseif interleaved == 11
          slice_ord  = [num_slices:-2:1 num_slices-1:-2:1];	% single descending interleaved. e.g. 5 slices would be: 5,3,1,4,2
        elseif interleaved == 2
          
          tmp  = [1:2:num_slices 2:2:num_slices];
          slice_ord = [1];
          for slc=2:num_slices
            slice_ord = [ slice_ord tmp(tmp(slc)) ]; % double ascending interleaved: 1, 5, 9, ..., 4, 8, 12, ..., 3, 7, 11, ..., 2, 6, 10,...
          end
        elseif interleaved == 3
          slice_ord  = [2:2:num_slices 1:2:num_slices];	% single ascending interleaved with second-bottom slice acquired first e.g. 5 slices would be: 2,4,1,3,5
        elseif interleaved == 13
          slice_ord  = [num_slices-1:-2:1 num_slices:-2:1]; % single descending interleaved with second-top slice acquired first e.g. 5 slices would be 4,2,5,3,1
        else
          error(sprintf('Unsupported value specified for iBT.what.pre.interleaved: %0i', interleaved ));
        end
        
        if size(slice_ord,2) ~= num_slices
          error(sprintf('Number of elements in slice_order array ( %0i ) differs from number of slices ( %0i )',size(slice_ord,2),num_slices));
        end
        
        if strcmp(target_time, 'middle')
          ref_slice = slice_ord( fix(length(slice_ord) /2) )
        elseif strcmp(target_time, 'first')
          ref_slice = slice_ord(1)
        else
          error(['Unsupported value specified for iBT.what.pre.target_time: ' target_time ' (only "first" and "middle" are presently supported)']);
        end
        
        disp(sprintf('Subject #%d, session #%d slice timing interpolation using the following slice order:',sub,ses))
        disp(slice_ord) % Print it so we can see
        disp(sprintf(' using %s to reference slice %i (iBrain %s interleaving scheme = %i).\n',iBT.what.SpmVersionReleaseString,ref_slice,iBT.iBTversionReleaseString,interleaved));
        disp(sprintf(' Time between start of slice acquisitions within a volume = %f seconds\n',slice_time));
        disp(sprintf(' Time between start of last acquired slice of a volume and first slice of the next = %f seconds\n',slice_gap));
        
        spm_slice_timing(P, slice_ord, ref_slice, [slice_time slice_gap]);
        
        flag_filename = fullfile(iBT.who.sub{sub}.sess{ses}.loc,['_slice_timed_to_',sprintf('%i',ref_slice),'.log']);
        flag_file =  fopen( flag_filename , 'w' );
        fprintf(flag_file,'Slice time corrected with %s to reference slice %i (iBrain %s interleaving = %i).\n',iBT.what.SpmVersionReleaseString,ref_slice,iBT.iBTversionReleaseString,interleaved);
        fprintf(flag_file,' i.e. assumed acquisition slice order = \n');
        fprintf(flag_file,'%i ',slice_ord);
        fprintf(flag_file,'\nTime (seconds) between start of slice acquisitions within a volume = \n');
        fprintf(flag_file,'%f \n',slice_time);
        fprintf(flag_file,'Time (seconds) between start of last acquired slice of a volume and first slice of the next = \n');
        fprintf(flag_file,'%f \n',slice_gap);
        fprintf(flag_file,'\n ');
        fclose( flag_file );
        
      end % for ses=ses_start
      
      disp(sprintf('Subject #%d, slice timing interpolation completed.',sub))
      
      diary off; diary on; %Flush the current diary logfile
      
    end % if do_SliceTiming
    
    if do_iBrainRealignTarget == 1
      
      if iBT.what.pre.sessions_together == 1
        cd (iBT.who.sub{sub}.sess{iBT.who.sub{sub}.MasterSession}.loc); % We create an alignment target for the master session only.
      else
        cd (iBT.who.sub{sub}.sess{ses}.loc); % This instance is only processing a single session: ses.
      end
      disp(['Determining realignment target for ' iBT.who.sub{sub}.sess{ses}.wild.raw ' in ' pwd]);
      diary off; diary on; %Flush the current diary logfile
      if do_SliceTiming == 0
        job=['-clobber -AlignTarget "' iBT.who.sub{sub}.sess{ses}.wild.raw '"'];
      else
        job=['-clobber -AlignTarget "' iBT.who.sub{sub}.sess{ses}.wild.slice '"'];
      end
      filespec_to_realign = slice_wild.sub{sub}.sess{ses}; % This was set earlier to a value dependent on do_SliceTiming
      iBrain( job, pwd, use_iBrainServer );
      cd (cwd)
      
    end % if do_iBrainRealignTarget
    
    
    disp(sprintf('Retrieving filenames in preparation for further processing...'))
    diary off; diary on; %Flush the current diary logfile
    P = cell(1,pnses); % Will contain elements for each session if iBT.what.pre.sessions_together
    for ses=ses_start:ses_end % Just does one session unless iBT.what.pre.sessions_together; Using the existing loop variable ses within another loop here is intentional,
      %                            see comments near start of first loop near "if iBT.what.pre.sessions_together" many lines above.
      if iBT.what.pre.sessions_together == 1, P_index = ses; else P_index = 1; end; % If processing sessions separately, P will always only have 1 element.
      switch (iBT.what.SpmVersion),
        case 'SPM2',
          P{P_index}   = spm_get('Files',iBT.who.sub{sub}.sess{ses}.loc, ref_wild);
          P{P_index}   = strvcat(P{P_index}, spm_get('Files', iBT.who.sub{sub}.sess{ses}.loc, slice_wild.sub{sub}.sess{ses}));
        otherwise
          P{P_index}   = spm_select('FPList',iBT.who.sub{sub}.sess{ses}.loc, ref_wild);
          P{P_index}   = strvcat(P{P_index}, spm_select('FPList', iBT.who.sub{sub}.sess{ses}.loc, slice_wild.sub{sub}.sess{ses}));
      end
      if isempty(P{P_index})
        disp(['No files found matching "' ref_wild ' or ' slice_wild.sub{sub}.sess{ses} '" in ' iBT.who.sub{sub}.sess{ses}.loc])
        throw(MException('iBT_preprocess:Retrieving_filenames',['Unable to continue without any image files.'] ));
      end
    end % for ses=ses_start
    
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % Spatial realignment and mean image generation
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    
    if iBT.what.pre.Realign > 0
      
      if iBT.what.pre.sessions_together == 1
        cd (iBT.who.sub{sub}.sess{iBT.who.sub{sub}.MasterSession}.loc) ; %Since we realign all sessions for same subject to the target
        %in the subject's master session, we'll output the spm*.ps there.
      else
        cd (iBT.who.sub{sub}.sess{ses}.loc) ; %output spm*.ps here
      end
      
      
      %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
      % Reload slice-time corrected image information
      
      disp(sprintf('Subject #%d, Loading alignment target and slice-time corrected image information...',sub))
      diary off; diary on; %Flush the current diary logfile
      
      V       = spm_vol(P);
      V	 = cat(1,V{:});
      
      %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
      % Realignment (coregister only)
      
      disp(sprintf('Subject #%d, Coregistering & Reslicing...',sub))
      diary off; diary on; %Flush the current diary logfile
      if iBT.what.pre.Realign == 2
        realign_flags.lkp = [1 2 6];% Fit only x & y translations and a rotation about z for "in-pane" realignment
        % (this is still a volumetric alignment, restricted to paremeters noted above, because
        %  an entire volume is fitted with a single translation & rotation. Ideally one should
        %  be able to do slice-by-slice realignment but SPM does not support this.)
        disp('Realignment options set to 2D (in-plane) only.');
      else
        realign_flags.lkp = 1:6; % Fit all parameters (conventional 3D realignment using 3 translationas and 3 rotations)
      end
      spm_realign(V,realign_flags);% SPM8 defaults to a two-pass realign-to-mean option, which is good.
      disp(sprintf('OK'))
      diary off; diary on; %Flush the current diary logfile
      
      %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
      % reslice and/or mean generation
      
      FlagsR = struct('interp',iBT.what.pre.realign.write.interp,...
        'wrap',iBT.what.pre.realign.write.wrap,...
        'mask',1,...
        'which',2*do_WriteRealign,'mean',1);
      % which=0...no reslice, only mean (if mean=1)
      % which=2 ...reslice all images
      
      disp(sprintf('Subject #%d, Reslicing & generating mean...',sub))
      diary off; diary on; %Flush the current diary logfile
      
      if iBT.what.pre.sessions_together == 1
        % In addition to writing out the realigned mean and (optionally) the individual realigned images
        % we also need to extract relevant realignment parameters from the rp_*  file in the master session
        % and save them to a separate rp file in ses's own folder for use later by stats analysis.
        % To avoid confusion later, we rename the original rp_*.txt file as master_rp_*.txt.
        
        switch (iBT.what.SpmVersion),
          case 'SPM2',
            rp_filename = spm_get('Files', iBT.who.sub{sub}.sess{iBT.who.sub{sub}.MasterSession}.loc, 'rp_*.txt'); % All parameters are in this one file.
          otherwise
            rp_filename = spm_select('FPList', iBT.who.sub{sub}.sess{iBT.who.sub{sub}.MasterSession}.loc, spm_wildconvert('rp_*.txt'));
        end % switch
        if  isempty(rp_filename)
          disp('Unable to find realignment parameters: Did something fail in the realignment of all the sessions together?');
          throw(MException('iBT_preprocess:realign',['Unable to find rp_*.txt in ' iBT.who.sub{sub}.sess{iBT.who.sub{sub}.MasterSession}.loc ] ));
        end
        if size(rp_filename,1) > 1
          disp(['WARNING: Found more that one file matching rp_*.txt in' iBT.who.sub{sub}.sess{iBT.who.sub{sub}.MasterSession}.loc]);
          throw(MException('iBT_preprocess:realign','Please remove extra file(s) before re-running processing.' ));
        end
        
        [rp_pathstr,rp_fname,rp_ext] = fileparts(rp_filename);
        master_rp_filename_unqualified = ['master_' rp_fname rp_ext];
        master_rp_filename = fullfile(rp_pathstr,master_rp_filename_unqualified) ;
        disp([' Renaming ' rp_fname rp_ext ' to master_' rp_fname rp_ext]);
        [SUCCESS,MESSAGE,MESSAGEID] = movefile(rp_filename,master_rp_filename,'f');
        if SUCCESS == 0
          throw(MException('iBT_preprocess:realign',[MESSAGE ':' MESSAGEID] ));
        end
        
        RPs = load(master_rp_filename,'-ascii');
        RP_start = 1; % Initialise for session 1.
        
        rp_name = rp_fname(4:end); % Remove the "rp_" prefix from the original name
        
        for ses=1:nses; % Using the existing loop variable ses within another loop here is intentional,
          %               see comments near start of first loop near "if iBT.what.pre.sessions_together" many lines above.
          
          spm_reslice(P{ses},FlagsR);
          
          num_files = size(P{ses},1); % Number of files realigned in this session
          if num_files > 0
            % extract relevant realignment parameters for this session
            % and save to an rp file in ses's own folder so can be used by stats analysis later.
            Q = RPs(RP_start:RP_start+num_files-1,:);
            new_rp_filename = fullfile(iBT.who.sub{sub}.sess{ses}.loc,sprintf('rp_toSession%03i_%s.txt',iBT.who.sub{sub}.MasterSession,rp_name));
            save(new_rp_filename,'Q','-ascii');
            disp(sprintf('Realignment parameters for session %i from %s in session %i saved to',ses,master_rp_filename_unqualified,iBT.who.sub{sub}.MasterSession));
            disp(sprintf('  %s', new_rp_filename));
            diary off; diary on; %Flush the current diary logfile
          end; % if ses > 1
          RP_start = RP_start + num_files;
          
        end % for
        
      else % Sessions are being realigned separately so just write out the present result:
        
        spm_reslice(P{1},FlagsR);
        
      end %if iBT.what.pre.sessions_together else
      
      disp(sprintf('OK'))
      
      switch(iBT.what.SpmVersion),
        case 'SPM2',
          global defaults; % Give us access to the SPM defaults variable.
          original_printstr = defaults.printstr;
          defaults.printstr = [spm_figure('DefPrintCmd'), 'spm2_realign.ps' ];
          spm_print;
          defaults.printstr = original_printstr;
        otherwise
          spm_print([iBT.what.SpmVersion '_realign.ps']);
      end
      
      cd (cwd)
      
    end % if iBT.what.pre.Realign
    
    
    if (iBT.what.pre.Realign <= 0 ) && (iBT.what.pre.mean > 0)
      % If we have't used SPM to realign, then we usually still need to geneate a mean image
      disp(sprintf('Subject #%d, Calculating mean...',sub))
      diary off; diary on; %Flush the current diary logfile
      cd (iBT.who.sub{sub}.sess{ses}.loc) ; %output spm*.ps here
      V = spm_vol(P);
      V = cat(1,V{:});
      Vo = ['mean' iBT.what.ext.SPM];
      Vo = spm_imcalc(V,Vo,'mean(X)',{1,0,0,4,'spm - algebra (mean)'});
      disp(sprintf('OK'))
      diary off; diary on; %Flush the current diary logfile
      cd (cwd)
    end % if iBT.what.pre.mean > 0
    
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % Create within-brain mask
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    
    if (iBT.what.pre.Create_r_iBrainWithinBrainMask == 1 ) || (iBT.what.pre.Create_r_BETwithinBrainMask == 1 )
      
      %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
      % Create within-brain mask of mean realigned image
      
      for ses=ses_start:ses_end % Using the existing loop variable ses within another loop here is intentional, see comments soon after start of the parent loop above.
        
        switch (iBT.what.SpmVersion),
          case 'SPM2',
            [realigned_mean,dir] = spm_get('Files', iBT.who.sub{sub}.sess{ses}.loc, realigned_mean_wild);
          otherwise,
            [dir,realigned_mean,ext] = fileparts(spm_select('FPList', iBT.who.sub{sub}.sess{ses}.loc, realigned_mean_wild));
            realigned_mean = strcat(realigned_mean,ext);
        end % switch
        
        disp(sprintf('Subject #%d, Calculating within-brain mask...',sub))
        diary off; diary on; %Flush the current diary logfile
        
        if size(realigned_mean,1) > 0
          realigned_mean = realigned_mean(1,:) % in case more than one - we'll take only the first
          iBrain_realigned_mean_mask = fullfile(dir,['iBrainMask_' realigned_mean])
          BET_realigned_mean_mask = fullfile(dir,['BETmask_' realigned_mean])
          
          cd(iBT.who.sub{sub}.sess{ses}.loc);
          if (iBT.what.pre.Create_r_iBrainWithinBrainMask == 1) && (iBT.who.sub{sub}.sess{ses}.iBrainWithinBrainMask.create ~= 0)
            job=['-clobber -wbt ' sprintf('%f',iBT.who.sub{sub}.sess{ses}.iBrainWithinBrainMask.create) ' -wbm "' realigned_mean '" "' iBrain_realigned_mean_mask '"'];
            iBrain(job,pwd,use_iBrainServer);
          end
          if iBT.what.pre.Create_r_BETwithinBrainMask == 1
            if ispc()
              disp('BET is not available on Windows PC platforms. Please use a GNU/Linux platform if you need BET.')
            else
              command=[iBT_BET  ' "' realigned_mean '" "' BET_realigned_mean_mask '" ' iBT.who.sub{sub}.sess{ses}.BETwithinBrainMask.create];
              [status,result] = system(command);
              disp(result)
              if status ~=0
                error(['Error: The following script returned an error - you may need to customise it for your site: ' iBT_BET ]);
              end
              % The command above generates two images - a masked brain with filename BET_normalised_mean_mask,
              % and the mask itself with a '_mask' appended to name. We only want the latter, so now delete the
              % former so wildcards in later scripts only get the _mask version:
              delete(BET_realigned_mean_mask);
              delete(strrep(BET_realigned_mean_mask,'.img','.hdr'));
            end %if ispc else
          end
          cd(cwd);
        end % if file found i.e. size(realigned_mean,1) > 0
        
        %if iBT.what.pre.sessions_together == 1
        
        % Here it would be nice if we also created a mean of all
        % images accross all sessions if iBT.what.pre.sessions_together == 1
        % or possibly  a mean of the session means (depending upon whether we
        % weight all scans equally, or all session equally irrespective of
        % differences in number of scans in each session). This would be useful for
        % overlays of multisession stats results. Code to do this should be
        % added above, along with appropriate adjustment to make use of it in display.
        %
        % Pseudo-code:
        % load mean from present session, add to a running total (possibly multiplied by number of scans in session if want each scan to contribute equally)
        % then, outside this loop, divide result by number of sessions (or scans) and save as mean_all_sessions (or scans).
        % Make sure running total data-type is long so don't get overflow.
        
        %end %if iBT.what.pre.sessions_together == 1
        
      end % for ses=ses_start:ses_end
      
    end % iBT.what.pre.Create_r_iBrainWithinBrainMask || (iBT.what.pre.Create_r_BETwithinBrainMask == 1 )
    
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % Smooth realigned images
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    
    if do_SmoothRealign == 1
      
      %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
      % Reload data into structure
      
      disp(sprintf('Subject #%d, Loading re-aligned image information...',sub))
      diary off; diary on; %Flush the current diary logfile
      
      cd(cwd);
      P = cell(1,pnses);
      if iBT.what.pre.sessions_together == 1
        for ses=1:nses % Using the existing loop variable ses within another loop here is intentional, see comments soon after start of the parent loop above.
          switch (iBT.what.SpmVersion),
            case 'SPM2',
              P{ses}   = spm_get('Files',iBT.who.sub{sub}.sess{ses}.loc, realign_wild.sub{sub}.sess{ses});
            otherwise
              P{ses}   = spm_select('FPList',iBT.who.sub{sub}.sess{ses}.loc, realign_wild.sub{sub}.sess{ses});
          end
        end
      else
        switch (iBT.what.SpmVersion),
          case 'SPM2',
            P{1}   = spm_get('Files',iBT.who.sub{sub}.sess{ses}.loc,realign_wild.sub{sub}.sess{ses});
          otherwise
            P{1}   = spm_select('FPList',iBT.who.sub{sub}.sess{ses}.loc, realign_wild.sub{sub}.sess{ses});
        end
      end
      
      %  We can't just load the slice-timed images here, because even though the realignment
      %  parameters are in the mat files, we want the transform applied to the image output (so other programs can
      %  use the data without needing to support .mat files. Therefore we need to re-load the written re-aligned
      %  image information here prior to smoothing. -dfa
      
      V       = spm_vol(P);
      V       = cat(1,V{:});
      
      %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
      % Smoothing
      
      disp(sprintf('Subject #%d, Smoothing...',sub))
      diary off; diary on; %Flush the current diary logfile
      
      for i=1:size(V,1),
        [pth,nam,ext] = fileparts(V(i).fname);
        fname         = fullfile(pth,['s' nam ext]);
        spm_smooth(V(i),fname,r_smth);
      end;
      disp(sprintf('OK'))
      
      if iBT.what.pre.sessions_together == 1
        for ses=1:nses % Using the existing loop variable ses within another loop here is intentional, see comments soon after start of the parent loop above.
          fclose( fopen( fullfile(iBT.who.sub{sub}.sess{ses}.loc,['_r_smoothed_to',sprintf('_%.2f',r_smth),'mm']) , 'W' ) );
        end
      else
        fclose( fopen( fullfile(iBT.who.sub{sub}.sess{ses}.loc,['_r_smoothed_to',sprintf('_%.2f',r_smth),'mm']) , 'W' ) );
      end
      
    end % if do_SmoothRealign
    
    cd(cwd)
    
    
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % Spatial normalisation
    % (without intra-subject inter-session norm)
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    
    if ((do_SpatialNorm == 1) && (do_norm_intra == 0)),
      
      %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
      % Reload data into structure
      
      disp(sprintf('Subject #%d, Re-loading non-smoothed slice-timed image information...',sub))
      
      if iBT.what.pre.sessions_together == 1
        cd (iBT.who.sub{sub}.sess{iBT.who.sub{sub}.MasterSession}.loc) ; %Since we realign all sessions for same subject to the target
        %in the subject's master session, we'll output the spm*.ps there.
      else
        cd (iBT.who.sub{sub}.sess{ses}.loc) ; %output spm*.ps here
      end
      
      % disp(sprintf('Retrieving filenames in preparation for spatial normalisation...'))
      diary off; diary on; %Flush the current diary logfile
      P = cell(1,pnses); % Will contain elements for each session if iBT.what.pre.sessions_together
      for ses=ses_start:ses_end % Just does one session unless iBT.what.pre.sessions_together; 
      %                           Using the existing loop variable ses within another loop here is intentional, see comments soon after start of the parent loop above.
        if iBT.what.pre.sessions_together == 1, P_index = ses; else P_index = 1; end; % If processing sessions separately, P will always only have 1 element.
        switch (iBT.what.SpmVersion),
          case 'SPM2',
            P{P_index}   = spm_get('Files',iBT.who.sub{sub}.sess{ses}.loc, slice_wild.sub{sub}.sess{ses});
          otherwise
            P{P_index}   = spm_select('FPList',iBT.who.sub{sub}.sess{ses}.loc, slice_wild.sub{sub}.sess{ses});
        end
        if isempty(P{P_index})
          throw(MException('iBT_preprocess:Retrieving_filenames',['No files found matching "' ref_wild ' or ' slice_wild.sub{sub}.sess{ses} '" in ' iBT.who.sub{sub}.sess{ses}.loc] ));
        end
      end % for ses=ses_start
      
      %  Now re-load the slice-timed data, along with the realignment
      %  transformation that is now in the .mat files:
      V       = spm_vol(P);
      V       = cat(1,V{:});
      
      switch (iBT.what.SpmVersion),
        case 'SPM2',
          if iBT.what.pre.sessions_together == 1
            if (iBT.what.pre.Realign == 0) && (iBT.what.pre.mean == 0), % no mean so use first image
              disp(sprintf('Subject #%d, Loading first image information from master session...',sub))
              meanf   = spm_get('Files',iBT.who.sub{sub}.sess{iBT.who.sub{sub}.MasterSession}.loc, slice_wild.sub{sub}.sess{iBT.who.sub{sub}.MasterSession});
              meanf   = spm_vol(meanf);
              meanf   = meanf(1).fname;
            else
              disp(sprintf('Subject #%d, Loading mean image information from master session...',sub))
              meanf   = spm_get('Files',iBT.who.sub{sub}.sess{iBT.who.sub{sub}.MasterSession}.loc, mean_wild);
            end
          else
            if (iBT.what.pre.Realign == 0) && (iBT.what.pre.mean == 0), % no mean so use first image
              disp(sprintf('Subject #%d, Loading first image information from current session...',sub))
              meanf   = spm_get('Files',iBT.who.sub{sub}.sess{ses}.loc, slice_wild.sub{sub}.sess{ses});
              meanf   = spm_vol(meanf);
              meanf   = meanf(1).fname;
            else
              disp(sprintf('Subject #%d, Loading mean image information from current session...',sub))
              meanf   = spm_get('Files',iBT.who.sub{sub}.sess{ses}.loc, mean_wild);
            end
          end
        otherwise
          if iBT.what.pre.sessions_together == 1
            if (iBT.what.pre.Realign == 0) && (iBT.what.pre.mean == 0), % no mean so use first image
              disp(sprintf('Subject #%d, Loading first image information from master session...',sub))
              meanf   = spm_select('FPList',iBT.who.sub{sub}.sess{iBT.who.sub{sub}.MasterSession}.loc, slice_wild.sub{sub}.sess{iBT.who.sub{sub}.MasterSession});
              meanf   = spm_vol(meanf);
              meanf   = meanf(1).fname;
            else
              disp(sprintf('Subject #%d, Loading mean image information from master session...',sub))
              meanf   = spm_select('FPList',iBT.who.sub{sub}.sess{iBT.who.sub{sub}.MasterSession}.loc, mean_wild);
            end
          else
            if (iBT.what.pre.Realign == 0) && (iBT.what.pre.mean == 0), % no mean so use first image
              disp(sprintf('Subject #%d, Loading first image information from current session...',sub))
              meanf   = spm_select('FPList',iBT.who.sub{sub}.sess{ses}.loc, slice_wild.sub{sub}.sess{ses});
              meanf   = spm_vol(meanf);
              meanf   = meanf(1).fname;
            else
              disp(sprintf('Subject #%d, Loading mean image information from current session...',sub))
              meanf   = spm_select('FPList',iBT.who.sub{sub}.sess{ses}.loc, mean_wild);
            end
          end
      end
      
      diary off; diary on; %Flush the current diary logfile
      Vm      = spm_vol(meanf);
      
      
      %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
      % Estimate unwarping parameters
      
      disp(sprintf('Subject #%d, Determining normalisation parameters...',sub))
      matname  = [spm_str_manip(Vm.fname,'sd') '_sn.mat'];
      disp(sprintf('  Normalisation template is %s',VG));
      
      try
        swi = iBT.who.sub{sub}.sess{ses}.swi;
      catch
        iBT.who.sub{sub}.sess{ses}.swi = 0; % Maintain backwards compatibility
      end % try
      
      if iBT.who.sub{sub}.sess{ses}.swi ~= 0
        Vswi = spm_vol(iBT.who.sub{sub}.sess{ses}.swi);
        disp(sprintf('  Source Weighted Image is %s',Vswi.fname));
      else
        Vswi = ''; % No SWI
      end
      
      diary off; diary on; %Flush the current diary logfile
      
      params   = spm_normalise(VG,Vm,matname,'',Vswi,to_template.normalise.estimate);
      
      % Now generate a mask for where there is data for all images,
      % so we can apply this when writing the mean
      msk      = spm_write_sn(V,params,pp_defaults.normalise.write,'mask');
      
      disp(sprintf('Subject #%d, Writing mean...',sub))
      diary off; diary on; %Flush the current diary logfile
      spm_write_sn(Vm,params,pp_defaults.normalise.write,msk);		% write nmean
      disp(sprintf(' Done.'))
      
      disp(sprintf('Subject #%d, Writing %s ...',sub,[iBT.what.SpmVersion '_normalise.ps']))
      diary off; diary on; %Flush the current diary logfile
      switch (iBT.what.SpmVersion),
        case 'SPM2'
          original_printstr = defaults.printstr;
          defaults.printstr = [spm_figure('DefPrintCmd'), 'spm2_normalise.ps' ];
          spm_print;
          defaults.printstr = original_printstr;
        otherwise
          spm_print([iBT.what.SpmVersion '_normalise.ps']);
      end
      disp(sprintf(' Done.'))
      diary off; diary on; %Flush the current diary logfile
      
      if do_WriteNorm == 1
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        % Write normalised
        
        disp(sprintf('Subject #%d, Calculating %d normalised images, smoothing, writing...',sub,size(V,1)))
        diary off; diary on; %Flush the current diary logfile
        for i=1:size(V,1)
          fprintf('\r    (currently applying warp to image %d)',i)
          VO = spm_write_sn(V(i),params,pp_defaults.normalise.write,msk);
          [pth,nam,ext] = fileparts(V(i).fname);
          fname         = fullfile(pth,['sw' nam ext]);
          spm_smooth(VO,fname,n_smth);
          clear VO; %Don't need it now that we've saved it
        end;
        fprintf('\rOK.                                             \n');
        diary off; diary on; %Flush the current diary logfile
        
        fclose( fopen( fullfile(iBT.who.sub{sub}.sess{ses}.loc,['_n_smoothed_to',sprintf('_%.2f',n_smth),'mm']) , 'W' ) );
      end; %if do_WriteNorm == 1
      
      cd(cwd);
      
    end % if do_SpatialNorm (without intra-subject inter-session norm)
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % Spatial normalisation
    % (including intra-subject inter-session norm)
    % Only done when we're processing the subject's final session
    % as need the pre-processed data from all other sessions
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    
    if ( (do_SpatialNorm == 1) && (do_norm_intra > 0 ) && ((ses == nses) || (iBT.what.pre.sessions_together == 1)) ),
      
      %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
      % Spatial normalisation
      % Read in the mean image from each session
      %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
      
      do_swi = 0; % Initialise
      
      disp(sprintf('*** Commencing spatial normalisation procedure for subject #%d...',sub))
      
      for ses2=1:nses,
        
        if iBT.what.pre.Realign == 0 % no realignment so use first image
          disp(sprintf('Subject #%d, Loading first image information from session %d...', ...
            sub, ses2))
          diary off; diary on; %Flush the current diary logfile
          switch (iBT.what.SpmVersion),
            case 'SPM2',
              meanf   = spm_get('Files',iBT.who.sub{sub}.sess{ses2}.loc, slice_wild.sub{sub}.sess{ses2});
            otherwise,
              meanf   = spm_select('FPList',iBT.who.sub{sub}.sess{ses2}.loc, slice_wild.sub{sub}.sess{ses2});
          end %case
          meanf   = spm_vol(meanf);
          meanf   = meanf(1).fname;
        else
          disp(sprintf('Subject #%d, Loading mean image information from session %d...', ...
            sub, ses2))
          diary off; diary on; %Flush the current diary logfile
          switch (iBT.what.SpmVersion),
            case 'SPM2',
              meanf   = spm_get('Files',iBT.who.sub{sub}.sess{ses2}.loc, mean_wild);
            otherwise,
              meanf   = spm_select('FPList',iBT.who.sub{sub}.sess{ses2}.loc, mean_wild);
          end %case
        end
        
        try
          swi = iBT.who.sub{sub}.sess{ses2}.swi;
        catch
          iBT.who.sub{sub}.sess{ses2}.swi = 0; % Maintain backwards compatibility
        end % try
        
        % Check for source weighted image and add it if it exists
        % Happy to accept number 0, string '0', empty string, or even
        % a string containing two single quotes as indicating none
        % These possibilities make it easier to specify nothing in a
        % studylist file.
        if (strcmp('',iBT.who.sub{sub}.sess{ses2}.swi) == 1)...
            || (strcmp('''''',iBT.who.sub{sub}.sess{ses2}.swi) == 1)...
            || (strcmp('0',iBT.who.sub{sub}.sess{ses2}.swi) == 1)...
            || max(iBT.who.sub{sub}.sess{ses2}.swi == 0)
          % Note that max in the line above works becasue '0' ~= 0 so any string (array of chars)
          % will always return an array of zeroes even if a char '0' is present in the string.
          Vswi{ses2} = ''; % No SWI
          iBT.who.sub{sub}.sess{ses2}.swi = 0; % Simplify logic for subsequent tests of this.
        else
          
          % Maintain compatibility if ImageSWIdrawnUpon set.
          try, tmp = iBT.who.sub{sub}.sess{ses2}.ImageSWIdrawnUpon; catch, iBT.who.sub{sub}.sess{ses2}.ImageSWIdrawnUpon = 0; end % Default if not specified
          if ~ ( (strcmp('',iBT.who.sub{sub}.sess{ses2}.ImageSWIdrawnUpon) == 1)...
              || (strcmp('''''',iBT.who.sub{sub}.sess{ses2}.ImageSWIdrawnUpon) == 1)...
              || (strcmp('0',iBT.who.sub{sub}.sess{ses2}.ImageSWIdrawnUpon) == 1)...
              || max(iBT.who.sub{sub}.sess{ses2}.ImageSWIdrawnUpon == 0) )
            
            display(sprintf('Co-registering mask %s...',iBT.who.sub{sub}.sess{ses2}.swi));
            clear matlabbatch;
            load iBT_spm8_coregister_job; % This is also compatible with SPM12
            matlabbatch{1}.spm.spatial.coreg.estwrite.ref{1} = meanf
            matlabbatch{1}.spm.spatial.coreg.estwrite.source{1} =  iBT.who.sub{sub}.sess{ses2}.ImageSWIdrawnUpon
            matlabbatch{1}.spm.spatial.coreg.estwrite.other{1} = iBT.who.sub{sub}.sess{ses2}.swi
            out = spm_jobman('run',matlabbatch);
            [mpth,mnam,mext] = fileparts(meanf);
            [spth,snam,sext] = fileparts(iBT.who.sub{sub}.sess{ses2}.swi);
            new_mask_name = fullfile(spth,['r' snam sext]);  %add r to iBT.who.sub{sub}.sess{ses2}.swi
            [SUCCESS,MESSAGE,MESSAGEID] = movefile(new_mask_name, mpth);
            if SUCCESS == 0
              throw(MException('iBT_preprocess:co-registering mask',[MESSAGE ':' MESSAGEID] ));
            end
            if strcmpi(sext,'.img') % then we need to move the .hdr also
              new_mask_name_hdr = fullfile(spth,['r' snam '.hdr']);  %add r to iBT.who.sub{sub}.sess{ses2}.swi
              [SUCCESS,MESSAGE,MESSAGEID] = movefile(new_mask_name_hdr, mpth);
            end
            if SUCCESS == 0
              throw(MException('iBT_preprocess:co-registering mask',[MESSAGE ':' MESSAGEID] ));
            end
            new_mask_name = fullfile(mpth,['r' snam sext]);  %We've moved it to mpth
            Vswi{ses2} = spm_vol(new_mask_name);
          else
            Vswi{ses2} = spm_vol(iBT.who.sub{sub}.sess{ses2}.swi);
          end
          do_swi = 1;
        end
        Vm(ses2) = spm_vol(meanf);
      end % for ses2=1:nses
      
      
      %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
      % Spatial normalisation
      %
      % perform an iterative process of registering the mean from
      % each session into a common space
      %
      % The code here for calculating the mean of a set of images
      % is taken from spm_mean_ui.m
      %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
      defaultPrefix = pp_defaults.normalise.write.prefix;
      if nses == 1
        plural = '';
      else
        plural = 's';
      end
      
      % First create our initial template (using bias corrected images if requested)
      
      if (do_BiasCorrect ~= 0)
        % Use SPM8's segment tool to create bias corrected versions of the mean image
        % in each session (in own-space) - these will then be used for normalisation.
        % Unfortunately SPM's segment tool requires the images to be roughly in the
        % right space, so we do an initial affine registration to achieve this and will
        % combine this with other transforms later.
        
        oVm = Vm;%  Keep the originals for use later
        if (do_BiasCorrect == 1) % Hasn't already been done
          display(sprintf('Re-orienting mean image%s approximately to %s...',plural,VG));
          diary off; diary on; %Flush the current diary logfile
          for ses2 = 1:nses,
            pp_defaults.normalise.write.prefix = '_w00_';
            [pth,nam,ext] = fileparts(Vm(ses2).fname);
            matname_affine = sprintf('%s_affine_to_template_sn.mat', spm_str_manip(Vm(ses2).fname, 'sd'));
            temp_estimate = intra_session.normalise.estimate;
            temp_estimate.nits = 0; % 0 because this 1st step always affine.
            norm_name = fullfile(pth, [pp_defaults.normalise.write.prefix nam ext]); % This matches the name that spm_write_sn will write
            display(sprintf('   session %d: re-orienting %s (save result as %s)',ses2,Vm(ses2).fname,norm_name));
            if iBT.what.pre.norm.intra.do_swi ~= 0
              if iBT.who.sub{sub}.sess{ses2}.swi ~= 0, disp(sprintf('  (using Source Weighted Image %s)',Vswi{ses2}.fname)); end
              params = spm_normalise(VG, Vm(ses2), matname_affine, '', Vswi{ses2},temp_estimate);
            else
              params = spm_normalise(VG, Vm(ses2), matname_affine, '', '',temp_estimate);
            end
            spm_write_sn(Vm(ses2),params,pp_defaults.normalise.write)
            Vm(ses2) = spm_vol(norm_name); %Use the affine normalised version from now on (originals still in oVm)
            
            if iBT.who.sub{sub}.sess{ses2}.swi ~= 0; % need to create a lesion mask in this space for template weighted image for the next iteration
              [pth_swi,nam_swi,ext_swi] = fileparts(Vswi{ses2}.fname);
              norm_name_swi = fullfile(pth_swi, [pp_defaults.normalise.write.prefix nam_swi ext_swi]); % matches the name spm_write_sn will write
              swi_defaults.normalise.write = pp_defaults.normalise.write;
              if swi_defaults.normalise.write.interp > 1, swi_defaults.normalise.write.interp = 1;end % Don't use any more than tri-linear interpolation for the mask
              spm_write_sn(Vswi{ses2},params,swi_defaults.normalise.write)
              Vtwi(ses2) = spm_vol(norm_name_swi);
              Vtwi_vol = spm_read_vols(Vtwi(ses2));
              Vtwi_vol = max(Vtwi_vol, 0); % Set any negative values to zero
              Vtwi(ses2) = spm_write_vol(Vtwi(ses2),Vtwi_vol); % Re-save this non-negative version of the mask.
              % When using bias-corrected images, we start iterating from the affine-registered verison, not the original image,
              % so need to also change the Vswi accordingly as follows:
              Vswi{ses2} = Vtwi(ses2);
            end; % if iBT.who.sub{sub}.sess{ses2}.swi ~= 0
            
          end %for ses2 = 1:nses
          pp_defaults.normalise.write.prefix = defaultPrefix;
          
          display(sprintf('Bias correcting image%s...',plural));
          diary off; diary on; %Flush the current diary logfile
          clear matlabbatch;
          load iBT_spm8_segment_job; %This is also OK for SPM12.
          switch(iBT.what.SpmVersion), %% Set location of the anatomical priors
            case {'SPM8'}
              templates_folder=fullfile(spm('Dir'),'toolbox','OldNorm');
              matlabbatch{1}.spm.spatial.preproc.opts.tpm = {
                fullfile(spm('Dir'),'tpm','grey.nii')
                fullfile(spm('Dir'),'tpm','white.nii')
                fullfile(spm('Dir'),'tpm','csf.nii')
                };
            otherwise %  SPM12 or later
              matlabbatch{1}.spm.spatial.preproc.opts.tpm = {
                fullfile(spm('Dir'),'toolbox','OldSeg','grey.nii')
                fullfile(spm('Dir'),'toolbox','OldSeg','white.nii')
                fullfile(spm('Dir'),'toolbox','OldSeg','csf.nii')
                };
          end % switch(iBT.what.SpmVersion)
          
          for ses2 = 1:nses,
            matlabbatch{1}.spm.spatial.preproc.data{ses2} = [Vm(ses2).fname ',1'];
          end
          %       out = spm_jobman('initcfg');
          out = spm_jobman('run',matlabbatch);
        end %if (do_BiasCorrect == 1) % Hasn't already been done
        
        for ses2 = 1:nses,
          [pth,nam,ext] = fileparts(oVm(ses2).fname);
          bias_corrected_norm_name = fullfile(pth, ['m' '_w00_' nam ext]);
          Vm(ses2) = spm_vol(bias_corrected_norm_name); %We'll use the bias corrected images from now on.
          if (do_BiasCorrect < 0) % Bias correction was done earlier so we won't yet have a loaded the TWI's.
            if iBT.who.sub{sub}.sess{ses2}.swi ~= 0
              [pth_swi,nam_swi,ext_swi] = fileparts(iBT.who.sub{sub}.sess{ses2}.swi);
              norm_name_swi = fullfile(pth_swi, [pp_defaults.normalise.write.prefix nam_swi ext_swi]); % matches the name spm_write_sn wrote earlier
              Vtwi(ses2) = spm_vol(norm_name_swi);
              Vswi{ses2} = Vtwi(ses2); % The SWI and TWI start off the same.
            end
          end
        end %for ses2 = 1:nses
        Vtarget = Vm(iBT.who.sub{sub}.MasterSession); % This will be the first template for multi-session normalisation
        
        if iBT.who.sub{sub}.sess{iBT.who.sub{sub}.MasterSession}.swi ~= 0
          Vtwi_m = Vtwi(iBT.who.sub{sub}.MasterSession); % The TWI for the master session
        else
          Vtwi_m = ''; % No master session TWI to begin with
        end
        
        display(sprintf('Initial template is %s',Vtarget.fname));
        diary off; diary on; %Flush the current diary logfile
        
      else % do_BiasCorrect must be 0
        
        oVm = Vm;%  Makes the code a bit simpler later as don't need to test for do_BiasCorrect to decide which to use. Wastes memory though as only really need the fname fields.
        Vtarget = Vm(iBT.who.sub{sub}.MasterSession);
        display(sprintf('Determining initial template (affine normalizing session %d mean image to %s)',iBT.who.sub{sub}.MasterSession,VG));
        diary off; diary on; %Flush the current diary logfile
        pp_defaults.normalise.write.prefix = '_w00_';
        [pth,nam,ext] = fileparts(Vtarget.fname);
        matname_affine = sprintf('%s_affine_to_template_sn.mat', spm_str_manip(Vm(iBT.who.sub{sub}.MasterSession).fname, 'sd')); % The spm_str_manip is used to remove the trailing suffix
        temp_estimate = intra_session.normalise.estimate;
        temp_estimate.nits = 0; % 0 because this 1st step always affine.
        if iBT.what.pre.norm.intra.do_swi ~= 0
          if iBT.who.sub{sub}.sess{iBT.who.sub{sub}.MasterSession}.swi ~= 0, disp(sprintf('  (using Source Weighted Image %s)',Vswi{iBT.who.sub{sub}.MasterSession}.fname)); end
          params = spm_normalise(VG, Vtarget, matname_affine, '', Vswi{iBT.who.sub{sub}.MasterSession},temp_estimate);
        else
          params = spm_normalise(VG, Vtarget, matname_affine, '', '',temp_estimate);
        end
        spm_write_sn(Vtarget,params,pp_defaults.normalise.write)
        norm_name = fullfile(pth, [pp_defaults.normalise.write.prefix nam ext]); % This is the name that spm_write_sn will have written
        Vtarget = spm_vol(norm_name);
        
        if iBT.who.sub{sub}.sess{iBT.who.sub{sub}.MasterSession}.swi ~= 0; % need to transform SWI into this space to use as a template weighted image for the next iteration
          [pth_swi,nam_swi,ext_swi] = fileparts( Vswi{iBT.who.sub{sub}.MasterSession}.fname);
          norm_name_swi = fullfile(pth_swi, [pp_defaults.normalise.write.prefix nam_swi ext_swi]); % matches the name spm_write_sn will write
          swi_defaults.normalise.write = pp_defaults.normalise.write;
          if swi_defaults.normalise.write.interp > 1, swi_defaults.normalise.write.interp = 1;end % Don't use any more than tri-linear interpolation for the mask
          spm_write_sn(Vswi{iBT.who.sub{sub}.MasterSession},params,swi_defaults.normalise.write)
          Vtwi_m(iBT.who.sub{sub}.MasterSession) = spm_vol(norm_name_swi);
          Vtwi_vol = spm_read_vols(Vtwi_m(iBT.who.sub{sub}.MasterSession));
          Vtwi_vol = max(Vtwi_vol, 0); % Set any negative values to zero
          Vtwi_m(iBT.who.sub{sub}.MasterSession) = spm_write_vol(Vtwi_m(iBT.who.sub{sub}.MasterSession),Vtwi_vol);% Re-save this non-negative version of the mask.
          % This is the TWI for the master session
        else
          Vtwi_m = ''; % No master session SWI to begin with
        end; % if iBT.who.sub{sub}.sess{iBT.who.sub{sub}.MasterSession}.swi ~= 0
        
      end %if (do_BiasCorrect ~= 0) else
      
      
      for coregistrations_i = 1:intra_niterations,
        
        sessions_used_for_new_Vtwi_m = 0; % Initialise flag
        for ses2 = 1:nses,
          display(sprintf('Determining next template (iteration %d, session=%d)',coregistrations_i,ses2));
          diary off; diary on; %Flush the current diary logfile
          matname = sprintf('%s_to_normtarget%02d_sn.mat', spm_str_manip(Vm(ses2).fname, 'sd'), coregistrations_i-1);
          
          if iBT.what.pre.norm.intra.do_swi ~= 0
            if iBT.who.sub{sub}.sess{iBT.who.sub{sub}.MasterSession}.swi ~= 0, disp(sprintf('  (using Template Weighted Image %s)',Vtwi_m.fname)); end;
            if iBT.who.sub{sub}.sess{ses2}.swi ~= 0, disp(sprintf('  (using Source Weighted Image %s)',Vswi{ses2}.fname)); end;
            params  = spm_normalise(Vtarget, Vm(ses2), matname, Vtwi_m, Vswi{ses2}, intra_session.normalise.estimate);
          else
            params  = spm_normalise(Vtarget, Vm(ses2), matname, '', '', intra_session.normalise.estimate);
          end
          
          pp_defaults.normalise.write.prefix = sprintf('_w%02d_',coregistrations_i);
          [pth,nam,ext] = fileparts(Vm(ses2).fname);
          norm_name = fullfile(pth, [pp_defaults.normalise.write.prefix nam ext]); % Name of file we are about to write;
          spm_write_sn(Vm(ses2),params,pp_defaults.normalise.write);		% write normalised image
          wVm(ses2) = spm_vol(norm_name); % Load the normalised image we just wrote
          
          if iBT.who.sub{sub}.sess{ses2}.swi ~= 0; % need to create a lesion mask in this space for template weighted image for the next iteration
            [pth_swi,nam_swi,ext_swi] = fileparts(Vswi{ses2}.fname);
            norm_name_swi = fullfile(pth_swi, [pp_defaults.normalise.write.prefix nam_swi ext_swi]); % matches the name spm_write_sn will write
            swi_defaults.normalise.write = pp_defaults.normalise.write;
            if swi_defaults.normalise.write.interp > 1, swi_defaults.normalise.write.interp = 1;end % Don't use any more than tri-linear interpolation for the mask
            spm_write_sn(Vswi{ses2},params,swi_defaults.normalise.write)
            Vtwi(ses2) = spm_vol(norm_name_swi);
            Vtwi_vol = spm_read_vols(Vtwi(ses2));
            Vtwi_vol = max(Vtwi_vol, 0); % Set any negative values to zero
            Vtwi(ses2) = spm_write_vol(Vtwi(ses2),Vtwi_vol); % Re-save this non-negative version of the mask.
            % The next template weighted image is a composite master TWI based on the
            % minimum of the values across all the sesson-specific TWI's:
            if sessions_used_for_new_Vtwi_m == 0
              new_Vtwi_m = Vtwi(ses2);
              new_Vtwi_m_vol = Vtwi_vol;
              sessions_used_for_new_Vtwi_m = 1;
            else
              if ses2 == iBT.who.sub{sub}.MasterSession, new_Vtwi_m = Vtwi(ses2); end % We prefer details from the master if present
              new_Vtwi_m_vol = min(new_Vtwi_m_vol,Vtwi_vol);
              sessions_used_for_new_Vtwi_m = sessions_used_for_new_Vtwi_m + 1;
            end
          end; % if iBT.who.sub{sub}.sess{ses2}.swi ~= 0
          
        end; % for ses2 = 1:nses
        
        pp_defaults.normalise.write.prefix = defaultPrefix;
        
        [pth,nam,ext] = fileparts(wVm(iBT.who.sub{sub}.MasterSession).fname);
        fname = fullfile(pth,sprintf(['normalisation_target_%02d' iBT.what.ext.SPM], coregistrations_i));
        Vtarget  = struct('fname',    fname,...
          'dim',      wVm(iBT.who.sub{sub}.MasterSession).dim(1:3),...
          'dt',       [8, spm_platform('bigend')],...
          'mat',      wVm(iBT.who.sub{sub}.MasterSession).mat,...
          'pinfo',    [1.0,0,0]',...
          'descrip',  'spm - mean image');
        
        %-Adjust scalefactors by 1/n to effect mean by summing
        for i=1:prod(size(wVm))
          wVm(i).pinfo(1:2,:) = wVm(i).pinfo(1:2,:)/nses;
        end;
        
        Vtarget            = spm_create_vol(Vtarget);
        Vtarget.pinfo(1,1) = spm_add(wVm, Vtarget);
        Vtarget            = spm_create_vol(Vtarget);
        display(sprintf('Template created for iteration %d (%s)',coregistrations_i, Vtarget.fname));
        diary off; diary on; %Flush the current diary logfile
        
        if sessions_used_for_new_Vtwi_m > 0
          new_Vtwi_m.fname = fullfile(pth,sprintf(['twi_for_normalisation_target_%02d' iBT.what.ext.SPM], coregistrations_i));
          Vtwi_m = spm_write_vol(new_Vtwi_m,new_Vtwi_m_vol); % Save the new mask.
          display(sprintf('Template weighted image created from %i source weighted image(s) for iteration %d (%s)',sessions_used_for_new_Vtwi_m,coregistrations_i, Vtwi_m.fname));
        end
        
      end % for coregistrations_i
      
      %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
      % Spatial normalisation
      % register the common space for that patient to a standard
      % template
      %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
      
      if iBT.what.pre.norm.intra.strategy == 2, % (Now the default strategy)
        
        % Here we first inverse-affine transform the target back into approximately the space of the master session
        % I (dfa) heuristically determined this yielded a superior registration for GE images at our site
        % despite this being a counter-intuative strategy. I suspect it may be due to field of view and
        % head position differences that the Spatial Normalisation routine only properly takes
        % into account when it is dealing with the images in its original bounding box.
        clear job;
        matname_affine = sprintf('%s_affine_to_template_sn.mat', spm_str_manip(oVm(iBT.who.sub{sub}.MasterSession).fname, 'sd')); % The spm_str_manip is used to remove the trailing suffix
        job.comp{1}.inv.comp{1}.sn2def.matname{1} = matname_affine;
        job.comp{1}.inv.comp{1}.sn2def.vox = [NaN NaN NaN]; % Use size as originally used to estimate template
        job.comp{1}.inv.comp{1}.sn2def.bb = [[NaN NaN NaN];[NaN NaN NaN]]; % Use bb as originally used to estimate template
        job.comp{1}.inv.space = {oVm(iBT.who.sub{sub}.MasterSession).fname};
        disp(sprintf('Creating inverse of affine transform %s', job.comp{1}.inv.comp{1}.sn2def.matname{1}))
        disp(sprintf('  image to base inverse on = %s', job.comp{1}.inv.space{1}))
        %fname of original from master session{'/brain6/data/dfa/test_pp/manual_affineToEPI_bc_WarpToEPI_masked/00001R.img,1'}
        
        diary off; diary on; %Flush the current diary logfile
        switch(iBT.what.SpmVersion),
          case {'SPM5','SPM8'}
            if do_swi ~= 0,
              job.fnames = {Vtarget.fname, Vtwi_m.fname}; else job.fnames = {Vtarget.fname}; end
            job.ofname = 'affine_inverse';
            job.savedir.savesrc=1
            job.interp = iBT.what.pre.norm.write.interp;
            out = spm_defs(job);
          otherwise %  SPM12 or later
            job.out{1}.savedef.ofname = 'affine_inverse';
            [PATHSTR,NAME,EXT] = fileparts(Vtarget.fname);
            job.out{1}.savedef.savedir.saveusr{1} = PATHSTR; % Where to save y_affine_inverse.nii deformation file
            if do_swi ~= 0, job.out{2}.pull.fnames = {Vtarget.fname, Vtwi_m.fname}; else job.out{2}.pull.fnames = {Vtarget.fname}; end
            job.out{2}.pull.savedir.savesrc = 1; % Save output images in same folder as source images
            job.out{2}.pull.interp = iBT.what.pre.norm.write.interp;
            job.out{2}.pull.mask = 0;
            job.out{2}.pull.fwhm = [0 0 0];
            out = spm_deformations(job);
        end % switch(iBT.what.SpmVersion)
        
        [pth_Vtarget,nam_Vtarget,ext_Vtarget] = fileparts(Vtarget.fname);
        new_Vtarget_fname = fullfile(pth_Vtarget, [pp_defaults.normalise.write.prefix nam_Vtarget ext_Vtarget]); % matches the name the inverse_affine wrote
        Vtarget = spm_vol(new_Vtarget_fname); % We now proceed from here using this back-transformed target
        
        if do_swi ~= 0
          [pth_swi,nam_swi,ext_swi] = fileparts(Vtwi_m.fname);
          nonfixed_new_swi_fname = fullfile(pth_swi, [pp_defaults.normalise.write.prefix nam_swi ext_swi]); % matches the name the inverse_affine wrote
          % The nonfixed_new_swi_fname image has been observed to have negative values and/or values
          % greater than 1 in some places outside the brain. This casues the resultant normalisation to fail.
          % We fix things up here
          nonfixed_Vtwi_m = spm_vol(nonfixed_new_swi_fname);
          new_swi_fname = fullfile(pth_swi, [pp_defaults.normalise.write.prefix nam_swi '_fixed' ext_swi]);
          Vtwi_m  = struct('fname',    new_swi_fname,...
            'dim',      nonfixed_Vtwi_m(1).dim(1:3),...
            'dt',       [4, spm_platform('bigend')],...
            'mat',      nonfixed_Vtwi_m(1).mat,...
            'pinfo',    [1.0,0,0]',...
            'descrip',  'fixed-up back transformed SWI');
          Vtwi_m = spm_imcalc([nonfixed_Vtwi_m],Vtwi_m,'(((i1.*(i1 > 0))-1).*(i1<1))+1',{0,0,0})
          % This zeroes all negative values and sets all values greater than 1 to 1.
          % It does this by first setting values below 0 to 0 via i1.(i1 > 0)
          % then subtracts 1 and sets all values originally more than 1 to 0.
          % then adds 1 to get original values back.
          
          % We now proceed from here using this back-transformed mask
        end %if do_swi ~= 0
        
      end % if iBT.what.pre.norm.intra.strategy == 2
      
      matname = sprintf('%s_to_template_sn.mat', spm_str_manip(Vtarget.fname, 'sd'));
      
      if do_swi ~= 0
        params   = spm_normalise(VG, Vtarget, matname, '', Vtwi_m, to_template.normalise.estimate); % The final TWI is the SWI for this step.
        
        % We now create a normalised lesion mask in the final template space for quality control purposes:
        [pth_swi,nam_swi,ext_swi] = fileparts(Vtwi_m.fname);
        norm_name_swi = fullfile(pth_swi, [pp_defaults.normalise.write.prefix nam_swi ext_swi]); % matches the name spm_write_sn wrote
        swi_defaults.normalise.write = pp_defaults.normalise.write;
        if swi_defaults.normalise.write.interp > 1, swi_defaults.normalise.write.interp = 1;end % Don't use any more than tri-linear interpolation for the mask
        spm_write_sn(Vtwi_m,params,swi_defaults.normalise.write)
        wVtwi_m = spm_vol(norm_name_swi);
        wVtwi_vol = spm_read_vols(wVtwi_m);
        wVtwi_vol = max(wVtwi_vol, 0); % Set any negative values to zero
        wVtwi_m = spm_write_vol(wVtwi_m,wVtwi_vol); % Re-save this non-negative version of the mask.
      else
        params   = spm_normalise(VG, Vtarget, matname, '', '', to_template.normalise.estimate);
      end % if do_swi ~= 0
      
      %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
      % Spatial normalisation
      % We also smooth the files in this section.
      %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
      
      %
      %
      % combine the own-space transformation with the transform to
      % the standard template into a single transformation and apply
      % that to every frame for the current patient
      %
      % This section of the code calls routines in spm_defs.m
      % In the nomenclature of the SPM documentation
      % (see the Deformations Toolbox GUI), we have two
      % transformations x:A->B and y:B->C where in our case
      % A is original image space, B is own-space template and
      % C is standard template space. Transform x (_iterN_sn.mat)
      % takes us from A to B, and y (_template_sn.mat) takes us
      % from B to C. The composition of transforms is denoted by
      % a circle, and in this case would be yox:A->C. (i.e. in this
      %  notation the right-most transform is applied first).
      % In the Deformations toolbox GUI the left-to-right order it is
      % claimed to be top to bottom (i.e. the left-most or top
      % transform which is the first selected should be applied last),
      % HOWEVER tests I undertook on 11/11/2010 revealed that this
      % is NOT the order in which they are applied.
      % I (dfa) have previously
      % verified that the top transform in the Deformations
      % toolbox GUI is passed as job.comp{1} to spm_defs.m, and
      % the next is job.comp{2} etc.
      %
      
      clear job;
      
      final_jobnum = 2; % Initialise
      if do_BiasCorrect ~= 0, final_jobnum = final_jobnum + 1; end
      if iBT.what.pre.norm.intra.strategy == 2, final_jobnum = final_jobnum + 1; end
      
      % The final transform (same for all sessions):
      job.comp{final_jobnum}.sn2def.bb = [[NaN NaN NaN];[NaN NaN NaN]]; % Use bb as originally used to estimate template
      job.comp{final_jobnum}.sn2def.vox = [NaN NaN NaN]; % Use size as originally used to estimate template
      job.comp{final_jobnum}.sn2def.matname{1} = matname; % Transform to standard template
      
      if iBT.what.pre.norm.intra.strategy == 2,
        % The penultimate transform (same for all sessions)
        [pth,nam,ext] = fileparts(oVm(iBT.who.sub{sub}.MasterSession).fname);
        job.comp{final_jobnum-1}.def{1} = fullfile( pth, 'y_affine_inverse.nii'); % Inverse-affine transform final target to approx. own-space of master session.
      end
      
      for ses2 = 1:nses,
        cd (iBT.who.sub{sub}.sess{ses2}.loc);
        jobnum = 1;
        if (do_BiasCorrect ~= 0), % Include the initial affine transform to the template as from
          % then on we used the bias corrected images in this space.
          job.comp{jobnum}.sn2def.bb = [[NaN NaN NaN];[NaN NaN NaN]]; % Use bb as originally used to estimate template
          job.comp{jobnum}.sn2def.vox = [NaN NaN NaN]; % Use size as originally used to estimate template
          job.comp{jobnum}.sn2def.matname{1} = sprintf('%s_affine_to_template_sn.mat', spm_str_manip(oVm(ses2).fname, 'sd'));
          jobnum = jobnum + 1
        end %if (do_BiasCorrect ~= 0)
        
        job.comp{jobnum}.sn2def.bb = [[NaN NaN NaN];[NaN NaN NaN]]; % Use bb as originally used to estimate template
        job.comp{jobnum}.sn2def.vox = [NaN NaN NaN]; % Use size as originally used to estimate template
        job.comp{jobnum}.sn2def.matname{1} = sprintf('%s_to_normtarget%02d_sn.mat', spm_str_manip(Vm(ses2).fname, 'sd'), ...
          coregistrations_i-1);
        
        if (do_BiasCorrect ~= 0), Vm(ses2) = oVm(ses2); end % Use the original mean image name from now on...
        
        if (do_WriteNorm == 1),
          [tmp,basename] = fileparts(Vm(ses2).fname);
          if (do_BiasCorrect == 1)
            job_ofname = sprintf('%s_to_template_biasc_to_normtarget%02d_to_template_sn.mat', basename, ...
              coregistrations_i-1);
          else
            job_ofname = sprintf('%s_to_normtarget%02d_to_template_sn.mat', basename, ...
              coregistrations_i-1);
          end
          job_fnames = spm_select('FPList',iBT.who.sub{sub}.sess{ses2}.loc, slice_wild.sub{sub}.sess{ses2});
        end %if (do_WriteNorm == 1)
        
        switch(iBT.what.SpmVersion),
          case {'SPM5','SPM8'}
            job.fnames=job_fnames;
            job.ofname=job_ofname;
            job.interp = iBT.what.pre.norm.write.interp;
            job.savedir.saveusr{1} = pwd; % Where to save deformation file and output images
            out = spm_defs(job);
          otherwise %  SPM12 or later
            job.out{1}.savedef.ofname = job_ofname;
            job.out{1}.savedef.savedir.saveusr{1} = pwd; %Where to save deformation file
            job.out{2}.pull.fnames = cellstr(job_fnames);
            job.out{2}.pull.savedir.saveusr{1} = pwd; % Where to save output images
            job.out{2}.pull.interp = iBT.what.pre.norm.write.interp;
            job.out{2}.pull.mask = 0;
            job.out{2}.pull.fwhm = [0 0 0];
            out = spm_deformations(job);
        end % switch(iBT.what.SpmVersion)
        
        P = spm_select('FPList',iBT.who.sub{sub}.sess{ses2}.loc, wslice_wild.sub{sub}.sess{ses2});
        
        V = spm_vol(P);
        
        VmAll = V;
        %%%%%% Create mean of non-smoothed fully normalised images (code duplicated from above)
        [pth,nam,ext] = fileparts(VmAll(1).fname);
        fname         = fullfile(pth,['wmean_all' iBT.what.ext.SPM]);
        Vtarget  = struct('fname',    fname,...
          'dim',      VmAll(1).dim(1:3),...
          'dt',       [4, spm_platform('bigend')],...
          'mat',      VmAll(1).mat,...
          'pinfo',    [1.0,0,0]',...
          'descrip',  'spm - mean image');
        
        %-Adjust scalefactors by 1/n to effect mean by summing
        for i=1:prod(size(VmAll))
          VmAll(i).pinfo(1:2,:) = VmAll(i).pinfo(1:2,:)/nses;
        end;
        
        Vtarget            = spm_create_vol(Vtarget);
        Vtarget.pinfo(1,1) = spm_add(VmAll, Vtarget);
        Vtarget            = spm_create_vol(Vtarget);
        %%%%%% End of create mean of non-smoothed fully normalised images (code duplicated from above)
        
        %%%%%% Smooth all fully normalised images
        
        for i = 1:length(V),
          [pth, nam, ext] = fileparts(V(i).fname);
          fname = fullfile(pth, ['s' nam ext]);
          spm_smooth(V(i), fname, n_smth);
        end
        fclose( fopen( fullfile(iBT.who.sub{sub}.sess{ses2}.loc,['_n_smoothed_to',sprintf('_%.2f',n_smth),'mm']) , 'W' ) );
        %      end % what was this end for?
      end% for ses2 = 1:nses,
      
      disp(sprintf('*** Completed spatial normalisation procedure for subject #%d...',sub))
      diary off; diary on; %Flush the current diary logfile
      
    end  %if do_SpatialNorm (including intra-subject inter-session norm)
    
    % Unless we don't yet have the required data for within-brain masking, then proceed...
    if ~( (do_SpatialNorm == 1) && (do_norm_intra > 0 ) && ((ses < nses) && (iBT.what.pre.sessions_together == 0)) ),
      
      %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
      % Create within-brain mask of mean normalised image
      %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
      if (( iBT.what.pre.Create_n_iBrainWithinBrainMask > 0 ) || ( iBT.what.pre.Create_n_BETwithinBrainMask == 1 )) && ( iBT.what.pre.norm.do ~= 0 )
        if do_norm_intra == 0
          if (iBT.what.pre.iBrainRealignTarget == 0) | (iBT.what.pre.Realign == 0) | ( (iBT.what.pre.Realign < 0) & (iBT.what.pre.fmriprep.done == 1) ),
            normalised_mean_wild = ['wmean' iBT.what.ext.SPM];
          else
            normalised_mean_wild = ['wmean' iBT.what.pre.wild.ref];
          end
        else
          normalised_mean_wild = ['wmean_all' iBT.what.ext.SPM];
        end; % if do_norm_intra
        
        if ( (do_SpatialNorm == 1) && (do_norm_intra > 0 ) ) || (iBT.what.pre.sessions_together == 1)
          % We only get here after all sessions have been done, so need to loop again
          ses_start2 = 1;
          ses_end2 = nses;
        else % Doing each session as we go.
          ses_start2 = ses_start;
          ses_end2 = ses_end;
        end
        
        for ses=ses_start2:ses_end2 % Just does one session unless iBT.what.pre.sessions_together or do_norm_intra > 0; 
        %                             Using the existing loop variable ses within another loop here is intentional, see comments soon after start of the parent loop above.
          
          switch (iBT.what.SpmVersion),
            case 'SPM2',
              [normalised_mean,dir] = spm_get('Files', iBT.who.sub{sub}.sess{ses}.loc, normalised_mean_wild);
            otherwise
              [dir,normalised_mean,ext] = fileparts(spm_select('FPList', iBT.who.sub{sub}.sess{ses}.loc, spm_wildconvert(normalised_mean_wild)));
              normalised_mean = strcat(normalised_mean,ext);
          end
          
          if size(normalised_mean,1) > 0
            normalised_mean = normalised_mean(1,:); % in case more than one - we'll take only the first
            if  strcmp(iBT.what.SpmVersion, 'SPM2')
              iBrain_normalised_mean_mask = fullfile(dir,['iBrainMask_' normalised_mean])
              BET_normalised_mean_mask = fullfile(dir,['BETmask_' normalised_mean])
            end
            cd(iBT.who.sub{sub}.sess{ses}.loc);
            if  ~strcmp(iBT.what.SpmVersion, 'SPM2')
              iBrain_normalised_mean_mask = fullfile(pwd,['iBrainMask_' normalised_mean])
              BET_normalised_mean_mask = fullfile(pwd,['BETmask_' normalised_mean])
            end
            if (iBT.what.pre.Create_n_iBrainWithinBrainMask == 1) && (iBT.who.sub{sub}.sess{ses}.iBrainWithinBrainMask.create ~= 0)
              job=['-clobber -wbt ' sprintf('%f',iBT.who.sub{sub}.sess{ses}.iBrainWithinBrainMask.create) ' -wbm "' normalised_mean '" "' iBrain_normalised_mean_mask '"'];
              iBrain( job, pwd, use_iBrainServer );
            end
            if iBT.what.pre.Create_n_BETwithinBrainMask == 1
              if ispc()
                disp('BET is not available on Windows PC platforms. Please use a GNU/Linux platform if you need BET.')
              else
                command=[iBT_BET  ' "' normalised_mean '" "' BET_normalised_mean_mask '" ' iBT.who.sub{sub}.sess{ses}.BETwithinBrainMask.create]
                [status,result] = system(command);
                disp(result)
                if status ~=0
                  error(['Error: The following script returned an error - you may need to customise it for your site: ' iBT_BET ]);
                end
                % The command above generates two images - a masked brain with filename BET_normalised_mean_mask,
                % and the mask itself with a '_mask' appended to name. We only want the latter, so now delete the
                % former so wildcards in later scripts only get the _mask version:
                delete(BET_normalised_mean_mask);
                delete(strrep(BET_normalised_mean_mask,'.img','.hdr'));
              end %if ispc else
            end
            cd(cwd);
          else % if size(normalised_mean,1) > 0
            disp(sprintf('Unable to find mean normalised image matching %s',normalised_mean_wild))
            % Only write a warning if this is unexpected.
            if ~ ( (iBT.what.pre.sessions_together == 1) && (ses ~= iBT.who.sub{sub}.MasterSession) && (iBT.what.pre.norm.intra.do == 0) )
              disp(sprintf('Warning: unable to create normalised mask image for subject %d, session %d.',sub, ses))
            else
              disp(sprintf('Within-brain mask applicable to subject %d, session %d is that of master session (session %d) ',sub, ses, iBT.who.sub{sub}.MasterSession))
            end
          end % if file found i.e. size(normalised_mean_mean,1) > 0
          
          if ( iBT.what.pre.Create_n_iBrainWithinBrainMask == 2 ) && (iBT.who.sub{sub}.sess{ses}.iBrainWithinBrainMask.create ~= 0),
            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            %%%%%%%%%%%%%%%%%%%%%%% iBrain Mask Ceanup %%%%%%%%%%%%%%%%%%%
            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            
            !\rm ./__Script2_CSIRO__
            script_filename = '__Script2_CSIRO__';
            script_file2 =  fopen( script_filename , 'w' );
            fprintf(script_file2,'#!/bin/sh\n');
            fprintf(script_file2,'export LD_LIBRARY_PATH=/usr/local/milxview/lib \n');
            fprintf(script_file2,'export PATH="/usr/local/milxview/bin:$PATH" \n');
            
            tmp = regexprep(iBrain_normalised_mean_mask, iBT.what.ext.iB, ['_backup' iBT.what.ext.iB]);
            fprintf(script_file2,'mv %s %s \n',iBrain_normalised_mean_mask, tmp);
            if strcmpi(iBT.what.ext.iB,'.img') % then we need to move the .hdr also
              iBrain_normalised_mean_mask_hdr = regexprep(iBrain_normalised_mean_mask, '.img', '.hdr');
              tmp = regexprep(iBrain_normalised_mean_mask_hdr, '.hdr', '_backup.hdr');
              fprintf(script_file2,'mv %s %s \n',iBrain_normalised_mean_mask_hdr, tmp);
            end
            
            fprintf(script_file2,'echo CSIRO WithinBrainMask Correction %s %s \n ',tmp, iBrain_normalised_mean_mask);
            fprintf(script_file2,'/usr/local/milxview/bin/milxBinaryMaskCorrect %s %s \n ', tmp, iBrain_normalised_mean_mask);
            fclose(script_file2 );
            !chmod 744 __Script2_CSIRO__
            !./__Script2_CSIRO__
            
            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            %%%%%%%%%%%%%%%%%%%%%%% End iBrain Mask Ceanup %%%%%%%%%%%%%%%
            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
          end % if ( iBT.what.pre.Create_n_iBrainWithinBrainMask == 2 )
          
        end % for ses=ses_start2:ses_end2
        
      end % if Create_n_iBrainWithinBrainMask
      
      %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
      % Clean up
      %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
      
      diary off; diary on; %Flush the current diary logfile
      if do_Cleanup == 1
        % remove no longer needed intermediate images (only need keep the smoothed images for most purposes)
        disp(sprintf('Cleaning up files for subject #%d',sub))
        diary off; diary on; %Flush the current diary logfile
        try, ses_start2; catch ses_start2 = ses_start; end
        try, ses_end2; catch ses_end2 = ses_end; end
        
        for ses=ses_start2:ses_end2 % Using the existing loop variable ses within another loop here is intentional, see comments soon after start of the parent loop above.
          % Unfortunately Matlab only supports '*' wildcards, not '?', in its delete function
          % so we can't just call delete on the wildcard. We therefore use the SPM machinery.
          delete_filespec{1} = iBT.who.sub{sub}.sess{ses}.wild.slice;
          delete_filespec{2} = iBT.who.sub{sub}.sess{ses}.wild.realign;
          for i=1:numel(delete_filespec)
            switch (iBT.what.SpmVersion),
              case 'SPM2',
                P = spm_get('Files', iBT.who.sub{sub}.sess{ses}.loc, delete_filespec{i});
              otherwise
                P = spm_select('FPList', iBT.who.sub{sub}.sess{ses}.loc, spm_wildconvert(delete_filespec{i}));
            end % switch
            P = cellstr(P); P = P(:);
            n = numel(P);
            if n==0 || (n==1 && isempty(P{1})),
              disp(['    no images matching "' delete_filespec{i} '" in ' iBT.who.sub{sub}.sess{ses}.loc]);
            else
              disp(sprintf('    deleting %i images matching "%s" in %s',n, delete_filespec{i}, iBT.who.sub{sub}.sess{ses}.loc));
              feval(@spm_unlink,P{:});
            end% if n==0 else
          end; % for filespec_no
        end %for ses=ses_start2:ses_end2
        
      end % if do_Cleanup
      
      if (iBT.what.pre.sessions_together == 1) && (ses_end > ses_start)
        disp(sprintf('*** Finished preprocessing subject #%d (all sessions %d - %d).',sub,ses_start,ses_end))
      else
        disp(sprintf('*** Finished preprocessing subject #%d session #%d.',sub,ses))
      end
      diary off; diary on; %Flush the current diary logfile
      cd(cwd);
      
    end % if ~( (do_SpatialNorm == 1) && (do_norm_intra > 0 ) && ((ses < nses) && (iBT.what.pre.sessions_together == 0)) )
    
  end % for ses=1:loop_ses %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  
end % For subject

disp(sprintf('****** Preprocessing complete.'))
diary off


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

function new_str = quoted_wildconvert(orig_str)
%
% Sub function to convert normal wildcards to quoted format that can be passed to a shell command
%
% Modification History
% ---------------------
%
% 2010-12-14: (dfa) Function creation

new_str =  strcat('"', orig_str, '"');
