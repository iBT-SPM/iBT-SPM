function iBT_display(iBT)
% Perform display of activation maps.
% This script is designed to be called by iBT_analysis
% FORMAT iBT_display(iBT)
%  where 'iBT.what' is a structure indicating what is to be done,
%    and 'iBT.who' is a structure indicating subject specific details.
%___________________________________________________________________________
% Copyright 2006-2018,2020,2024,2025 The Florey Institute 
%                          of Neuroscience and Mental Health
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
% 2025-12-30: (dfa) Support iBT.what.disp.output_subfolder
% 2025-01-16: (dfa) Correct mean filename for fmriprep 'norm' analysis
% 2024-09-07: (dfa) Support iBT.what.disp.show.procVoxelSize
% 2024-08-19: (dfa) If the deprecated iBT.what.disp.revflip is defined, then 
%    display the requested sagittal slices from iBT.what.disp.orient.slices 
%    in reverse order, to result in a display order for legacy runinfos  
%    simlar to that experienced with iBT version 3.x (the change in iBT 4.x
%    then being just the displayed sign of the sagittal slice labels
%    (now consistent with MNI), without a change to slice order.
%    i.e. Now no need to change iBT.what.disp.orient.slices in legacy
%    runinfos to preserve previous slice order behaviour.
% 2024-08-04: (dfa) Implemented iBT.what.disp.orient.AnteriorOnRight
%                   (default 'no') and fixed a sagittal slice clipping
%                   problem introduced with changes on 2024-06-22.
% 2024-08-02: (bt) Prevent SPM welcome screen obscuring output when headless
% 2024-06-22: (dfa) 
%			Bugfix: "Left|Right is top" labels were incorrect on sagittal
%             display when step sign of iBT.what.disp.orient.slices
%             was changed to a positive value (label was correct
%             for default negative step sizes). Users of prior
%             iBT versions who generated sagittal montages with a
%             positive step size are urged to regenerate these.
%			  The text is now "Left|Right is 1st top" so can tell
%			  new montages from old.
%			Bugfix: The signs of the Sagittal display x-coordinate were
%             not following MNI convention if iBT.what.disp.revflip was 'yes'.
%             Now they do. A consequence of this change is that the displayed
%             slice ordering will differ. One can retain previous
%             ordering by altering direction in iBT.what.disp.orient.slices.
%			Replaced iBT.what.disp.revflip with more flexible:
%				iBT.what.disp.orient.LeftOnRight (implemented) &
%				iBT.what.disp.orient.AnteriorOnRight (coming soon)
%			    (for old runinfo files with revflip, this is auto-translated
%				 to these new variables for backwards compatibility)
%			Added iBT.what.disp.orient.showText so that users can
%			  disable orientation text if their input images 
%			  do not conform to Talairach / MNI convention, where
%			  positive is subject right, anterior or superior.
%			  (if users need labels for non-conformant orientations
% 			   then please let us know and we might implement override options)
%           Added warning if iBT.what.disp.revOrder is defined,
%             as it was long deprecated and ignored (slice ordering
%             is governed by iBT.what.disp.orient.slices)
%			Improved handling of iBT.what.disp.output_parent_folder
%			  so top level common elements of source and destination path
%			  are not duplicated within the destination.
%			Now use case-insensitive matching for most string comparisons
% 2020-12-02: (dfa) Do not show smooth if not present in _header_info.txt
% 2020-11-09: (dfa) Support iBT.what.disp.structuralLabel
% 2020-11-07: (dfa) Support iBT.what.disp.threshold.structural.low & high
%                     to limit the allowed adaptive grey-range of structural
%                     image display.
%                   Support iBT.what.disp.threshold.structural.min & max
%                     to force fixed (non-adaptive) grey-range.
%                   Support optional n_templateLabel in header info file.
% 2020-11-06: (dfa) Support info file version 3.0, which may include n_regtype
%                      (the normalisation constraint).
%                   Halt with error if unfamiliar major info file version.
% 2020-07-07: (dfa) Indicate 'liT' in png filenames of an li T threshold
% 2020-06-08: (dfa) Support iBT.what.disp.show.LR
%             	    Tweaked header linespace to stop crowding
% 2020-05-21: (dfa) Support iBT.what.pre.fmriprep.done
% 2020-04-26: (dfa) Support iBT.what.disp.show.voxelSize
% 2020-04-19: (dfa) Bugfix: support case where realignment target made but
%                           realignment to it not performed.
% 2018-09-19: (dfa) Support iBT.what.ext.SPM, iBT.what.ext.FSL and
%			iBT.what.ext.iB (image filename extensions)
% 2016-05-04: (dfa) Log Matlab Version
% 2014-03-19: (dfa) Fixed iBT.what.disp.contrast = 'auto' not omitting -effects
% 2014-01-13: (dfa) Support iBT.what.disp.ROI for semi-transparent ROI overlay.
% 2014-01-10: (dfa) Support auto- and +effect options in iBT.what.disp.contrast
% 		    and alter logic order so this filter is also applied when
%		    a contrast subset is specified via iBT.what.disp.con
%		    Support iBT.what.disp.invert
% 2013-12-11: (dfa) Honour iBT.what.disp.copyStructural & copyHeaderInfo
% 2013-12-09: (dfa) Support iBT.what.disp.mc.rejLog
% 2013-11-15: (dfa) Support iBT.what.disp.together
% 		    Support iBT.what.disp.contrast='-effect'
% 2013-11-04: (dfa) Fixed bug in last update that crashed if no contrast
%                   of interest in SPM.mat.
% 2013-10-22: (dfa) Support iBT.what.disp.show.mrej
%		    Bugfix (bug introduced on 2013-10-01): If displaying
%		     FWEc or FDRc and theshold was Inf, correct zero-length 
%		     indicator file was written but then incorrectly 
% 		     generated an image using a display title for the current
%		     image but image data, threhsold and filename
%		     of the previously displayed contrast, overwriting it.
% 2013-10-18: (dfa) Support iBT.what.disp.OTHER.do and
% 		     iBT.who.sub{sub}.sess{ses}.disp.loc 
%		     ( replaces iBT.what.disp.RFX.do and
%		      iBT.who.sub{sub}.sess{ses}.RFX.loc )
%		    Allow missing _header_info.txt for other display types
% 2013-10-04: (dfa) New iBT.what.disp.fontsize options & iBT.what.disp.con
% 2013-10-03: (dfa) New iBT.what.disp.show.conName
% 2013-10-01: (dfa) Support generation of FWEc and FDRc thresholded images.
%		    Bugfix: use correct li threshold pertaining to relevant
%                   contrast (previously that from 1st contrast only)
% 2013-09-30: (dfa) Bugfix: li file not found if only display at li threshold
% 2013-09-24: (dfa) Support fallback to display onto the copy of structural
% 		    if it had previously been copied to the stats folder.
%		    Moved loop over orientations to make execution more
%		    efficient.
%		    Support wildcards in iBT.what.disp.contrast
%		    Support iBT.what.disp.contrastMissingOK
%		    Tweak text display so more likely to fit on the page
% 2013-08-09: (dfa) Support iBT.what.disp.copyHeaderInfo
% 2013-04-28: (dfa) Support iBT.what.disp.copyStructural (default enabled
%                   except for RFX and ICA display).
% 2013-04-15: (dfa) If iBT.what.disp.threshold.type is set, use its number
%		    of elements to determine what is displayed (i.e.
%		    ignore extra elements in iBT.what.disp.threshold.thresh)
%		    Removed obsolete code including iBT.what.disp.fewerSlices.
% 2013-01-30: (dfa) Began support of iBT.what.disp.threshold.type.
%		    (presently supports 't', 'p', 'FWE', 'li').
%		    Support iBT.what.disp.threshold.write_t_value_file
% 2013-01-29: (dfa) Exit gracefully if no SPM.mat or no contrasts configured.
% 2012-09-11: (dfa) Rename iBT.what.disp.hold to iBT.what.disp.resample.
% 2012-09-05: (dfa) Support iBT.what.disp.MELODIC.SOCK
% 2012-06-08: (dfa) Support iBT.what.disp.output_parent_folder 
%		    and iBT.what.disp.hold
% 2012-03-12: (dfa) Display version info from iBT_version()
% 2012-02-27: (dfa) Work around bug in SPM8 r4667 without user intervention.
% 2012-02-27: (dfa) Support RFX folder specification via
%		   iBT.who.sub{}.sess{}.RFX.loc to allow multiple RFX 
%		  (in this context each sub and/or sess is actually a
%		   different group random effects analysis, e.g. each
%		   might be a different seed voxel in an FC analysis).
% 2012-02-13: (dfa) Use master session mean when appropriate
% 2011-12-05: (dfa) Corrected display of path when absolute path 
%                   specified in iBT.what.disp.structural
% 2011-11-09: (dfa) Support operation without a header info file.
% 2011-10-17: (dfa) Honour wildcards iBT.what.disp.structural
% 2011-10-12: (dfa) Honour absolute path in iBT.what.disp.structural
% 2011-09-05: (dfa) Updated for public release
%

% 2004-2011: iBT analysis scripts were developed by David Abbott based
%            upon the who & what structure concept of Tony Waites, with 
%            substantial contributions from Matt Harvey, Kaushik Bhaganagarapu,
%            and Richard Masterton.
%            In 2006 the script that ultimately became iBT_display.m
%            was created by Matt Harvey and included code adapted 
%            from display_slices.m v0.4 by R. Henson.
%

% Note; iBT.what.disp.analysis (set in iBT_start.m and possibly modifed in iBT_analysis.m) can have
% the following values:
% 1 - Single-subject single session analysisn_templateLabel
% 2 - Single-subject multi-session analysis
% 3 - Multi-subject random effects analysis
% 4 - Multi-subject single-session fixed-effects analysis
% 5 - ICA
%

disp(sprintf('Entered display script...'));

if isempty(iBT), error('Error: "iBT" structure is not defined.'), end

[iBTversion,iBTrelease,iBTversionReleaseString] = iBT_version(1); % Display version so it gets logged.
[iBT.what.SpmVersion, iBT.what.SpmRelease, iBT.what.SpmVersionReleaseString, MatlabVersion, MatlabDate, MatlabArch] = iBT_SPMversion(2);

% check orientation attributes match in length
if length(iBT.what.disp.orient.transform) > length(iBT.what.disp.orient.slices), error('Error: insufficient number of slice ranges specified (require a range for each orientation).'), end

% check threshold attributes match in length
if length(iBT.what.disp.threshold.mini) ~= length(iBT.what.disp.threshold.maxi), error('Error: threshold lengths do not match.'), end

% check we have enough threshold ranges for the number of thresholds specified (it is OK if we have more):
if length(iBT.what.disp.threshold.thresh) > length(iBT.what.disp.threshold.maxi), error('Error: Not enough entries in iBT.what.disp.threshold.mini and .maxi arrays.'), end

spm('defaults','fmri'); % Ensure we are starting with standard defaults

spm_jobman('initcfg'); % (Prevents welcome screen obscuring output if running headless)

OriginalBackgroundColour=get(spm_figure('GetWin', 'Graphics'), 'Color'); %So can restore later
OriginalTextColor = get(gcf,'DefaultTextColor'); %So can restore later
set(gcf,'DefaultTextColor','white');


try, 
	iBT.what.disp.revOrder; 
	disp(['Warning: iBT.what.disp.revOrder is long deprecated and ignored; Slice ordering is governed by iBT.what.disp.orient.slices']);
catch 
; % Do nothing if it is not defined - that is good
end

try, iBT.what.disp.positive; catch iBT.what.disp.positive=1; end
try, iBT.what.disp.negative; catch iBT.what.disp.negative=1; end

try, iBT.what.disp.show.name; catch iBT.what.disp.show.name=1; end
try, iBT.what.disp.show.ID; catch iBT.what.disp.show.ID=1; end
try, iBT.what.disp.show.date; catch iBT.what.disp.show.date=1; end
try, iBT.what.disp.show.task; catch iBT.what.disp.show.task=1; end
try, iBT.what.disp.show.conName; catch iBT.what.disp.show.conName=0; end
try, iBT.what.disp.show.thresh; catch iBT.what.disp.show.thresh=1; end
try, iBT.what.disp.show.warp; catch iBT.what.disp.show.warp=1; end
try, iBT.what.disp.show.smooth; catch iBT.what.disp.show.smooth=1; end
try, iBT.what.disp.show.voxelSize; catch iBT.what.disp.show.voxelSize=1; end
try, iBT.what.disp.show.procVoxelSize; catch iBT.what.disp.show.procVoxelSize=1; end
try, iBT.what.disp.transparency; catch iBT.what.disp.transparency=0; end
try, iBT.what.disp.show.colourbar; catch iBT.what.disp.show.colourbar=1; end
try, iBT.what.disp.pngname.include.contrastnum; catch iBT.what.disp.pngname.include.contrastnum=1; end
try, iBT.what.disp.statistic; catch iBT.what.disp.statistic='T'; end
try, iBT.what.li.save.path; catch iBT.what.li.save.path=''; end
try, iBT.what.disp.fix_for_SPM8r4667_applied; catch iBT.what.disp.fix_for_SPM8r4667_applied = 0; end
try, iBT.what.disp.output_parent_folder; catch iBT.what.disp.output_parent_folder = ''; end
try, iBT.what.disp.output_subfolder; catch iBT.what.disp.output_subfolder = ''; end
try, iBT.what.disp.resample; catch iBT.what.disp.resample = 1; end
try, iBT.what.disp.MELODIC.SOCK; catch iBT.what.disp.MELODIC.SOCK = 0; end
try, iBT.what.disp.threshold.write_t_value_file; catch iBT.what.disp.threshold.write_t_value_file = 0; end
try, iBT.what.disp.con; catch iBT.what.disp.con = 0; end
try, iBT.what.disp.fontsize1; catch iBT.what.disp.fontsize1 = 15; end
try, iBT.what.disp.fontsize2; catch iBT.what.disp.fontsize2 = 11; end
try, iBT.what.disp.show.mrej; catch iBT.what.disp.show.mrej = 1; end
try, iBT.what.disp.mc.rejLog; catch iBT.what.disp.mc.rejLog = 0; end % 1 = write to logfile number of task and rest scans affected by rejection regressors, 2 = do this this even if no rejection.
try, copyStructural=iBT.what.disp.copyStructural; catch copyStructural = 2; end
try, copyHeaderInfo=iBT.what.disp.copyHeaderInfo; catch copyHeaderInfo = 2; end
try, iBT.what.disp.invert; catch iBT.what.disp.invert = 0; end
try, iBT.what.disp.ROI.img; catch iBT.what.disp.ROI.img = ''; end
try, iBT.what.disp.ROI.transparency; catch iBT.what.disp.ROI.transparency = 0.2; end
try, iBT.what.disp.ROI.colour; catch iBT.what.disp.ROI.colour = [0.19 0.0 0.38]; end % [0.19 0 0.38] is purple.
try, iBT.what.ext.SPM; catch iBT.what.ext.SPM = '.img'; end % Used for images generated by SPM
try, iBT.what.ext.FSL; catch iBT.what.ext.FSL = '.img'; end % Used for images generated by MELODIC
try, iBT.what.ext.iB;  catch iBT.what.ext.iB  = '.img'; end % Used for images generated by iBrain
try, iBT.what.disp.MELODIC.wild; catch iBT.what.disp.MELODIC.wild=['thresh_zstat*' iBT.what.ext.FSL]; end
try, iBT.what.pre.fmriprep.done; catch iBT.what.pre.fmriprep.done = 0; end % Input data is from fmriprep 
try, iBT.what.disp.show.LR; catch iBT.what.disp.show.LR = 0; end % Display a Left ot Right indicator to the left of axial and coronal images
try, iBT.what.disp.fontsizeLR; catch iBT.what.disp.fontsizeLR = 11; end % Not mentioned in runinfo, because changing this will not change line spaces at present
try, iBT.what.disp.orient.showText; catch iBT.what.disp.orient.showText = 1; end % If 0, no orientation strings are generated to label the images. 

try, iBT.what.disp.OTHER.do;
catch 
	try, iBT.what.disp.OTHER.do = 2 * (iBT.what.disp.RFX.do > 0); % use legacy variable if present
	catch;
		iBT.what.disp.OTHER.do = 0;
	end
end

try, iBT.what.disp.together; % Indicates what to do if both iBT.what.disp.positive and .negative are 1
	display_together = iBT.what.disp.together;
	put_sign_in_filename = 2; % Put sign in filename for separate and together maps, so we can generate together and separate maps without a naming conflict
catch
	display_together = 1; % Legacy was one map with both + and - together if both requested
	put_sign_in_filename = 0; % Preserve legacy behaviour for + and - separate files for old runinfo files.
end

reorder_sagittal = 0; % Initialise
try, iBT.what.disp.revflip; % deprecated, but interpreted here for backwards compatibility
	disp(['Warning: iBT.what.disp.revflip is deprecated: Please consider replacing with ']);
	disp(['  iBT.what.disp.orient.LeftOnRight, and if you do, you will probably want to change ']);
	disp(['  the order of sagittal slices in iBT.what.disp.orient.slices, since the x-sign in ']);
	disp(['  iBT 4.x is right positive (consistent with MNI), whereas was left positive in iBT 3.x.']);
	try iBT.what.disp.orient.LeftOnRight; 
		error(["Conflict: Old iBT.what.disp.revflip should not be used when iBT.what.disp.orient.LeftOnRight is set"]);
	catch; 
		reorder_sagittal = 1; 	% Display the requested sagittal slices from iBT.what.disp.orient.slices in reverse order,
								% to result in a display order for legacy runinfos (containing revflip) simlar to that 
								% experienced with iBT version 3.x
	end;
	if strcmpi(iBT.what.disp.revflip,'yes')
		iBT.what.disp.orient.LeftOnRight = 'yes';
		iBT.what.disp.orient.AnteriorOnRight = 'no';
	elseif strcmpi(iBT.what.disp.revflip,'no')
		iBT.what.disp.orient.LeftOnRight = 'no';
		iBT.what.disp.orient.AnteriorOnRight = 'no';
	else
		error(["Error: iBT.what.disp.revflip has unexpected value (should be 'yes' or 'no', but see above re: deprecation)" ]);
	end
catch
; % OK (good) if not defined
end
try, iBT.what.disp.orient.LeftOnRight; catch iBT.what.disp.orient.LeftOnRight = 'no'; end % = SPM default neurological orientation
try, iBT.what.disp.orient.AnteriorOnRight; catch iBT.what.disp.orient.AnteriorOnRight = 'no'; end 

config_num = 0; % initialise

% The following prevents legacy use of negative values in iBT.what.disp.positive and iBT.what.disp.negative
% that was only permitted in the (non-public) iBT versions 3.5 and v3.5a.
% Should now use iBT.what.disp.invert instead.
if ((iBT.what.disp.positive < 0) || (iBT.what.disp.negative < 0))
		disp('WARNING: Obsolete display options: please check iBT.what.disp.positive and iBT.what.disp.negative');
		disp('         They should be 0 or 1');
		disp('         Consider use of iBT.what.disp.invert');
		disp('         Display aborted.');
		return
end
 
if (iBT.what.disp.positive ~= 0) && (iBT.what.disp.negative ~=0)
	if display_together > 0,
		config_num = config_num + 1;
		disp_config{config_num}.positive = iBT.what.disp.positive;
		disp_config{config_num}.negative = iBT.what.disp.negative;
		disp_config{config_num}.invert = iBT.what.disp.invert;
	end
else
	display_together = 0;
end
if (display_together == 0) || (display_together == 2) % then (also) want to display separately
	if iBT.what.disp.positive ~= 0
		config_num = config_num + 1;
		disp_config{config_num}.positive = iBT.what.disp.positive;
		disp_config{config_num}.negative = 0;
		disp_config{config_num}.invert = iBT.what.disp.invert;
	end
	if iBT.what.disp.negative ~= 0
		config_num = config_num + 1;
		disp_config{config_num}.positive = 0;
		disp_config{config_num}.negative = iBT.what.disp.negative;
		disp_config{config_num}.invert = iBT.what.disp.invert;
	end
end


try, iBT.what.disp.threshold.thresh; 
	catch  %This emulates the old behaviour:
		iBT.what.disp.threshold.thresh=iBT.what.disp.threshold.mini
	end

% If iBT.what.disp.threshold.type is not defined, the following creates it as a cell array with as 
% many elements as iBT.what.disp.threshold.thresh, each set to 'li' if thresh == -1, otherwise 't'.
try, iBT.what.disp.threshold.type; catch 
	tmp = {'t','li'};
	iBT.what.disp.threshold.type = cellfun(@(x) tmp{1+(x==-1)}, iBT.what.disp.threshold.thresh, 'UniformOutput', false);
end


awd = iBT.what.where.awd;
sub = iBT.what.do_sub;    % The subject from which to obtain study detail information
sess = iBT.what.disp.ses; % The session from which to obtain study detail information

last_Structural = '';
last_spmT_loc = '';

% Determine spmT.loc and the correct spmT contrasts to use
%===========================================================================
try iBT.who.sub{sub}.sess{sess}.disp.loc;
; % OK if it exists
catch; 
	% If it does not exist, perhaps using legacy RFX.loc:
	try, iBT.who.sub{sub}.sess{sess}.disp.loc = iBT.who.sub{sub}.sess{sess}.RFX.loc; 
	catch; end % 
end

if iBT.what.disp.analysis == 3 
	try
		spmT.loc = iBT.who.sub{sub}.sess{sess}.disp.loc; % Set in runinfo
	catch
		try, spmT.loc = iBT.what.disp.RFX.loc; % Perhaps using this even older legacy variable.
		catch;
			spmT.loc = '.'; % If not set, stay where we are.
		end % 
	end
elseif iBT.what.disp.analysis == 5
	spmT.loc = fullfile(iBT.what.MELODIC.loc,'filtered_func_data.ica','stats'); % Set in iBT_start.m (not really spmT in this context, but location of thresholded stats)
else	
	spmT.loc = iBT.what.stats.dir; % Set in iBT_start.m and possibly adjusted for offset in iBT_analysis.m.		
end

% Determine where we want to write our output 
if strcmp(iBT.what.disp.output_parent_folder, '') == 0
	out_tmp = spmT.loc;
	new_tmp = iBT.what.disp.output_parent_folder;
	if ispc
		if out_tmp(2) == ':' % then need to remove drive letter
			out_tmp = out_tmp(3:end);
		end
		if new_tmp(2) == ':' % create version with removed drive letter for comparison
			new_tmp = new_tmp(3:end);
		end
	end
	% We want to preserve the folder structure of the stats files
	% within the new display location, but do not need to duplicate 
	% top level folders that are the same for both:
	str = pad([string(out_tmp);string(new_tmp)]); % Create strings of same length so can compare
	cmp = char(str(1)) == char(str(2)); % Compare each character
	pos = find(~cmp,1); % Find the position of the first disparate character
	out_tmp = out_tmp(pos:end) ; % remove the common starting elements
	if strcmp(iBT.what.disp.output_subfolder, '') == 0 % i.e. not empty string
		output_folder=fullfile(iBT.what.disp.output_parent_folder,out_tmp,iBT.what.disp.output_subfolder);
	else
		output_folder=fullfile(iBT.what.disp.output_parent_folder,out_tmp);					
	end
else % No new output_parent_folder, but still might be a subfolder:
	if strcmp(iBT.what.disp.output_subfolder, '') == 0 % i.e. not empty string
		output_folder=fullfile(spmT.loc,iBT.what.disp.output_subfolder);
	else
		output_folder=spmT.loc;
	end
end

% Create the output folder if new:
output_folder_exist = exist(output_folder,'dir');
if output_folder_exist ~= 7
	[SUCCESS,MESSAGE,MESSAGEID] = mkdir (output_folder);
	if ~SUCCESS
		disp(['Failed to create display output folder ' output_folder]);
		error(['Please check permissions and try again.']);
	end
end

% We no longer loop over orientations here - moved loop start to much later on to avoid unnecessary re-initialisation - dfa 2013-09-09. 
% Should now re-indent the code to reflect the new loop structure.

	if iBT.what.disp.analysis == 3 % RFX
        	dir_remove = 0; % number of directories to remove from data location for display purposes
	elseif iBT.what.disp.analysis == 5 % ICA
        	dir_remove = 2; % number of directories to remove from data location for display purposes
	else 
        	dir_remove = 4; % number of directories to remove from data location for display purposes
		if (copyStructural == 2), copyStructural = 1; end % Default to enabled if not RFX or ICA.
		if (copyHeaderInfo == 2), copyHeaderInfo = 1; end % Default to enabled if not RFX or ICA.
	end
        cd(awd);

	    % Set some defaults in case the fields we use are not present in the infoFile file:
	    mapHeader{sess}.name = ''; 
	    mapHeader{sess}.ID = '';
	    mapHeader{sess}.sdate = '';
	    mapHeader{sess}.TR = '?';
	    mapHeader{sess}.r_smooth = '?';
	    mapHeader{sess}.n_smooth = '?';
	    mapHeader{sess}.slicethickness = '?';
	    mapHeader{sess}.spacing = '?';
	    mapHeader{sess}.pixelx = '?';
	    mapHeader{sess}.pixely = '?';
	    mapHeader{sess}.StudyNo = '?';
	    mapHeader{sess}.StudyDesc = '?';
	    mapHeader{sess}.SeriesDesc = '?';
	    mapHeader{sess}.n_templateFile = '?';
	    mapHeader{sess}.n_templateLabel = '?';
	    mapHeader{sess}.n_iterations = '?';
	    mapHeader{sess}.n_regtype = '?';

	    % Read infoFile for subject details to be used in map header
	    % We use the one in the spmT.loc if it is present there, otherwise we get it from the 
	    % pre-processing folder (and optionally copy it into the spmT.loc for later use)
	    infoFile = fullfile(spmT.loc,'_header_info.txt');
	    if ~exist(infoFile,'file')
	       try infoFile = iBT.who.sub{sub}.sess{sess}.infoFile; catch infoFile = ''; end
	       if  ~strcmp( infoFile, '') 
	       		if (copyHeaderInfo == 1) && exist(infoFile,'file') && exist(spmT.loc,'dir'),
				[SUCCESS,MESSAGE,MESSAGEID] = copyfile(infoFile,spmT.loc);
				if (SUCCESS ~= 1 ) && strcmp(iBT.what.disp.output_parent_folder, '')
					% We failed to copy, even though our output directory is here 
					% so we expect to be able to write to this folder.
					disp(['Warning: Unable to copy ' infoFile ' into ' spmT.loc]);
					disp(MESSAGE);
				end % Alternatively, if output_parent_folder is specified, then not unexpected that cannot write to spmT.loc so ignore the error. Maybe we should copy it to output_folder in this case?
			end
	       end 
	    end %if ~exist(infoFile,'file')
	    
	    if  ~strcmp( infoFile, '')
	    	fsettings = fopen(infoFile); % open settings file and read the data
		if fsettings == -1
			if (iBT.what.pre.header_info == 0) || (iBT.what.disp.analysis == 3) || (iBT.what.disp.analysis == 5) 
				%That's OK, we probably deliberately did not create one.
				mapHeader{sess}.info_version = -1;
			else
				disp(['Unable to open header information file ' infoFile]);
				disp('If you want to proceed without header information, set iBT.what.pre.header_info = 0.');
				error(['Please fix the error noted above and try again']);
			end
		else
	    	    mapHeader{sess}.buffer = fgetl(fsettings); % version header text

	            if ~strcmpi(mapHeader{sess}.buffer, 'Info_version=')
			disp(sprintf('Error: Unfamiliar header info file format in: %s',infoFile));
			error('Now terminating script execution prematurely - please fix problem and try again.');
	            end
	            mapHeader{sess}.buffer = fgetl(fsettings); % version header number
	            mapHeader{sess}.info_version = str2num( mapHeader{sess}.buffer );
		end
	    else % no info file
	    	mapHeader{sess}.info_version = -1
	    end % if ~strcmpi( infoFile, '')
		
	    if mapHeader{sess}.info_version == -1 % No info file
		disp(sprintf('Warning: No header info file specified.'));
	    elseif mapHeader{sess}.info_version == 1.0
		disp(sprintf('Found old-style header info file.')); % No labels, so we need to rely on the position of fields.
		mapHeader{sess}.name = fgetl(fsettings); 
	        mapHeader{sess}.ID = fgetl(fsettings);
	        mapHeader{sess}.sdate = fgetl(fsettings);
		mapHeader{sess}.TR = fgetl(fsettings);
        	mapHeader{sess}.r_smooth = fgetl(fsettings);
        	mapHeader{sess}.n_smooth = fgetl(fsettings);
		mapHeader{sess}.slicethickness = fgetl(fsettings);
        	mapHeader{sess}.spacing = fgetl(fsettings);
	       	mapHeader{sess}.pixelx = fgetl(fsettings);
	        mapHeader{sess}.pixely = fgetl(fsettings);
	        % The remaining items we might like are not present in this type of header.
	    else 				
		if  ( (mapHeader{sess}.info_version) >= 4.0 )
			disp(sprintf('ERROR: Unfamiliar version of header info found (v%s).',mapHeader{sess}.buffer));
			error(sprintf('      Please update your iBT_display.m '));
		end
		if  ( (mapHeader{sess}.info_version) < 2.0 ) || ( (mapHeader{sess}.info_version) > 3.0 )
			disp(sprintf('Warning: Unfamiliar version of header info found (v%s).',mapHeader{sess}.buffer));
			disp(sprintf('         ( consider updating your iBT_display.m ).'));
		end

	    	  while 1 == 1 ; % We will continue until we get an error or EOF
	    	    buffer = fgetl(fsettings); % read next line from file
		    if buffer == -1 break; end;
	            [token,value] = strtok(buffer,':'); % Split it into token and value
		    if strcmpi(token,'iBrain(TM) standard header output version') 
        	 	mapHeader{sess}.iBrainHeaderVersionString = strrep(value,':',' ');
        	 	mapHeader{sess}.iBrainHeaderVersion = str2num( mapHeader{sess}.iBrainHeaderVersionString  );
			if (mapHeader{sess}.iBrainHeaderVersion) ~= 1.1 
			  disp(sprintf('Warning: Unfamiliar version of iBrain header output (v%s).',mapHeader{sess}.iBrainHeaderVersionString));
			  disp(sprintf('         Consider updating your iBT_display.m ).'));
		 	end	
		    elseif strcmpi(token,'Subject name') mapHeader{sess}.name = strrep(value,':',' '); 
		    elseif strcmpi(token,'Subject ID') mapHeader{sess}.ID = strrep(value,':',' '); 
		    elseif strcmpi(token,'Study date (yyyy-mm-dd)') mapHeader{sess}.sdate = strrep(value,':',' '); 
		    % Note that currently we do not actually use the following three fields, but if we ever do in future here they are...
		    elseif strcmpi(token,'Study Number') mapHeader{sess}.StudyNo = strrep(value,':',' '); 
		    elseif strcmpi(token,'Study Description') mapHeader{sess}.StudyDesc = strrep(value,':',' '); 
		    elseif strcmpi(token,'Series Description') mapHeader{sess}.SeriesDesc = strrep(value,':',' '); 
		    elseif strcmpi(token,'TR(s)') mapHeader{sess}.TR = strrep(value,':',' ');
		    elseif strcmpi(token,'Voxel size X (mm)') mapHeader{sess}.pixelx = strrep(value,':',' '); 
		    elseif strcmpi(token,'Voxel size Y (mm)') mapHeader{sess}.pixely = strrep(value,':',' '); 
		    elseif strcmpi(token,'Slice thickness') mapHeader{sess}.slicethickness = strrep(value,':',' ');
		    elseif strcmpi(token,'Slice gap') mapHeader{sess}.spacing = strrep(value,':',' ');
		    elseif strcmpi(token,'r_smth') mapHeader{sess}.r_smooth = strrep(value,':',' ');
		    elseif strcmpi(token,'n_smth') mapHeader{sess}.n_smooth = strrep(value,':',' ');
		    elseif strcmpi(token,'n_template') mapHeader{sess}.n_templateFile = deblank(strrep(value,':',' '));
		    elseif strcmpi(token,'n_templateLabel') mapHeader{sess}.n_templateLabel = deblank(strrep(value,':',' '));
		    elseif strcmpi(token,'n_iterations') mapHeader{sess}.n_iterations = str2num( strrep(value,':',' ') );
		    elseif strcmpi(token,'n_regtype') mapHeader{sess}.n_regtype = strrep(value,':','');
		    end; %if strcmp
		  end; % while.  Ignore everything else in file. 
	    end; % if mapHeader{sess}.info_version
	    if mapHeader{sess}.info_version ~= -1, fclose(fsettings); end

	    % In the following, the test for '?' (i.e. indicating not present, as set above) ensures backwards compatibility:
	    if strcmp(mapHeader{sess}.n_regtype,'?') || strcmpi(mapHeader{sess}.n_regtype,'mni') ,
	    	if mapHeader{sess}.n_iterations == 0,
			mapHeader{sess}.n_constraint = 'affine';
		else
			mapHeader{sess}.n_constraint = 'warp';		
		end
	    else
	    	mapHeader{sess}.n_constraint = mapHeader{sess}.n_regtype ; % Likely 'rigid'.
	    end
	    
	    % We will display the spatial normalisation template descriptive label if present, otherwise the filename:
	    if ~strcmp(mapHeader{sess}.n_templateLabel,'?') ,
	    	mapHeader{sess}.n_template = mapHeader{sess}.n_templateLabel;
	    else
	    	mapHeader{sess}.n_template = mapHeader{sess}.n_templateFile;
	    end
	    	    	    
	    % Check for header overrides
	    try,
	      fsettings = fopen(fullfile(iBT.who.sub{sub}.sess{sess}.loc, '_header_override.txt')); % open override settings file and read the data
	    catch
	      fsettings = -1;
	    end
	    if fsettings ~= -1
	    	try
	    	  while 1 == 1 ; % We continue until get an error or EOF
	   		buffer = fgetl(fsettings);
		        if buffer == -1 break; end;
	         	[token,value] = strtok(buffer,':');
			if strcmpi(token,'Subject name') mapHeader{sess}.name = strrep(value,':',' '); 
			elseif strcmpi(token,'Subject ID') mapHeader{sess}.ID = strrep(value,':',' '); 
			% Note that currently we don not actually use the following, but if 
			% we ever do in future here they are (also should set up a default above)...
			elseif strcmpi(token,'Study Number') mapHeader{sess}.StudyNo = strrep(value,':',' '); 
			elseif strcmpi(token,'Study Description') mapHeader{sess}.StudyDesc = strrep(value,':',' '); 
			elseif strcmpi(token,'Series Description') mapHeader{sess}.SeriesDesc = strrep(value,':',' '); 
			end; % Ignore everything else in file.
		  end % while
	    	catch
		  fclose(fsettings);
		end
	    end; % if fsettings ~= -1

      	    if iBT.what.disp.show.mrej
      		mrej = iBT_mrej(fullfile(spmT.loc,'SPM.mat'), iBT.what.disp.mc.rejLog, mapHeader{sess}.ID);
      	    end

	for current_threshold = 1:1:length(iBT.what.disp.threshold.type)
	    stats_k_threshold = 0; %initialise
	    xmini = iBT.what.disp.threshold.mini{current_threshold};
	    xmaxi = iBT.what.disp.threshold.maxi{current_threshold};

	    cd_OK = 1; % initialise
	    if iBT.what.disp.analysis == 3
        	cd(awd); % The user can create the mean for overlay and store it here
		try, mapHeader{sess}.name = [iBT.what.disp.OTHER.title ': '];
		catch
			if iBT.what.disp.OTHER.do == 2,
				mapHeader{sess}.name = 'RFX group analysis: ';
			else
				mapHeader{sess}.name = '';
			end
		end
		mapHeader{sess}.ID = iBT.who.sub{sub}.ID;
		mapHeader{sess}.sdate = '';
		if strcmpi(iBT.what.disp.structural,'auto')
		        structural_filename = ['group_mean_wmean' iBT.what.ext.SPM];
		else
			structural_filename = iBT.what.disp.structural;
		end
		mapHeader{sess}.smooth = mapHeader{sess}.n_smooth;		
	    elseif strcmpi(iBT.what.disp.type,'norm') == 1
                if(strcmpi(iBT.what.disp.structural,'auto'))
			if iBT.what.disp.analysis == 5;
				cd(iBT.who.sub{sub}.sess{1}.MELODIC.loc)
		        	structural_filename = ['mean_func' iBT.what.ext.SPM];			
			elseif (iBT.what.pre.norm.intra.do ~= 0) || ( iBT.what.pre.norm.BiasCorrect ~= 0 )
				cd(iBT.who.sub{sub}.sess{sess}.loc);% The mean on which to overlay is stored here
				structural_filename = ['wmean_all' iBT.what.ext.SPM];
			else
		       		if (iBT.what.pre.sessions_together == 1) 
					cd(iBT.who.sub{sub}.sess{iBT.who.sub{sub}.MasterSession}.loc); % The mean on which to overlay is stored here
				else
					cd(iBT.who.sub{sub}.sess{sess}.loc);% The mean on which to overlay is stored here
				end
				if (iBT.what.pre.fmriprep.done == 1),
					if (iBT.what.pre.iBrainRealignTarget == 0) | (iBT.what.pre.Realign == 0),
					  structural_filename = ['mean*' iBT.what.ext.SPM];
					else
					  % by dfa 2025-01-16: Untested: Note we would get here if do multi-session realignment in iBT, which has not been tested with fmriprep input.
					  structural_filename = ['mean_median*' iBT.what.ext.SPM]; 
					end
				else
					if (iBT.what.pre.iBrainRealignTarget == 0) | (iBT.what.pre.Realign == 0),
					  structural_filename = ['wmean*' iBT.what.ext.SPM];
					else
					  structural_filename = ['wmean_median*' iBT.what.ext.SPM]; 
					end
				end % if (iBT.what.pre.fmriprep.done == 1) else
			end % if iBT.what.disp.analysis == 5 else
		else
			if iBT.what.disp.analysis == 5;
				cd(iBT.who.sub{sub}.sess{1}.MELODIC.loc)
			else
       				cd(iBT.who.sub{sub}.sess{sess}.loc);% The mean on which to overlay is stored here
			end
                        structural_filename = iBT.what.disp.structural;
		end
		mapHeader{sess}.smooth = mapHeader{sess}.n_smooth;
	    else
       		try, cd(iBT.who.sub{sub}.sess{sess}.loc);% The mean on which to overlay is stored here
		catch, cd_OK = 0; end % catch
		if (strcmpi(iBT.what.disp.structural,'auto'))
			if iBT.what.disp.analysis == 5;
       			    try, cd(iBT.who.sub{sub}.sess{1}.MELODIC.loc); catch, cd_OK = 0; end %catch
		            structural_filename = ['mean_func' iBT.what.ext.SPM];			
			elseif (iBT.what.pre.iBrainRealignTarget == 0) | (iBT.what.pre.Realign == 0) | ( (iBT.what.pre.Realign < 0) & (iBT.what.pre.fmriprep.done == 1) ),
			  structural_filename = ['mean*' iBT.what.ext.SPM];
			else
			  structural_filename = ['mean_median*' iBT.what.ext.SPM];
			end
		else
			if iBT.what.disp.analysis == 5;
				cd(iBT.who.sub{sub}.sess{1}.MELODIC.loc)
			end			
                        structural_filename = iBT.what.disp.structural;
		end
		mapHeader{sess}.smooth = mapHeader{sess}.r_smooth;
	    end

	    if iBT.what.disp.hideName == 1 && iBT.what.disp.analysis ~= 3
                mapHeader{sess}.name = iBT.who.sub{sub}.ID;
	    end
	    
	    
	    if iBT.what.disp.analysis == 5
	    
	      try tmp = iBT.what.MELODIC.smooth; % Only use this if set, otherwise ignore.
	      	mapHeader{sess}.smooth = sprintf('%g',tmp);
	      catch
	      end % try

	    end %if iBT.what.disp.analysis == 5
	    
	try
	    if iBT.what.disp.hideHeaderID == 1 && iBT.what.disp.analysis ~= 3
                mapHeader{sess}.ID = '';
	    end
	catch
	end; % try iBT.what.disp.hideHeaderID 
	
		[PATHSTR_structural,NAME_structural, EXT_structural] = fileparts(structural_filename);
		NAME_structural = [NAME_structural EXT_structural];
                if strcmp(PATHSTR_structural,'') % If no path specified, assume in pwd, otherwise leave as-is...
			mapHeader{sess}.imgs{1} = fullfile(pwd,NAME_structural);
		else
			mapHeader{sess}.imgs{1} = fullfile(structural_filename); 
		end; %  if strcmp()
		if cd_OK == 1,
			wildcard_expanded_names = dir(mapHeader{sess}.imgs{1}); % Wildcard-expanded filename without path
		else
			wildcard_expanded_names = []; % The desired directory could not be found, so do not set any names.	
		end
		if length(wildcard_expanded_names) == 0,
			% Not found in pre-processing directory, so we hope it was previously copied to stats folder.
			original_structural_found = 0;
			wildcard_expanded_names = dir( fullfile(spmT.loc,NAME_structural) );
			if length(wildcard_expanded_names) == 0,
		  	  % Not found in stats directory either, so warn user
		  	  %dbclear if error; %Do not enter debug at this error
		  	  error(sprintf('Error: unable to find file matching %s or %s',fullfile(mapHeader{sess}.imgs{1}),fullfile(spmT.loc,NAME_structural)));
			  mapHeader{sess}.imgs{1} = '';
			else
		  	  % We support wildcards and take the first match:
			  PATHSTR_structural = spmT.loc;
		  	  mapHeader{sess}.imgs{1} = fullfile(PATHSTR_structural,wildcard_expanded_names(1).name);
			end
		else
		  % We support wildcards and take the first match:
		  [PATHSTR_structural,NAME_structural, EXT_structural] = fileparts(mapHeader{sess}.imgs{1});
		  mapHeader{sess}.imgs{1} = fullfile(PATHSTR_structural,wildcard_expanded_names(1).name);
		  original_structural_found = 1;
		end		
		
		cd(awd);

		if ~strcmp(mapHeader{sess}.imgs{1},'')
			disp(['Overlay on to ' mapHeader{sess}.imgs{1}]);
		else
			disp(['Overlay on to nothing.']);
		end

		% shorten header by removing unnecesssary directory info
		
		if dir_remove > 0 
			[tok,rem] = strtok(PATHSTR_structural,filesep());
		else
			rem = PATHSTR_structural;
		end
		for i=2:dir_remove
	        	[tok,rem] = strtok(rem,filesep());
		end
		rem = strcat(rem,filesep());
		mapHeader{sess}.himgs{1} = char(strcat(rem, NAME_structural));
		stringsize = size(mapHeader{sess}.himgs{1},2);
		if (stringsize > 55), % Limit size of string to the final 55 characters
			mapHeader{sess}.himgs{1} = [ '...' mapHeader{sess}.himgs{1}(stringsize-54:end) ];
		end 

		if (~strcmp(mapHeader{sess}.imgs{1},'')) && (copyStructural == 1) && original_structural_found && ( (~strcmpi(last_Structural,mapHeader{sess}.imgs{1})) || (~strcmpi(last_spmT_loc,spmT.loc)) )
			  last_Structural = mapHeader{sess}.imgs{1};
			  last_spmT_loc = spmT.loc;
			  if exist(spmT.loc,'dir')
			  	[SUCCESS,MESSAGE,MESSAGEID] = copyfile(mapHeader{sess}.imgs{1},spmT.loc);
			  else
			  	SUCCESS = 0;
				MESSAGE = ['Directory not found.'];
			  end
          		  if SUCCESS ~=1,
            			disp(['WARNING: Unable to copy ' mapHeader{sess}.imgs{1} ' into ' spmT.loc]);            	
				disp(MESSAGE);
			  else
			    disp(['Copied ' mapHeader{sess}.imgs{1} ' into directory ' spmT.loc]);
			    if strcmpi(EXT_structural, iBT.what.ext.SPM),
			      [PATHSTR_structural,NAME_structural, EXT_structural] = fileparts(mapHeader{sess}.imgs{1});
			      tmp = fullfile(PATHSTR_structural,[NAME_structural '.hdr']);
			      if exist(tmp,'file') == 2,
				[SUCCESS2,MESSAGE,MESSAGEID] = copyfile(tmp,spmT.loc);
          			if SUCCESS2 ~=1,
            				disp(['WARNING: Unable to copy ' tmp ' into ' spmT.loc]);            	
					disp(MESSAGE);
			   	else
					disp(['Copied ' tmp ' into directory ' spmT.loc]);
				end %if SUCCESS2 ~=1,
			      end % if exist(tmp)
			    end % if strcmpi(EXT_structural, iBT.what.ext.SPM)
			  end %if SUCCESS ~=1,
		end %if (copyStructural == 1)
		
		offset_text =  iBT.what.stats.offset_text;  % Set in iBT_analysis.m
		
 		mapHeader{sess}.task = offset_text;
		
		if iBT.what.disp.analysis == 1
		  
		  mapHeader{sess}.task = strcat(iBT.what.stats.tasks,mapHeader{sess}.task);	
		  
		end %if iBT.what.disp.analysis == 1
		cd(spmT.loc);
		manual_con = 0; %Initialise
		manual_filenames = ~ ( strcmpi(iBT.what.disp.contrast,'auto') || strcmpi(iBT.what.disp.contrast,'+effect') || strcmpi(iBT.what.disp.contrast,'auto-') || strcmpi(iBT.what.disp.contrast,'-effect') );
		
		if ~manual_filenames,
		  if iBT.what.disp.analysis == 5
		    melodic_files = dir(iBT.what.disp.MELODIC.wild); 
		    total_contrasts = numel(melodic_files);
		    T_start = 1
		    if iBT.what.disp.MELODIC.SOCK > 0
			SOCK_results = fullfile(iBT.what.MELODIC.loc,'IC.mat') % iBT.what.MELODIC.loc set in iBT_start.m
			% Load info for ICs for which display is required (SOCK classifies components as 1,2 or 3, with 3 being definite artefact)
			ICstruct = load(SOCK_results);
			cat1_ICs = find(ICstruct.IC.cat(:)==1);
			
			cat2_ICs = find(ICstruct.IC.cat(:)==2);
			
			cat3_ICs = find(ICstruct.IC.cat(:)==3);
			
			switch iBT.what.disp.MELODIC.SOCK % Up to what category do we want to display?
				case 1
					observe_ICs = [cat1_ICs];
				case 2
					observe_ICs = [cat1_ICs;cat2_ICs];
				otherwise
					observe_ICs = [cat1_ICs;cat2_ICs;cat3_ICs];
			end % switch
			
		  	total_contrasts = length(observe_ICs);
		  	T_start = 1;
		  	spmm.SPM.xCon(1).name='';
		  	spmm.SPM.xCon(1).Vspm.fname = iBT.what.disp.contrast;
		    else
	              for T_count = 1:1:total_contrasts
		      	[pathstr,name,ext] = fileparts(melodic_files(T_count).name);
		      	spmm.SPM.xCon(T_count).name=name;
		      	spmm.SPM.xCon(T_count).Vspm.fname = [name ext];
		      end % for
		    end %if iBT.what.disp.MELODIC.SOCK
		  else % Must be SPM analysis
		    try 
                    	spmm = load('SPM.mat'); % load in SPM.mat
		    catch
		    	disp(['WARNING: Unable to load' fullfile(pwd,'SPM.mat')]);
			disp('          Skipping display of this analysis.');
			return
		    end
		    try 
                    	total_contrasts = length(spmm.SPM.xCon); % find the total number of contrasts
		    catch
		    	disp(['WARNING: No contrasts configured in ' fullfile(pwd,'SPM.mat')]);
			disp('          Skipping display of this analysis.');		    	
			return
		    end
		    if sum(abs(iBT.what.disp.con)) ~= 0, % then have requested specific contrast number(s)
		    	con_range = []; % initialise
			for rc = 1:numel(iBT.what.disp.con)
				if iBT.what.disp.con(rc) < 0, 
					if (-iBT.what.disp.con(rc)) > total_contrasts
						if total_contrasts > 1, plural = 's'; else plural = ''; end
		    	    			disp(sprintf('WARNING: Skipping contrast %i counting back from the end, since only %i contrast%s configured in %s',-iBT.what.disp.con(rc),total_contrasts,plural,fullfile(pwd,'SPM.mat')));
					else
						con_range = [con_range (total_contrasts + 1 + iBT.what.disp.con(rc))]; % Negative means we count down with the last contrast indicated by -1.	    	
			  		end
		      		elseif  iBT.what.disp.con(rc) > total_contrasts,
					if total_contrasts > 1, plural = 's'; else plural = ''; end
		    			disp(sprintf('WARNING: Skipping contrast %i as only %i contrast%s configured in %s',T_start,total_contrasts,plural,fullfile(pwd,'SPM.mat')));		    	
				else
		      			con_range = [con_range iBT.what.disp.con(rc)];
		      		end
			end %for rc
			if sum(abs(con_range)) == 0,
		    		disp(sprintf('WARNING: No valid contrasts selected for display.'));		    	
				return % All selected contrasts must have been out of range, so nothing to do!
			else
				manual_con = 1; % Flag that we have explicitly selected valid contrasts.
			end
		    else % a particular iBT.what.disp.con not specified so display all T contrasts (except -effects which will be skipped later on)
	              for T_start = 1:1:total_contrasts
         		if spmm.SPM.xCon(T_start).STAT ~= 'F' % exit once the first T contrast is found
                		break;
                	end
                      end
		    end % if iBT.what.disp.con ~= 0,
		  end % if iBT.what.disp.analysis
		else % if ~manual_filenames else 
		   % Contrast name(s) specified explicitly
		   [pathstr,name,ext]=fileparts(iBT.what.disp.contrast);
		   if isempty(pathstr), pathstr = pwd; end
		   filespec = [name,ext];
		   fullfilespec = fullfile(pathstr,filespec);
		   switch (iBT.what.SpmVersion),
        	   case 'SPM2',
            	  	contrast_filenames   = spm_get('Files',pathstr, filespec);
            	   otherwise,
          	  	contrast_filenames   = spm_select('FPList',pathstr, iBT_spm_wildconvert(filespec));
     	    	   end %case	    
		   total_contrasts = size(contrast_filenames,1);		 
		   if  (total_contrasts < 1),
		   	if (iBT.what.disp.contrastMissingOK ~= 1)
				error(['File not found matching ' fullfilespec]);
		   	else
				disp(['Warning: File not found matching ' fullfilespec]);
				return % Nothing to do.
		   	end %if (iBT.what.disp.contrastMissingOK ~= 1)
		   end % if  (total_contrasts < 1),
		   T_start = 1;
	           for T_count = 1:1:total_contrasts
		      	[pathstr,name,ext] = fileparts(contrast_filenames(T_count,:));
		      	spmm.SPM.xCon(T_count).name=name;
		      	spmm.SPM.xCon(T_count).Vspm.fname = contrast_filenames(T_count,:);
		   end % for	   		   
		end % if ~manual_filenames else 
		
                if ~manual_con,
			con_range = T_start:1:total_contrasts;  % Unless con_range already specified, loop for each T contrast in this analysis
		end
                for T_count = con_range % for each selected contrast, generate a map (except -effects which will be skipped later on)
			contrast_columns = []; % Initialise
			if manual_filenames && (~strcmpi(iBT.what.disp.threshold.type{current_threshold},'t'))
				display_this = 0; % We only support 't' for raw display of whatever filename is manually specified as we have no access to SPM contrast information.
				disp('Warning: iBT.what.disp.contrast specifies filenames manually, however iBT.what.disp.threshold.type is not ''t''.');
				disp('         We only support ''t'' for iBT.what.disp.threshold.type to effect raw display of data in manually specified filenames');
				disp('         as we have no access to contrast information from SPM.mat. Please consider using iBT.what.disp.contrast = ''auto'',')
				disp('         and perhaps also iBT.what.disp.con')
			elseif iBT.what.disp.MELODIC.SOCK
				switch ICstruct.IC.cat(observe_ICs(T_count))
				  case 1
					SOCK_label = 'Accept';
				  case 2
					SOCK_label = 'Maybe';
				  case 3
					SOCK_label = 'Reject';
				end % switch observe_ICs

				if ICstruct.IC.edge_flag(observe_ICs(T_count)) == 1
					High_edge_label = '_High_edge';
				else
					High_edge_label = '';
				end

				if ICstruct.IC.CSF_flag(observe_ICs(T_count)) == 1
					High_CSF_label = '_High_CSF';
				else
					High_CSF_label = '';
				end
				SOCK_string = ['_' SOCK_label High_edge_label High_CSF_label];
				display_this = 1; % We have already selected the observations, so OK to display
			else
				SOCK_string = '';
				if strcmpi(iBT.what.disp.contrast, 'auto')
					display_this = ~length(findstr(spmm.SPM.xCon(T_count).name,'-effect')); % Do not display if it is a negative contrast
				elseif strcmpi(iBT.what.disp.contrast,'auto-')
					display_this = 1; % With this option we display all contrasts of interest
					display_this = display_this - ( 2 * (0<length(findstr(spmm.SPM.xCon(T_count).name,'-effect'))) ); % Change sign of flag if it is a negative contrast
				elseif strcmpi(iBT.what.disp.contrast,'+effect')
					display_this = 0 < length(findstr(spmm.SPM.xCon(T_count).name,'+effect')); % Only display if it is labelled with +effect
				elseif strcmpi(iBT.what.disp.contrast,'-effect')
					display_this = -1 * (0 < length(findstr(spmm.SPM.xCon(T_count).name,'-effect'))); % Only display if it is labelled with -effect
				else
					display_this = 1;
				end % if strcmpi(iBT.what.disp.contrast... else ...
				
				try,
					contrast_columns = find( spmm.SPM.xCon(T_count).c ) ; % Which design matrix columns are involved in this contrast? 
				catch
				end
			end
						
			if ( display_this ~= 0  )

			  if strcmpi(iBT.what.disp.threshold.type{current_threshold},'li')   % Load threshold details determined by laterality script
				% Get filename (code copied from iBT_laterality.m):
				active.subject_ID = iBT.who.sub{iBT.what.do_sub}.ID;
				active.session = iBT.what.disp.ses; % Session number
				active.tasks = iBT.what.stats.tasks; 
				active.offset_text = iBT.what.stats.offset_text; %If we have shifted the HRF in time
				fileIDstring = sprintf('%s_%02d_%s%s',active.subject_ID, active.session, active.tasks, active.offset_text); 
				if strcmp(iBT.what.li.save.path,''),
					li_path = spmT.loc;
				else
					li_path = iBT.what.li.save.path;
				end
		 		mat_filename=fullfile(li_path, sprintf('li_display_thresholds_%s_%.4i.mat',[fileIDstring '_nvox'],T_count) );
				try
					load(mat_filename);
					li_load_ok = 1;
				catch
					li_load_ok = 0;
				end
				if li_load_ok == 1
				  stats_threshold = li_display_thresholds.actThres_at_min_valid_Nvox;
				  li_string = sprintf('Activation LI = %.2f @ N=%i, t=%.2f' , li_display_thresholds.actLI_at_min_valid_Nvox, li_display_thresholds.controls_min_valid_Nvox, li_display_thresholds.actThres_at_min_valid_Nvox);
				  disp(li_string);
				else 
				  disp(sprintf('Warning: unable to load adaptive threshold file: %s',mat_filename));
				  stats_threshold = NaN;
				end
			  elseif strcmpi(iBT.what.disp.threshold.type{current_threshold},'t')   
				stats_threshold = iBT.what.disp.threshold.thresh{current_threshold};
				disp(sprintf('Specified threshold = %f', stats_threshold));
			  else % Need to work out the t-value and k-value to use from SPM.mat for each contrast

				% Determine t-threshold for this contrast from SPM.mat. FWE and uncorrected p thresholds need only be calculated 
				% once per session; FWEc and FDRc requires re-calculation for each contrast.
				new.swd = spmT.loc; % directory containing current SPM.mat
				new.Ic=T_count; % indices of contrasts (in SPM.xCon)
				new.title= spmm.SPM.xCon(T_count).name;
				new.k=0; % initialise extent threshold {voxels}
				if strcmpi(iBT.what.disp.threshold.type{current_threshold},'FWE')
					new.thresDesc = 'FWE';
				else
					new.thresDesc = 'none';
				end
				if strcmpi(iBT.what.disp.threshold.type{current_threshold},'FWEc') || strcmpi(iBT.what.disp.threshold.type{current_threshold},'FDRc')
					new.u=iBT.what.disp.threshold.thresh{current_threshold}(1); % Feature-selecting uncorrected p value				
				elseif strcmpi(iBT.what.disp.threshold.type{current_threshold},'p') || strcmpi(iBT.what.disp.threshold.type{current_threshold},'FWE') 
					new.u=iBT.what.disp.threshold.thresh{current_threshold}; %  p value				
				else
		  			error(sprintf('Error: Unknown threshold type: %s',iBT.what.disp.threshold.thresh{current_threshold}));
				end
				new.pm=[]; % p-value for masking (uncorrected)
				new.Im=[]; % indices of masking contrasts (in xCon)  
				new.Ex=[]; % flag for exclusive or inclusive masking
				new.n=1;  % conjunction number <= number of contrasts
				% TODO: This can be potentially computationally inefficient because spm_getSPM could be called multiple times
				% for the same info (for example, if we generate both FWE and FDRc for the same contrast). Ideally we
				% shouldn't call again if the previous calculation remains relevant (or does spm_getSPM already realise this?)...
				[nSPM,nxSPM] = spm_getSPM(new); % Determine thresholds for FWE/FDR voxel and cluster for this contrast.		   		
				stats_threshold = nxSPM.u;
	
				if (iBT.what.disp.threshold.write_t_value_file) && strcmpi(iBT.what.disp.threshold.type{current_threshold},'p'),
					indicator_filename=sprintf('threshFor_spmT_%04i_p=%g_is_t=%g',T_count,iBT.what.disp.threshold.thresh{current_threshold},stats_threshold);
					disp(indicator_filename);
				    	fclose( fopen(  indicator_filename, 'W' ) );
				end	

				if (iBT.what.disp.threshold.write_t_value_file) && strcmpi(iBT.what.disp.threshold.type{current_threshold},'FWE'),
					indicator_filename = sprintf('threshFor_spmT_%04i_FWE=%g_is_t=%g',T_count,iBT.what.disp.threshold.thresh{current_threshold},stats_threshold);
					disp(indicator_filename);
				    	fclose( fopen(  indicator_filename, 'W' ) );
				end	

				if strcmpi(iBT.what.disp.threshold.type{current_threshold},'FWEc'),
					if iBT.what.disp.threshold.thresh{current_threshold}(2) ~= 0.05,
						disp('WARNING: iBT presently only supports calculation of FWEc<0.05; ignoring different user-supplied threshold');
						% TODO: This is a limitation of spm_getSPM that could easily be fixed.
						iBT.what.disp.threshold.thresh{current_threshold}(2) = 0.05;
					end
				 	stats_k_threshold = nxSPM.uc(3);
					if isinf(stats_k_threshold),
						disp(sprintf('k threshold required to achieve FWEc < %f for contrast #%i is Inf',iBT.what.disp.threshold.thresh{current_threshold}(2),T_count));
						if (iBT.what.disp.threshold.write_t_value_file)
							indicator_filename=sprintf('kthreshFor_spmT_%04i_p=%g_t=%g_FWEc=0.05_is_Inf',T_count,iBT.what.disp.threshold.thresh{current_threshold}(1),stats_threshold);
				    			fclose( fopen( indicator_filename, 'W' ) );
						end
					else
						disp(sprintf('k threshold required to achieve FWEc < %f for contrast #%i is %i',iBT.what.disp.threshold.thresh{current_threshold}(2),T_count,stats_k_threshold));
						new.k=stats_k_threshold;
						[nSPM,nxSPM] = spm_getSPM(new); %Obtain an image thresholded in this way
						activation = zeros(nxSPM.DIM');
						ind = sub2ind(nxSPM.DIM', nxSPM.XYZ(1,:), nxSPM.XYZ(2,:), nxSPM.XYZ(3,:));
						activation(ind) = nxSPM.Z;
						%  Write out the volumes:
						Vout = nSPM.xCon(T_count).Vspm;
						Vout.fname = [ sprintf('kthresh_p=%g_t=%g_FWEc=0.05_k=%i_' ,iBT.what.disp.threshold.thresh{current_threshold}(1),stats_threshold,stats_k_threshold), Vout.fname];
						spm_write_vol(Vout, activation);
	                        		mapHeader{sess}.imgs{2} = fullfile(spmT.loc, Vout.fname);
						mapHeader{sess}.imgs{3} = mapHeader{sess}.imgs{2};
						mapHeader{sess}.himgs{2} = Vout.fname;
						mapHeader{sess}.himgs{3} = Vout.fname;
					end

				elseif strcmpi(iBT.what.disp.threshold.type{current_threshold},'FDRc')
					if iBT.what.disp.threshold.thresh{current_threshold}(2) ~= 0.05,
						disp('WARNING: iBT presently only supports calculation of FDRc<0.05; ignoring different user-supplied threshold');
						% TODO: This is a limitation of spm_getSPM that could easily be fixed.
						iBT.what.disp.threshold.thresh{current_threshold}(2) = 0.05;
					end
				 	stats_k_threshold = nxSPM.uc(4);
					if isinf(stats_k_threshold),
						disp(sprintf('k threshold required to achieve FDRc < %f (at a feature-selecting uncorrected threshold of p<%f) for contrast #%i is Inf',iBT.what.disp.threshold.thresh{current_threshold}(2),iBT.what.disp.threshold.thresh{current_threshold}(1),T_count));
						if (iBT.what.disp.threshold.write_t_value_file)
							indicator_filename=sprintf('kthreshFor_spmT_%04i_p=%g_t=%g_FDRc=0.05_is_Inf',T_count,iBT.what.disp.threshold.thresh{current_threshold}(1),stats_threshold);
				    			fclose( fopen( indicator_filename, 'W' ) );
						end
					else
						disp(sprintf('k threshold required to achieve FDRc < %f (at a feature-selecting uncorrected threshold of p<%f) for contrast #%i is %i',iBT.what.disp.threshold.thresh{current_threshold}(2),iBT.what.disp.threshold.thresh{current_threshold}(1),T_count,stats_k_threshold));
						new.k=stats_k_threshold;
						[nSPM,nxSPM] = spm_getSPM(new); %Obtain an image thresholded in this way
						activation = zeros(nxSPM.DIM');
						ind = sub2ind(nxSPM.DIM', nxSPM.XYZ(1,:), nxSPM.XYZ(2,:), nxSPM.XYZ(3,:));
						activation(ind) = nxSPM.Z;
						%  Write out the volumes:
						Vout = nSPM.xCon(T_count).Vspm;
						Vout.fname = [ sprintf('kthresh_p=%g_t=%g_FDRc=0.05_k=%i_' ,iBT.what.disp.threshold.thresh{current_threshold}(1),stats_threshold,stats_k_threshold), Vout.fname];
						spm_write_vol(Vout, activation);
	                        		mapHeader{sess}.imgs{2} = fullfile(spmT.loc, Vout.fname);
						mapHeader{sess}.imgs{3} = mapHeader{sess}.imgs{2};
						mapHeader{sess}.himgs{2} = Vout.fname;
						mapHeader{sess}.himgs{3} = Vout.fname;
					end
				else 
					stats_k_threshold = nxSPM.k; % should be 0 unless we set new.k to something else.
				end									  
		   	  end %if strcmpi(iBT.what.disp.threshold.type{current_threshold},'li') else ....		
			
			  if (~isnan(stats_threshold)) && (~isinf(stats_k_threshold)) 
			    if stats_k_threshold == 0 % (if k ~= 0 we have already determined the name)
			    	mapHeader{sess}.imgs{2} = spmT.loc;
			    	mapHeader{sess}.imgs{3} = mapHeader{sess}.imgs{2};
			    end
			    
			    if iBT.what.disp.MELODIC.SOCK
			        if (strcmpi(iBT.what.disp.contrast, 'auto') == 0) && (strcmpi(iBT.what.disp.contrast,'-effect') == 0)
			    		mapHeader{sess}.himgs{2} =  strcat('thresh_zstat', num2str(observe_ICs(T_count)), iBT.what.ext.FSL);
			    		mapHeader{sess}.imgs{2} = fullfile(mapHeader{sess}.imgs{2}, strcat('thresh_zstat', num2str(observe_ICs(T_count)), iBT.what.ext.FSL));
			    		mapHeader{sess}.himgs{3} = mapHeader{sess}.himgs{2};
			    		mapHeader{sess}.imgs{3} = mapHeader{sess}.imgs{2};
			        else
                        		mapHeader{sess}.himgs{2} = strcat('thresh_zstat', num2str(observe_ICs(T_count)), iBT.what.ext.FSL);
	                        	mapHeader{sess}.imgs{2} = fullfile(mapHeader{sess}.imgs{2}, strcat('thresh_zstat', num2str(observe_ICs(T_count)), iBT.what.ext.FSL));
	        	                mapHeader{sess}.himgs{3} = strcat('thresh_zstat', num2str(observe_ICs(T_count)), iBT.what.ext.FSL);
        	        	        mapHeader{sess}.imgs{3} = fullfile(mapHeader{sess}.imgs{3}, strcat('thresh_zstat', num2str(observe_ICs(T_count)), iBT.what.ext.FSL));
			        end
			    elseif stats_k_threshold == 0 % (if k ~= 0 we have already determined the name)
			        if (strcmpi(iBT.what.disp.contrast, 'auto') == 0) && (strcmpi(iBT.what.disp.contrast,'-effect') == 0)
			        	% We have already fully qualified the name
			    		[pathstr,name,ext]=fileparts(spmm.SPM.xCon(T_count).Vspm.fname);
			    		mapHeader{sess}.imgs{2} = pathstr;
		                	mapHeader{sess}.imgs{3} = pathstr;
                        		mapHeader{sess}.himgs{2} = name;
	                        	mapHeader{sess}.imgs{2} = spmm.SPM.xCon(T_count).Vspm.fname;
	        	                mapHeader{sess}.himgs{3} = name;
        	        	        mapHeader{sess}.imgs{3} = spmm.SPM.xCon(T_count).Vspm.fname;
			        else
                        		mapHeader{sess}.himgs{2} = spmm.SPM.xCon(T_count).Vspm.fname;
	                        	mapHeader{sess}.imgs{2} = fullfile(mapHeader{sess}.imgs{2}, spmm.SPM.xCon(T_count).Vspm.fname);
	        	                mapHeader{sess}.himgs{3} = spmm.SPM.xCon(T_count).Vspm.fname;
        	        	        mapHeader{sess}.imgs{3} = fullfile(mapHeader{sess}.imgs{3}, spmm.SPM.xCon(T_count).Vspm.fname);
			        end
			    end %if iBT.what.disp.MELODIC.SOCK else

			    mrej_header=''; %Initialise
			    if numel(contrast_columns) > 0
			    % We know which design matrix columns we are using, so we can display motion rejection information if we have it.
			    	try, mrej;
					num_cons = numel(contrast_columns);
					for col_count = 1:num_cons
						%if mrej{contrast_columns(col_count)}.num_rejected > 0
							N_active = numel(find(mrej{contrast_columns(col_count)}.active));
							N_inactive = numel(find(mrej{contrast_columns(col_count)}.inactive));
							num_active_rejected = mrej{contrast_columns(col_count)}.num_active_rejected;
							num_inactive_rejected = mrej{contrast_columns(col_count)}.num_inactive_rejected;
							pc_active = round( 100 * (N_active-num_active_rejected) / N_active );
							pc_inactive=round( 100 * (N_inactive-num_inactive_rejected) / N_inactive ); 
							mrej_header = [mrej_header sprintf('Col%i:%i=%i%% of tasks,%i=%i%% of rests. ',contrast_columns(col_count),N_active-num_active_rejected,pc_active,N_inactive-num_inactive_rejected,pc_inactive)];
						%end
					end
				catch;
				end
			    end			    
			    
			    % Loop over each desired combination of postive and/or negative map requested...
			    %===============================================================================
			    done_pos_only = 0 ; % initialise
			    for config_index = 1:1:config_num

			     if strcmpi(iBT.what.disp.threshold.type{current_threshold},'FDRc') || strcmpi(iBT.what.disp.threshold.type{current_threshold},'FWEc') 
			      	do_neg = 0; % Cluster-corrected spmTs are presently one-tailed only so there will be no negative values present.				
			     else
			      	do_neg = disp_config{config_index}.negative;
			     end
			     if (do_neg == 0) && (done_pos_only == 1)
			     	do_pos = 0; % Do not re-do something that we have already done.
			     else
			     	do_pos = disp_config{config_index}.positive;
			     end
			     
			     if ( do_pos ~= 0 )	|| ( do_neg ~= 0 )
			      if do_neg == 0,done_pos_only = 1; end % flag that we will have done a positive-only map in this loop	
			      positive_image_index = 2; % unless we change it below
		       	      if (disp_config{config_index}.invert == 1), 
			      	invert = 1; % Show positive values in cool colours and negative in warm colours
		       	      elseif (disp_config{config_index}.invert == 2), 
			      	if display_this < 0 
					invert = 1; % Only invert a "negative" contrast so its positive values are shown in cool colours and negative in warm colours
				else
					invert = 0;
				end
			      else 
			      	invert = 0; 
			      end 
 
		       	      if do_pos && do_neg
			    	images = strvcat(mapHeader{sess}.imgs{1}, mapHeader{sess}.imgs{2}, mapHeader{sess}.imgs{3});
			      elseif do_pos
		       			images = strvcat(mapHeader{sess}.imgs{1}, mapHeader{sess}.imgs{2});
			      elseif do_neg
		       			images = strvcat(mapHeader{sess}.imgs{1}, mapHeader{sess}.imgs{2});
			    		positive_image_index = -1; % There isn't a positive image.
			      else
		       			images = strvcat(mapHeader{sess}.imgs{1});
			      end % if do_pos && do_neg

			      if ~strcmp(iBT.what.disp.ROI.img,'')
			      		ROI_image_index = 1 + (do_pos~=0) + (do_neg~=0) + 1;
					mapHeader{sess}.imgs{ROI_image_index} = fullfile(iBT.what.disp.ROI.img); 
					mapHeader{sess}.himgs{ROI_image_index} = fullfile(iBT.what.disp.ROI.img);
					images = strvcat(images,mapHeader{sess}.imgs{ROI_image_index});
			      else
			      		ROI_image_index = -1;  % There isn't a ROI image.	
			      end
				
		              images = cellstr(images);
			    

			      % Now loop over each desired orientation...
			      %==========================================
			    
			      for current_transform = 1:1:length(iBT.what.disp.orient.transform)
			    
				xtransform = iBT.what.disp.orient.transform{current_transform};
				xslices = iBT.what.disp.orient.slices{current_transform};

				if iBT.what.disp.analysis == 5 % MELODIC
				  if iBT.what.disp.MELODIC.SOCK
				    disp(sprintf('Generating %s activation map at threshold %f for %s...',xtransform, stats_threshold, strcat('thresh_zstat', num2str(observe_ICs(T_count)), iBT.what.ext.FSL)))
				  else
				    disp(sprintf('Generating %s activation map at threshold %f for thresh_zstat #%d ...',xtransform, stats_threshold, T_count))
				  end
				else % Not MELODIC
				    disp(sprintf('Generating %s activation map at threshold %f using offset of %s from contrast #%d ...',xtransform, stats_threshold, num2str(iBT.what.stats.offset), T_count))
	                	end

				%% display_slices code begins here
				clear global SO;
				global SO;

				spm_input('!SetNextPos', 1);

                if (strcmpi(iBT.what.SpmVersion, 'SPM2')) 
				    original_flip = spm_flip_analyze_images; % preserve original flip value
				    flip = original_flip;
					% Note: The following attempts to maintain backwards compatibility but has not been tested since revflip deprecated.
					if strcmpi(iBT.what.disp.orient.LeftOnRight,'yes')  % equivilant to old revflip
					        flip = flip -1;
					end
				end %SPM2

				if ~strcmpi(xtransform,'sagittal')
 					if ( strcmpi(iBT.what.disp.orient.LeftOnRight,'yes') )
 						mapHeader{sess}.format = 'radiological'; % set for later display as title
					else
						mapHeader{sess}.format = 'neurological';
					end % if LeftOnRight,'yes' else
				else % must be sagittal:
					if reorder_sagittal == 1,
						xslices = xslices(end:-1:1); %reverse order displayed
						disp(['Will reverse display order of sagittal slices because legacy runinfo detected (iBT.what.disp.revflip is defined).']);

					end
				end % if ~strcmpi(xtransform,

				%
				% Note: The following label logic depends on the input images following MNI sign convention. If not, 
				%       then the user should consider disabling orientation text as it may be incorrect:
				if iBT.what.disp.orient.showText == 1,
					if strcmpi(xtransform,'axial')
						if ( (xslices(2) - xslices(1)) > 0 ) %if step is positive then most negative (inferior according to MNI) slice is first displayed (top)
							mapHeader{sess}.format = [mapHeader{sess}.format ', inferior 1st top']; % set for later display as title
						else
							mapHeader{sess}.format = [mapHeader{sess}.format ', superior 1st top']; % set for later display as title
						end
					elseif strcmpi(xtransform,'coronal')
						if ( (xslices(2) - xslices(1)) > 0 ) %if step is positive then most negative (posterior according to MNI) slice is first displayed (top)
							mapHeader{sess}.format = [mapHeader{sess}.format ', posterior 1st top']; % set for later display as title
						else
							mapHeader{sess}.format = [mapHeader{sess}.format ', anterior 1st top']; % set for later display as title
						end
					elseif strcmpi(xtransform,'sagittal')
						if ( (xslices(2) - xslices(1)) > 0 ) %if step is positive then most negative (Left according to MNI) slice is first displayed (top)
							mapHeader{sess}.format = 'left 1st top'; % set for later display as title
						else
							mapHeader{sess}.format = 'right 1st top'; % set for later display as title
						end
					else % perhaps tilted, which we do not support with an orientation description
							mapHeader{sess}.format = ''; % we do not know the orientation
					end
				else
							mapHeader{sess}.format = ''; % we do not know the orientation 
				end % if iBT.what.disp.orient.showText, else

				if ~iBT.what.disp.fix_for_SPM8r4667_applied
				  if ( strcmpi(iBT.what.SpmVersion, 'SPM8') && ((iBT.what.SpmRelease >= 4341) && (iBT.what.SpmRelease < 4669)) && strcmpi(xtransform,'sagittal') )
				  ; % in the SPM8 public release 4667 slover was broken for sagittal 
				  ; % by an innapropriate change to fill_defaults 
				  ; % (the change was actually made in private release 4341) 
				  ; % As a work-around, we revert to the previous version. 
				  ; % I submitted a bug report a few days after r4667 was released, 
				  ; % and as a result the bug was fixed in private r4669. The first
				  ; % public release containing the fix was r5236.
					iBT.what.disp.fix_for_SPM8r4667_applied = 1; %So we know that we do not need to do this again this session
					[pathstr,name,ext]=fileparts(which('iBT_display'));
					newpath = fullfile(pathstr,'fix_for_SPM8r4667');
  					addpath(newpath);
  					disp(['Added ' newpath ]);
  					disp(['  to Matlab path to avoid bug introduced in SPM8_r4667.']);
				  end
				end %if ~iBT.what.disp.fix_for_SPM8r4667_applied
				
                                if (strcmpi(iBT.what.SpmVersion, 'SPM2')) 
				  global defaults;  % Give us access to the SPM defaults variable.
				end
 
				% load images
				nimgs = size(images,1); 

				% process names
				nchars = 20;
				imgns = spm_str_manip(images, ['rck' num2str(nchars)]);
				% identify image types
				cscale = [];
				deftype = 1;
				SO.cbar = [];
				for i = 1:nimgs
                                        if (strcmpi(iBT.what.SpmVersion, 'SPM2')) 
					        defaults.analyze.flip = flip; % flip (if selected, otherwise has no effect)
                                        end

					SO.img(i).vol = spm_vol(images{i});
					SO.img(i).hold = iBT.what.disp.resample;

                                        if (strcmpi(iBT.what.SpmVersion, 'SPM2')) 
					        defaults.analyze.flip = original_flip; % set global back to original value as soon as possible to avoid ctrl+c problems
                                        end

					if i == 1 % set the appropriate image type
					    itype{1} = 'Structural';
					elseif i == ROI_image_index 
					    itype{1} = 'ROI';
					else
					    itype{1} = 'Blobs';
					end

					imgns(i) = {sprintf('Img %d (%s)',i,itype{1})};

                                        if (strcmpi(iBT.what.SpmVersion, 'SPM2')) 
					        [mx mn] = iBT_slice_overlay_spm2('volmaxmin', SO.img(i).vol);
                                        else
                                                [mx mn] = slover('volmaxmin', SO.img(i).vol);
                                        end

					if ~isempty(strmatch('Structural', itype))
						try,
							if (mn < iBT.what.disp.threshold.structural.low); mn = iBT.what.disp.threshold.structural.low; end
						catch; end
						try,
							if (mx > iBT.what.disp.threshold.structural.high); mx = iBT.what.disp.threshold.structural.high; end
						catch; end
						try,
							mn = iBT.what.disp.threshold.structural.min; 
						catch; end
						try,
							mx = iBT.what.disp.threshold.structural.max; 
						catch; end

			    			SO.img(i).cmap = gray;
						SO.img(i).range = [mn mx];
						deftype = 2;
    						cscale = [cscale i];
	    					if strcmpi(itype,'Structural with SPM blobs')
                                                        if (strcmpi(iBT.what.SpmVersion, 'SPM2')) 
							        iBT_slice_overlay_spm2('addspm',[],0);
                                                        else
                                                                slover('addspm',[],0);
                                                        end
						end
					else
			    			if iBT.what.disp.show.colourbar && isempty(strmatch('ROI', itype))
						  SO.cbar = [SO.cbar i];
						  cprompt = ['Colormap: ' imgns{i}];
						end

	    					switch itype{1}
					     		case 'Truecolour'
      								dcmap = 'actc';
      								drange = [mn mx];
	      							cscale = [cscale i];
					     		case 'ROI'
      								drange = [0 1]; % Range for the ROI image.
								new_cmap(1,:) = [0 0 0];
								new_cmap(2,:) = iBT.what.disp.ROI.colour;
								SO.img(i).cmap = new_cmap;
								SO.img(i).range = drange;
								SO.img(i).type='truecolour';
      								SO.img(i).prop = 1.0 - iBT.what.disp.ROI.transparency;
								ROI_func = 'i1(i1~=0)=1;';;
								SO.img(i).func =  ROI_func;
	     						case 'Blobs'
								drange = [xmini xmaxi];
      								if i == positive_image_index && ~invert
									dcmap = 'hot'; %  Display positive in hot scale
									SO.img(i).cmap = hot;
									SO.img(i).range = drange;
									if iBT.what.disp.transparency > 0
									  SO.img(i).type='truecolour';
									  transparent_i = i;
									end
									if isfinite(stats_threshold) && (stats_threshold ~= 0) && (stats_threshold ~= iBT.what.disp.threshold.mini{current_threshold})
									  thresh_func = sprintf('i1(i1<%f)=0;',stats_threshold);
									  SO.img(i).func =  thresh_func;
									end
      								elseif i == positive_image_index && invert
	 								dcmap = 'winter'; %  Display positive in winter scale
									SO.img(i).cmap = winter;
									SO.img(i).range = drange;
									if iBT.what.disp.transparency > 0
									  SO.img(i).type='truecolour';
									  transparent_i = i;
									end
									if isfinite(stats_threshold) && (stats_threshold ~= 0) && (stats_threshold ~= iBT.what.disp.threshold.mini{current_threshold})
									  thresh_func = sprintf('i1(i1<%f)=0;',stats_threshold);
									  SO.img(i).func =  thresh_func;
									end
      								elseif ~invert
	 								dcmap = 'winter'; %  Display negative in winter scale
									SO.img(i).cmap = winter;
									SO.img(i).range = -drange;
									if iBT.what.disp.transparency > 0
									  SO.img(i).type='truecolour';
									  transparent_i = i;
									end
									if isfinite(stats_threshold) && (stats_threshold ~= 0) && (stats_threshold ~= iBT.what.disp.threshold.mini{current_threshold})
									  thresh_func = sprintf('i1(i1>-%f)=0;',stats_threshold);
									  SO.img(i).func =  thresh_func;
									end
      								else
									dcmap = 'hot'; %  Display negative in hot scale
									SO.img(i).cmap = hot;
									SO.img(i).range = -drange;
									if iBT.what.disp.transparency > 0
									  SO.img(i).type='truecolour';
									  transparent_i = i;
									end
									if isfinite(stats_threshold) && (stats_threshold ~= 0) && (stats_threshold ~= iBT.what.disp.threshold.mini{current_threshold})
									  thresh_func = sprintf('i1(i1>-%f)=0;',stats_threshold);
									  SO.img(i).func =  thresh_func;
									end
      								end
      								SO.img(i).prop = Inf;
     							case 'Negative blobs'
      								dcmap = 'winter';
	      							drange = [0 mx];
      								SO.img(i).prop = Inf;
	    					end
  					end
				end

				ncmaps=length(cscale);
				if ncmaps == 1
	  				SO.img(cscale).prop = 1;
				else
  					remcol=1;
			  		for i = 1:ncmaps
					    	ino = cscale(i);
    						SO.img(ino).prop = spm_input(sprintf('%s intensity',imgns{ino}),...
								'+1', 'e', ...
								remcol/(ncmaps-i+1),1);
	    					remcol = remcol - SO.img(ino).prop;
				 	end
				end
				if (iBT.what.disp.transparency > 0) && (~isempty(strmatch('Blobs', itype)))
					% notwithstanding above, we can force a particular transparency here
					SO.img(transparent_i).prop = 1.0-iBT.what.disp.transparency;
				end

				SO.transform = xtransform;
				mapHeader{sess}.view = SO.transform;

				if strcmpi(SO.transform, 'tilted')
					%  transform the image
					angle = spm_input('Slice angle (degrees)', '+1', 'r', '0', NaN, [-180 180]) ;
					transform = [0 0 0 (angle/180)*pi 0 0 1 1 1 0 0 0];
					SO.transform = spm_matrix(transform);
				end

				% Construct transformation matrix here for a sagittal variation not directly
				% supported by slover, so when slover is called it can deal with slice positioning optimally.
				if ( strcmpi(xtransform,'sagittal') ~= 0 )
					if strcmpi(iBT.what.disp.orient.AnteriorOnRight,'yes')
						transform = [0 0 0 pi/2 0 pi/2 1 1 1]; % Rotate in opposite direction to default behaviour and do not flip; yields AnteriorOnRight, retaining MNI-convention signs
					else
						transform = [0 0 0 pi/2 0 -pi/2 -1 1 1]; % This is what slover normally does and is typical Anterior towards left behaviour (see fill_defaults)
					end
					SO.transform = spm_matrix(transform);
				end

				% use SPM figure window
				SO.figure = spm_figure('GetWin', 'Graphics');
				set(SO.figure, 'Color', [0 0 0])

                                if (strcmpi(iBT.what.SpmVersion, 'SPM2'))
				        % slices for display
        				iBT_slice_overlay_spm2('checkso');

        				SO.slices = xslices;
        				SO.labels.colour = 'white';

        				% and do the display
        				iBT_slice_overlay_spm2

        				temp = defaults.printstr;
        				defaults.printstr =  [spm_figure('DefPrintCmd'),SO.printfile]; 
                                else
                                        SO.slices = xslices;
                                        SO.labels.colour = 'white';
                                        SO = slover(SO);
                                        SO.figure_struct.Position = [518 15 604 871];
                                        SO.figure_struct.Units = 'pixels';

                                         % flip if necessary
										if ( ( strcmpi(xtransform,'sagittal') == 0 ) && strcmpi(iBT.what.disp.orient.LeftOnRight,'yes') )
 									   			% Display left of brain on right of image (Radiological convention; opposite to SPM's Neurologocal convention)
                                            	SO.transform(1,1:3) = -SO.transform(1,1:3); % Flip the horizonally displayed axis;
											end
                                        % and do the display
										SO = paint(SO);
                                end

                                F=spm_figure('FindWin','Graphics');

				% Create the various titles and footnotes
				HNote = '';
				if iBT.what.disp.show.name
					HNote = strcat(HNote, sprintf('%s', mapHeader{sess}.name));				
				end
				if iBT.what.disp.show.ID
					HNote = strcat(HNote, sprintf('%s', mapHeader{sess}.ID));
				end
				if iBT.what.disp.show.date
					HNote = strcat(HNote, sprintf('%s', mapHeader{sess}.sdate));
				end
				if ~strcmp(HNote,'')
					HNote = strcat(HNote, sprintf('; '));
				end
				if iBT.what.disp.show.task
					HNote = strcat(HNote, sprintf('%s', mapHeader{sess}.task));
				end
				if iBT.what.disp.show.conName
					HNote = strcat(HNote, sprintf('%s', spmm.SPM.xCon(T_count).name));
				end
				mapHeader{sess}.thresh = num2str(stats_threshold);
				if iBT.what.disp.show.thresh
					if strcmpi(iBT.what.disp.threshold.type{current_threshold},'p') || strcmpi(iBT.what.disp.threshold.type{current_threshold},'FWE')
					  HNote = strcat(HNote, sprintf( ' %s; ', nxSPM.thresDesc));
					end
					if ( (stats_threshold == 0) && (iBT.what.disp.threshold.mini{current_threshold} > 0)) G_or_GE = '>'; else G_or_GE = '>='; end
					if do_neg
					  HNote = strcat(HNote, sprintf( ' |%s|%s%s', iBT.what.disp.statistic, G_or_GE, mapHeader{sess}.thresh));
					else
					  HNote = strcat(HNote, sprintf( ' %s%s%s',   iBT.what.disp.statistic, G_or_GE, mapHeader{sess}.thresh));
					end
				end
				HNote1a = sprintf('%s', mapHeader{sess}.view);
				% Capitalise the first letter:
				if (HNote1a(1) >= 'a'), HNote1a(1)=char(HNote1a(1)-'a'+'A'); end
				HNote1a = strcat(HNote1a, sprintf(' %s', mapHeader{sess}.format)); 
				if iBT.what.disp.show.voxelSize == 1, % Display original scanner voxel sizes (before any spatial normalisation)
					HNote2 = strcat(sprintf('Scanner voxel size: %s', mapHeader{sess}.pixelx), sprintf(' x %s', mapHeader{sess}.pixely), sprintf(' x (%s +', mapHeader{sess}.slicethickness), sprintf(' %s). ', mapHeader{sess}.spacing)); 	% setup the 2nd title line
				else
					HNote2 = ['']; % setup the 2nd title line
				end % if iBT.what.disp.hideVoxelSize
				if iBT.what.disp.show.procVoxelSize == 1; % Display processed volume voxel sizes
					Vol=spm_vol(mapHeader{sess}.imgs{2}); % mapHeader{sess}.imgs{2} is the first spmT (mapHeader{sess}.imgs{1} is the anatomy)
					voxel_size = sqrt(diag(Vol.mat(1:3, 1:3)).^2); % Voxel sizes in mm are on the diagonal of the affine transformation from voxels to "real world" coordinates. Calculate absolute in case negative.
					HNote2 = [HNote2, sprintf('Processed voxel size: %g x %g x %g. ', voxel_size(1), voxel_size(2), voxel_size(3))]; 
				end % iBT.what.disp.show.procVoxelSize
				HNote2a = '';
				if iBT.what.disp.show.warp
				  if strcmpi(iBT.what.disp.type,'norm') == 1
				    if ~strcmp(mapHeader{sess}.n_template,'')
					  if strcmpi(mapHeader{sess}.n_constraint,'affine'),
						HNote2a = [HNote2a, sprintf('Affine transform to%s. ', mapHeader{sess}.n_template)]; 
					  elseif strcmpi(mapHeader{sess}.n_constraint,'rigid'),
						HNote2a = [HNote2a, sprintf('Realigned (~ rigid body) to%s. ', mapHeader{sess}.n_template)]; 
					  elseif strcmpi(mapHeader{sess}.n_constraint,'warp'),
						HNote2a = [HNote2a, sprintf('Warped to%s. ', mapHeader{sess}.n_template)]; 
					  else
						HNote2a = [HNote2a, sprintf('Unknown transform to%s. ', mapHeader{sess}.n_template)]; 
			  		  end
				    end % if ~ ''
				  elseif strcmpi(iBT.what.disp.type,'own') == 1
						HNote2a = [HNote2a, sprintf('Own space. ')]; 
				  end
				end % iBT.what.disp.show.warp 

				if iBT.what.disp.show.smooth && (~strcmp(mapHeader{sess}.smooth,'?')),
				 HNote2a = [HNote2a, sprintf('Smoothing: %smm, ',mapHeader{sess}.smooth), sprintf('TR: %s', mapHeader{sess}.TR)];
				end
				
				% setup the 3rd title line (HNote3)
				if isfield(iBT.what.disp,'structuralLabel'),
				    if ~strcmp(iBT.what.disp.structuralLabel,''),
					HNote3 = strcat(sprintf('Anatomy: %s', iBT.what.disp.structuralLabel)); 
				    else
					HNote3 = '';
				    end % if strcmp(iBT.what.disp.structuralLabel,'');			
				elseif iBT.what.disp.hideName == 0,		
					HNote3 = strcat(sprintf('Anatomy: %s', mapHeader{sess}.himgs{1}));
				else
					HNote3 = strcat(sprintf('Anatomy: %s', NAME_structural));
				end

				% FNote = strcat(sprintf('Activation & Deactivation contrast: %s', mapHeader{sess}.himgs{2}));
				switch 1
				  case iBT.what.disp.resample == 0, 
				  	resampling = 'nearest neighbour resampling';
				  case iBT.what.disp.resample == 1, 
				  	resampling = 'trilinear resampling';
				  case iBT.what.disp.resample > 1, 
				  	resampling = sprintf('Lagrange (polynomial) resampling, order %i',iBT.what.disp.resample);
				  case iBT.what.disp.resample < 1, 
				  	resampling = sprintf('sinc resampling, order %i',-iBT.what.disp.resample);
				end
				
				FNote = strcat(sprintf('Overlay: %s (%s)',mapHeader{sess}.himgs{2},resampling));
				if strcmpi(iBT.what.disp.threshold.type{current_threshold},'li')
					FNote = [li_string '.  ' FNote];
				end
				axes('Position',[0.005,0.005,0.1,0.1],'Visible','off','Tag','SPMprintTitle')
				set(gcf,'DefaultTextColor','white');

				% Header display
				vposition = 9.82;
				line_space = 0.18;
				text(0,vposition,HNote,'FontSize',iBT.what.disp.fontsize1);
				vposition = vposition - (line_space * 1.2); % 1.2 as a bit more space for the first linespace is required (otherwise underscores on first line hit the second line text)
				text(0.1,vposition,[HNote1a '. ' HNote2],'FontSize',iBT.what.disp.fontsize2);				
				vposition = vposition - line_space;
				if iBT.what.disp.show.smooth && (~strcmp(mapHeader{sess}.smooth,'?')),
				  text(0.1,vposition,[HNote2a ', ' HNote3],'FontSize',iBT.what.disp.fontsize2);
				  vposition = vposition - line_space;
				else
				  text(0.1,vposition,HNote3,'FontSize',iBT.what.disp.fontsize2);
				  vposition = vposition - line_space;
				end
				
%				if strcmpi(iBT.what.disp.type,'norm') == 1
%					% Footer display if normalised analysis (no room on ownspace)
%					text(0,line_space,FNote,'FontSize',10);
%				end
				text(0.1,vposition,FNote,'FontSize',iBT.what.disp.fontsize2);
				vposition = vposition - line_space;

				if iBT.what.disp.show.mrej,
					text(max([10.0 - (((1+iBT.what.disp.fontsize2)/4)*numel(mrej_header)/28.0), 0.1]),line_space,mrej_header,'FontSize',iBT.what.disp.fontsize2);
				end

				if iBT.what.disp.show.LR,
					% Put an L or R label near the left of the axial and coronal images to indicate orientation:
					% We put it at half the current vposition, near the left margin of the page:
					if strcmpi(xtransform,'axial') || strcmpi(xtransform,'coronal'),
					    if strcmpi(iBT.what.disp.orient.LeftOnRight,'yes'),
						% For fontize = 18, vertical spacing of 0.22 works
						% 	and for RIGHT, xpos = 0.03,0.08,0.02,.0.2 and 0.025 centres the letters
						% For fontize = 11 (our default), vertical spacing of 0.14 works
						vspace = 0.14;
						%vstart = 2*vposition/3;
						vstart = vposition - vspace;
						text(0.005 ,vstart+(vspace*0.4),char(8592),'FontSize',iBT.what.disp.fontsizeLR);
						text(0.03 ,vstart-(vspace*1),'R','FontSize',iBT.what.disp.fontsizeLR);
						text(0.08 ,vstart-(vspace*2),'I','FontSize',iBT.what.disp.fontsizeLR);
						text(0.02 ,vstart-(vspace*3),'G','FontSize',iBT.what.disp.fontsizeLR);
						text(0.02 ,vstart-(vspace*4),'H','FontSize',iBT.what.disp.fontsizeLR);
						text(0.025,vstart-(vspace*5),'T','FontSize',iBT.what.disp.fontsizeLR);
						% text(0.005 ,vstart-(vspace*6.2),char(8592),'FontSize',iBT.what.disp.fontsizeLR);
					    elseif strcmpi(iBT.what.disp.orient.LeftOnRight,'no'),
						vspace = 0.14;
						vstart = vposition - vspace;
						text(0.005 ,vstart+(vspace*0.4),char(8592),'FontSize',iBT.what.disp.fontsizeLR);
						text(0.03 ,vstart-(vspace*1),'L','FontSize',iBT.what.disp.fontsizeLR);
						text(0.02 ,vstart-(vspace*2),'E','FontSize',iBT.what.disp.fontsizeLR);
						text(0.02 ,vstart-(vspace*3),'F','FontSize',iBT.what.disp.fontsizeLR);
						text(0.025,vstart-(vspace*4),'T','FontSize',iBT.what.disp.fontsizeLR);
						% text(0.005 ,vstart-(vspace*5.2),char(8592),'FontSize',iBT.what.disp.fontsizeLR);					    
					    end
					end
				end

				tmp = strfind( mapHeader{sess}.himgs{2}, iBT.what.ext.SPM ); 
				if length(tmp) > 0 % Remove trailing file extension .img if present:
					contrastname= mapHeader{sess}.himgs{2}(1:tmp(length(tmp))-1);
				else
					contrastname= mapHeader{sess}.himgs{2};
				end
				if strcmpi(iBT.what.disp.threshold.type{current_threshold},'p') || strcmpi(iBT.what.disp.threshold.type{current_threshold},'FWE')
					thresDesc = fixname([nxSPM.thresDesc]);
				elseif strcmpi(iBT.what.disp.threshold.type{current_threshold},'li')
					thresDesc = 'li';	 
				else 
					thresDesc = '';
				end
 
				if iBT.what.disp.pngname.include.contrastnum
				 pngname = sprintf('slices_%s_%s_%.2i_%s_%s%s%s_%s%s',mapHeader{sess}.view, iBT.who.sub{sub}.ID,sess,mapHeader{sess}.task,thresDesc,iBT.what.disp.statistic,mapHeader{sess}.thresh,contrastname,SOCK_string);
				else
				 pngname = sprintf('slices_%s_%s_%s_%s%s%s_%s%s',mapHeader{sess}.view, iBT.who.sub{sub}.ID,mapHeader{sess}.task,thresDesc,iBT.what.disp.statistic,mapHeader{sess}.thresh,contrastname,SOCK_string);
				end

				if (put_sign_in_filename == 2) || (put_sign_in_filename && ~( do_pos && do_neg && (invert == 0) )) % || is for legacy runinfos where we maintain a lack of +-
					if do_pos && do_neg,
						 pngname = [pngname '+-']; 
					elseif do_neg, pngname = [pngname '_-'];
					elseif do_pos, pngname = [pngname '_+'];
					end
				end
				










				pngname = fullfile(output_folder,[pngname '.png']);
				disp(['Saving to ' pngname]);
				set(F,'PaperPositionMode','auto');
				print(F,'-dpng',pngname); % save to png format
				set(F,'PaperPositionMode','manual');
				
				tmp_pwd = pwd;
				cd(output_folder);
				if iBT.what.disp.postscript > 0
                                        if (strcmpi(iBT.what.SpmVersion, 'SPM2')) 
	        	        	        spm_print; % save to postscript file
                                        else
                                                slover('printfile'); % save to postscript file
                                        end
				end
				cd(tmp_pwd);

                                if (strcmpi(iBT.what.SpmVersion, 'SPM2')) 
				        defaults.printstr =  temp; % return to previous flip settings
        				defaults.analyze.flip = original_flip;
                                end
			      end % for current_transform = 1:1:length(iBT.what.disp.orient.transform)
			     end % if ( do_pos ~= 0 )	|| ( do_neg ~= 0 )	
			    end % for config_num
			  end % if ~isnan(stats_threshold)
			end % if display_this
		end % for T_count
	end %for current_threshold
        cd(awd)
%end % for current_transform = 1:1:length(iBT.what.disp.orient.transform); NOTE: This loop has now been moved up to avoid much unnecessary duplicate processing for each orientation.

% Reset colours to those in place on entry:
set(gcf,'DefaultTextColor',OriginalTextColor);
try,set(SO.figure, 'Color', OriginalBackgroundColour);catch;end % If SO does not exist then that is OK - means we did not change anything.

end %of function iBT_display(iBT)

function cmap = return_cmap(prompt,defmapn)
 cmap = [];
 while isempty(cmap)
        if (strcmpi(iBT.what.SpmVersion, 'SPM2')) 
                cmap = iBT_slice_overlay_spm2('getcmap', spm_input(prompt,'+1','s', defmapn));
        else
                cmap = slover('getcmap', spm_input(prompt,'+1','s', defmapn));
        end
 end
end %of function cmap

function varargout = iBT_slice_overlay_spm2(action, varargin);
	varargout = slice_overlay(action, varargin);
end %of function iBT_slice_overlay_spm2

function new_str = fixname(orig_str)
	new_str = strrep(orig_str,'<','');
	new_str = strrep(new_str,' ','');
end %of function fixname
