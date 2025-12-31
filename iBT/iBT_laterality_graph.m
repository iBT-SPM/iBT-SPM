function  thresh_detail = iBT_laterality_graph(iBT,con,...
               con_deactivation,li_group,li_individual,step,...
	       xaxis,maxX,xaxis_label,fileIDstring,fileBaseName,fignums)
%
% Assist with analysis and plotting of regional laterality
% over a range of thresholds
% This script is designed to be called by iBT_laterality.
%___________________________________________________________________________
% Copyright 2009-2015,2018,2020,2024 The Florey Institute 
%                               of Neuroscience and Mental Health
%
% This file is part of iBT (the Integrated Brain Analysis Toolbox for SPM).
% See iBT_laterality.m for more information.

% Please refer to the following paper for a description of the
% methods implemented in this toolbox. If you publish results using
% these methods or variations of them, please remember to cite this paper: 
%     Abbott DF, Waites AB, Lillywhite LM, Jackson GD.
%     fMRI assessment of language lateralization: An objective approach. 
%     Neuroimage 50(4):1446-1455 (2010).
%
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
% Known bugs: LI index value at max n for normal controls is not calculated or
%             displayed on summary plot unless deactivation
%             contrast is also analysed ( see con_deactivation > 0 )
%             This is an unnecessary limitation of the present code.
%
% Recent Version History:
% 2024-10-30: (dfa) Support iBT.what.li.lowNpercentValid - minimum % of 
%                   controls (default = 95) that must not have exhausted
%                   positive voxels in the ROI for N to be in the valid range.
%                   Bugfix: Write correct percentage to log 
%                          (previously said half when was actualy 95%)
% 2024-10-18: (dfa) Bugfix: Prevent crash if there are no control data 
%                           files inside the control data folder.
% 2024-05-28 (dfa) Bugfix: Don't try to plot P value calculation point when 
%                          there is not a valid range to calculate a P value
% 2024-05-26 (dfa) Bugfix: Support iBT.what.li.FullPlot.AllY when no controls
% 2020-06-08 (dfa) Limit P-value display to 2 significant figures.
%                  Support iBT.what.li.LineWidth.Default
%                  Support iBT.what.li.LineWidth.Individual
%                  Support iBT.what.li.SummaryPlot.LineWidth.Default
%                  Replaced *fullPlotAllY and *summaryPlotAllY
%                    with *.FullPlot.* & *.SummaryPlot.* versions
% 2020-05-24: (dfa) Halt and provide instructions for anyone who may have
%   uncommented and used the now retired fields:
%   iBT.what.li.GroupMeanColour, .GroupControlsColour and/or .IndividualColour
% 2020-04-30: (dfa) Support iBT.what.li.SummaryPlot.show100 to display LI %
% 2020-04-28: (dfa) Allow full RGB specification of shading colors via new
%                   iBT.what.li.RGB.* iBT.what.li.LineStyle.* arrays
%                   Support iBT.what.li.SummaryPlot.RGB.GroupMean and
%		      iBT.what.li.SummaryPlot.LineStyle.GroupMean
%                   Replaced legacy *.summary_plot_* variables with 
%		    *.SummaryPlot.* for better naming consistency.
% 2020-04-19: (dfa) Support iBT.what.li.SummaryPlot.[a]typicalString
% 2018-02-02: (dfa) Fixed sumamry plot x-axis label problem that 
%                   occurred in Matlab versions R2014b and later.
% 2015-10-29: (dfa) Replaced nanmean, nanmedian and nanstd with ibt_ functions 
%                   to regain control of which functions actually are called, 
%		    avoiding potentially buggy replacements in SPM or elsewhere.
% 2014-07-02: (dfa) More flexible options for defining valid summary region
%		    and associated summary LI value
%                 (iBT.what.li.lowNpercentOffCeiling, ceiling and jb_alpha )
% 2013-09-30: (dfa) More informative error message when too few contrasts
%		    in control group.
% 2013-06-01: (dfa) Added some extra plotting options.
%		Extend y axis a bit so axes don't cover-up data 
%		(bug: only works properly at y=-1 for summary plot)
%		Stop annoying window-focus-stealing by Matlab.
% 2012-03-12: (dfa) Display version info from iBT_version()
% 2011-11-11 (dfa) Add option for iBT.what.li.fullPlotAllY = 2 to do a half-range plot.
%		   and same for iBT.what.li.summaryPlotAllY = 2.
%		   Added what.li.save.control_group_stats_filename and documented a semi-automated
%                  method to compare two groups. Search for "to compare two groups" in this file.
% 2011-11-03: (dfa) Work around a bug in some older versions of Matlab that had jbtest of
%		a flat distribution incorrectly returning 0 or crashing.
% 2011-08-26: (dfa) Updated for public release in the 
%			Integrated Brain Analysis Toolbox for SPM
% 2010-11-22: (dfa) Ensure saved plots have tick marks & labels consistent with those displayed on screen.
% 2010-09-21: (dfa) Added what.li.conclusion.do (default to 1 if not specified)
% 2010-08-09: (dfa) Stop Matlab using 10 to power notation for number of voxel plot axes
% 2010-08-05: (dfa) Fixed crash if what.li.save.thresh is 1 and either act or deact thresh is not applicable.
% 2010-06-12: (dfa) Support what.li.save.thresh
% 2010-03-25: (dfa) Added more options to summary_plot_in_range_only.
%                   Added what.li.summary_plot_lower_range_fraction and _upper_range_fraction
% 2010-03-23: First public release in the iBrain Laterality Toolbox.
%

% 2004-2010: iBT analysis scripts were developed by David Abbott based
%            upon the who & what structure concept of Tony Waites, with 
%            substantial contributions from Matt Harvey, Kaushik Bhaganagarapu,
%            and Richard Masterton.
%            In 2009 the script that became iBT_laterality_graph
%            was created by David Abbott.
%

disp(sprintf('Entered Laterality Index Graphing script...'));

if isempty(iBT), error('Error: "iBT" structure is not defined.'), end

[iBTversion,iBTrelease,iBTversionReleaseString] = iBT_version(1); % Display version so it gets logged.

% These  values are good for saved plots
axis_label_fontsize=28;
axis_fontsize=20;
title_fontsize=12; % Smaller so have room to display everything

% These values are better with on-screen display in smaller windows:
axis_label_fontsize=14;
axis_fontsize=14;

% Set defaults for some parameters if not already set:
try iBT.what.li.errorbars; catch iBT.what.li.errorbars=0; end

% Allow full RGB specification for shading colours, whilst preserving backwards compatibility with prior monochrome "shading" options:
try default_CI95shading=iBT.what.li.CI95shading; catch default_CI95shading=0.85; end % iBT.what.li.CI95shading is deprecated: See iBT.what.li.RGB.CI95
try default_CI99shading=iBT.what.li.CI99shading; catch default_CI99shading=0.80; end % iBT.what.li.CI95shading is deprecated: See iBT.what.li.RGB.CI95
try default_InvalidShading = iBT.what.li.InvalidShading; catch default_InvalidShading=0.9; end % Deprecated: See iBT.what.li.RGB.Invalid
try iBT.what.li.RGB.CI95; catch iBT.what.li.RGB.CI95=[default_CI95shading default_CI95shading default_CI95shading]; end
try iBT.what.li.RGB.CI99; catch iBT.what.li.RGB.CI99=[default_CI99shading default_CI99shading default_CI99shading]; end
try iBT.what.li.RGB.Invalid; catch iBT.what.li.RGB.Invalid=[default_InvalidShading default_InvalidShading default_InvalidShading]; end %Shading for invalid region (non-normal or too few control subjects contribute)

try iBT.what.li.PlotGroupOnly = iBT.what.li.PlotGroupOnly; catch iBT.what.li.PlotGroupOnly=0; end
try iBT.what.li.PlotJB = iBT.what.li.PlotJB; catch iBT.what.li.PlotJB=0; end

% The following 3 items have been retired:
%   try iBT.what.li.GroupMeanColour = iBT.what.li.GroupMeanColour; catch iBT.what.li.GroupMeanColour='b'; end % OBSOLETED BY iBT.what.li.RGB.GroupMean & iBT.what.li.LineStyle.GroupMean
%   try iBT.what.li.GroupControlsColour = iBT.what.li.GroupControlsColour; catch iBT.what.li.GroupControlsColour='c--'; end % OBSOLETED BY iBT.what.li.RGB.GroupControls & iBT.what.li.LineStyle.GroupControls
%   try iBT.what.li.IndividualColour = iBT.what.li.IndividualColour; catch iBT.what.li.IndividualColour='r'; end % OBSOLETED BY iBT.what.li.RGB.Individual & iBT.what.li.LineStyle.Individual
% It is not worth providing backward compatibility, because these fields were always commented out in distributed runinfos 
% I suspect they were only ever used by the developer, but just in case, for anyone who might have uncommented and used
% these fields, we detect this, halt and provide instructions here for the more flexible replacements.
if (isfield(iBT.what.li,'GroupMeanColour')) | (isfield(iBT.what.li,'GroupControlsColour')) | (isfield(iBT.what.li,'IndividualColour')),
    disp('ERROR: The following fields in iBT.what.li have been retired:');
    disp('iBT.what.li.GroupMeanColour');
    disp('   - Replaced by iBT.what.li.RGB.GroupMean     & iBT.what.li.LineStyle.GroupMean');
    disp('iBT.what.li.GroupControlsColour');
    disp('   - Replaced by iBT.what.li.RGB.GroupControls & iBT.what.li.LineStyle.GroupControls');
    disp('iBT.what.li.IndividualColour');
    disp('   - Replaced by iBT.what.li.RGB.Individual    & iBT.what.li.LineStyle.Individual');
    error('iBT will now end. Please update your runinfo file and try again.');	    
end 

try iBT.what.li.RGB.GroupMean; catch iBT.what.li.RGB.GroupMean=[0 1 1]; end        % Colour of the control group mean on the detailed plot [0 1 1] is cyan
try iBT.what.li.LineStyle.GroupMean; catch iBT.what.li.LineStyle.GroupMean='-'; end % Style of the control group mean on the detailed plot 
try iBT.what.li.RGB.GroupControls; catch iBT.what.li.RGB.GroupControls=[0.4 1 1]; end        % Colour of each control subject on the detailed plot ( [0.4 1 1] is light cyan)
try iBT.what.li.LineStyle.GroupControls; catch iBT.what.li.LineStyle.GroupControls='--'; end % Style of each control subject on the detailed plot 
try iBT.what.li.RGB.Individual; catch iBT.what.li.RGB.Individual=[1 0 0]; end        % Colour of the individual subject of interest on the detailed plot ( [0 1 1] is red )
try iBT.what.li.LineStyle.Individual; catch iBT.what.li.LineStyle.Individual='-'; end % Style of the individual subject of interest on the detailed plot 
try iBT.what.li.LineWidth.Default; catch iBT.what.li.LineWidth.Default = 1; end % 1 is MATLAB's usual default.
try iBT.what.li.LineWidth.Individual; catch iBT.what.li.LineWidth.Individual = iBT.what.li.LineWidth.Default; end % 1 is MATLAB's usual default.

try iBT.what.li.SummaryPlot.RGB.GroupMean; catch iBT.what.li.SummaryPlot.RGB.GroupMean=iBT.what.li.RGB.GroupMean; end        % Colour of the control group mean on summary plot 
try iBT.what.li.SummaryPlot.LineStyle.GroupMean; catch iBT.what.li.SummaryPlot.LineStyle.GroupMean=iBT.what.li.LineStyle.GroupMean; end % Style of the control group mean on summary plot
try iBT.what.li.SummaryPlot.RGB.LIpoint; catch iBT.what.li.SummaryPlot.RGB.LIpoint=iBT.what.li.RGB.Individual; end        % Colour of the LI point on the summary plot 
try iBT.what.li.SummaryPlot.LineWidth.Default; catch iBT.what.li.SummaryPlot.LineWidth.Default = iBT.what.li.LineWidth.Default; end % 1 is MATLAB's usual default.
try iBT.what.li.SummaryPlot.LineWidth.Individual; catch iBT.what.li.SummaryPlot.LineWidth.Individual = iBT.what.li.SummaryPlot.LineWidth.Default; end % 1 is MATLAB's usual default.
try iBT.what.li.SummaryPlot.LineWidth.LIpoint; catch iBT.what.li.SummaryPlot.LineWidth.LIpoint = 1; end % 
try iBT.what.li.SummaryPlot.MarkerSize.LIpoint; catch iBT.what.li.SummaryPlot.MarkerSize.LIpoint = 5 ; end 

% Use iBT.what.li.[Full | SummaryPlot].* syntax wherever relevant, whilst preserving backwards compatibility with legacy variables:
 try default_fullPlotAllY = iBT.what.li.fullPlotAllY; catch default_fullPlotAllY = 0; end
 try iBT.what.li.FullPlot.AllY; catch iBT.what.li.FullPlot.AllY = default_fullPlotAllY; end
 try summaryPlotAllY = iBT.what.li.summaryPlotAllY; catch summaryPlotAllY = 0; end % We use an internal variable for this one because we sometimes change it later in this function.
 try summaryPlotAllY = iBT.what.li.SummaryPlot.AllY; catch iBT.what.li.SummaryPlot.AllY = summaryPlotAllY; end % Prefer the newer form variable if it exists.
 try default_summary_plot_in_range_only=iBT.what.li.summary_plot_in_range_only; catch default_summary_plot_in_range_only=0; end
 try default_summary_plot_lower_range_fraction=iBT.what.li.summary_plot_lower_range_fraction; catch default_summary_plot_lower_range_fraction=2/3; end
 try default_summary_plot_upper_range_fraction=iBT.what.li.summary_plot_upper_range_fraction; catch default_summary_plot_upper_range_fraction=2/3; end
 try iBT.what.li.SummaryPlot.in_range_only; catch iBT.what.li.SummaryPlot.in_range_only=default_summary_plot_in_range_only; end
 try iBT.what.li.SummaryPlot.lower_range_fraction; catch iBT.what.li.SummaryPlot.lower_range_fraction=default_summary_plot_lower_range_fraction; end
 try iBT.what.li.SummaryPlot.upper_range_fraction; catch iBT.what.li.SummaryPlot.upper_range_fraction=default_summary_plot_upper_range_fraction; end

try iBT.what.li.sanity.do = iBT.what.li.sanity.do; catch iBT.what.li.sanity.do = 0; end
try iBT.what.li.conclusion.do = iBT.what.li.conclusion.do; catch iBT.what.li.conclusion.do = 1; end
try iBT.what.li.save.thresh = iBT.what.li.save.thresh; catch iBT.what.li.save.thresh = 0; end
try iBT.what.li.save.control_group_stats_filename = iBT.what.li.save.control_group_stats_filename; catch iBT.what.li.save.control_group_stats_filename = ''; end
try iBT.what.li.FullPlot.Box=iBT.what.li.FullPlot.Box; catch iBT.what.li.FullPlot.Box = 1; end
try iBT.what.li.SummaryPlot.Box=iBT.what.li.SummaryPlot.Box; catch iBT.what.li.SummaryPlot.Box = 1; end
try iBT.what.li.jb_alpha = iBT.what.li.jb_alpha; catch iBT.what.li.jb_alpha=0.05; end % Alpha level to use for the Jarque-Bera normality test. Set to 0 for no JB test.
try iBT.what.li.lowNpercentValid; catch iBT.what.li.lowNpercentValid=95; end %  Valid range and limit for Jarque-Bera normality test is subject to this % of controls having not exhausted their positive voxels.
try iBT.what.li.ceiling = iBT.what.li.ceiling; catch iBT.what.li.ceiling=1.0; end % 1 for maximum.
try iBT.what.li.lowNpercentOffCeiling = iBT.what.li.lowNpercentOffCeiling; catch iBT.what.li.lowNpercentOffCeiling=0; end % Set to 0 for no percentOffCeiling constraint.
default_typicalString = 'Conclusion: typical';
default_atypicalString = 'Conclusion: atypical';
try iBT.what.li.SummaryPlot.typicalString; catch iBT.what.li.SummaryPlot.typicalString=default_typicalString; end
try iBT.what.li.SummaryPlot.atypicalString; catch iBT.what.li.SummaryPlot.atypicalString=default_atypicalString; end
try iBT.what.li.SummaryPlot.show100; catch iBT.what.li.SummaryPlot.show100 = 0; end
try iBT.what.li.SummaryPlot.LIpoint; catch iBT.what.li.SummaryPlot.LIpoint = 0; end

% We use the "sanity" section to caclulate our conclusion, 
% so if we don't want the sanity check but do want conclusion,
% we'll just avoid displaying the sanity check result (flag here with do_sanity=-1)
do_sanity = iBT.what.li.sanity.do;
if (iBT.what.li.conclusion.do == 1)
	 if (do_sanity == 0), do_sanity = -1; end
	 summaryPlotAllY = 1; % We override any previous setting to ensure room for text. 
end

sanity_done = 0; % Initialise
Pstring_neg = ''; % Initialise
li_string01=''; 		% Initialise
li_string01space='                ';% Initialise
li_string02='';			% Initialise
typical_string='';% Initialise
controls_min_valid_Nvox = -1; % Initialise
li_at_controls_min_valid_Nvox = -999; % Initialise

warningjbtest = warning('off', 'stats:jbtest:OutOfRangeP'); % This doesn't matter.
warningjbtest = warning('off', 'stats:jbtest:PTooSmall');   % This doesn't matter either.

	num_controls = size(li_group,2); % Will be zero if li_group is an empty structure {}
	control_length = 0; % Initialise
    ccon = iBT.what.li.control_index_mapping(con); % Contrast in the control group corresponding to that which we are currently analysing

    	if num_controls ~= 0
		if ccon > numel(li_group{1})
			disp('Error (iBT_latrerality_graph): Control data contains fewer contrasts than current analysis')
			error('Error (iBT_latrerality_graph): Unable to continue. Please consider use of iBT.what.li.control_index_mapping');
		end
		for i=1:num_controls, % Work out how many points we have for each subject/contrast combination
		  li_controls{ccon}.length{i} = size(li_group{i}{ccon}.li,2);
		  control_length = max(li_controls{ccon}.length{i}, control_length); % maximum length we need to accomodate
		end
		
		% Create two fixed length arrays: one will contain the li's padded by zeroes
		% the other will contain 1 if data in that element is present for this subject, 
		% or 0 if that element is padding. This facilitates summing and calculation of
		% a softmean.
		li_controls{ccon}.li = zeros(num_controls,control_length); %Initialise
		li_controls{ccon}.valid = zeros(num_controls,control_length); %Initialise
 		li_controls{ccon}.sum_li= zeros(1,control_length); %Initialise
		li_controls{ccon}.sum_valid= zeros(1,control_length); %Initialise
		for i=1:num_controls
		  li_controls{ccon}.li(i,1:li_controls{ccon}.length{i})  = li_group{i}{ccon}.li;
		  li_controls{ccon}.valid(i,1:li_controls{ccon}.length{i}) = 1;
		  li_controls{ccon}.sum_li= li_controls{ccon}.sum_li + li_controls{ccon}.li(i,:);
		  li_controls{ccon}.sum_valid= li_controls{ccon}.sum_valid+ li_controls{ccon}.valid(i,:);
		  li_controls{ccon}.li(i,li_controls{ccon}.length{i}+1:control_length) = NaN; % Set any padding to NaN's.
		end; %for i=1:num_controls

		li_controls{ccon}.mean_li= ibt_nanmean(li_controls{ccon}.li); % Mean ignoring NaN's
		li_controls{ccon}.stdev_li= ibt_nanstd(li_controls{ccon}.li); % StDev ignoring NaN's
		[ ttest_li{ccon} p_li{ccon} CI_li_mean{ccon} stats{ccon} ] = ttest( bsxfun(@minus, li_controls{ccon}.li, li_controls{ccon}.mean_li), 0, 0.05 * 2 ); 
		% CI_li_mean is the single tailed 95% confidence interval for the population mean, 
		% but we really want the 95% CI for individual data points rather than the mean, for our question is
		% what is the probability that a test subject comes from a similar population.
		% That is, what is the StDev multiplied by the critical t-value for P=0.05 single tailed.
		% Can use crit_t = tinv(single_tailed_prob,df)
		li_controls{ccon}.critical_t95 = tinv(0.95,stats{ccon}.df);
		li_controls{ccon}.critical_t99 = tinv(0.99,stats{ccon}.df);
		li_controls{ccon}.CI95_li = stats{ccon}.sd .* li_controls{ccon}.critical_t95;
		li_controls{ccon}.CI99_li = stats{ccon}.sd .* li_controls{ccon}.critical_t99;
		
		
		% Where should we start considering N's (i.e lowest N that provides a reasonable dynamic range accross subjects,
		% without instability (N too low) or ceiling effect (everyone |LI| == 1). Original method we published was to 
		% use Jarque-Bera test to test for normality, but this might not work well for initial determination of groups 
		% where we need to first identify and remove atypical (right and bilateral) cases. So we now also implement another
		% option: lowest N where at least a specified % of subjects have |LI| < ceiling. We also continue to demand at
		% least percentage iBT.what.li.lowNpercentValid (default 95) of subjects have not exhausted all of their positive voxels
		data = li_controls{ccon}.li;
	        valid_percentOffCeiling_found = 0; % initialise
	        max_valid_found = 0; % initialise
		for n = 1:control_length
		    if stats{ccon}.df(n) >= (iBT.what.li.lowNpercentValid / 100.0 * num_controls) - 1.5 ; % -1 for df -> N and -0.5 for rounding. Tried num_controls/2 but gave too high N so curve too near 0
		 	li_controls{ccon}.OffCeiling(n) = sum( abs(data(:,n)) < iBT.what.li.ceiling ); % This is the number of subjects whose |LI| at this n is less than Ceiling
			if li_controls{ccon}.OffCeiling(n) >= ( num_controls * iBT.what.li.lowNpercentOffCeiling / 100 ); % We've enough subjects away from the ceiling
				if ~valid_percentOffCeiling_found,
					valid_percentOffCeiling_found = 1;
					li_controls{ccon}.min_valid_index = n;
					li_controls{ccon}.min_valid = n * step;
				end
			end
		    else % if stats{ccon}.df(n) >= (iBT.what.li.lowNpercentValid...
		    		if ~max_valid_found
					max_valid_found = 1;
					li_controls{ccon}.max_valid_index = max([0,n-1]);
					li_controls{ccon}.max_valid = li_controls{ccon}.max_valid_index * step;
				end
		    end % if stats{ccon}.df(n) >= (iBT.what.li.lowNpercentValid...
		end % for n = 1:control_length
		
		if ~valid_percentOffCeiling_found, % no valid region - ceiling effect always present (this would be strange!)
					li_controls{ccon}.min_valid_index = 0;
					li_controls{ccon}.max_valid_index = 0;
					li_controls{ccon}.min_valid = 0;
					li_controls{ccon}.max_valid = 0;
		elseif ~max_valid_found, % There were no problems right up to the end
					li_controls{ccon}.max_valid_index = control_length;
					li_controls{ccon}.max_valid = control_length * step;
		end

		if iBT.what.li.SummaryPlot.in_range_only==2
			xhigh_index{ccon} = find(stats{ccon}.df >= (iBT.what.li.SummaryPlot.upper_range_fraction * num_controls) - 1.5,1,'last');
		end
		
		if iBT.what.li.jb_alpha ~= 0,
		 % For which N's is the data normally distributed? Use the Jarque-Bera test 
		 % of null hypotheisis that data is normally distributed, with alpha = jb_alpha
		 % (you might want to use a lenient alpha because the test is for non-normality,
		 % whereas we're more interested in normality so need the null to be a reasonable
		 % possibility - in practice however there was not much difference when I tried 0.05 and 0.2).
		 % Unfortunately jbtest doesn't appear to support matrices, so need to do in a loop.
		 % Further, there is not much point testing for non-normality when the degrees of
		 % freedom fall too low (as will happen for higher N values because not all subjects
		 % will have that many active voxels). Also, as N increases t decreases and
		 % eventually noise begins to dominate, so probably not very helpful to look at
		 % very high N. Therefore we demand at least percentage iBT.what.li.lowNpercentValid
		 % (default 95) of control subjects have not exhausted their positive voxels at the maximum N.
		 for n = 1:control_length
			if stats{ccon}.df(n) >= (iBT.what.li.lowNpercentValid / 100.0 * num_controls) - 1.5 ; % -1 for df -> N and -0.5 for rounding. Tried num_controls/2 but gave too high N so curve too near 0
				% work around a bug in some older versions of Matlab that had jbtest of a flat distribution returning 0 or crashing.
				if sum( data(:,n) ~= data(1,n) ) == 0
					% All elements are the same, so reject the NULL that data is normally distributed:
					li_controls{ccon}.jbH(n) = 1; li_controls{ccon}.jbP(n) = NaN;
				else
					[li_controls{ccon}.jbH(n) li_controls{ccon}.jbP(n)] = jbtest(data(:,n),iBT.what.li.jb_alpha);
				end
			else
				li_controls{ccon}.jbH(n) = 1; li_controls{ccon}.jbP(n) = NaN;
			end
		 end

		 % Find the valid range. We'll use the largest contiguous region
		 % of valid values as the valid range. If more than one valid range
		 % and the longest have the same maximum length, use the first of the
		 % longest. We use a run-length encoding type scheme to find this.
		 x = li_controls{ccon}.jbH; % Makes code simpler
		 % The following three elegant lines of code are adapted from rle.m Version 1.0 submitted
		 %  to matlabcentral by Stefan Eireiner:
		 %	http://www.mathworks.com/matlabcentral/fileexchange/4955-rle-deencoding	
		  i = [ find(x(1:end-1) ~= x(2:end)) length(x) ]; % The transitions between 0 & 1, and the final index.
		  xlengths = diff([ 0 i ]); % Length of each run of contiguous values
		  xvals = x(i); % Value of each run
		 v = find(xvals==0); % Where there are runs of zeroes
		 if isempty(v) % then we have no valid region.
			li_controls{ccon}.min_valid_index = 0;
			li_controls{ccon}.max_valid_index = 0;
			li_controls{ccon}.min_valid = 0;
			li_controls{ccon}.max_valid = 0;
		 else % find longest contiguous valid region
			vlengths = xlengths(v); % Length of each run of contiguous zeroes
		 	vlength_max = max(vlengths); % length of longest run of contiguous zeroes
			max_v = find( vlengths==vlength_max, 1 ); % First run of maximum length of contigous zeroes
			% v(max_v) is the index into xvals of the start of desired valid run
			min_valid_index = 1 + sum( xlengths(1:v(max_v)-1) ); % Start of desired valid run.
			% Note if v(max_v) is 1 then the valid run is at the start of x; Matlab evaluates sum( xlengths(1:0) ) =0 so above still OK.
			li_controls{ccon}.min_valid_index = max([li_controls{ccon}.min_valid_index, min_valid_index]); % further restrict extisting limit if necessary
			li_controls{ccon}.max_valid_index = min([li_controls{ccon}.max_valid_index, li_controls{ccon}.min_valid_index + vlength_max - 1]); % further restrict extisting limit if necessary
			li_controls{ccon}.min_valid = li_controls{ccon}.min_valid_index * step;
			li_controls{ccon}.max_valid = li_controls{ccon}.max_valid_index * step;
		 end % if isempty(v)	
		end %if iBT.what.li.jb_alpha ~= 0
			
    	end; %num_controls ~= 0

%	maxXis = max(control_length,size(li_individual{con}.li)); Problem with this is we can't plot error bars beyond the control length 
	maxXis = control_length;
        if maxX == 0, Xend  = maxXis; % Size plot to fit in all control data
	else Xend = min(maxX/step,maxXis); end % Manually limit max of plot
	
    	x=(1:Xend)*step; % Size plot to fit in all control and subject data
   
	full_fignum = fignums(1);
    	fignum = fignums(2);
	
    	if num_controls ~= 0 
   		sfigure(full_fignum); %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% The full plot %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
		if iBT.what.li.LineWidth.Default > 0,  set(full_fignum, 'DefaultLineLineWidth', iBT.what.li.LineWidth.Default); end
		% Determine whether we're dealing with a group result that is mainly positive or negative.
		% (we plot atypical as the probability of deviation towards the opposite sign)
		if ibt_nanmedian(li_controls{ccon}.mean_li) >= 0 
			GroupSign = 1;
			
			miny_controls = min(min(li_controls{ccon}.li));
			miny_subject = min(li_individual{con}.li);
			miny=min([0 miny_controls miny_subject]);
			if iBT.what.li.FullPlot.AllY == 1, 
				yfill_to = -GroupSign; 
			elseif iBT.what.li.FullPlot.AllY == 2
				yfill_to = 0; 
			else 
				yfill_to = miny; 
			end			
		else 
			GroupSign = -1;

			maxy_controls = max(max(li_controls{ccon}.li));
			maxy_subject = max(li_individual{con}.li);
			maxy=max([0 maxy_controls maxy_subject]);			
			if iBT.what.li.FullPlot.AllY == 1, 
				yfill_to = -GroupSign; 
			elseif iBT.what.li.FullPlot.AllY == 2
				yfill_to = 0; 
			else 
				yfill_to = maxy;
			end			
		end
		
		y95 = li_controls{ccon}.mean_li- ( GroupSign * li_controls{ccon}.CI95_li(1,:) ); % The 95% confidence line
		y99 = li_controls{ccon}.mean_li- ( GroupSign * li_controls{ccon}.CI99_li(1,:) ); % The 99% confidence line
		% Now create versions with NaN's set to zero for fill operations:
		yy95 = y95(1:Xend); yy99 = y99(1:Xend);
		y95NaN = find(isnan(yy95));
		y99NaN = find(isnan(yy99));
		if length(y95NaN) ~= 0 yy95(y95NaN) = 0;end
		if length(y99NaN) ~= 0 yy99(y99NaN) = 0;end
		% fill([0 0 x x(end)], [0 yy95(1) yy95 0 ], [0.9 0.9 0.9], 'EdgeColor','w')
		hold off; % In case we've already plotted something, for example the previous contrast
		fill([0 x(end) x(end) 0], [1 1 -1 -1 ], iBT.what.li.RGB.Invalid, 'EdgeColor',iBT.what.li.RGB.Invalid); %Fill entire plot
		hold on;
 		fill([li_controls{ccon}.min_valid li_controls{ccon}.min_valid li_controls{ccon}.max_valid li_controls{ccon}.max_valid], [-1 1 1 -1 ], 'w', 'EdgeColor','w');
		fill([0 0 x x(end)], [yfill_to yy95(1) yy95 yfill_to ], iBT.what.li.RGB.CI95, 'EdgeColor',iBT.what.li.RGB.CI95);
		fill([0 0 x x(end)], [yfill_to yy99(1) yy99 yfill_to ], iBT.what.li.RGB.CI99, 'EdgeColor',iBT.what.li.RGB.CI99);	
		if iBT.what.li.errorbars == 1, errorbar(x, li_controls{ccon}.mean_li(1:Xend),li_controls{ccon}.mean_li(1:Xend) - yy95,0*x, iBT.what.li.LineStyle.GroupMean, 'color', iBT.what.li.RGB.GroupMean); end
		tmp = li_controls{ccon}.li' ;
		plot(x,tmp(1:Xend,:), iBT.what.li.LineStyle.GroupControls, 'color', iBT.what.li.RGB.GroupControls);
		plot(x, li_controls{ccon}.mean_li(1:Xend), iBT.what.li.LineStyle.GroupMean, 'color', iBT.what.li.RGB.GroupMean);		
 		if iBT.what.li.PlotGroupOnly == 0, plot(xaxis, li_individual{con}.li, iBT.what.li.LineStyle.Individual, 'color', iBT.what.li.RGB.Individual, 'LineWidth', iBT.what.li.LineWidth.Individual); end
                if iBT.what.li.PlotJB, plot(x,(GroupSign * li_controls{ccon}.jbP(1:Xend)),'g'); end %If we want to look at the Jarque-Bera P values.
		plot([0 x(end)], [0 0],'k:'); % The li=0 (Bilateral) line
		axis tight
		set(full_fignum,'Color',[1 1 1]); % Set plot background to white
		try full_plot_titleString = iBT.what.li.FullPlot.Title; % If title is being forced
		catch, full_plot_titleString = [fileIDstring ' ' fileBaseName{con}]; % Normal title
		end; % catch
    	        title_handle = title(full_plot_titleString,'Interpreter','none','fontsize',title_fontsize,'FontWeight','Bold'); 
	        % The 'Interpreter','none' above stops underscores being interpreted as a next-character-lowercase instruction.
		v = axis;
		v(3) = max(-1,v(3)); % Don't plot y lower than -1  (could happen if CI's extend out past this)
		v(4) = min(1,v(4));  % Don't plot y higher than +1 (could happen if CI's extend out past this)
		if GroupSign == 1, v(3) = max(yfill_to,v(3)); end;
		if GroupSign == -1, v(4) = min(yfill_to,v(4)); end;
		axis(v)
		xlabel(xaxis_label,'fontsize',axis_label_fontsize,'FontWeight','Bold') %label the xaxis
		ylabel('Laterality index','fontsize',axis_label_fontsize,'FontWeight','Bold') %label the yaxis
		set(gca,'fontsize',axis_fontsize,'layer','top');
		set(title_handle,'fontsize',title_fontsize);
		if iBT.what.li.FullPlot.Box==0, box off; else box on; end
		% Now do a clean plot within the valid range
		% First we find the minimum y of this subject within the valid range of the group:
		sizey = size(li_individual{con}.li,2);  %Number of valid elements in this subject.

		xrange_strict_controls = li_controls{ccon}.min_valid_index:li_controls{ccon}.max_valid_index;
		xrange_strict_individual = li_controls{ccon}.min_valid_index:min([li_controls{ccon}.max_valid_index,sizey]);

		if (iBT.what.li.lowNpercentOffCeiling == 0) && (iBT.what.li.jb_alpha ~= 0)
			disp(sprintf(['Valid range over which control li''s for ' fileBaseName{con} ' are normally distributed and at least %g%% contribute: %i:%i'],...
		        	iBT.what.li.lowNpercentValid, li_controls{ccon}.min_valid, li_controls{ccon}.max_valid ));
		elseif (iBT.what.li.lowNpercentOffCeiling ~= 0) && (iBT.what.li.jb_alpha == 0)
			disp(sprintf(['Valid range over which %i%% of control li''s for ' fileBaseName{con} ' have magnitude less than %g and at least %g%% contribute: %i:%i'],...
		        	iBT.what.li.lowNpercentOffCeiling, iBT.what.li.ceiling,...
				iBT.what.li.lowNpercentValid, li_controls{ccon}.min_valid, li_controls{ccon}.max_valid ));
		elseif (iBT.what.li.lowNpercentOffCeiling ~= 0) && (iBT.what.li.jb_alpha ~= 0)
			disp(sprintf(['Valid range over which control li''s for ' fileBaseName{con} ' are normally distributed and %i%% have magnitude less than %g and at least %g%% contribute: %i:%i'],...
		        	iBT.what.li.lowNpercentOffCeiling, iBT.what.li.ceiling,...
				iBT.what.li.lowNpercentValid, li_controls{ccon}.min_valid, li_controls{ccon}.max_valid ));
		end
		
		if strcmp(xaxis_label,'Number of voxels above threshold') && ( (do_sanity ~= 0) || (iBT.what.li.save.thresh == 1) ) && (con_deactivation > 0) && (li_controls{ccon}.max_valid_index > 0)
			actThres = NaN  ; % Initilise - if remains unchanged then we'll know there was no valid actThres
			deactThres = NaN; % Initilise - if remains unchanged then we'll know there was no valid deactThres
			sizey_deactivation = size(li_individual{con_deactivation}.li,2);  %Number of valid elements in this subject.
			if ( li_controls{ccon}.min_valid_index <= sizey ) && ( li_controls{ccon}.min_valid_index <= sizey_deactivation )			
				actThres= li_individual{con}.thresh(li_controls{ccon}.min_valid_index)
				deactThres=li_individual{con_deactivation}.thresh(li_controls{ccon}.min_valid_index)
				sanity_string = { sprintf('Activation at N=%i, t=%.2f: LI=%+.2f' ,li_controls{ccon}.min_valid , actThres,li_individual{con}.li(li_controls{ccon}.min_valid_index)) };
				if iBT.what.li.SummaryPlot.show100 == 1,
% li_individual{con}.li(li_controls{ccon}.min_valid_index)=0.995; %testing
					li_string01 = sprintf('LI = %+.0f%%',100*li_individual{con}.li(li_controls{ccon}.min_valid_index));
					if ( 100*li_individual{con}.li(li_controls{ccon}.min_valid_index)) >= 95.5, li_string01space = [li_string01space '  ']; end; % Make more space if three digits
					if (-100*li_individual{con}.li(li_controls{ccon}.min_valid_index)) >= 95.5, li_string01space = [li_string01space '  ']; end; % Make more space if three digits
				else
					li_string01 = sprintf('LI = %+.2f',li_individual{con}.li(li_controls{ccon}.min_valid_index)); % This is the value we display in big writing on the summary plot
				end
				if (li_individual{con}.li(li_controls{ccon}.min_valid_index)) >= 0, li_string01space = [li_string01space '  ']; end; % Make more space for plus sign (it is wider than the minus?) if positive
				if iBT.what.li.SummaryPlot.LIpoint == 1,
					li_string01 = ['*' li_string01];
					li_string01space = [li_string01space '  '];
				end
				if iBT.what.li.SummaryPlot.lower_range_fraction == 0, % then offset the text to keep it a couple of spaces away from the axis
					li_string01 = ['  ' li_string01];
					li_string01space = [li_string01space '  '];
				end
				
				li_string02 = sprintf('(at N=%i)',li_controls{ccon}.min_valid);
				controls_min_valid_Nvox = li_controls{ccon}.min_valid
				li_at_controls_min_valid_Nvox = li_individual{con}.li(li_controls{ccon}.min_valid_index);
				if actThres > deactThres
					sanity_string{2} = 'Check passed: activation > deactivation'
				else
					sanity_string{2} = 'WARNING: Deactivation is more significant than activation!'
					sanity_string{3} = sprintf('Deactivation @ N=%i, t= %.2f: LI=%+.2f' ,li_controls{ccon}.min_valid , deactThres, li_individual{con_deactivation}.li(li_controls{ccon}.min_valid_index));
				end	
			else
				%sanity_string = {sprintf('Unable to determine activation:deactivation significance at N = %i voxels.', li_controls{ccon}.min_valid)};
				if  li_controls{ccon}.min_valid_index <= sizey 
					actThres= li_individual{con}.thresh(li_controls{ccon}.min_valid_index)
					sanity_string = { sprintf('Activation @ N=%i, t=%.2f: LI=%+.2f' ,li_controls{ccon}.min_valid , actThres,li_individual{con}.li(li_controls{ccon}.min_valid_index)) };
					sanity_string{2} = sprintf('Check passed (Deactivation maxN = %i).',sizey_deactivation * step);
				elseif li_controls{ccon}.min_valid_index <= sizey_deactivation
					deactThres=li_individual{con_deactivation}.thresh(li_controls{ccon}.min_valid_index)
					sanity_string = {'WARNING: Deactivation is more significant than activation!'}
					sanity_string{2} = sprintf('Deactivation @ N=%i, t=%.2f: LI=%+.2f' ,li_controls{ccon}.min_valid , deactThres, li_individual{con_deactivation}.li(li_controls{ccon}.min_valid_index));
				else
					sanity_string = {'WARNING: Unable to perform check (no thresholds within valid range).'}
				end
			end
			if (iBT.what.li.sanity.do == 1) sanity_done = 1; end; % If we did calculation solely due to iBT.what.li.save.thresh then we don't report the sanity result.
			if (iBT.what.li.save.thresh == 1) 
				li_display_thresholds.structure_version = 1.0; %Version number of the 'li_display_thresholds' structure.
				li_display_thresholds.fileBaseName = fileBaseName{con};
				li_display_thresholds.controls_min_valid_Nvox = li_controls{ccon}.min_valid;
				li_display_thresholds.controls_max_valid_Nvox = li_controls{ccon}.max_valid;
				li_display_thresholds.actThres_at_min_valid_Nvox = actThres;
				li_display_thresholds.deactThres_at_min_valid_Nvox = deactThres;
				if isnan(actThres)
					li_display_thresholds.actLI_at_min_valid_Nvox = NaN;
				else
					li_display_thresholds.actLI_at_min_valid_Nvox =   li_individual{con}.li(li_controls{ccon}.min_valid_index);
				end; % if isnan(actThres)
				if isnan(deactThres)
					li_display_thresholds.deactLI_at_min_valid_Nvox = NaN;
				else
					li_display_thresholds.deactLI_at_min_valid_Nvox = li_individual{con_deactivation}.li(li_controls{ccon}.min_valid_index);
				end
	 			mat_filename=fullfile( iBT.what.li.save.path, sprintf('li_display_thresholds_%s_%.4i.mat',fileIDstring,con) );
        			save (mat_filename, 'li_display_thresholds');
				disp(['Saved ' mat_filename]);
			end
		end

		if (li_controls{ccon}.max_valid_index > 0) && ( li_controls{ccon}.min_valid_index <= sizey )
		    xrange = li_controls{ccon}.min_valid_index:min(li_controls{ccon}.max_valid_index,sizey);
		    P_available = 1;
		    if sizey < li_controls{ccon}.max_valid_index 
		    	disp('Warning: Number of voxels at most lenient threshold is below the maximum of the valid control range')
		    end
		    min_x_to_plot = li_controls{ccon}.min_valid;
		else
		    % In this case there are no valid values of x to plot 
		    % so we plot from the start rather than min_valid .
		    xrange = 1:min(li_controls{ccon}.max_valid_index,sizey);
		    P_available = 0;
		    if (li_controls{ccon}.max_valid_index > 0) && ( li_controls{ccon}.min_valid_index > sizey )
		    	disp('Warning: Number of voxels at most lenient threshold is below the minimum of the valid control range,')
		        disp('         so will not be able to calculate a P value.')
		    	Pstring = 'P unavailable; data not in valid control range.'; % Used in title later.
		    elseif (li_controls{ccon}.max_valid_index == 0)
		    	disp('Warning: No valid control range')
		        disp('         so will not be able to calculate a P value.')
		        Pstring = 'P unavailable; no valid control range.'; % Used in title later.
		    end
		    min_x_to_plot = 0;
		end

		xrange_strict_controls = li_controls{ccon}.min_valid_index:li_controls{ccon}.max_valid_index;
		xrange_strict_individual = li_controls{ccon}.min_valid_index:min([li_controls{ccon}.max_valid_index,sizey]);

		if strcmp(iBT.what.li.save.control_group_stats_filename,'') == 0 , % A quick hack so that we can later compare two groups (the toolbox is presently designed for individual vs group comparisons)
			control_group_stats_filename = sprintf('%s_con%0.2i',iBT.what.li.save.control_group_stats_filename,ccon)
			save(control_group_stats_filename,'ccon','li_controls','mean_li_diff_over_valid_n_range_for_each_control', 'xrange_strict_controls');            
			disp(['Control_group_stats saved to' control_group_stats_filename]);
			%% To use this to compare two groups, for example subjects performing a right and a left-handed task, 
			%% can do the following from the Matlab command line after using the toolbox to
			%% generate summary stats for each of the two groups and the groups combined. Use a distinct
			%% name (set in iBT.what.li.save.control_group_stats_filename) for each of the three LI analyses, 
			%% for example 'Rt_group_stats', 'Lt_group_stats' and both_group_stats. Then from MATLAB:
			% R = load('Rt_group_stats_con01.mat') % Stats saved for first group
			% L = load('Lt_group_stats_con01.mat') % Stats saved for second group
			% B = load('both_group_stats_con01.mat') % Stats saved for both groups (i.e. need to run a plot with all subjects from both groups included)
			% li_controls{B.ccon}.mean_li = ibt_nanmean( [B.li_controls{B.ccon}.li] ) % Using the mean accross both groups establishes an unbiased adjustment
			%% Calculate scores for each individual, adjusted for the whole group mean, and average accross nvox
			% mean_li_diff_over_valid_n_range_for_each_L = ibt_nanmean( bsxfun(@minus, L.li_controls{L.ccon}.li(:,B.xrange_strict_controls), B.li_controls{B.ccon}.mean_li(B.xrange_strict_controls))' )
			% mean_li_diff_over_valid_n_range_for_each_R = ibt_nanmean( bsxfun(@minus, R.li_controls{R.ccon}.li(:,B.xrange_strict_controls), B.li_controls{B.ccon}.mean_li(B.xrange_strict_controls))' )
			%% Now compare the adjusted summary scores accross the two groups (the following example does 
			%% a two-tailed paired t-test; use ttest2 function for unpaired; change 'both' to 'left' or 'right'
			%% for a single-tailed test):
			% [ ttest_li p_li ] = ttest( mean_li_diff_over_valid_n_range_for_each_L, mean_li_diff_over_valid_n_range_for_each_R, 0.05, 'both' )
		end

		% Notwithstanding the default range for the summary plot selected above, we might want to do something different:	
		if ( iBT.what.li.SummaryPlot.in_range_only == -1 ) || (li_controls{ccon}.min_valid == li_controls{ccon}.max_valid) % plot full range
			xSummaryRange_controls = 1:control_length;
			xSummaryRange_individual = 1:sizey;
		    	min_x_to_plot = 0;
			max_x_to_plot = x(end);
		elseif iBT.what.li.SummaryPlot.in_range_only == 0 % plot full range at lower end only
			xSummaryRange_controls = 1:min(li_controls{ccon}.max_valid_index,sizey);
			xSummaryRange_individual = xSummaryRange_controls;
			min_x_to_plot = 0;
			max_x_to_plot = li_controls{ccon}.max_valid;
		elseif iBT.what.li.SummaryPlot.in_range_only == 1 % plot in valid range only
		 	xSummaryRange_controls = li_controls{ccon}.min_valid_index:min(li_controls{ccon}.max_valid_index,sizey);
			xSummaryRange_individual = xSummaryRange_controls;
		 	min_x_to_plot = li_controls{ccon}.min_valid;
			max_x_to_plot = li_controls{ccon}.max_valid;
		elseif iBT.what.li.SummaryPlot.in_range_only > 1 % include some fraction of invalid range in plot
			if min_x_to_plot == 0
				xlow_index = 1;
			else 
		    		xlow_index = max(1, floor(li_controls{ccon}.min_valid_index * (1-iBT.what.li.SummaryPlot.lower_range_fraction)));
			end
			%try,disp(sprintf('xhigh_index{%i} = %i', ccon, xhigh_index{ccon}));,catch,end % Display this if it was calculated
		        xSummaryRange_controls = (xlow_index:xhigh_index{ccon});
			xSummaryRange_individual = (xlow_index:min(xhigh_index{ccon},sizey));
			min_x_to_plot = x(xlow_index);
			max_x_to_plot = x( min( [xhigh_index{ccon} numel(x)]) ); % numel(x) might be lower if we specified a fixed maxNvox
		end
			
		if (li_controls{ccon}.max_valid_index > 0)
			 whatsign = ibt_nanmedian(li_controls{ccon}.mean_li(xrange_strict_controls));
		else
			 whatsign = ibt_nanmedian(li_controls{ccon}.mean_li(xSummaryRange_controls));
		end
					
		if whatsign >= 0 
			GroupSign = 1;  % We'll t-test againt the null hypothesis that the li is higher than
			tail = 'right'; % the group mean, so a significant P indicates bilateral or right laterality
			tail_meaning = 'null: li >= controls'; 
			neg_tail_meaning = 'null: -li >= controls'; 
			
			miny_controls = min(y95(xSummaryRange_controls));
			miny_subject = min(li_individual{con}.li(xSummaryRange_individual));
			miny=min([0 miny_controls miny_subject]);
			
			if summaryPlotAllY == 1, 
				yfill_to = -GroupSign; 
			elseif summaryPlotAllY == 2, 
				yfill_to = 0; 
			else
				yfill_to = miny;
			end			
		else 
			GroupSign = -1; % We'll t-test againt the null hypothesis that the li is lower than
			tail = 'left';  % the group mean, so a significant P indicates bilateral or left laterality
			tail_meaning = 'null: li <= controls'; 
			neg_tail_meaning = 'null: -li <= controls'; 

			maxy_controls = max(y95(xSummaryRange_controls));
			maxy_subject = max(li_individual{con}.li(xSummaryRange_individual));
			maxy=max([0 maxy_controls maxy_subject]);			

			if summaryPlotAllY == 1, 
				yfill_to = -GroupSign; 
			elseif summaryPlotAllY == 2, 
				yfill_to = 0; 
			else
				yfill_to = maxy;
			end			
		end
		
		if P_available == 1,
		    % Determine probability that subject is typical by comparing the average, over valid xrange, of the zero-mean-adjusted group li's to the similarly adjusted individual:
		    mean_li_diff_over_valid_n_range_for_each_control{ccon} = ibt_nanmean( bsxfun(@minus, li_controls{ccon}.li(:,xrange_strict_individual), li_controls{ccon}.mean_li(xrange_strict_individual))' );
		    mean_li_diff_over_valid_n_range_for_individual{con} = ibt_nanmean(li_individual{con}.li(xrange_strict_individual) - li_controls{ccon}.mean_li(:,xrange_strict_individual));		    
		    [ ttest2_li{ccon} p2_li{ccon} ] = ttest2( mean_li_diff_over_valid_n_range_for_each_control{ccon}, mean_li_diff_over_valid_n_range_for_individual{con}, 0.05, tail );
		    Pvalue = p2_li{ccon};
		    Pstring=sprintf('%.2g',Pvalue);
		    disp(['P (' tail '-tailed t-test, ' tail_meaning ') for ' fileIDstring '_' fileBaseName{con} ' = ' Pstring]);
		    Pstring = ['P(' tail_meaning ') = ' Pstring ]; % Used in title later.
		    if Pvalue < 0.05
		     ; % The following test is conditional on the first test, so strictly the probability needs to be corrected for this, but we don't do this
		     ; % because we want to see what it would have been had we been performing the exact same test in the other direction. Need to be careful when
		     ; % reporting this to acknowledge what has been done so can interpret appropriately.
		     typical_string = iBT.what.li.SummaryPlot.atypicalString;
		     disp([fileIDstring ' ' fileBaseName{con} ' ' default_atypicalString]); %Irrespective of what is displayed on the plot, we write this to the log file
		     neg_mean_li_diff_over_valid_n_range_for_individual{con} = ibt_nanmean(-li_individual{con}.li(xrange_strict_individual) - li_controls{ccon}.mean_li(:,xrange_strict_individual));		    
		     [ ttest2_li_neg{ccon} p2_li_neg{ccon} ] = ttest2( mean_li_diff_over_valid_n_range_for_each_control{ccon}, neg_mean_li_diff_over_valid_n_range_for_individual{con}, 0.05, tail );
		     if p2_li_neg{ccon} > 0.01 % We won't bother displaying the neg test result unless this is true (trend to the opposite laterality) 
		    	Pstring_neg=sprintf('%.2g',p2_li_neg{ccon});
		    	disp(['P (' tail '-tailed t-test, ' neg_tail_meaning ') for ' fileIDstring '_' fileBaseName{con} ' = ' Pstring_neg]);
		    	Pstring_neg = [' ,  P(' neg_tail_meaning ') = ' Pstring_neg ]; % Used in title later.
		     end
		    else
		     typical_string = iBT.what.li.SummaryPlot.typicalString;		    
		     disp([fileIDstring ' ' fileBaseName{con} ' ' default_typicalString]); %Irrespective of what is displayed on the plot, we write this to the log file
		    end
		end 


   		sfigure(fignum); %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% The summary plot when there are controls %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
		if iBT.what.li.SummaryPlot.LineWidth.Default > 0,  set(fignum, 'DefaultLineLineWidth', iBT.what.li.SummaryPlot.LineWidth.Default); end
		hold off; % In case we've already plotted something, for example the previous contrast
		fill([0 x(end) x(end) 0], [1 1 -1 -1 ], iBT.what.li.RGB.Invalid, 'EdgeColor',iBT.what.li.RGB.Invalid);  %Fill entire plot
		hold on;		
 		fill([li_controls{ccon}.min_valid li_controls{ccon}.min_valid li_controls{ccon}.max_valid li_controls{ccon}.max_valid], [-1 1 1 -1 ], 'w', 'EdgeColor','w');
		fill([0 0 x x(end)], [yfill_to yy95(1) yy95 yfill_to ], iBT.what.li.RGB.CI95, 'EdgeColor',iBT.what.li.RGB.CI95);
		if iBT.what.li.errorbars == 1, errorbar(x, li_controls{ccon}.mean_li(1:Xend),li_controls{ccon}.mean_li(1:Xend) - yy95,0*x, 'color', iBT.what.li.RGB.GroupMean); end
%		plot(x, li_controls{ccon}.mean_li(1:Xend), 'color', iBT.what.li.RGB.GroupMean);
		plot(x, li_controls{ccon}.mean_li(1:Xend), iBT.what.li.SummaryPlot.LineStyle.GroupMean, 'color', iBT.what.li.SummaryPlot.RGB.GroupMean);

		if iBT.what.li.PlotGroupOnly == 0, 
			% Plot the individual subject
			plot(xaxis, li_individual{con}.li, iBT.what.li.LineStyle.Individual, 'color', iBT.what.li.RGB.Individual, 'LineWidth', iBT.what.li.SummaryPlot.LineWidth.Individual);			
		else 
			% Plot the group on the Summary plot (was used to prepare a figure for the 2010 paper)
			tmp = li_controls{ccon}.li' ;
			plot(x,tmp(1:Xend,:), iBT.what.li.LineStyle.GroupControls, 'color', iBT.what.li.RGB.GroupControls);
		end; % if iBT.what.li.PlotGroupOnly == 0, else
			
		plot([0 x(end)], [0 0],'k:'); % The li=0 (Bilateral) line
		
		if (iBT.what.li.SummaryPlot.LIpoint == 1) && (P_available == 1), % Plot the point where the single LI is determined (i.e. at highest N where controls are normaly distributed)
		    plot(li_controls{ccon}.min_valid, li_individual{con}.li(li_controls{ccon}.min_valid_index), '*', 'Color', iBT.what.li.SummaryPlot.RGB.LIpoint,...
		    	'MarkerSize',iBT.what.li.SummaryPlot.MarkerSize.LIpoint, 'LineWidth', iBT.what.li.SummaryPlot.LineWidth.LIpoint);
%		    plot(li_controls{ccon}.min_valid, li_individual{con}.li(li_controls{ccon}.min_valid_index), 'o', 'Color', iBT.what.li.SummaryPlot.RGB.LIpoint,...
%		    	'MarkerSize',iBT.what.li.SummaryPlot.MarkerSize.LIpoint, 'LineWidth', iBT.what.li.SummaryPlot.LineWidth.LIpoint);
%		    plot(li_controls{ccon}.min_valid, li_individual{con}.li(li_controls{ccon}.min_valid_index), '+', 'Color', iBT.what.li.SummaryPlot.RGB.LIpoint,...
%		        'MarkerSize',4*iBT.what.li.SummaryPlot.MarkerSize.LIpoint, 'LineWidth', iBT.what.li.SummaryPlot.LineWidth.LIpoint);
%		    %, 'MarkerSize',10);
		end% if iBT.what.li.SummaryPlot.LIpoint
		
		axis( [min_x_to_plot max_x_to_plot min(yfill_to, GroupSign) max(yfill_to, GroupSign)] );  
		v = axis;
		v(3) = max(-1,v(3))-0.001; % Don't plot y lower than -1  (could happen if CI's extend out past this); -extra 0.001 so that axes don't obliterate data
		v(4) = min(1,v(4))+0.001;  % Don't plot y higher than +1 (could happen if CI's extend out past this); +extra 0.001 so that axes don't obliterate data
		axis(v)
		set(fignum,'Color',[1 1 1]); % Set plot background to white
		xlabel(xaxis_label,'fontsize',axis_label_fontsize,'FontWeight','Bold') %label the xaxis
		ylabel('Laterality index','fontsize',axis_label_fontsize,'FontWeight','Bold') %label the yaxis
		set(gca,'fontsize',axis_fontsize,'layer','top');
   	else  % num_controls must be 0 
   		sfigure(fignum); %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% The summary plot when there are no controls %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
		if iBT.what.li.SummaryPlot.LineWidth.Default > 0,  set(fignum, 'DefaultLineLineWidth', iBT.what.li.SummaryPlot.LineWidth.Default); end
		hold off; % In case we've already plotted something, for example the previous contrast
		plot(xaxis, li_individual{con}.li, iBT.what.li.LineStyle.Individual, 'color', iBT.what.li.RGB.Individual, 'LineWidth', iBT.what.li.SummaryPlot.LineWidth.Individual); 
		if iBT.what.li.FullPlot.AllY > 0, 
			currentAxis = axis;
			axis( [currentAxis(1) currentAxis(2) -1 1] );  
		end % if iBT.what.li.FullPlot.AllY == 1
		Pstring = '';
		xlabel(xaxis_label,'fontsize',axis_label_fontsize,'FontWeight','Bold') %label the xaxis
		ylabel('Laterality index','fontsize',axis_label_fontsize,'FontWeight','Bold') %label the yaxis
		set(gca,'fontsize',axis_fontsize,'layer','top');
    	end; % if num_controls ~= 0 
	
	if verLessThan('matlab', '8.4')
	  % The  bit didn't work as expected from MatlabR2014b (when new graphics system introduced)
	  XTicks=get(gca,'XTick');
	  if XTicks(end) > 999 % Then stop Matlab using 10 to power notation.
		% XTickLabels = get(gca,'XTickLabel')
		set(gca,'XTickLabel',sprintf('%i|',XTicks));
		% We've set these labels based on screen-displayed ticks;
		% now we need to ensure the printed output uses the same ticks,
		% otherwise the labels won't match the ticks in saved plots!
		set(gca,'XTickMode',sprintf('manual')); % Ensure printed output doesn't change ticks.
	  end
	end %if verLessThan('matlab', '8.4')
	
	try titleString = iBT.what.li.SummaryPlot.Title; % manual override for title 
	catch
    	  if iBT.what.li.PlotGroupOnly == 1, 
		try titleString = iBT.what.li.PlotGroupOnlyTitle,
		catch
			titleString = [fileIDstring '_' fileBaseName{con} ];
		end
	  else
		titleString = {[fileIDstring '_' fileBaseName{con}]};
		if ~strcmp(Pstring_neg,'') 
			titleString = [titleString {[Pstring Pstring_neg]}];
		else  
			titleString = [titleString {Pstring}]; 
		end
		if (sanity_done == 1) && (do_sanity > 0)
			titleString = [titleString sanity_string]; % sanity_string is also a cell	
		end
	  end 
	end; %catch
	title_handle = title(titleString,'Interpreter','none','fontsize',title_fontsize,'FontWeight','Bold'); 
	set(title_handle,'fontsize',title_fontsize);
	% The 'Interpreter','none' above stops underscores being interpreted as a next-character-lowercase instruction.
 	if (iBT.what.li.conclusion.do == 1) && (controls_min_valid_Nvox ~= -1) && (li_at_controls_min_valid_Nvox ~= -999),
		if li_at_controls_min_valid_Nvox < 0 
			text_placement1 = 0.7;
			text_placement2 = 0.5;
			text_placement2sub = 0.48;
		else
			text_placement1 = -0.5;
			text_placement2 = -0.7;		
			text_placement2sub = -0.72;
		end	
		annotation1 = text(controls_min_valid_Nvox,text_placement1,typical_string,'fontsize',1.5 * title_fontsize, 'FontWeight','Bold'); 
		%annotation2 = text(controls_min_valid_Nvox,text_placement2,li_string,'fontsize',1.5 * title_fontsize, 'FontWeight','Bold'); 
		annotation2a = text(controls_min_valid_Nvox,text_placement2,[li_string01],'fontsize',1.5 * title_fontsize, 'FontWeight','Bold'); % LI = ...
		annotation2b = text(controls_min_valid_Nvox,text_placement2sub,[li_string01space li_string01space li_string02],'fontsize',0.75 * title_fontsize, 'FontWeight','Bold'); % (@N=...)
	end
	if iBT.what.li.SummaryPlot.Box==0, box off; else box on; end
	
    	if iBT.what.li.save.graph == 1 
		try full_fileIDstring = iBT.what.li.FullPlot.fileIDstring; % use override value if specified.
		catch, full_fileIDstring = fileIDstring; end
		try summary_fileIDstring = iBT.what.li.SummaryPlot.fileIDstring; % use override value if specified.
		catch, summary_fileIDstring = fileIDstring; end
	   	% print(fignum,'-dpng',['li_plot-' active.tasks])
		if num_controls ~= 0, saveas(full_fignum, fullfile(iBT.what.li.save.path, ['li_full_plot_' full_fileIDstring '_' fileBaseName{con}]), 'png'); end;
	   	saveas(fignum, fullfile(iBT.what.li.save.path, ['li_plot_' summary_fileIDstring '_' fileBaseName{con}]), 'png');
    	end; %if iBT.what.li.save.graph == 1 

% Restore warning state:  
warning(warningjbtest.state, warningjbtest.identifier);

return
% End of iBT_laterality_graph

function h = sfigure(h)
% SFIGURE  Create figure window (minus annoying focus-theft).
%
% Usage is identical to figure.
%
% Daniel Eaton, 2005
%% See also "help figure"
if nargin>=1
if ishandle(h)
    set(0, 'CurrentFigure', h);
else
    h = figure(h);
end
else
    h = figure;
end
% End of function sfigure
