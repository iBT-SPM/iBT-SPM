function iBT_analysis(iBT)
% Called by iBT_start to manage calls to other scripts.
% FORMAT iBT_analysis(iBT)
%  where 'iBT.what' is a structure indicating what is to be done,
%    and 'iBT.who' is a structure indicating subject specific details.
%
%___________________________________________________________________________
% Copyright 2009-2016,2020,2023 The Florey Institute 
%                               of Neuroscience and Mental Health
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
% 2023-11-24: (dfa) Ensure laterality is called if iBT.what.li.li_only == 1
%                   regardless of iBT.what.stats.norm (to fixed change on
%                   2020-04-19 which failed for runinfo_li_template.m) 
% 2020-04-19: (dfa) Support iBT.what.li.own.do & iBT.what.li.norm.do
% 2016-05-04: (dfa) Log Matlab Version
% 2013-10-18: (dfa) Support iBT.what.disp.OTHER.do (replaces .RFX.do)
% 2013-08-12: (dfa) Support iBT.what.custom_stats.do
% 2012-04-28: (dfa) Added a comment.
% 2012-03-12: (dfa) Display version info from iBT_version()
% 2012-02-27: (dfa) Separated RFX section from 1st level analysis section.
% 2011-08-26: (dfa) Updated for public release
%

% 2004-2011: iBT analysis scripts were developed by David Abbott based
%            upon the who & what structure concept of Tony Waites, with 
%            substantial contributions from Matt Harvey, Kaushik Bhaganagarapu,
%            and Richard Masterton.
%            In 2009 the script that became iBT_analysis was 
%            created by David Abbott.
%

disp(sprintf('Entered analysis script...'));

if isempty(iBT), error('Error: "iBT" structure is not defined.'), end

[iBTversion,iBTrelease,iBTversionReleaseString] = iBT_version(1); % Display version so it gets logged.
[iBT.what.SpmVersion, iBT.what.SpmRelease, iBT.what.SpmVersionReleaseString, MatlabVersion, MatlabDate, MatlabArch] = iBT_SPMversion(2);

spm('defaults','fmri'); % Ensure we are starting with standard defaults

format compact

% Maintain compatibility with old runinfo files that may not have 
% some new fields that are used in this function.
try, iBT.what.stats.offsets; catch iBT.what.stats.offsets = 0; end 
try, iBT.what.stats.fc_flag; catch iBT.what.stats.fc_flag = 0; end 
try, iBT.what.custom_stats.do; catch iBT.what.custom_stats.do = 0; end 
try, iBT.what.stats.norm; catch iBT.what.stats.norm = 0; end
try, iBT.what.stats.own; catch iBT.what.stats.own = 0; end

try, global_li_flag=iBT.what.li.do; catch global_li_flag = 0; end % We previoulsy could not distinguish between own and norm space when requesting LI analysis.
try, iBT.what.li.own.do; catch iBT.what.li.own.do = global_li_flag; end % Maintain backward compatibility
try, iBT.what.li.norm.do; catch iBT.what.li.norm.do = global_li_flag; end % Maintain backward compatibility
iBT.what.li.do = iBT.what.li.norm.do | iBT.what.li.own.do; 
try, iBT.what.disp.MELODIC.do; catch iBT.what.disp.MELODIC.do = 0; end 
try, iBT.what.disp.OTHER.do;
catch 
	try, iBT.what.disp.OTHER.do = 2 * (iBT.what.disp.RFX.do > 0); % use legacy variable if present
	catch;
		iBT.what.disp.OTHER.do = 0;
	end
end
try, iBT.what.disp.SPM.do; catch iBT.what.disp.SPM.do = (iBT.what.disp.OTHER.do == 0); end  

try, iBT.what.li.li_only; catch, 
	% Work out if not specified
	iBT.what.li.li_only = (iBT.what.li.do~=0) && ( (iBT.what.pre.do ==0) && (iBT.what.stats.do ==0) && (iBT.what.cons.do ==0) && (iBT.what.custom_stats.do ==0)  && (iBT.what.disp.do ==0) )
end 

noffsets = length(iBT.what.stats.offsets);
if (iBT.what.stats.fc_flag == 0),
  Nanalysis = noffsets;
else
  Nanalysis = 1;
end
 
stats_dir = iBT.what.stats.dir;

for(j = 1:1:Nanalysis) % Create a separate analysis for each event timing adjustment (offset)

 iBT.what.stats.offset = iBT.what.stats.offsets(j); % Analyse with timings adjusted by this offset

 % Create an appropriate name for the analysis folder
 iBT.what.stats.offset_text = ''; % Initialise
 if (iBT.what.stats.fc_flag == 0),
  if iBT.what.stats.offset > 0
    iBT.what.stats.offset_text = strcat('_+', num2str(iBT.what.stats.offset));
  else 
    if iBT.what.stats.offset < 0
      iBT.what.stats.offset_text = strcat('_', num2str(iBT.what.stats.offset));
    else
      iBT.what.stats.dir = stats_dir;
      try
        if iBT.what.stats.subdirZeroOffset == 1
          iBT.what.stats.offset_text = strcat('_+', num2str(iBT.what.stats.offset));
        end
      catch
        iBT.what.stats.dir = stats_dir;
      end; % try     
    end
  end
 end
 if ~strcmp(iBT.what.stats.offset_text, '') 
 	iBT.what.stats.dir = strcat(stats_dir,iBT.what.stats.offset_text);
 end

 if (iBT.what.stats.do == 1)
    iBT_stats(iBT)
 end

 if (iBT.what.cons.do == 1)
   iBT_contrasts(iBT)
 end

 if (iBT.what.custom_stats.do == 1)
    iBT_custom_stats(iBT)
 end

 if ( (iBT.what.stats.norm == 0) | (iBT.what.li.li_only == 1) ) && (iBT.what.li.own.do == 1),
   iBT_laterality(iBT)
 end

 if ((iBT.what.stats.norm == 1) | (iBT.what.li.li_only == 1) ) && (iBT.what.li.norm.do == 1),
   iBT_laterality(iBT)
 end

 if (iBT.what.disp.do == 1)
   % iBT.what.disp.type used to be set in runinfo, but now we determine it automatically
   % here which permits multiple types to be procesed from a single runinfo.
   if (iBT.what.disp.SPM.do == 1)  
     if iBT.what.stats.norm == 0; iBT.what.disp.type = 'own'; else iBT.what.disp.type = 'norm'; end; % iBT.what.stats.norm is set in iBT_start.m
      iBT_display(iBT)
   end
   if (iBT.what.disp.OTHER.do > 0)  
     if iBT.what.stats.norm == 0; iBT.what.disp.type = 'own'; else iBT.what.disp.type = 'norm'; end % iBT.what.stats.norm is set in iBT_start.m
     % iBT.what.disp.analysis is set in iBT_start.m according to the analyis being undertaken,
     % however if we want to display a manual analysis then we override this.
     iBT.what.disp.analysis = 3;
     iBT_display(iBT)
   end 
   if (iBT.what.disp.MELODIC.do == 1) 
     if iBT.what.MELODIC.norm_space.do == 0; iBT.what.disp.type = 'own'; else iBT.what.disp.type = 'norm'; end
     iBT.what.disp.analysis = 5 
     iBT_display(iBT)
   end 
 end

end % for (j = 1:1:Nanalysis)

iBT.what.stats.dir = stats_dir; % Reset to original state

return
