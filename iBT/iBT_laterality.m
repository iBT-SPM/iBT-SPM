function  iBT_laterality(iBT)
% Perform analysis of regional laterality over a range of thresholds
% This script is designed to be called by iBT_analysis 
% FORMAT iBT_laterality(iBT)
%  where 'iBT.what' is a structure indicating what is to be done,
%    and 'iBT.who' is a structure indicating subject specific details.
%___________________________________________________________________________
% Copyright 2009-2013,2024 The Florey Institute 
%                          of Neuroscience and Mental Health
%
% This file is part of iBT (the Integrated Brain Analysis Toolbox for SPM).
% See iBT_SPM.m for more information.
%
% Please refer to the following paper for a description of the
% methods implemented in this toolbox. If you publish results using
% these methods or variations of them, please remember to cite this paper: 
%     Abbott DF, Waites AB, Lillywhite LM, Jackson GD.
%     fMRI assessment of language lateralization: An objective approach. 
%     Neuroimage 50(4):1446-1455 (2010).

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
% Version History:
% 2024-10-18: (dfa) Bugfix: Prevent crash if there are no control data 
%                           files inside the control data folder.
% 2024-04-28: (dfa) Bugfix: Prevent crash if saving to CSV file if loaded
%                           previously saved results from a matfile.
% 2013-08-12: (dfa) Don't try to make destination folder if it already exists.
% 2013-08-07: (dfa) Abort if input images and ROIs not all in same orientation
% 2013-06-01: (dfa) Change x-axis label to 'Number of voxels above threshold'
% 2012-03-12: (dfa) Display version info from iBT_version()
% 2011-08-26: (dfa) Updated for public release in the 
%			Integrated Brain Analysis Toolbox for SPM
% 2011-01-13: (dfa) Removed obsolete 'versn' from call to fileparts
% 2011-01-12: (dfa) Use new calling procedure to set SPM defaults as
%                   recommended in spm8_r4010 update.
% 2010-08-09: (dfa) Added support for what.li.reuse_existing
% 2010-03-23: First public release in the iBrain Laterality Toolbox.
%

%
% 2004-2011: iBT analysis scripts were developed by David Abbott based
%            upon the who & what structure concept of Tony Waites, with 
%            substantial contributions from Matt Harvey, Kaushik Bhaganagarapu,
%            and Richard Masterton.
%            In 2009 the script that became iBT_laterality was 
%            created by David Abbott.
%

disp(sprintf('Entered Laterality Index script...'));

[iBTversion,iBTrelease,iBTversionReleaseString] = iBT_version(1); % Display version so it gets logged.

nvox_step = iBT.what.li.nvox_step;
thresh_step = iBT.what.li.thresh_step;
try iBT.what.li.save.path = iBT.what.li.save.path; catch iBT.what.li.save.path=''; end
try iBT.what.li.maxNvox = iBT.what.li.maxNvox; catch iBT.what.li.maxNvox=0; end % 0 = unlimited.
try iBT.what.li.maxT = iBT.what.li.maxT; catch iBT.what.li.maxT=0; end % 0 = unlimited.
try iBT.what.li.reuse_existing = iBT.what.li.reuse_existing; catch iBT.what.li.reuse_existing=0; end % 0 = Default is to not re-use existing li_data.

if (~strcmp(iBT.what.li.save.path, '')) & (  (iBT.what.li.save.mat == 1) | (iBT.what.li.save.text == 1) | (iBT.what.li.save.graph == 1) )
	if ~exist(iBT.what.li.save.path,'dir'), mkdir(iBT.what.li.save.path); end
end


pwd_orig=pwd;

try
    [iBT.what.SpmVersion,iBT.what.SpmRelease]   = spm('Ver','',1); iBT.what.SpmRelease = str2double(iBT.what.SpmRelease);
catch
    error('SPM cannot be found in MATLAB path.');
end

if  isnan(iBT.what.SpmRelease)
	iBT.what.SpmVersionReleaseString = iBT.what.SpmVersion;
else
	iBT.what.SpmVersionReleaseString = sprintf('%s (release %i)',iBT.what.SpmVersion,iBT.what.SpmRelease);
end
disp(sprintf('SPM version %s',iBT.what.SpmVersionReleaseString))

if isempty(who), error('Error: "who" structure is not defined.'), end

spm('defaults','fmri'); % Ensure we are starting with standard defaults

format compact

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
%  loading ROIs
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


wild1=iBT.what.li.wild;
roi_dir=iBT.what.li.roi_dir;
cd(roi_dir);
try iBT.what.li.do_tgraph = iBT.what.li.do_tgraph; catch iBT.what.li.do_tgraph = 0; end
	  
%load rois  (individual, left and right)

switch (iBT.what.SpmVersion),
          case 'SPM2',
		l_allw    = spm_get('files',iBT.what.li.roi_dir,iBT.what.li.left_roi);  
		r_allw    = spm_get('files',iBT.what.li.roi_dir,iBT.what.li.right_roi);  
          otherwise
            	l_allw = spm_select('FPList',iBT.what.li.roi_dir, ...
             		spm_wildconvert(iBT.what.li.left_roi))
            	r_allw = spm_select('FPList',iBT.what.li.roi_dir, ...
             		spm_wildconvert(iBT.what.li.right_roi))
end


VL=spm_vol(l_allw);
        l_all      = spm_read_vols(VL); 
        
VR=spm_vol(r_allw);
        r_all      = spm_read_vols(VR); 

             
all=l_all+r_all;

%for each subject, load act maps

%calc LI, thresh and number of vox 


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
% calc and write to data structure
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    cd(pwd_orig); % Back to where we started in case iBT.what.stats.dir is a relative path

    switch (iBT.what.SpmVersion),
          case 'SPM2',
		v_act = spm_get('files',iBT.what.stats.dir , wild1);
          otherwise
            	v_act = spm_select('FPList',iBT.what.stats.dir, ...
             		spm_wildconvert(wild1))
    end

ncon=size(v_act,1);

if ncon == 0
   disp(['No contrasts found matching ' fullfile(iBT.what.stats.dir,wild1) ]);
   return
end
 
cd(iBT.what.stats.dir); % Set in iBT_start.m and possibly adjusted for offset in iBT_analysis.m.

%String to embed in output filenames:
fileIDstring = sprintf('%s_%02d_%s%s',iBT.who.sub{iBT.what.do_sub}.ID, iBT.what.disp.ses, iBT.what.stats.tasks, iBT.what.stats.offset_text); 
mat_filename=fullfile( iBT.what.li.save.path, sprintf('li_active_%s.mat',fileIDstring) );

%We create a structure 'active' in which results are saved.
active.structure_version = 1.0; %Version number of the active structure. This field was added on 2009-10-20

if iBT.what.li.reuse_existing == 1
    try
	previous = load(mat_filename);
	if active.structure_version > previous.active.structure_version
		calculate_new = 1;
		clear previous
	else
		calculate_new = 0;
		clear active
		active = previous.active;
	end
    catch
	calculate_new = 1;
    end
else
    calculate_new = 1;
end % if iBT.what.li.reuse_existing ==1

if calculate_new == 1	
	active.subject_ID = iBT.who.sub{iBT.what.do_sub}.ID;
	active.session = iBT.what.disp.ses; % Session number
	active.tasks = iBT.what.stats.tasks; % This is not ideal for some situations, but at 
	   % least gives us a clue as to what task was analysed. Ideally should obtain 
	   % contrast label.
	active.offset_text = iBT.what.stats.offset_text;%If we've shifted the HRF in time
	active.stats_dir = iBT.what.stats.dir; % Location of the analysed data
	active.left_ROI_fname = l_allw;
	active.right_ROI_fname = r_allw;
	active.nvox_step = nvox_step;
	active.thresh_step = thresh_step;
else
	disp(['Successfully read LI data from ',mat_filename]);
end %if calculate_new == 1

active  % Display the current contents of the active structure

for con=1:ncon,
    
    V=spm_vol(v_act(con,:));
    %Check that the con image is in exactly the same orientation as the ROIs
    spm_check_orientations([V VL VR]);
    % If get to here, orientation check must have passed.
    act = spm_read_vols(V); 
    active.fname{con} = V.fname;
    disp(['Processing ' V.fname ]);
    [PATHSTR,NAME,EXT] = fileparts( active.fname{con} );
    fileBaseName{con} = NAME;
    	
    if calculate_new == 1	
      %calc li as a function of nvox
      [l,r,li,nvox,thresh,maxn]= iBT_LI_nvox(l_all, r_all, act, nvox_step);
      active.nvox{con}.l=l;
      active.nvox{con}.r=r;
      active.nvox{con}.nvox=nvox;
      active.nvox{con}.thresh=thresh;
      active.nvox{con}.li=li;
      active.nvox{con}.maxn=maxn;
    else
      l=active.nvox{con}.l;
      r=active.nvox{con}.r;
      nvox=active.nvox{con}.nvox;
      thresh=active.nvox{con}.thresh;
      li=active.nvox{con}.li;
      maxn=active.nvox{con}.maxn;
    end %if calculate_new == 1

    if iBT.what.li.save.text==1
	%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
	%
	% write data structure to excel readable file
	%
	%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
 	script_filename=fullfile( iBT.what.li.save.path, sprintf('li_nvox_%s_%s.csv', fileIDstring, fileBaseName{con}));
    	script_file =  fopen( script_filename , 'w' );
    	fprintf(script_file,'Laterality index analysis by number of active voxels.\n');
    	fprintf(script_file,'nvox \t LI \t thresh \t l \t r \n');
   	for i=1:length(l)
            fprintf(script_file, '%8.4f \t %8.4f \t %8.4f \t %8.4f \t %8.4f \n',nvox(i), li(i), thresh(i), l(i), r(i));
   	end % for i
  	fclose(script_file );
    end %if iBT.what.li.save.text==1
  
    if calculate_new == 1
      %calc li as a function of thresh
      [l,r,li,nvox,thresh,maxt]= iBT_LI_thresh(l_all, r_all, act, thresh_step); 
      active.thresh{con}.l=l;
      active.thresh{con}.r=r;
      active.thresh{con}.nvox=nvox;
      active.thresh{con}.thresh=thresh;    
      active.thresh{con}.li=li;   
      active.thresh{con}.maxt=maxt;   
    else
      l=active.thresh{con}.l;
      r=active.thresh{con}.r;
      nvox=active.thresh{con}.nvox;
      thresh=active.thresh{con}.thresh;
      li=active.thresh{con}.li;
      maxt=active.thresh{con}.maxt;
    end %if calculate_new == 1

    if iBT.what.li.save.text==1
	%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
	%
	% write data structure to excel readable file
	%
	%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
 	script_filename=fullfile( iBT.what.li.save.path, sprintf('li_thresh_%s_%s.csv', fileIDstring, fileBaseName{con}));
    	script_file =  fopen( script_filename , 'w' );
    	fprintf(script_file,'Laterality index analysis by threshold.\n');
    	fprintf(script_file,'nvox \t LI \t thresh \t l \t r \n');
   	for i=1:length(l)
            fprintf(script_file, '%8.4f \t %8.4f \t %8.4f \t %8.4f \t %8.4f \n',nvox(i), li(i), thresh(i), l(i), r(i));
   	end % for i
  	fclose(script_file );
    end %if iBT.what.li.save.text==1
	  
end; %for con = 1:ncon


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
% save data structure with all contrasts
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%save relevant data to be reloaded for subsequent stats        
if ( iBT.what.li.save.mat==1 ) && ( calculate_new == 1 )
 	mat_filename=fullfile( iBT.what.li.save.path, sprintf('li_active_%s.mat',fileIDstring) );
        save (mat_filename, 'active');
	disp(['Saved ' mat_filename]);
end %if iBT.what.li.save_mat==1



if iBT.what.li.do_graph == 1
   li_c_nvox = {}; % initialise so variable exists even if no control data
   li_c_thresh = {};  % initialise so variable exists even if no control data
   if ~strcmp(iBT.what.li.control_dir,'')
	    control_files = dir( fullfile(iBT.what.li.control_dir, iBT.what.li.control_wild) );
	    num_controls = size(control_files,1);
      	for i=1:num_controls % load the li control data
		  li_c{i}  = load(fullfile(iBT.what.li.control_dir,control_files(i).name));
		  li_c_nvox{i} = li_c{i}.active.nvox; % Rejig so can pass to iBT_laterality_graph
		  li_c_thresh{i} = li_c{i}.active.thresh; % Rejig so can pass to iBT_laterality_graph
     	end %for
   end % if ~strcmp
    
    try
    	  control_index_mapping = iBT.what.li.control_index_mapping; % Index into control data for this contrast
    catch
	  iBT.what.li.control_index_mapping=[1:ncon]; % Default if not specified is 1 to 1 mapping;
    end

    for con=1:ncon,
	if (ncon > 1) 
          if (con == 1) con_deactivation = 2; elseif (con == 2) con_deactivation = 1; end 
	else 
          con_deactivation = 0; 
	end % con_deactivation used for sanity check
        fignums = [1:2] + (100*con);
        iBT_laterality_graph(iBT,con,con_deactivation,li_c_nvox, active.nvox, active.nvox_step, active.nvox{con}.nvox, iBT.what.li.maxNvox, 'Number of voxels above threshold', [fileIDstring '_nvox'], fileBaseName,fignums);
      	if iBT.what.li.do_tgraph == 1 
	  if (ncon > 1) && (con == 1) con_deactivation = 2; else con_deactivation = 0; end % con_deactivation used for sanity check
	  fignums = [3:4] + (100*con);
	  iBT_laterality_graph(iBT,con,con_deactivation,li_c_thresh, active.thresh, active.thresh_step, active.thresh{con}.thresh, iBT.what.li.maxT, 't-threshold', [fileIDstring '_t-thresh'], fileBaseName,fignums);
	end
      
    end; % for con
    
 
end % if iBT.what.li.do_graph == 1


% back to initial directory
cd(pwd_orig);

return
% End of iBT_laterality



function new_str = spm_wildconvert(orig_str)
%
% Sub function to convert normal wildcards to a regexp format for spm*
%
% Modification History
% ---------------------
%
% 2008-08-21 - matt: Function creation

new_str =  strcat('^', orig_str, '$');
new_str = strrep(new_str,'.','\.');
new_str = strrep(new_str,'?','.');
new_str = strrep(new_str,'*','(.*)');
