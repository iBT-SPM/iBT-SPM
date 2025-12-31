function [l,r,li,nvox,thresh,maxt]= iBT_LI_thresh(left, right, act, step) 
% Calculate li (and nvox etc.) for a range of thresholds with ganularity = step
%
% FORMAT [l,r,li,nvox,thresh,maxt]= iBT_LI_thresh(left, right, act, step)
% where
%   left = matrix of voxels for the left ROI
%  right = matrix of voxels for the right ROI
%   act = matrix of voxels of the activation map (e.g. loded from spmT*.img)
%  step = the granularity of thresh for which to generate results
% 
% This script is designed to be called by iBT_analysis 
% FORMAT iBT_laterality(iBT)
%  where 'iBT.what' is a structure indicating what is to be done,
%    and 'iBT.who' is a structure indicating subject specific details.
%___________________________________________________________________________
% Copyright 2006,2010-2011 The Florey Institute of Neuroscience and Mental Health
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
% Version History:
% 2011-08-26: (dfa) Updated for public release in the 
%			Integrated Brain Analysis Toolbox for SPM
% 2010-08-09: (dfa) supported images where there are no positive voxels 
%			at any threshold..
% 2009-11-09: (dfa) added maxt return.
% 2006: Original version by Tony Waites.
%
%
%first load the currently used ROIs for left and right activation
pwdorig = pwd;
 
 if nargin < 1,
    left     = spm_get(1,'*.img','load left roi');
    V2=spm_vol(left);
    left       = spm_read_vols(V2);
 end
 if nargin < 2,
    right     = spm_get(1,'*.img','load right roi');
    V3=spm_vol(right);
    right       = spm_read_vols(V3);
 end
%now load the spmT files for each analysis
if nargin < 3,
    
        P     = spm_get(1,'*.img','load spmT volume');
        V=spm_vol(P);
        act      = spm_read_vols(V) ;   
end 
if nargin < 4,
    
      
 step=spm_input('step size in # voxels?',1, 'e', 10);   
end 

all=(left+right)>0;
nl=sum(sum(sum(left)));
nr=sum(sum(sum(right)));
nall=sum(sum(sum(all)));


 maxt=max(max(max(act(all>0))))   % this is the max t that will be activated

if maxt > 0
 act_roi=act(find(all > 0));
% sort_act_roi=sort(act_roi,'descend');
 
 nsteps=ceil(maxt/step);
 for n=1:nsteps,
    thresh(n)=n*step;
    
    [l(n),r(n),li(n)]= iBT_LI(left, right, act,thresh(n));
    nvox(n)=l(n)+r(n);
 end
else
	l(1)=0;
	r(1)=0;
	li(1)=0;
	nvox(1)=0;
	thresh(1)=0;
end % if maxt > 0
%
return;
 
 
 
