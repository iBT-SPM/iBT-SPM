function [l,r,li]= iBT_LI(left, right, spmt, thresh)
% Determine number of activated voxels in left and right ROIs and ...
% the laterality index at a single specified threshold.
% FORMAT [l,r,li]= iBT_LI(left, right, spmt, thresh)
% where
%   left = matrix of voxels for the left ROI
%  right = matrix of voxels for the right ROI
%   spmt = matrix of voxels of the activation map (e.g. loded from spmT*.img)
%  thresh = the overlaping voxels between activation maps and the 
%           combined ROI (i.e. left + right)
%
%___________________________________________________________________________
% Copyright 2006,2010-2011 The Florey Institute of Neuroscience and Mental Health
%
% This file is part of iBT (the Integrated Brain Analysis Toolbox for SPM).
% See iBT_laterality.m for more information.

% Please refer to the following paper for a description of the laterality
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
% 2010-03-23: First public release in the iBrain Laterality Toolbox.
%

%
% 2004-2011: iBT analysis scripts were developed by David Abbott based
%            upon the who & what structure concept of Tony Waites, with 
%            substantial contributions from Matt Harvey, Kaushik Bhaganagarapu,
%            and Richard Masterton.
%            In 2006 the script that became iBT_LI was 
%            created by Tony Waites.
%


 
%first load the currently used ROIs for left and right activation
pwdorig = pwd;
 
if nargin < 1
    left = spm_get(1,'*.img','load left roi');
    V2 = spm_vol(left);
    left = spm_read_vols(V2);
end

if nargin < 2
    right = spm_get(1,'*.img','load right roi');
    V3 = spm_vol(right);
    right = spm_read_vols(V3);
end

%now load the spmT files for each analysis
if nargin < 3
        P = spm_get(1,'*.img','load spmT volume');
        V = spm_vol(P);
        act = spm_read_vols(V)    
else    
        act = spmt;    
end 

li=0;

thr = thresh;
ltemp = ((left>0).*(act>thr));
l = sum(sum(sum(ltemp)));
rtemp = ((right>0).*(act>thr));
r = sum(sum(sum(rtemp)));
lr=l+r;
if lr == 0, li = NaN; else li = (l-r)/lr; end
           
    
cd(pwdorig);
  
end
