function K = iBT_spm_filter(K)
%
% Creates a band-pass filter for connectivity analyses in SPM
%
% FORMAT K  = iBT_spm_filter(K)
%
%
% INPUT:
% 
%   K(s)        - struct array containing partition-specific specifications
%
%       K(s).RT     - observation interval in seconds
%       K(s).row    - row of Y constituting block/partition s
%       K(s).HParam - high-pass cut-off period in seconds
%       K(s).LParam - low-pass cut-off period in seconds
%
%
% OUTPUT: The following field of the structure is filled in:
%
%   K(s).X0     - low frequencies to be removed (DCT)
%
% USAGE:
%
%   This structure should be assigned to the SPM.xX.K field of the SPM structure.
%
%   Example:
%
%    K.RT = SPM.xY.RT;
%    K.HParam = 128;         % high-pass cut-off (in seconds)
%    K.LParam = 20;          % low-pass cut-off (in seconds)
%    K.row = [1:SPM.nscan];  % the partition of the design matrix 
%
%    K = iBT_spm_filter(K);
%    SPM.xX.K =  K;
%
% NOTE:
%
%   SPM also populates this structure, so the following code
%   __must__ be inserted __after__ the following line in an SPM analysis script:
%
%   SPM = spm_fmri_spm_ui(SPM);
%
%___________________________________________________________________________
% Copyright 2009,2011 The Florey Institute of Neuroscience and Mental Health
% Richard Masterton
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
% Version History:
% 2011-08-26: (dfa) Updated comments for public release
% Code by Richard Masterton 2009
%
% 2004-2011: iBT analysis scripts were developed by David Abbott based
%            upon the who & what structure concept of Tony Waites, with 
%            substantial contributions from Matt Harvey, Kaushik Bhaganagarapu,
%            and Richard Masterton. The script that became
%            iBT_spm_filter was written by Richard Masterton in 2009, 
%            based upon spm_filter.m v2.10 by Karl Friston 03/03/04 (in spm2).
%

	% set K.X0
	%-------------------------------------------------------------------
	for s = 1:length(K)

        try
            K(s).LParam;
        catch
            K(s).LParam = K(s).RT;
        end
        
		k       = length(K(s).row);
		n       = fix(2*(k*K(s).RT)/K(s).HParam + 1);
		N       = ceil(2*(k*K(s).RT)/K(s).LParam + 1);
		X0      = spm_dctmtx(k,k+1);
		K(s).X0 = X0(:,2:n);
        if N < (k+1)
            K(s).X0 = [K(s).X0 X0(:,N:end)];
        end
	end

end
