function iBT_contrasts(iBT)
% Add extra contrasts to statistical design.
% This script is designed to be called by iBT_analysis 
% FORMAT iBT_cons(iBT)
%  where 'iBT.what' is a structure indicating what is to be done,
%    and 'iBT.who' is a structure indicating subject specific details.
%___________________________________________________________________________
% Copyright 2008-2016 The Florey Institute of Neuroscience and Mental Health
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
% 2016-05-04: (dfa) Log Matlab Version
% 2012-03-12: (dfa) Display version info from iBT_version()
% 2011-08-26: (dfa) Updated for public release
%

% 2004-2011: iBT analysis scripts were developed by David Abbott based
%            upon the who & what structure concept of Tony Waites, with 
%            substantial contributions from Matt Harvey, Kaushik Bhaganagarapu,
%            and Richard Masterton.
%            In 2008 the script that became iBT_contrasts was created
%            by David Abbott.

disp(sprintf('Entered cons script...'));

if isempty(iBT), error('Error: "iBT" structure is not defined.'), end

[iBTversion,iBTrelease,iBTversionReleaseString] = iBT_version(1); % Display version so it gets logged.
[iBT.what.SpmVersion, iBT.what.SpmRelease, iBT.what.SpmVersionReleaseString, MatlabVersion, MatlabDate, MatlabArch] = iBT_SPMversion(2);

spm('defaults','fmri'); % Ensure we are starting with standard defaults

set(gcf,'DefaultTextColor','white');


% Add contrasts to an existing analysis and evaluate.

start = iBT.what.cons.start; 	% Number of the first new contrast we are defining. If this
		% does not match the next available contrast number in the SPM.mat file,
		% then an error is printed and nothing is done.
		% (Useful to prevent repeated executions of this script adding the same contrasts)

n_ev=iBT.what.cons.n_ev;   % number of events (columns) of interest
n_conf=iBT.what.cons.n_conf; % Number of confounds (events of no interest e.g. motion correction parameters)


% Define new contrast(s) (will pad automatically with zeros for confounds and constant)
cnam = iBT.what.cons.cnam;
contrast = iBT.what.cons.contrast;
ctyp = iBT.what.cons.ctyp;



nsub = length(iBT.who.sub);
awd = iBT.what.where.awd;

for j = 1:1:length(iBT.what.stats.offsets) 

for sub = 1:nsub
	nsess = length(iBT.who.sub{sub}.sess);
        cd(awd)
	for sess = 1:nsess

cd(iBT.what.stats.dir); % Set in iBT_start.m and possibly adjusted for offset in iBT_analysis.m.


load SPM;
ocon = length(SPM.xCon);			% number of contrasts already present

if ( ocon + 1 ) == start

  % add new contrasts to SPM.xCon
  for c = 1:length(cnam)
	cwgt{c} = [contrast{c} zeros(1, n_conf+1)];
	cw = [cwgt{c}]';
	SPM.xCon(c+ocon)     = spm_FcUtil('Set',cnam{c},ctyp{c},'c',cw,SPM.xX.xKXs);
  end
 
  % And evaluate
  %===========================================================================

  switch (iBT.what.SpmVersion),
    case 'SPM2'
      % Make sure we've already estimated the model.
      check_files = spm_get('Files',cd, 'beta_0001.hdr'); 
    otherwise
      % Make sure we've already estimated the model.
      check_files = spm_select('List',cd, spm_wildconvert('beta_0001.hdr')); 
  end

  if length(check_files) > 0
    switch (iBT.what.SpmVersion),
      case 'SPM2'
  	spm_contrasts(SPM);
      otherwise
        Ic = length(SPM.xCon);
        spm_contrasts(SPM, [ocon:Ic]);
    end
  else
	  disp ' ' 
	  disp '***********************************************************************' 
 	  disp(sprintf('WARNING: No betas present - it looks like the model has not'))
 	  disp(sprintf('         yet been estimated. Unable to configure contrasts. '))
	  disp '***********************************************************************' 
	  disp ' ' 

  end; % if length(check_files)

else
	  disp ' ' 
	  disp '***********************************************************************' 
 	  disp(sprintf('WARNING: Specified contrast number (start = %i) does not match ', start))
 	  disp(sprintf('         next available contrast (%i). No contrasts added. ',ocon + 1))
 	  disp(sprintf('         Please check existing contrasts, and start value. '))
	  disp '***********************************************************************' 
	  disp ' ' 


end; %if ocon + 1 ~= start


	end; %for sess = 1:nsess
        cd(awd)
end
end


function new_str = spm_wildconvert(orig_str)
%
% Sub function to convert normal wildcards to a regexp format for spm8
%
% Modification History
% ---------------------
%
% 2008-08-18 - matt: Function creation

new_str =  strcat('^', orig_str, '$');
new_str = strrep(new_str,'.','\.');
new_str = strrep(new_str,'?','.');
new_str = strrep(new_str,'*','(.*)');

