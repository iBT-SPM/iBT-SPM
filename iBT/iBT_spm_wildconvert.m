function new_str = iBT_spm_wildconvert(orig_str)
%
% Function to convert normal wildcards to a regexp format for spm
%
% Modification History
% ---------------------
%
% 2013-09-20 - dfa: Moved to separate file
% 2008-08-18 - matt: Function creation

new_str =  strcat('^', orig_str, '$');
new_str = strrep(new_str,'.','\.');
new_str = strrep(new_str,'?','.');
new_str = strrep(new_str,'*','(.*)');

end % of function iBT_spm_wildconvert
