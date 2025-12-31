%
% Start the GUI for the Integrated Brain Analysis Toolbox for SPM
%
% FORMAT iBT_SPM
%
% This is an SPM toolbox to automate processing, including
% optionally accessing features of the iBrain software suite.
%
% The file README_iBT.txt is displayed when starting the toolbox
% from the Toolbox menu of the SPM user interface. Please see this for
% instructions, acknowledgements and citation details.

%___________________________________________________________________________
% Copyright 2006, 2010-2014 The Florey Institute of Neuroscience and Mental Health
%
% This file is part of iBT (the Integrated Brain Analysis Toolbox for SPM).
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
% Recent Version History:
% 2014-07-32: (dfa) Increase window size to mitigate menu focus problems on
%                   systems with a buggy focus-follows-mouse implementation.
% 2013-10-18: (dfa) Support runinfo_other_display
% 2012-03-12: (dfa) Obtain version info from iBT_version()
% 2011-09-01: (dfa) Updated for public release
%

function varargout = iBT_SPM(varargin)

    [iBTversion,iBTrelease, iBTversionReleaseString,iBTversionReleaseSentence,iBTmenuTitle] = iBT_version(0);

    iBTtag = 'spm_iBT';
    
    %
    %  Handle closing call-back 
    %
    
    if (nargin >= 1) & strcmp(varargin{1}, 'CLOSE')
        FiBT = findobj('Tag',iBTtag);
        close(FiBT);
        
        Finter=spm_figure('GetWin','Interactive');
        spm_input('!DeleteInputObj');
        spm_figure('Clear',Finter)

        Fgraph =spm_figure('GetWin','Graphics');
        spm_figure('Clear',Fgraph)
        
        return;
    end
    
    %
    %  Don't do anything if the window is already open
    %
    
    if length(findobj('Tag',iBTtag))
    
        return;
    
    end
    
    % 
    %  Display user guide in the Graphics window
    %
    
    spm_help('!Disp','README_iBT.txt','','Graphics');

    %
    %  Clear the interactive window
    %
    
    Finter=spm_figure('GetWin','Interactive');
    spm_input('!DeleteInputObj');
    spm_figure('Clear',Finter)
    
    
    %
    %  Create the toolbox window
    %
    
    Pos = get(Finter,'Position');
    
%    FiBT = figure('Position',[Pos(1) Pos(2)+Pos(4) Pos(3) 1], ...
% The above should have been fine (1 pixel figure window height) but on some systems
%  with focus follows mouse the window focus can be lost when moving on to the menu
%  if there is not the parent window under it (this can be a problem for all apps
%  including even selecting the toolbox from the SPM toolbox menu). We work around
%  it for the iBT_SPM.menu by using a larger figure window):
    FiBT = figure('Position',[Pos(1) Pos(2)+Pos(4) Pos(3) 0.8*Pos(4)], ...
                  'MenuBar','none', ...
                  'Resize','off', ...
                  'Name',iBTmenuTitle, ...
                  'NumberTitle','off', ...
                  'Tag',iBTtag);

    %
    %  The file menu
    %
    
    MFile = uimenu(FiBT,...
                  'Label','File',...
                  'HandleVisibility','off');

    iBT_runinfo = dir(fullfile(spm('dir'),'toolbox','iBT','runinfo*template.m'));
    lir_j = 1;
    odr_j = 1;
    for i = 1:numel(iBT_runinfo)
    	if strncmpi(iBT_runinfo(i).name,'runinfo_li',10)
	  iBT_runinfo_li(lir_j) = iBT_runinfo(i);
	  lir_j=lir_j+1;
    	elseif strncmpi(iBT_runinfo(i).name,'runinfo_other_display',21)
	  iBT_runinfo_od(odr_j) = iBT_runinfo(i);
	  odr_j=odr_j+1;
	else
    	  junk = uimenu(MFile,...
                  'Label',['Create new runinfo.m from ' iBT_runinfo(i).name] ,...
                  'Callback',['iBT_copy_config(''' iBT_runinfo(i).name ''', ''runinfo.m'' )'], ...
                  'HandleVisibility','off');
	end %if ~strncmp else
    end %for


    iBT_studylist = dir(fullfile(spm('dir'),'toolbox','iBT','studylist*template.txt'));
    lis_j = 1;
    ods_j = 1;
    for i = 1:numel(iBT_studylist)
    	if strncmpi(iBT_studylist(i).name,'studylist_li',12)
	  iBT_studylist_li(lis_j) = iBT_studylist(i);
	  lis_j=lis_j+1;
    	elseif strncmpi(iBT_studylist(i).name,'studylist_other_display',23)
	  iBT_studylist_od(ods_j) = iBT_studylist(i);
	  ods_j=ods_j+1;
	else
    	  junk = uimenu(MFile,...
                  'Label',['Create new studylist.txt from ' iBT_studylist(i).name] ,...
                  'Callback',['iBT_copy_config(''' iBT_studylist(i).name ''', ''studylist.txt'' )'], ...
                  'HandleVisibility','off');
	end %if ~strncmp else
    end %for


    if lir_j > 1
	for i = 1:lir_j-1
    	  junk = uimenu(MFile,...
                  'Label',['Create new runinfo_li.m (for laterality only) from ' iBT_runinfo_li(i).name] ,...
                  'Callback',['iBT_copy_config(''' iBT_runinfo_li(i).name ''', ''runinfo_li.m'' )'], ...
                  'HandleVisibility','off');
	end % for i
    end% if lir_j > 1
    

    if lis_j > 1
	for i = 1:lis_j-1
    	  junk = uimenu(MFile,...
                  'Label',['Create new studylist_li.txt (for laterality only) from ' iBT_studylist_li(i).name] ,...
                  'Callback',['iBT_copy_config(''' iBT_studylist_li(i).name ''', ''studylist_li.txt'' )'], ...
                  'HandleVisibility','off');
	end % for i
    end% if lis_j > 1

    if odr_j > 1
	for i = 1:odr_j-1
    	  junk = uimenu(MFile,...
                  'Label',['Create new runinfo_other_display.m (for display of other images) from ' iBT_runinfo_od(i).name] ,...
                  'Callback',['iBT_copy_config(''' iBT_runinfo_od(i).name ''', ''runinfo_other_display.m'' )'], ...
                  'HandleVisibility','off');
	end % for i
    end% if odr_j > 1
    

    if ods_j > 1
	for i = 1:ods_j-1
    	  junk = uimenu(MFile,...
                  'Label',['Create new studylist_other_display.txt (for display of other images) from ' iBT_studylist_od(i).name] ,...
                  'Callback',['iBT_copy_config(''' iBT_studylist_od(i).name ''', ''studylist_other_display.txt'' )'], ...
                  'HandleVisibility','off');
	end % for i
    end% if ods_j > 1
    

    junk = uimenu(MFile,...
                  'Label','Generate a studylist.txt from an old runinfo.m',...
                  'Callback','iBT_make_studylist', ...
                  'HandleVisibility','on');
    
    junk = uimenu(MFile,...
                  'Label','Quit',...
                  'Callback','iBT_SPM(''CLOSE'')', ...
                  'HandleVisibility','on');

    %
    %  The edit menu
    %
    
    MEdit = uimenu(FiBT,...
                  'Label','Edit',...
                  'HandleVisibility','off');
    
    junk = uimenu(MEdit,...
                  'Label','Edit runinfo.m',...
                  'Callback','edit runinfo.m', ...
                  'HandleVisibility','off');

    junk = uimenu(MEdit,...
                  'Label','Edit studylist.txt',...
                  'Callback','edit studylist.txt', ...
                  'HandleVisibility','off');

    junk = uimenu(MEdit,...
                  'Label','Edit runinfo_li.m',...
                  'Callback','edit runinfo_li.m', ...
                  'HandleVisibility','off');

    junk = uimenu(MEdit,...
                  'Label','Edit studylist_li.txt',...
                  'Callback','edit studylist_li.txt', ...
                  'HandleVisibility','off');

    junk = uimenu(MEdit,...
                  'Label','Edit runinfo_other_display.m',...
                  'Callback','edit runinfo_other_display.m', ...
                  'HandleVisibility','off');

    junk = uimenu(MEdit,...
                  'Label','Edit studylist_other_display.txt',...
                  'Callback','edit studylist_other_display.txt', ...
                  'HandleVisibility','off');

    junk = uimenu(MEdit,...
                  'Label','Edit another file',...
                  'Callback','edit', ...
                  'HandleVisibility','off');

    %
    %  The iBt menu
    %
    
    MLaterality  = uimenu(FiBT,...
                  'Label','Run',...
                  'HandleVisibility','on');

    junk  = uimenu(MLaterality,...
                  'Label', 'Run a job',...
                  'CallBack', 'iBT_start;',...
                  'HandleVisibility','on');


