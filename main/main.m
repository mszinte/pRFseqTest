function main(const)
% ----------------------------------------------------------------------
% main(const)
% ----------------------------------------------------------------------
% Goal of the function :
% Launch all function of the experiment
% ----------------------------------------------------------------------
% Input(s) :
% const : struct containing a lot of constant configuration
% ----------------------------------------------------------------------
% Output(s):
% none
% ----------------------------------------------------------------------
% Function created by Martin SZINTE (martin.szinte@gmail.com)
% Last update : 13 / 01 / 2020
% Project :     pRFseqTest
% Version :     1.0
% ----------------------------------------------------------------------

tic;
% Connection
% ----------
if ismac && const.tracker
    [~,~] = system('echo invibe | sudo -S networksetup -switchtolocation "Eyelink"');
else
    [~,~] = system('echo invibe | sudo -S networksetup -switchtolocation "Web"');
end
pause(2);

% File director
% -------------
[const]                 =   dirSaveFile(const);

% Screen configurations
% ---------------------
[scr]                   =   scrConfig(const);

% Triggers and button configurations
% ----------------------------------
[my_key]                =   keyConfig(const);

% Experimental constant
% ---------------------
[const]                 =   constConfig(scr,const);

% Experimental design
% -------------------
[expDes]                =   designConfig(const);

% Open screen window
% ------------------
[scr.main,scr.rect]     =   Screen('OpenWindow',scr.scr_num,const.background_color,[], scr.clr_depth,2);
[~]                     =   Screen('BlendFunction', scr.main, GL_SRC_ALPHA, GL_ONE_MINUS_SRC_ALPHA);
priorityLevel           =   MaxPriority(scr.main);Priority(priorityLevel);

% Initialise eye tracker
% ----------------------
el                      =   [];
if const.tracker
    [el]                    =   initEyeLink(scr,const);
end

% Draw and save all textures
% --------------------------
if const.genStimuli
    makeTextures(scr,const,expDes);
end

% Trial runner
% ------------
[const]                 = runExp(scr,const,expDes,el,my_key);

% End
% ---
overDone(const,my_key)

end