function [const]=sbjConfig(const)
% ----------------------------------------------------------------------
% [const]=sbjConfig(const)
% ----------------------------------------------------------------------
% Goal of the function :
% Define subject configurations (initials, gender...)
% ----------------------------------------------------------------------
% Input(s) :
% const : struct containing constant configurations
% ----------------------------------------------------------------------
% Output(s):
% const : struct containing constant configurations
% ----------------------------------------------------------------------
% Function created by Martin SZINTE (martin.szinte@gmail.com)
% Last update : 13 / 08 / 2019
% Project :     pRFseqTest
% Version :     1.0
% ----------------------------------------------------------------------

if const.expStart
    const.sjctNum           =  input(sprintf('\n\tParticipant number: '));
    if isempty(const.sjctNum)
        error('Incorrect participant number');
    end
    if const.sjctNum > 9
        const.sjct              =  sprintf('sub-%i',const.sjctNum);
    else
        const.sjct              =  sprintf('sub-0%i',const.sjctNum);
    end
end

const.runNum            =   input(sprintf('\n\tRun number (1 to 10): '));
if isempty(const.runNum)
    error('Incorrect run number');
end
if const.runNum > 10
    error('Only 10 runs');
end

if const.expStart == 0
    const.cond1             =   1;
    const.cond2             =   input(sprintf('\n\tSequence 1 (1), Sequence 2 (2) : '));
    if ~(const.cond2 == 1 || const.cond2 == 2)
        error('Type either (1) or (2)')
    end
else
    const.cond1         =   const.cond_run_order(const.runNum,1);
    const.cond2         =   const.cond_run_order(const.runNum,2);    
end

const.cond1_txt          =  'AttendStim';

if const.cond2 == 1
    const.cond2_txt         =  '_acq-2p5mm';
elseif const.cond2 == 2
    const.cond2_txt         =  '_acq-2mm';
end

fprintf(1,'\n\tTask: %s%s\n',const.cond1_txt,const.cond2_txt);

if const.expStart
    if const.tracker
        const.sjct_DomEye   =   'L';  % for all subjects
        const.recEye        =   1;
    else
        const.sjct_DomEye   =   'DM';
        const.recEye        =   1;
    end
    if const.runNum == 1
        const.sjctName      =   upper(strtrim(input(sprintf('\n\tParticipant identity: '),'s')));
        const.sjct_age      =   input(sprintf('\n\tParticipant age: '));
        if isempty(const.sjct_age)
            error('Incorrect participant age');
        end
        const.sjct_gender   =   upper(strtrim(input(sprintf('\n\tParticipant gender (M or F): '),'s')));
    end
else
    const.sjct          =   'sub-00X';
    const.sjct_age      =   '00';
    const.sjct_gender   =   'X';
    const.sjct_DomEye   =   'DM';
    const.recEye        =   1;
end

end