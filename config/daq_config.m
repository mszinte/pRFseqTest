%% NI parameters
AcqRate = 1000; % Hz. Rate at which analog NI data will be acquired. Monitoring sampling rate will be AcqRate / AcqBufSize, because NI background function (NI_OnDataAvailable) runs every AcqBufSize scans acquired
AcqBufSize = 100;
deviceID = 'Dev1'; % if not correct, check with 'daq.getDevices'
LeftPressCh = 7; % check with 'panneaux de test' of NI driver
RightPressCh = 3; % check with 'panneaux de test' of NI driver
Sens_Thresh = [1.4;1.6]; % Threshold in volt at which a hand will be considered present

NI_Press = daq.createSession('ni'); % create session for pression sensors

%% Add digital and analog channels
PressIn(1) = addAnalogInputChannel(NI_Press,deviceID,LeftPressCh,'Voltage');
PressIn(2) = addAnalogInputChannel(NI_Press,deviceID,RightPressCh,'Voltage');
% PowerIn = addAnalogInputChannel(NI_Press,deviceID,PowerCh,'Voltage');

for i = 1:length(PressIn)
    PressIn(i).TerminalConfig = 'SingleEnded'; % set pressure sensor channels to SingleEnded (referencee)
end
NI_Press.Rate = AcqRate;
NI_Press.NotifyWhenDataAvailableExceeds = AcqBufSize;

NI_Press.IsContinuous = 1;

fprintf('done\n')


l_NI_Press = addlistener(NI_Press,'DataAvailable',@NI_OnDataAvailable);


startBackground(NI_Press);




function NI_OnDataAvailable(~,event)
event.Data


stop(NI_Press);
