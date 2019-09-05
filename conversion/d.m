% Use the daq.getDevices command to display a list of devices available to your machine and MATLAB
daq.getDevices

% To learn more about an individual device, click the name of the device in the list in the Command window, or access the device in the array returned by daq.getDevices command.
d = daq.getDevices;

s = daq.createSession('ni');
s.Rate = 5000;
s.DurationInSeconds = 2;

addAnalogInputChannel(s,'cDAQ1Mod1',0,'Voltage');
% https://ch.mathworks.com/help/daq/ref/adddigitalchannel.html

addDigitalChannel(s,'dev1','Port0/Line0:1','InputOnly');
ch = addDigitalChannel(s,'dev1','Port0/Line2:3','OutputOnly');
[ch,idx] = addDigitalChannel(s,'dev1','Port2/Line0:1','Bidirectional')

data = s.startForeground();


% https://ch.mathworks.com/help/daq/examples/software-analog-triggered-data-capture.html