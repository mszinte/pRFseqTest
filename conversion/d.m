daq.reset;
ni_devices = daq.getDevices;
ni_session = daq.createSession(ni_devices.Vendor.ID);
ni_channel_button = ni_session.addDigitalChannel('Dev2','port0/line0:7','InputOnly');
ni_channel_band = ni_session.addDigitalChannel('Dev2','port1/line0','InputOnly');
ni_session.inputSingleScan;
% return = (0     0     0     0     0     0     0     0     1) order is
% given in order of 
mb_factor = 4;
slice_number = 40;
tr_dur = 1.6;
% 10 trigger per TR from 0 to 1 until the end (toggle mode)