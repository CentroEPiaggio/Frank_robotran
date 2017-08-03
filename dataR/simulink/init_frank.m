% init_Frank
addpath('/home/parallels/.robotran/mbsysc/MBsysC/build/lib')

mbs_load('Frank_segway');

prjname = 'Frank_segway';

% I use this in controller function
rob_str.dpt = MBS_data.dpt;
rob_str.m = MBS_data.m;

% step of simulation
Ts = 0.001;

