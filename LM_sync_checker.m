%% ------------------------------------------------------------------------
% This script is intended to evaluate how well the MT pose data and the
% LM gate tag data match. In theory, there should be just as many MT poses
% as there are gate-tags in the LM data (e.g. a pose/gate every 33ms if the 
% signal generator [which triggers camera & gate-tag insertion] is set to 
% 30Hz). Giorgos also increases the signal gener. frequency every 10 minutes 
% (meaning the time interval between poses/gates will be <33ms). This helps
% check how well MT and LM are synced.
% For use with mPET_sync_listmode_v1 or v2 (new software)
%
% Author: Nicolas Darmanthe
%% ------------------------------------------------------------------------

%Data Import
notSynced = []; %from _before_sync.log
mcSynced = PTSDRat2Scan1v1aftersync; %import data from _after_sync.log 

%Helps estimate location of control events in MT/LM data
TrackerStart = 615; %input tracker start time
SignalGenControlTime = [970 1570 2170 2770 3370 3970]; %input SignalGen control timepoints in seconds as vector

%Converts SingalGenControl timepoints to posepoints 
SignalGenControlPose = [];
for i=1:6
SignalGenControlPose = [SignalGenControlPose ((SignalGenControlTime(i)-TrackerStart)*33)];
end

%Plot Signal Gen Controls pre_sync
figure('Name', 'Signal generator control 1')
hold on
ylim([-20 80]); xlim([SignalGenControlPose(1)-3000 SignalGenControlPose(1)+3000])
plot(notSynced)


figure('Name', 'Signal generator control 2')
hold on
ylim([-20 80]); xlim([SignalGenControlPose(2)-3000 SignalGenControlPose(2)+3000])
plot(notSynced)

figure('Name', 'Signal generator control 3')
hold on
ylim([-20 80]); xlim([SignalGenControlPose(3)-3000 SignalGenControlPose(3)+3000])
plot(notSynced)

figure('Name', 'Signal generator control 4')
hold on
ylim([-20 80]); xlim([SignalGenControlPose(4)-3000 SignalGenControlPose(4)+3000])
plot(notSynced)

figure('Name', 'Signal generator control 5')
hold on
ylim([-20 80]); xlim([SignalGenControlPose(5)-3000 SignalGenControlPose(5)+3000])
plot(notSynced)

figure('Name', 'Signal generator control 6')
hold on
ylim([-20 80]); xlim([SignalGenControlPose(6)-3000 SignalGenControlPose(6)+3000])
plot(notSynced)

%Plot mc/synced data, corrected for dropped poses and gate8 glitches
figure('Name', 'Signal generator control 1 synced')
hold on
ylim([-20 80]); xlim([SignalGenControlPose(1)-3000 SignalGenControlPose(1)+3000])
plot(mcSynced);

figure('Name', 'Signal generator control 2 synced')
hold on
ylim([-20 80]); xlim([SignalGenControlPose(2)-3000 SignalGenControlPose(2)+3000])
plot(mcSynced);

figure('Name', 'Signal generator control 3 synced')
hold on
ylim([-20 80]); xlim([SignalGenControlPose(3)-3000 SignalGenControlPose(3)+3000])
plot(mcSynced);

figure('Name', 'Signal generator control 4 synced')
hold on
ylim([-20 80]); xlim([SignalGenControlPose(4)-3000 SignalGenControlPose(4)+3000])
plot(mcSynced);

figure('Name', 'Signal generator control 5 synced')
hold on
ylim([-20 80]); xlim([SignalGenControlPose(5)-3000 SignalGenControlPose(5)+3000])
plot(mcSynced);

figure('Name', 'Signal generator control 6 synced')
hold on
ylim([-20 80]); xlim([SignalGenControlPose(6)-3000 SignalGenControlPose(6)+3000])
plot(mcSynced);


