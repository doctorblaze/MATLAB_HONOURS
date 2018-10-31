%%replace char with num function, run seperately
%tocorrect = 'MT_scan5.txt';
%Text = fileread(tocorrect);
%Text = strrep(Text, ',', '.');
%Text = strrep(Text, '*', '-1');                                        
%Text = strrep(Text, 'None', '0');                                    
%fid = fopen(tocorrect, 'w');                       
%fwrite(fid, Text, 'char');                                          
%fclose(fid);                                                           
%clear Text       
%The user should import the time (t) column, Test X Y and Z and the Hazard column
% NEED INTRO!
%
% Author: Nicolas Darmanthe
%% ------------------------------------
experimentalmode = 'off';
studyname = 'Rat3Scan2';
ScanL = 3800; %scan time in seconds (sets axis length)
InjecEffectScan_TimeInterval = 12; %specify in seconds
data = sortrows(Rat8Scan1MTData,1); %sortrows is used to make sure that the time stamps are in chronological order

%data = data(3696:end,:);

nans = find(isnan(data)); %find NaNs
%data(nans(1,1),:) = [];   %delete NaN - this step is not necessary

ifstmsr = find(data(:,5)==0);  %find location of first measured pose
data(1:ifstmsr(1,1),:) = [];   %delete all poses before first measured pose                            

datatemp = data(:,1:4);                         %this
[n, bin] = histc(datatemp, unique(datatemp));   %deletes
multiple = find(n > 1);                         %replicate
index    = find(ismember(bin, multiple));       %values
data(index,:)=[];                               %in data


fstmsr = data(1,1);
data(:,1) = data(:,1)-fstmsr;    %makes start time 0
figure('Color','w','name',studyname) %Change this lol, not always rat1scan1

%% --dropped poses--------------------------------------
t=data(:,1);
s = size(data,1);
tdiff=zeros(s,1);

for i=2:s %start at cell 2 cell 1 is first measurement
    tdiff(i,1)=(t(i)-t(i-1)); %calculates the time difference between two poses
end

Hz = round(1/median(tdiff,1)); %Finds frequency of MT in Hz
PoseT = median(tdiff,1);       %Gives interval between poses

droppedpose=tdiff(tdiff>(PoseT));   %time differences between the dropped poses
locdroppedpose=t(tdiff>(PoseT));    %location of dropped poses
nlostpose=droppedpose./(PoseT);     %number of lost poses
pctdrop=sum(nlostpose)/size(data(:,1),1); %calculates pct of dropped poses

subplot (20,2,2:2:8); bar(locdroppedpose+(InjecEffectScan_TimeInterval),nlostpose,'LineWidth',1),
    %xlabel('time (s)'),
    ylabel('# of dropped poses'),
    xlim([0 ScanL]),
    ylim([0 4]),
    set(gca, 'XTickLabel', [])
    set(gca,'XMinorTick','off','YMinorTick','on'),
    hold on
    pd = plot(NaN,NaN,'>r');
    legend(pd, ['poses dropped = ' num2str(pctdrop) '%'],'location','best')

    
%% --lost pose replace---------------------------

logic = data(:,5)==-1;
locelim = t(logic);         %location in time (t) of lost poses  
numlost = size(locelim,1);  %number of lost poses
pctlost = (numlost/s)*100;  %percent of lost poses   

nbins = 30;
subplot (20,2,1:2:7); histogram(locelim+InjecEffectScan_TimeInterval,nbins),
    %xlabel('time in scan (s)'),
    ylabel('number of lost poses'),
    xlim([0 ScanL]),
    set(gca, 'XTickLabel', [])
    set(gca,'XMinorTick','off','YMinorTick','on'),
    hold on
    pl = plot(NaN,NaN,'>r');
    legend(pl, ['poses lost=' num2str(pctlost) '%'],'location','best');
    

nidx = data(:,5)==-1; %when a pose is lost, the corresponding Tx Ty Tz takes value 0.... holes in data...
data(nidx,2) = interp1(data(~nidx,1),data(~nidx,2),data(nidx,1),'linear','extrap'); %replace zeros in Tx by linear extrapolation
data(nidx,3) = interp1(data(~nidx,1),data(~nidx,3),data(nidx,1),'linear','extrap'); %replace zeros in Ty 
data(nidx,4) = interp1(data(~nidx,1),data(~nidx,4),data(nidx,1),'linear','extrap'); %replace zeros in Tz
%note: linear extrapolation underestimates the amout of actual motion.

%% --dropped pose replace-------------------------- to be used with caution (section can be commented out)
% --This 'resamples' the data and fills in the dropped poses. It removes noise from the time domain, as pose intervals
% are all exactly PoseT now
ti = min(data(:,1)):PoseT:max(data(:,1));
data = [ti' interp1(data(:,1), data(:,2:4),ti,'linear','extrap')];
%--------------------------------------------------------------------------
s = size(data,1);
t=data(:,1);

x=data(:,2);
y=data(:,3);
z=data(:,4);

    
%% --xaxis mouvement---------------------------------------
subplot(20,2,9:2:15); plot(t+InjecEffectScan_TimeInterval,x,'r-'), 
    %xlabel('time in scan (s)'),
    %ylabel('displacement (mm)'), 
    xlim([0 ScanL]),
    ylim([-100 100]),
    set(gca, 'XTickLabel', []),
    legend('X','location','best');
    
%% --yaxis mouvement --------------------------------------
subplot(20,2,17:2:23); plot(t+InjecEffectScan_TimeInterval,y,'b-'), 
   %xlabel('time in scan (s)'),
   ylabel('displacement (mm)'),
   xlim([0 ScanL]),
   ylim([-100 100])
   set(gca, 'XTickLabel', []),
   legend('Y','location','best');
   
%% --zaxis mouvement --------------------------------------
subplot(20,2,25:2:31); plot(t+InjecEffectScan_TimeInterval,z,'g-'), 
   %xlabel('time in scan (s)'),
   %ylabel('displacement (mm)'), 
   xlim([0 ScanL]),
   ylim([400 800])
   set(gca, 'XTickLabel', []),
   legend('Z','location','best');

%% --distance per pose-----------------------------------
posedist=zeros(s,1);
for i=3:s
   posedist(i,1)=sqrt((x(i)-x(i-1)).^2+(y(i)-y(i-1)).^2+(z(i)-z(i-1)).^2); %calculates euclidian distance between frames
end
medposedist = median(posedist);
linemedposedist = linspace(posedist(1,1), posedist(1,1),ScanL);

subplot(20,2,10:2:24); scatter(t+InjecEffectScan_TimeInterval,posedist,1)
    xlabel('time (s)'),
    ylabel('mm/pose'),
    xlim([0 ScanL]),
    ylim([0 10])
    hold on
    pm = plot(linemedposedist,'--','Color','r');
    legend(pm,['Global median speed (mm/pose) =' num2str(medposedist)],'Location','best');
    

%% --distance per second (at Hz)------------------------------------
posedist1 = posedist(1 : fix(numel(posedist)/Hz)*Hz);        % First x elements (multiples of Hz) to fit in boxes of 1s = Hz poses(30) 
posedist2 = posedist(fix(numel(posedist)/Hz)*Hz+1 : end);    % Last x Elements to be placed in box<Hz (30)
A = reshape(posedist1(:), Hz , []);                            % Each column is 1 second; in each column all the mm/pose measurements for that given second (to be summed)
secdist = [sum(A) sum(posedist2)];                             % each column is 1 second, takes the distance travelled in that given second

t=data(:,1);

t1 = t(1:fix(numel(posedist)/Hz)*Hz);                            
B = reshape(t1(:), Hz , []);                    %the first row is the time every 30 poses
B(1,end+1) = B(1,end)+1;                        %glitch prevention - adds a 1 second buffer at the end of the time series
 
tscan = B(1,:);                                 %the first row is the time every 30 poses
medsecdist = median(secdist);
linemedsecdist = linspace(medsecdist(1,1), medsecdist(1,1),ScanL);

subplot(20,2,33:2:39); scatter(tscan+InjecEffectScan_TimeInterval,secdist,1)
    xlabel('time in scan (s)'),
    ylabel('speed (mm/s)'),
    xlim([0 ScanL]),
    ylim([0 50]),
    set(gca,'XMinorTick','off','YMinorTick','on'),
    hold on
    pm = plot(linemedsecdist,'--','Color','r');
    legend(pm,['Global median speed (mm/s) =' num2str(medsecdist)],'Location','best');
    
subplot(20,2,26:2:40); histogram(secdist), 
    xlabel('mm/s'),
    ylabel('# occurences'),
    xlim([0 200]),
    set(gca,'XMinorTick','off','YMinorTick','on')

%Distance per frame (e.g. per 2mn dynamic frame)
framelength = 120;
framedist = [0]; %0mm/s at 0s
%calculate mm/s starting at 120-240s post-injection (will break for , where there was 130s?)
%this discards any motion data located between 0-120s post-injection
%for d = 110:framelength:(length(secdist)-framelength)
     %framedist = [framedist median(secdist(d:d+framelength))];
%end
for d = 1+(framelength-InjecEffectScan_TimeInterval):framelength:(length(secdist)-framelength)
    framedist = [framedist median(secdist(d:d+framelength))];
end
%above does not calculate distance between 0-120s post-injection (because missing data), add nan
%framedist = [NaN NaN framedist];
framedist = [framedist(1,1) nan framedist(1,2:end)];

framedisttime = 0:framelength:length(secdist)+13;
framedisttime = [framedisttime(1,1) framedisttime(1,2:end)-framelength/2]; %plot the median distance of a 2mn frame in the middle of the 2mn frame
subplot(20,2,33:2:39); plot(framedisttime,framedist,'LineWidth',3,'DisplayName','Frame median speed (mm/s)');

%"Distance" variables above actually represent median velocity per second (mm/second) or velocity per frame (mm/frame)!
%Bellow calculates the actual distance travelled by the head marker

MOT2 = [(tscan+InjecEffectScan_TimeInterval)' secdist']; %median motion every second
croploc = find(floor(MOT2(:,1))==120); %deletes the data between effective scan start (e.g 16s) and 120s
MOT2 = MOT2(croploc:end,:);
entirescandist = trapz(MOT2(:,2));
if size(MOT2)<3600
   toadd = 3600-size(MOT2,1);
   MOT2 = [MOT2; NaN([toadd 2])];
end
MOT2short = MOT2(1:(floor(size(MOT2,1)/120)*120),2);
motbinner = reshape(MOT2short,[120,size(MOT2short,1)/120]); %put into 2mn bins, columnwise
realframedistance = [[120+60:120:3600+60]' trapz(motbinner)']; %calculates AUC under velocity-time graph (is distance)
realframedistance = [NaN NaN;NaN NaN;realframedistance];
subplot(20,2,33:2:39); yyaxis right, plot(realframedistance(:,1),realframedistance(:,2),'DisplayName','Frame distance (mm)','Marker','s','LineStyle','none','MarkerSize',5,'MarkerEdgeColor',[0.75, 0, 0.75],'MarkerFaceColor',[0.75, 0, 0.75]), ylabel('Distance (mm)')

%% --savetopdf--------
%print(gcf,'rat4_2','-fillpage','-dpdf')

%% --include rotation matrix data-----
%currently, rotation and pitch, yaw, roll data is neglected from the distance calculations
%% --freezing episodes-----------------------EXPERIMENTAL-----------------
%freezing = complete cessation of movement with the exception of that required for respiration
%preselection of microfreeze timepoints where mvmnt <5mm over 3s * (arbitrary value..)
%Machine learning would be cool to do this
%freezelength calculated by estimating the number of consecutive
%microfreezes (sliding window type method)

if strcmp('on',experimentalmode)
    freezeloc = [];
    freezelength = [];
    freezelocrl = [];
    for i=1:(length(secdist)-2)
        if sum(secdist(i:2+i))<5 
            freezeloc = [freezeloc i]; %this makes time series integers, only broadly accurate freeze time locus (necessary for phase two of algorithm)
            freezelocrl = [freezelocrl tscan(i)]; %real loc of freeze
        end
    end

    freezelocwrk = freezeloc;
    while length(freezelocwrk)>31 
        for o=30:-1:0
            if (freezelocwrk(1,o+1))-freezelocwrk(1,1)==o
                freezelength = [freezelength o+1];
                freezelocwrk = freezelocwrk(1,o+2:end); %crops the data that's just been used and proceeds along
            else
                continue
            end
        end
    end

    binedges = [0:1.00001:30];

    figure('Name', strcat(studyname,' Freeze'))
    histogram(freezelength,binedges); xlim([0 30])
    xlabel('freezelength (s)')
    ylabel('# occurences')

    locelim_tst = round(locelim,2); 
    freezelocrl_tst =round(freezelocrl,2);
    causedby = (numel(intersect(locelim_tst,freezelocrl_tst))/numel(freezeloc))*100;
    fprintf('%d\n percent of detected microfreezes caused caused by interpolation of lost poses\n',causedby)
    %modify algorithm to remove artificial microfreezes from freezing analysis
end
%% --bursts of activity-----------------------EXPERIMENTAL-----------------
%microburst defined as activity above 60mm/3s
%burstlength calculated by estimating the number of consecutive microbursts
%
if strcmp('on',experimentalmode)
    burstloc=[];
    burstlocrl=[];
    burstlength=[];
    for i=1:(length(secdist)-2)
        if sum(secdist(i:2+i))>60
            burstloc = [burstloc i]; 
            burstlocrl = [burstlocrl tscan(i)]; 
        end
    end

    burstlocwrk = burstloc;
    while length(burstlocwrk)>31;
        for o=30:-1:0
            if burstlocwrk(1,1+o)-burstlocwrk(1,1)==o
                burstlength = [burstlength o+1];
                burstlocwrk = burstlocwrk(o+2:end);
            else
                continue
            end
        end
    end

    binedges = [0:1.00001:30];

    figure('Name', strcat(studyname,' Burst'))
    histogram(burstlength,binedges); xlim([0 30])
    xlabel('burstlength (s)')
    ylabel('# occurences')

    locelim_tst = round(locelim,2); 
    burstlocrl_tst =round(burstlocrl,2);
    causedby = (numel(intersect(locelim_tst,burstlocrl_tst))/numel(burstloc))*100;
    fprintf('%d\n percent of detected bursts caused caused by interpolation of lost poses\n',causedby)
end
