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

%% ------------------------------------
data = sortrows(rat6_1,1); %sortrows is used to make sure that the time stamps are in chronological order

nans = find(isnan(data)); %find NaNs
data(nans(1,1),:) = [];   %delete NaN - this is not necessary

ifstmsr = find(data(:,5)==0);  %find location of first measured pose
data(1:ifstmsr(1,1),:) = [];   %delete all poses before first measured pose
                               %because for our purposes we don't care about these

datatemp = data(:,1:4);                         %this
[n, bin] = histc(datatemp, unique(datatemp));   %deletes
multiple = find(n > 1);                         %replicate
index    = find(ismember(bin, multiple));       %values
data(index,:)=[];                               %in data


fstmsr = data(1,1);
data(:,1) = data(:,1)-fstmsr;    %makes start time 0
figure('Color','w','name','Rat1Scan1 motion data')

ScanL = 3700; %scan time in seconds 

%% --dropped poses--------------------------------------
t=data(:,1);
s = size(data,1);
tdiff=zeros(s,1);

for i=2:s %start at cell 2 cell 1 is first measurement
    tdiff(i,1)=(t(i)-t(i-1));
end

Hz = round(1/median(tdiff,1)); %Finds frequency of MT in Hz
PoseT = median(tdiff,1);

droppedpose=tdiff(tdiff>(PoseT));
locdroppedpose=t(tdiff>(PoseT));
nlostpose=droppedpose./(PoseT);
pctdrop=sum(nlostpose)/size(data(:,1),1);

subplot (20,2,2:2:8); bar(locdroppedpose,nlostpose,'LineWidth',1), 
    %xlabel('time (s)'),
    ylabel('# of dropped poses'),
    xlim([0 ScanL]),
    ylim([0 4]),
    set(gca, 'XTickLabel', [])
    set(gca,'XMinorTick','off','YMinorTick','on'),
    hold on
    pd = plot(NaN,NaN,'>r');
    legend(pd, ['poses dropped=' num2str(pctdrop) '%'],'location','best')

    
%% --lost pose replace---------------------------

logic = data(:,5)==-1;
locelim = t(logic); %location (t) of lost poses  
numlost = size(locelim,1);
pctlost = (numlost/s)*100;      

nbins = 34;
subplot (20,2,1:2:7); histogram(locelim,nbins), 
    %xlabel('time in scan (s)'),
    ylabel('number of lost poses'),
    xlim([0 ScanL]),
    set(gca, 'XTickLabel', [])
    set(gca,'XMinorTick','off','YMinorTick','on'),
    hold on
    pl = plot(NaN,NaN,'>r');
    legend(pl, ['poses lost=' num2str(pctlost) '%'],'location','best');
    

nidx = data(:,5)==-1;
data(nidx,2) = interp1(data(~nidx,1),data(~nidx,2),data(nidx,1),'linear','extrap'); %replace zeros by linear extrapolation
data(nidx,3) = interp1(data(~nidx,1),data(~nidx,3),data(nidx,1),'linear','extrap');
data(nidx,4) = interp1(data(~nidx,1),data(~nidx,4),data(nidx,1),'linear','extrap');

% --dropped pose replace-------------------------- to be used with caution 
ti = min(data(:,1)):PoseT:max(data(:,1));
data = [ti' interp1(data(:,1), data(:,2:4),ti,'linear','extrap')];
%----
s = size(data,1);
t=data(:,1);

x=data(:,2);
y=data(:,3);
z=data(:,4);

    
%% --xaxis mouvement---------------------------------------
subplot(20,2,9:2:15); plot(t,x,'r-'), 
    %xlabel('time in scan (s)'),
    %ylabel('displacement (mm)'), 
    xlim([0 ScanL]),
    ylim([-100 100]),
    set(gca, 'XTickLabel', []),
    legend('X','location','best');
    
%% --yaxis mouvement --------------------------------------
subplot(20,2,17:2:23); plot(t,y,'b-'), 
   %xlabel('time in scan (s)'),
   ylabel('displacement (mm)'),
   xlim([0 ScanL]),
   ylim([-100 100])
   set(gca, 'XTickLabel', []),
   legend('Y','location','best');
   
%% --zaxis mouvement --------------------------------------
subplot(20,2,25:2:31); plot(t,z,'g-'), 
   %xlabel('time in scan (s)'),
   %ylabel('displacement (mm)'), 
   xlim([0 ScanL]),
   ylim([350 650])
   set(gca, 'XTickLabel', []),
   legend('Z','location','best');

%% --distance per frame-----------------------------------
framedist=zeros(s,1);
for i=3:s
   framedist(i,1)=sqrt((x(i)-x(i-1)).^2+(y(i)-y(i-1)).^2+(z(i)-z(i-1)).^2); %euclidian distance
end
subplot(20,2,10:2:24); scatter(t,framedist,1)
    xlabel('time (s)'),
    ylabel('mm/pose'),
    xlim([0 ScanL])

%% --distance per second (at Hz)------------------------------------
framedist1 = framedist(1 : fix(numel(framedist)/Hz)*Hz); % First x elements (multiples of Hz) to fit in boxes of 1s = Hz poses(30) 
framedist2 = framedist(fix(numel(framedist)/Hz)*Hz+1 : end); % Last x Elements to be placed in box<Hz (30)
A = reshape(framedist1(:), Hz , []); 
secdist = [sum(A) sum(framedist2)];

t=data(:,1);

t1 = t(1:fix(numel(framedist)/Hz)*Hz);
B = reshape(t1(:), Hz , []);
B(1,end+1) = B(1,end)+1;
 
tscan = B(1,:);
medsecdist = median(secdist);
linemedsecdist = linspace(medsecdist(1,1), medsecdist(1,1),ScanL);

subplot(20,2,33:2:39); scatter(tscan,secdist,1)
    xlabel('time in scan (s)'),
    ylabel('speed (mm/s)'),
    xlim([0 ScanL]),
    %ylim([0 100]),
    set(gca,'XMinorTick','off','YMinorTick','on'),
    hold on
    pm = plot(linemedsecdist,'--','Color','r');
    legend(pm,['median=' num2str(medsecdist)]);
    
subplot(20,2,26:2:40); histogram(secdist), 
    xlabel('mm/s'),
    ylabel('# occurences'),
    xlim([0 200]),
    set(gca,'XMinorTick','off','YMinorTick','on')


%% --savetopdf--------
%print(gcf,'rat4_2','-fillpage','-dpdf')

%% --include rotation matrix data-----

%% --freezing episodes-----------------------
%freezing = complete cessation of movement with the exception of that required for respiration
%preselection of microfreeze timepoints where mvmnt <5mm over 3s
%freezelength calculated by estimating the number of consecutive
%microfreezes

freezeloc = [];
freezelength = [];
freezelocrl = [];
for i=1:(length(secdist)-2)
    if sum(secdist(i:2+i))<5
        freezeloc = [freezeloc i]; %this makes time series integers, only broadly accurate freeze time locus (necessary for phase two of analysis)
        freezelocrl = [freezelocrl tscan(i)]; %real loc of freeze
    end
end

freezelocwrk = freezeloc;
while length(freezelocwrk)>31;
    for o=30:-1:0
        if freezelocwrk(1,1+o)-freezelocwrk(1,1)==o
            freezelength = [freezelength o+1];
            freezelocwrk = freezelocwrk(o+2:end);
        else
            continue
        end
    end
end

binedges = [0:1.00001:30];

figure('Name', 'Freeze')
histogram(freezelength,binedges); xlim([0 30])
xlabel('freezelength (s)')
ylabel('# occurences')

locelim_tst = round(locelim,2); 
freezelocrl_tst =round(freezelocrl,2);
causedby = (numel(intersect(locelim_tst,freezelocrl_tst))/numel(freezeloc))*100;
fprintf('%d\n percent of detected microfreezes caused caused by interpolation of lost poses\n',causedby)
%modify algorithm to remove artificial microfreezes from freezing analysis
%% --bursts of activity
%microburst defined as activity above 60mm/3s
%burstlength calculated by estimating the number of consecutive microbursts
%

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

figure('Name', 'Burst')
histogram(burstlength,binedges); xlim([0 30])
xlabel('burstlength (s)')
ylabel('# occurences')

locelim_tst = round(locelim,2); 
burstlocrl_tst =round(burstlocrl,2);
causedby = (numel(intersect(locelim_tst,burstlocrl_tst))/numel(burstloc))*100;
fprintf('%d\n percent of detected microfreezes caused caused by interpolation of lost poses\n',causedby)

