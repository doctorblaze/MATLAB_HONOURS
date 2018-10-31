%% PRIME
%clear;
Studies = {'Rat1Scan1';'Rat2Scan1';'Rat4Scan1';'Rat5Scan1'; 'Rat6Scan1';...
           'Rat7Scan1';'Rat8Scan1';'Rat9Scan1'; 'Rat1Scan2';'Rat2Scan2';...
           'Rat3Scan2';'Rat5Scan2';'Rat6Scan2';'Rat7Scan2';...
           'Rat8Scan2';'Rat9Scan2'};
 %Rat3Scan1 removed because very strange scan and possibly misregistered
 %Rat6Scan1 has no dynamic data
 %Rat4Scan2 has no NORM120
       
pth='/Users/nicolas/Desktop/baseline_FDG_workdir';

%% IMPORT DATA
%import the individual .mat VOI files from SAMIT and strip out data 
%place into ordered structure called DATA

for s = 1:4 %Rat1Scan1 to Rat5Scan1
    
    filename = sprintf('ID/%s_VOI.mat',Studies{s});
    infile = fullfile(pth,Studies{s},filename);
    DATA(s).studyname = Studies{s};
    DATA(s).RAWdata = load(infile,'T');
    DATA(s).staticVOIs = DATA(s).RAWdata.T(1:51,3:8); %get VOI number, name, mean, SD, min, max

    for f=1:30
        frames = 52:51:height(DATA(s).RAWdata.T);
        DATA(s).dynamicVOIs(f).VOIs = DATA(s).RAWdata.T(frames(f):frames(f)+50,3:8); %get VOI number, name, mean, SD, min, max
    end
    
end
for s = 5 %Rat6Scan1 (has no dynamic data)
    filename = sprintf('ID/%s_VOI.mat',Studies{s});
    infile = fullfile(pth,Studies{s},filename);
    DATA(s).studyname = Studies{s};
    DATA(s).RAWdata = load(infile,'T');
    DATA(s).staticVOIs = DATA(s).RAWdata.T(1:51,3:8);
end
for s = 6:13 %6:8 is Rat7Scan1 to Rat9Scan1; 9:13 is Rat1Scan2 to Rat6Scan2 (excluding Rat4Scan2, no NORM120) 
    filename = sprintf('ID/%s_VOI.mat',Studies{s});
    infile = fullfile(pth,Studies{s},filename);
    DATA(s).studyname = Studies{s};
    DATA(s).RAWdata = load(infile,'T');
    DATA(s).staticVOIs = DATA(s).RAWdata.T(1:51,3:8); %get VOI number (:,1), name (:,2), mean (:,3), SD (:,4), min (:,5), max (:,6)

    for f=1:30
        frames = 52:51:height(DATA(s).RAWdata.T);
        DATA(s).dynamicVOIs(f).VOIs = DATA(s).RAWdata.T(frames(f):frames(f)+50,3:8); %get VOI number, name, mean, SD, min, max
    end
end
for s = 14:16 %Rat7Scan2 to Rat9Scan2
    filename = sprintf('ID/%s_VOI.mat',Studies{s});
    infile = fullfile(pth,Studies{s},filename);
    DATA(s).studyname = Studies{s};
    DATA(s).RAWdata = load(infile,'T');
    DATA(s).staticVOIs = DATA(s).RAWdata.T(1:51,3:8); %get VOI number (:,1), name (:,2), mean (:,3), SD (:,4), min (:,5), max (:,6)

    for f=1:29  %only have 29 dynamic frames for these guys
        frames = 52:51:height(DATA(s).RAWdata.T);
        DATA(s).dynamicVOIs(f).VOIs = DATA(s).RAWdata.T(frames(f):frames(f)+50,3:8); %get VOI number, name, mean, SD, min, max
    end
end

%% Static Image VOI Plot

studylook = [14 16];         %specify first and last study to look at as vector
studynum = 3;              %specify how many studies are being plotted
switcher = 'on';          %when on, the data is sorted according to index variable rather than sorting the mean in descending order
%plotcolor = 'b';
%plotcolor = [0.3 0.75 0.93];
plotcolor = [0, 0.5, 0];        %color of datapoints
plotlegend = 'on';      %activate pre- vs post- legend
pointjoiner = 'off';     %join rat datapoints with a grey line
noSEMpool = 'on';       %plots errorbar using nonPooled SEM
plotareasize = 'off';    %plots the areas size in voxels

VOInum_array = [];
mean_array = [1:51]'; %VOI number
SD_array = [1:51]';
min_array = [1:51]';
max_array = [1:51]';
whole_array = [1:51]';
data_array = [];
%place the data into arrays, which are easier to manipulate
%rows are areas 1 to 51; columns (:,2:end) are studies{s} (except for whole_array)

for s = studylook(1):studylook(2)
    VOInum_array = [VOInum_array  table2array(DATA(s).staticVOIs(:,1))]; %check that VOIs are in chronological order from 1:51 (they should)
    mean_array = [mean_array table2array(DATA(s).staticVOIs(:,3))];
    SD_array = [SD_array table2array(DATA(s).staticVOIs(:,4))];
    min_array = [min_array table2array(DATA(s).staticVOIs(:,5))];
    max_array = [max_array table2array(DATA(s).staticVOIs(:,6))];
    whole_array = [whole_array table2array(DATA(s).staticVOIs(:,3:end))]; %rows are VOIs, columns 2:5 is Studies{1} data (mean, SD, min, max), 6:9 is Studies{2} etc
end

whole_table = DATA(1).staticVOIs(:,1:2);                        %holds VOI numbers (:,1) and names (:,2)
whole_table = [whole_table array2table(whole_array(:,2:end))];  %all animals data in one table

data_array = table2array(whole_table(:,3:end));                 %this is exactly the same as whole array, except first column doesn't have VOI numbers... consider deleting this variable

overall_mean = array2table(mean(data_array(:,1:4:end),2)); %calculates overall (or pooled) mean
overall_mean.Properties.VariableNames = {'overall_mean'};

pooled_SD = sqrt(mean(data_array(:,2:4:end).^2,2));
overall_SD = array2table(pooled_SD); %determine composite/pooled SD for each area; https://en.wikipedia.org/wiki/Pooled_variance.
                                     %because the number of voxels used to calculate SD in e.g area 1 is the same across all Studies{s}
                                     %an estimate for the pooled SD is square root of the average of the squared SD values
                                     %(https://www.researchgate.net/post/how_to_compute_composite_standard_deviation_for_multiple_samples)
overall_SD.Properties.VariableNames = {'overall_SD'};

overall_SEM = array2table(pooled_SD./sqrt(studynum));
overall_SEM.Properties.VariableNames = {'overall_SEM'};

overall_min = array2table(mean(data_array(:,3:4:end),2));
overall_min.Properties.VariableNames = {'overall_min'};

overall_max = array2table(mean(data_array(:,4:4:end),2));
overall_max.Properties.VariableNames = {'overall_max'};

whole_table = [whole_table overall_mean overall_SD overall_min overall_max overall_SEM];

if strcmp('off',switcher)
%USE FOR PLOTTING THE INITIAL SET OF DATA
    [sorted_whole_table, index] = sortrows(whole_table, find(strcmpi(whole_table.Properties.VariableNames,'overall_mean')), 'descend'); %uses pooled mean column of whole_table
    sorted_whole_array = table2array(sorted_whole_table(:,3:end));
elseif strcmp('on',switcher)
%USE FOR PLOTTING COMPARISON (e.g. post-stress) TO FIRST PLOTTED DATA 
%(e.g. pre-stress), AVOIDS DIFFERENTIAL SORTING OF PRE-STRESS AND POST-STRESS
%by sortrows
    sorted_whole_array = table2array(whole_table(:,3:end));
    sorted_whole_array = sorted_whole_array(index,:);           %can't use indexing to sort a table unfortunately
end

labels = table2cell(sorted_whole_table(:,2));
figure(1); hold on; xticks([1:51]), xticklabels(labels), xtickangle(90),

if strcmp('on',plotlegend)  %used when comparing pre and post stress, makes individual data points color blue or red
    for s = 1:4:(size(sorted_whole_array(1,:),2)-5)
        hold all
        scatter(1:51,sorted_whole_array(:,s),10,'filled','MarkerFaceColor',plotcolor)
    end
else        %used when plotting a single group (e.g. pre-stress), individual points/rats take different colors
     for s = 1:4:(size(sorted_whole_array(1,:),2)-5)
        hold all
        scatter(1:51,sorted_whole_array(:,s),10,'filled')
    end   
end

if strcmp('on',pointjoiner) %joins the points with a grey line.
    greylines = plot(1:51,sorted_whole_array(:,1:4:(size(sorted_whole_array(1,:),2)-5)),'Color','k');
    for i = 1:studynum
        greylines(i).Color(4) = 0.1;
    end
end

legend(Studies{studylook(1):studylook(2)})

if strcmp('off',noSEMpool)  %plots pooled mean and SEM
    figure(1);
    pooled_mean = sorted_whole_array(:,end-4);
    pooled_SEM = sorted_whole_array(:,end);
    errorbar(1:51,pooled_mean,pooled_SEM,'s','MarkerFaceColor',plotcolor,'MarkerSize',5,'Color',plotcolor);

elseif strcmp('on',noSEMpool) %used to plot the non pooled SEM bars
    SD_of_mean_notpooled = [];
    for i = 1:51
        SD_of_mean_notpooled = [SD_of_mean_notpooled; std(data_array(i,1:4:end))];
    end
    pooled_mean = sorted_whole_array(:,end-4);
    SD_of_mean_notpooled = SD_of_mean_notpooled(index);
    SEM_of_mean_notpooled = SD_of_mean_notpooled./sqrt(studynum);
    figure(1);
    errorbar(1:51,pooled_mean,SEM_of_mean_notpooled,'s','MarkerFaceColor',plotcolor,'MarkerSize',5,'Color',plotcolor,'LineWidth',2)
end

if strcmp('on',plotlegend)
    L{1} = 'pre-stress';
    L{2} = 'post-stress';
    L{3} = 'post-stress control';
    LH(1) = plot(nan,nan,'s','Color','b','MarkerFaceColor','b');
    LH(2) = plot(nan,nan,'s','Color','r','MarkerFaceColor','r');
    LH(3) = plot(nan,nan,'s','Color',[0, 0.5, 0],'MarkerFaceColor',[0, 0.5, 0]);
    legend(LH, L)
end

if strcmp('on',plotareasize)
    areastoextract = [1:50]; 
    % Specify paths
    samit_mask = '/Users/nicolas/Desktop/workdir/SAMIT-Merged.nii';
    
    % Import samit mask
    samit_V = spm_vol(samit_mask);       %imports nifti into structure
    samit_Vm = spm_read_vols(samit_V);   %converts structure into array
    
    % Extract index of all voxels in samit_mask which satisfy areastoextract
    areastoextract_size = [];
    for i = 1:size(areastoextract,2)
        areastoextract_size = [areastoextract_size; size(find(samit_Vm == areastoextract(i)),1)];
    end
    areastoextract_size = [sum(areastoextract_size); areastoextract_size];
    
    %rearrange the values to be sorted like in DATA_plotter
    areastoextract_size = areastoextract_size(index);
    
    yyaxis right
    plot(1:51, areastoextract_size,'Color','r')
    hand = gca; hand.YColor = 'r';
    
    clear samit_V
    clear samit_Vm
    yyaxis left
end

%% Dynamic Image VOI TAC plot


studylook = [7 7];         %specify first and last study to look at as vector
studynum = 3;
switcher = 'off';          %can specify studylook as [1 4] and will do studies 6:8, skipping Study{5} (rat 6_1 which has no dynamic data)
plotcolor = 'g';
plotlegend = 'off';
noSEMpool = 'on';       %do not pool SDs of all the means, just calculate SD/SEM of the mean values as if they did not have their own SD
for_cointigration = [];
%plot on top of motion data
%gcf; subplot(20,2,33:2:39); hold on; plot([dynamicdatamean_array(1,1); dynamicdatamean_array(2:end,1)-60],overall_uptake*20,'Color','k');

for i = 19           %open SAMIT-Merged.txt or DATA(x).RAWdata.T and find index of the area you would like to plot
    area_index = i+1;    %if area is 24 in SAMIT.txt, it is stored in row 25 here (because first row is whole brain)
    dynamicdatamean_array = [[0:120:3600]'];
    dynamicdataSD_array = [[0:120:3600]'];
    for s = studylook(1):studylook(2)
        for f=1:30
            dynamicdatamean_array(f+1,s+1) = table2array(DATA(s).dynamicVOIs(f).VOIs(area_index,3));
            dynamicdataSD_array(f+1,s+1) = table2array(DATA(s).dynamicVOIs(f).VOIs(area_index,4));
        end
    end
    if strcmp('on',switcher)
        for s = 6:8
            for f=1:30
                dynamicdatamean_array(f+1,s+1) = table2array(DATA(s).dynamicVOIs(f).VOIs(area_index,3));
                dynamicdataSD_array(f+1,s+1) = table2array(DATA(s).dynamicVOIs(f).VOIs(area_index,4));
            end
        end
    end
    dynamicdatamean_array(dynamicdatamean_array==0) = nan; %in case there are missing values equal to zero
    dynamicdatamean_array(1,:) = 0;                        %makes first row = 0
    dynamicdataSD_array(dynamicdataSD_array==0) = nan;
    dynamicdataSD_array(1,:) = 0; 
    
    overall_uptake = mean(dynamicdatamean_array(:,2:end),2,'omitnan'); %mean of observations
    
    if strcmp('off',noSEMpool)
        if studylook(2)-studylook(1)>0
            pooled_uptake_SD = sqrt(mean(dynamicdataSD_array(:,2:end).^2,2,'omitnan'));
            overall_uptake_SEM = pooled_uptake_SD./sqrt(studynum);
        else
            voxelnum = areastoextract_size(i+1);
            overall_uptake_SEM = dynamicdataSD_array(:,studylook(1)+1)./sqrt(voxelnum); %determine voxelnum in this area.. use mask?
        end
    end
    if strcmp('on',noSEMpool)
        if studylook(2)-studylook(1)>0
            pooled_uptake_SD = std(dynamicdatamean_array(:,2:end),0,2,'omitnan');
            overall_uptake_SEM = pooled_uptake_SD./sqrt(studynum);
        else
            voxelnum = areastoextract_size(i+1);
            overall_uptake_SEM = dynamicdataSD_array(:,studylook(1)+1)./sqrt(voxelnum); %determine voxelnum in this area.. use mask?
        end
    end

    figure(5);
    %plot(0:120:3600,overall_uptake);
    %shadedErrorBar(0:120:3600,overall_uptake,overall_uptake_SEM,'r')
    hold on;
    %for s = studylook(1):studylook(2)
        %scatter(dynamicdata_array(:,1),dynamicdata_array(:,s),2,'filled','k')
    %end
    %plot(dynamicdata_array(:,1),dynamicdata_array(:,studylook(1):studylook(2)),'LineWidth',1,'Color','k')
    plot(dynamicdatamean_array(:,1)/60,overall_uptake,'Color',plotcolor);
    errorbar(dynamicdatamean_array(:,1)/60,overall_uptake,overall_uptake_SEM,'s','MarkerFaceColor',plotcolor,'MarkerSize',5,'Color',plotcolor)

for_cointigration = [for_cointigration overall_uptake];
end

if strcmp('on',plotlegend)
    L{1} = 'pre-stress';
    L{2} = 'post-stress';
    L{3} = 'post-stress control';
    LH(1) = plot(nan,nan,'s','Color','b','MarkerFaceColor','b');
    LH(2) = plot(nan,nan,'s','Color','r','MarkerFaceColor','r');
    LH(3) = plot(nan,nan,'s','Color','m','MarkerFaceColor','m');
    legend(LH, L)
    xlabel(table2cell(DATA(1).dynamicVOIs(1).VOIs(area_index,2)))
    ylabel('FDG uptake (standardized)')
    
end

%figure(1); subplot(20,2,33:2:39); hold on; yyaxis left, plot([dynamicdatamean_array(1,1); dynamicdatamean_array(2:end,1)-60],overall_uptake*20,'Color','k');

%figure(2); hold on; xlabel('Time (s)'); yyaxis left; pl1 = plot([dynamicdatamean_array(1,1); dynamicdatamean_array(2:end,1)-60], overall_uptake,'Marker','o','MarkerSize',2); ylabel('FDG uptake (standardized)'); yyaxis right; pl2 = plot(realframedistance(:,1),realframedistance(:,2),'Marker','o','MarkerSize',2), ylabel('Distance travelled (mm)'), xlim([0 3600])
%figure(2); yyaxis left; pl3 = plot([dynamicdatamean_array(1,1); dynamicdatamean_array(2:end,1)-60], overall_uptake)