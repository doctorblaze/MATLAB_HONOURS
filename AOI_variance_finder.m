%% --Description----------------------------------------------------------------------------
%   This script is intended to find the area (A) with least fluctuation relative
%   to other brain areas across differerent studies/images. The area with least fluctuation is the best area to
%   use for normalisation
%
%   A is found by:
%
%   (i) within a subject/image, find the absolute difference between a given area A1 
%   and all other areas (e.g. A2-A4) (how to avoid artificial inflation of variance)
%
%   (ii) reiterate for each area (combinatorial; diff between A2 and all other areas, 
%   then diff between A3 and all other areas etc
%   
%   (iii) do this in all subjects/images
%
%   (iii) the area which has the least overall variance relative to other areas across all subjects
%   is A and is used to normalise the data
%           %   We use a penalty function such that areas with fewer data
%               are penalised and given higher variance values
%               Areas that have >= 75% of data points get their
%               variance scores divided by 3
%               Areas that have >=50% of data points get their variance
%               divided by 2
%           
%   Authors: Georgios Angelis & Nicolas Darmanthé
%%--------------------------------------------------------------------------------------

%%
Labels = {'F','P','Te','O','Cd','Pt', 'Hipp','Th','Cb', 'Ins','Cing','Amy','Hyp'};  % Labels for regions
Values = cmrglchum1;    % Subjects x SUVs
         
      
%% --Calculate Difference comparison matrix ---------------------------------      
%C = zeros(size(Values,1),size(Values,2),size(Values,2)); % Subjects x Regions x Regions
C = zeros(size(Values,2),size(Values,2),size(Values,1));  % Regions x Regions x Subjects
for s=1:size(Values,1)                     % run along subjects    
    start=1;                               % start from region 1
    for i=start:size(Values,2)             % go along rows
        for j=start+1:size(Values,2)       % go along columns
            C(i,j,s) = (Values(s,i)-Values(s,j))/(Values(s,i)); % calculate absolute difference
            C(j,i,s) = (Values(s,j)-Values(s,i))/(Values(s,j)); % calculate absolute difference 
        end
        start=start+1;                     % move one up to avoid double comparison computation
    end
end

%% --Find overall variability of an area across all subjects-----------------
%C = permute(C, [2 3 1]);                %make C (area x area x subject)

for i=1:length(Values(1,:))             %remove self-comparison values in C
    C(i,i,:) = nan;
end



vars1 = [];                                                     %FvsF (col1), FvsP(col2), FvsTe etc in all subjects -> Fvar (row 1 of vars1 variable)
w8 = [];
for p = 1:length(Values(1,:))                                   %PvsF(col1), PvsP(col2), PvsTe in all  subjects -> Pvar (row 2 of vars1 variable)
    vars1 = [vars1; var((squeeze(C(p,:,:))'),'omitnan')];       %then TevsF, TevsP, TevsTe -> Tvar -- etc until all areas have been compared to every other area
    
    for i=1:length(Values(1,:))                                 %calcultes number of comparisons uses to calculate individual vars1 values
        w8 = [w8 numel(find(~isnan(squeeze(C(p,i,:)))))];      
    end
    
end
w8 = reshape(w8(:),length(Values(1,:)),[]);



%Overall variability calculation
ovvariance = [];
for i = 1:length(Values(1,:))
    if numel(find(~isnan(Values(:,i))))/numel(Values(:,1))>= 0.75
         ovvariance = [ovvariance mean(vars1(i,:),'omitnan')/2]; 
    elseif numel(find(~isnan(Values(:,i))))/numel(Values(:,1))>= 0.5
        ovvariance = [ovvariance mean(vars1(i,:),'omitnan')];
    else
        ovvariance = [ovvariance 0];                %regions that have less than 50% of values are excluded because                                                   %
                                                    %(i) would mean using <50% of review papers for normalisation 
                                                    %(ii) our overall variability calculation is not weighted per se,
                                                    %meaning that variability from areas with few datapoints
                                                    %contributes the same to the overall variability
    end
end

%dont use synthetics to synthetics comparisons in overall var calc, only synthetics to regions
for i=1:16  
    if numel(find(~isnan(Values(:,i))))/numel(Values(:,1))>= 0.75                                                    %16 is the index where synthetics finish currelty
        ovvariance(1,i) = mean(vars1(i,17:end),'omitnan')/2;  %only mean SYNvF, SYNvP, etc. not SYNvSYN
    elseif numel(find(~isnan(Values(:,i))))/numel(Values(:,1))>= 0.5
        ovvariance(1,i) = mean(vars1(i,17:end),'omitnan');
    else
        ovvariance(1,i) = 0;
    end
end

maxnum = max(ovvariance);
ovvariance = ovvariance/maxnum;



%% Find value of interest -------------------------------------------------
%   e.g.
%   subject = 2;
%   region1 = 'WB';
%   region2 = 'C';
%   idx1 = find(strcmp(Labels,region1));
%   idx2 = find(strcmp(Labels,region2));

%   mv=mean(C(:,idx1,idx2),1);
%   sv=std(C(:,idx1,idx2),0,1);

%   fprintf(1,'\n');
%   fprintf(1,'%svs%s for subject %d is: %f\n',Labels{idx1},Labels{idx2},subject,C(subject,idx1,idx2))
%   fprintf(1,'Mean and std across subjects: %.3f (%.3f)\n',mv,sv)
%%---

%threedee = [ovvariancerat' ovvariancehum' ovvariance02'];
%bar3(threedee)
      