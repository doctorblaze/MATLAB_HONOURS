%%-------------------------------------------------------------------------
% This scipt uses the PET_SclFactor_finder factor/or IDxweight to scale the 
% voxel intensities of associated images. Conceptually analogous to SUV 
% scaling.
%
% Author: ND
%-------------------------------------------------------------------------
%% Uses PET_SclFactor_finder scale factor, norm images are unitless

Studies = {'Rat1Scan1';'Rat2Scan1';'Rat3Scan1';'Rat4Scan1';'Rat5Scan1';...
    'Rat6Scan1';'Rat7Scan1';'Rat8Scan1';'Rat9Scan1';'Rat1Scan2';...
    'Rat2Scan2';'Rat3Scan2';'Rat5Scan2';'Rat6Scan2';...
    'Rat7Scan2';'Rat8Scan2';'Rat9Scan2'}; 
    %no Rat4Scan2 because injection was in 3 gos
    %Rat6Scan1 has no dynamic
    %Rat3Scan1 is a strange scan (registration is also possibly faulty)

pth='/Users/nicolas/Desktop/baseline_FDG_workdir';
mask = '/Users/nicolas/Desktop/baseline_FDG_workdir/brain_mask.nii';
datatype = 16;

% Uses PET_SclFactor_finder.m and determines the average voxel intensity
% in NORM120 within the whole brain mask (also plots intensity distribution, 
% which is typically normal)
figure(1); title('Distribution of Vox Intensities')
for s = 1:14
    filename = sprintf('%s_NORM120_PaxinosSpace.nii',Studies{s});
    infile = fullfile(pth,Studies{s},filename);
    SclFactors(s).name = Studies{s};
    SclFactors(s).SclF = PET_SclFactor_finder(infile, mask,s);
end
fprintf(1,'\n');

%saves the scaled images
for s = 7:14
    filename = sprintf('%s_PaxinosSpace.nii',Studies{s});
    infile = fullfile(pth,Studies{s},filename);
    newdir = 'NORM120';
    mkdir(fullfile(pth,Studies{s}),newdir);
    newfilename = sprintf('%s_PaxinosSpace_NORM120.nii',Studies{s});
    outfile = fullfile(pth,Studies{s},newdir,newfilename);
    V = spm_vol(infile); %loads nifti image as structure
    Vm = spm_read_vols(V); %reads array
    Vm = Vm./(SclFactors(s).SclF); %applies scaling
    ScldVol = make_nii(Vm, [0.2 0.2 0.2], [48 85 71], datatype); %make nifti volume
    nii.hdr.dime.datatype = datatype; %sets datatype for save_nii function
    nii.hdr.dime.bitpix = datatype; %sets datatype for save_nii function
    save_nii(ScldVol, outfile)
    fprintf(1, 'Processing frame: xx');
    for f = 1:30
        fprintf(1, '\b\b%2d', f);
        filename = sprintf('%s_dynamic%d_PaxinosSpace.nii',Studies{s},f);
        infile = fullfile(pth,Studies{s},filename);
        newfilename = sprintf('%s_dynamic%d_PaxinosSpace_NORM120.nii',Studies{s},f);
        outfile = fullfile(pth,Studies{s},newdir,newfilename);
        V = spm_vol(infile); %loads nifti image as structure
        Vm = spm_read_vols(V); %reads array
        Vm = Vm./(SclFactors(s).SclF); %apply scaling
        ScldVol = make_nii(Vm, [0.2 0.2 0.2], [48 85 71], datatype); %make nifti volume
        nii.hdr.dime.datatype = datatype; %sets datatype for save_nii function
        nii.hdr.dime.bitpix = datatype; %sets datatype for save_nii function
        save_nii(ScldVol, outfile) 
    end
end

%% Uses ID/weight as scaling factor (image becomes SUV image or %ID/mL)
Studies = {'Rat1Scan1';'Rat2Scan1';'Rat3Scan1';'Rat4Scan1';'Rat5Scan1';...
    'Rat6Scan1';'Rat7Scan1';'Rat8Scan1';'Rat9Scan1';'Rat1Scan2';...
    'Rat2Scan2';'Rat3Scan2'; 'Rat4Scan2';'Rat5Scan2';'Rat6Scan2';...
    'Rat7Scan2';'Rat8Scan2';'Rat9Scan2'};
%note: ID are "known" to be much cleaner for Rats7,8,9 (and perhaps 5_2 and 1_2 - Gita magic hand)
%i.e. can be used for anaesthetic comparison

datatype = 16;
pth='/Users/nicolas/Desktop/baseline_FDG_workdir';
mask = '/Users/nicolas/Desktop/baseline_FDG_workdir/brain_mask.nii';

%import injected doses and weights as vector in same order as Studies list
%convert injected dose from MBq to Bq (because our images are in Bq/mL)
IDinBq = InjectedDose*10^6;
IDdividedbyweight = IDinBq./weight;


for s = 16:18
    filename = sprintf('%s_PaxinosSpace.nii',Studies{s});
    infile = fullfile(pth,Studies{s},filename);
    %newdir = 'ID';
    newdir = 'SUV';
    mkdir(fullfile(pth,Studies{s}),newdir);
    %newfilename = sprintf('%s_PaxinosSpace_ID.nii',Studies{s});
    newfilename = sprintf('%s_PaxinosSpace_SUV.nii',Studies{s});
    outfile = fullfile(pth,Studies{s},newdir,newfilename);
    V = spm_vol(infile); %loads nifti image as structure
    Vm = spm_read_vols(V); %reads array
    %Vm = (Vm./IDinBq(s))*100;
    Vm = (Vm./IDdividedbyweight(s))*100; %applies scaling
    ScldVol = make_nii(Vm, [0.2 0.2 0.2], [48 85 71], datatype); %make nifti volume
    nii.hdr.dime.datatype = datatype; %sets datatype for save_nii function
    nii.hdr.dime.bitpix = datatype; %sets datatype for save_nii function
    save_nii(ScldVol, outfile)
    fprintf(1, '\n Processing frame: xx');
    for f = 1:29
        fprintf(1, '\b\b%2d', f);
        filename = sprintf('%s_dynamic%d_PaxinosSpace.nii',Studies{s},f);
        infile = fullfile(pth,Studies{s},filename);
        newfilename = sprintf('%s_dynamic%d_PaxinosSpace_ID.nii',Studies{s},f);
        %newfilename = sprintf('%s_dynamic%d_PaxinosSpace_SUV.nii',Studies{s},f);
        outfile = fullfile(pth,Studies{s},newdir,newfilename);
        V = spm_vol(infile); %loads nifti image as structure
        Vm = spm_read_vols(V); %reads array
        Vm = Vm./(IDinBq(s))*100;
        %Vm = Vm./(IDdividedbyweight(s))*100; %apply scaling
        ScldVol = make_nii(Vm, [0.2 0.2 0.2], [48 85 71], datatype); %make nifti volume
        nii.hdr.dime.datatype = datatype; %sets datatype for save_nii function
        nii.hdr.dime.bitpix = datatype; %sets datatype for save_nii function
        save_nii(ScldVol, outfile) 
    end
end

%% Uses the whole brain static scan activity as scaling factor (regional values will be RUV_wholebrain)
% We have shown that behavior during the scan does not affect FDG kinetic rates and does not cause an increase in FDG metabolism
% in areas such as the motor/somatosensory cortex/..
% It is thus "justifiable" for us to make the claim "activity in the "whole brain"* should be the same between all animals"
% regardless of their behavior during the scan
% We * exclude the hippocampus, mPFC, amygdala when determining the "whole brain" scaling factor, because these regions
% are certainly NOT expected to be the same between animals (esp. pre vs post stress). see FDG_CMRglc_calc paper

Studies = {'Rat1Scan1';'Rat2Scan1';'Rat3Scan1';'Rat4Scan1';'Rat5Scan1';...
    'Rat6Scan1';'Rat7Scan1';'Rat8Scan1';'Rat9Scan1';'Rat1Scan2';...
    'Rat2Scan2';'Rat3Scan2';'Rat4Scan2';'Rat5Scan2';'Rat6Scan2';...
    'Rat7Scan2';'Rat8Scan2';'Rat9Scan2'};

pth='/Users/nicolas/Desktop/baseline_FDG_workdir';
mask = '/Users/nicolas/Desktop/workdir/WB.nii';
datatype = 16;
dynamic = 'off'

% Uses PET_SclFactor_finder.m and determines the average voxel intensity
% in STATIC within the WB_minus_affected_mask 

for s = 1:18
    filename = sprintf('%s_PaxinosSpace.nii',Studies{s});
    infile = fullfile(pth,Studies{s},filename);
    SclFactors(s).name = Studies{s};
    SclFactors(s).SclF = PET_SclFactor_finder(infile, mask);
end
fprintf(1,'\n');

%saves the scaled images
for s = 1:18 %start by scaling the static frames
    filename = sprintf('%s_PaxinosSpace.nii',Studies{s});
    infile = fullfile(pth,Studies{s},filename);
    newdir = 'WB';
    mkdir(fullfile(pth,Studies{s}),newdir);
    newfilename = sprintf('%s_PaxinosSpace_WB.nii',Studies{s});
    outfile = fullfile(pth,Studies{s},newdir,newfilename);
    V = spm_vol(infile); %loads nifti image as structure
    Vm = spm_read_vols(V); %reads array
    Vm = Vm./(SclFactors(s).SclF); %applies scaling
    ScldVol = make_nii(Vm, [0.2 0.2 0.2], [48 85 71], datatype); %make nifti volume, with [48 35 71] bregma origin
    nii.hdr.dime.datatype = datatype; %sets datatype for save_nii function
    nii.hdr.dime.bitpix = datatype; %sets datatype for save_nii function
    save_nii(ScldVol, outfile) %saves the new scaled static image
    if strcmp(dynamic,'on')
        fprintf(1, 'Processing frame: xx'); 
        for f = 1:30   %then scale the dynamic frames
            fprintf(1, '\b\b%2d', f);
            filename = sprintf('%s_dynamic%d_PaxinosSpace.nii',Studies{s},f);
            infile = fullfile(pth,Studies{s},filename);
            newfilename = sprintf('%s_dynamic%d_PaxinosSpace_WB.nii',Studies{s},f);
            outfile = fullfile(pth,Studies{s},newdir,newfilename);
            V = spm_vol(infile); %loads nifti image as structure
            Vm = spm_read_vols(V); %reads array
            Vm = Vm./(SclFactors(s).SclF); %apply scaling
            ScldVol = make_nii(Vm, [0.2 0.2 0.2], [48 85 71], datatype); %make nifti volume
            nii.hdr.dime.datatype = datatype; %sets datatype for save_nii function
            nii.hdr.dime.bitpix = datatype; %sets datatype for save_nii function
            save_nii(ScldVol, outfile) %saves the new scaled dynamic image
        end
    end
end
