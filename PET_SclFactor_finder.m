%%-------------------------------------------------------------------------
% This scipt finds the average intensity of an image within the boundaries
% of a mask.
%
% [SclFactor] = PET_SclFactor_finder((image,mask,i) - i is used by PET_intensity_norm.m
% for diagnostic plotting and is an optional term
%
% e.g. SclFactor =
% PET_SclFactor_finder('/Users/nicolas/Desktop/baseline_FDG_workdir/Rat1Scan1/Rat1Scan1_PaxinosSpace.nii', ...
%   '/Users/nicolas/Desktop/baseline_FDG_workdir/brain_mask.nii')
%    
%
% Author: ND
%-------------------------------------------------------------------------
function [SclFactor] = PET_SclFactor_finder(image,mask,i)
                

V = spm_vol(image); %loads nifti image as structure
Vm = spm_read_vols(V); %reads array

M = spm_vol(mask); %loads nifti mask as structure
Mm = spm_read_vols(M); %reads array
Mm = logical(Mm); %converts to logical

brainvoxels = Vm(Mm); %extracts voxel values that are within mask
if nargin > 2
    subplot(ceil(14/2),2,i); histfit(brainvoxels); xlabel('Voxel Intensity'); ylabel('Frequency'); title(sprintf('Rat%d',i));
end
%[SclFactor] = mean(brainvoxels)/100000; 
[SclFactor] = mean(brainvoxels); 
end

