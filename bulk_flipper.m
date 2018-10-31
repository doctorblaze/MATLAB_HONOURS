Studies = {'Rat1Scan1';'Rat2Scan1';'Rat3Scan1';'Rat4Scan1';'Rat5Scan1';'Rat6Scan1'};
pth='/Users/nicolas/OneDrive/Documents/Honours/X/3.2baseline_MRI/PET_PTSD_baseline_nifti';

for s = 1:6
    filename = sprintf('%s.nii',Studies{s});
    infile = fullfile(pth,filename);
    flippedfilename = sprintf('%s_flipped.nii',Studies{s});
    outfile = fullfile(pth,flippedfilename);
    flipper(infile, outfile, [0.15625 0.5972190 0.15625], [0 0 0], 16)
    fprintf(1, 'Processing frame: xx');
end
fprintf(1,'\n');