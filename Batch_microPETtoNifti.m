%------------------------------------------------------------------------
%   This script bulk-converts binary microPET images to nifti images      
%
%   Authors: ND
%------------------------------------------------------------------------
%%
Studies = {'Rat1Scan1';'Rat2Scan1';'Rat3Scan1';'Rat4Scan1';'Rat5Scan1';'Rat7Scan1';'Rat8Scan1';'Rat9Scan1';...
    'Rat1Scan2';'Rat2Scan2';'Rat3Scan2';'Rat4Scan2';'Rat5Scan2';'Rat6Scan2';'Rat7Scan2';'Rat8Scan2';'Rat9Scan2'};


pth='/Users/nicolas/Desktop/baseline_FDG_workdir';
%pth = '/Users/nicolas/Desktop/debug';

for s = 1:13
    filename = sprintf('/static/%s_static_recon3D_psf_frame2_it10.img',Studies{s});
    infile = fullfile(pth,Studies{s},filename);
    outfile = fullfile(pth,Studies{s},sprintf('%s.nii',Studies{s}));
    microPETtoNifti(infile, outfile, [128 128 95], [0.949 0.949 0.796], [0 0 0], 16)
    fprintf(1, '\n %s: Processing frame: xx',Studies{s});
    for f = 1:30
        fprintf(1, '\b\b%2d', f);
        filename = sprintf('120x30/%s_dynamic_recon3D_psf_frame%d_it10.img',Studies{s},f);
        infile = fullfile(pth,Studies{s},filename);
        outfile = fullfile(pth,Studies{s},sprintf('%s_dynamic%d.nii',Studies{s},f));
        microPETtoNifti(infile, outfile, [128 128 95], [0.949 0.949 0.796], [0 0 0], 16)
    end
end
fprintf(1,'\n');