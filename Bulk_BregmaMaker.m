%%------------------------------------------------------------------------%
%   This scirpt bulk-sets image origins
%
%   Authors: ND and GA
%%------------------------------------------------------------------------%

Studies = {'Rat1Scan1';'Rat2Scan1';'Rat3Scan1';'Rat4Scan1';'Rat5Scan1';...
    'Rat6Scan1';'Rat7Scan1';'Rat8Scan1';'Rat9Scan1';'Rat1Scan2';...
    'Rat2Scan2';'Rat3Scan2';'Rat4Scan2';'Rat5Scan2';'Rat6Scan2';...
    'Rat7Scan2';'Rat8Scan2';'Rat9Scan2'}; 
pth='/Users/nicolas/Desktop/baseline_FDG_workdir';

for s = 16:18
    filename = sprintf('WB_minus/%s_PaxinosSpace_WB_minus.nii',Studies{s});
    infile = fullfile(pth,Studies{s},filename);
    makebregma(infile, infile, [0.2 0.2 0.2], [48 35 71], 16)
    fprintf(1, 'Processing frame: xx');
    for f = 1:29
        fprintf(1, '\b\b%2d', f);
        filename = sprintf('WB_minus/%s_dynamic%d_PaxinosSpace_WB_minus.nii',Studies{s},f);
        infile = fullfile(pth,Studies{s},filename);
        makebregma(infile, infile, [0.2 0.2 0.2], [48 35 71], 16)
    end
end
fprintf(1,'\n');