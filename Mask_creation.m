%=====================================================
%=====================================================
% Example 1:  Ripping out BA 10 from a Brodmann Area Atlas

% Step 1:  Load in the Atlas
BrodmannAtlas = spm_read_vols(spm_vol('/Users/mvlombardo/Documents/BrainAtlases/brodmann.img'));  % Reads in the brodmann atlas into a variable called BrodmannAtlas

% Step 2:  Make a mask that extracts BA 10.  I know that the Brodmann atlas has its voxels labeled based on the Brodmann Area number. I simply have to extract all voxels whose value is 10
VOI = 10;
BA10 = ismember(BrodmannAtlas,VOI);  
% ismember is a crucial function. Give it a matrix in the first argument and tell it the value you are looking for in that matrix within the second argument. 
% Note that ismember takes BrodmannAtlas, assigns a logical value of "true" (or 1) to voxels whose value is 10, and assigns all other voxels a logical value of false (or 0). 

% Step 3:  Write out BA 10 as its own image.
% Load the header information of another file of similar dimensions and voxel sizes
V = spm_vol('/Users/mvlombardo/Documents/BrainAtlases/brodmann.img');

% Change the name in the header
V.fname = 'BA10.nii';
V.private.dat.fname = V.fname;

% Write out the new header and data
spm_write_vol(V,BA10);
%=====================================================
%=====================================================


% Example 2:  Extracting all voxels whose value is above a specified t-threshold

% Step 1:  Load spmT*.img which is a map containing t-values from a statistical contrast
TMap = spm_read_vols(spm_vol('/Users/mvlombardo/Documents/fMRI/Contrast1/spmT_0001.img'));
tthresh = 3.1768423;  % This is the corresponding t-threshold I'll need to extract all voxels whose value is greater than this.

% Step 2:  Extract all voxels in TMap that are greater than tthresh
TMap_suprathresh = TMap>=tthresh;  
% Note that TMap_suprathresh is a matrix of the same size as TMap filled with 0 or 1, where voxels greater than tthresh are 1, otherwise they are 0

% Step 3:  Write out TMap_suprathresh as its own image.
V = spm_vol('/Users/mvlombardo/Documents/fMRI/Contrast1/spmT_0001.img');
V.fname = 'tmap_suprathresh.nii';
V.private.dat.fname = V.fname;
spm_write_vol(V,TMap_suprathresh);
%=====================================================
%=====================================================


% Example 3:  Making a conjunction map between two suprathreshold t-maps

% Read in both t-maps
TMap1 = spm_read_vols(spm_vol('/Users/mvlombardo/Documents/fMRI/Contrast1/spmT_0001.img'));
TMap2 = spm_read_vols(spm_vol('/Users/mvlombardo/Documents/fMRI/Contrast1/spmT_0002.img'));

% Define t-threshold
tthresh = 3.1768423;

% Mask out suprathreshold voxels from both t-maps
TMap1_suprathresh = TMap1>=tthresh;
TMap2_suprathresh = TMap2>=tthresh;

% Make conjunction map as the logical AND of both suprathreshold t-maps
ConjunctionMap = TMap1_suprathresh & TMap2_suprathresh;

% Write out the new conjunction map
V = spm_vol('/Users/mvlombardo/Documents/fMRI/Contrast1/spmT_0001.img');
V.fname = 'ConjunctionMap_TMap1_AND_TMap2.nii';
V.private.dat.fname = V.fname;
spm_write_vol(V,ConjunctionMap);
%=====================================================
%=====================================================