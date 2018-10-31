%% ------------------------------------------------------------------------
% (First, DICOM MRS images are saved as nii using VINCI)
% This script permutes/reorients the y and z axes of the the MRI images acquired 
% on the 3T MR Solutions scanner using axial acquisition settings (transverse 
% rat slices) to match orientation the SAMIT/Schwartz Atlas in Paxinos Space.
% This step is necessary in order to proceed to the next step; MRI image
% (moving image) co-registration/warping with the SAMIT atlas (fixed image)
% Ideally, this step will be performed automatically in elastix... 
% 
% Author: Nicolas Darmanthe
%% ------------------------------------------------------------------------

function flipper(infile, outfile, voxelsize, origin, datatype)

% Inputs:
%   infile - is a string for input filename
%
%   outfile - is a string for output filename
%
%   voxelsize - is a vector containing the size of voxels in each dimension
%               x, y and z e.g. [0.2 0.2 0.2]
%
%   origin - is a vector containing the matrix coordinates of the origin (Bregma)
%   x, y and z e.g. [60 64 62]
%
%   datatype - is data format of outputfile. Takes values
%             2 - uint8,  4 - int16,  8 - int32,  16 - float32,
%             64 - float64,  128 - RGB24
%
% e.g. flipper('/Users/nicolas/Desktop/test.nii', '/Users/nicolas/Desktop/processedtest.nii', [0.2 0.6 0.2], [50 42 32], 16)
%%-------------------------------------------------------------------------

%%-Flip matrix-------------------------------------------------------------
V = spm_vol(infile); %imports nifti as structure
Vm = spm_read_vols(V); %makes matrix 


Vmp = permute(Vm, [1 3 2]); %permutes y and z dims to match template
Vmr = fliplr(Vmp); %when the dimensions are permuted, the brain is upside down
                   %This step puts the brain back up right

%%-Make new nifti structure------------------------------------------------
nii = make_nii(Vmr, voxelsize, origin, datatype); %makes nifti structure

%%-Export images-------------------------------------------------------------
nii.hdr.dime.datatype = datatype; %sets datatype for save_nii function
nii.hdr.dime.bitpix = datatype; %sets datatype for save_nii function
save_nii(nii, outfile); %saves image

