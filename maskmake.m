%%------------------------------------------------------------------
% This script is intended to create masks of our PET images. While this
% process is relatively easy in humans using toolboxes such as FSL BET
% there does not appear to be a standard way of doing this in rats.
%           see https://www.nature.com/articles/s41598-017-00148-1
% The user is required to specify a threshold: voxel values in the image
% above the threshold are selected and their coordinates are used to create
% the mask.
%
% The mask will have the same dimensions (and be coregistered to voxelwise)
% as the original image. Voxels in the PET image above the threshold will
% be given a value of 1, and voxels bellow the threshold a value of 0. This
% data will be saved in a matrix.
%
% Finally, we clean up the matrix such that voxels which reached the
% threshold but are clearly outside of the brain are not retained for 
% mask creation (see xbound, ybound and zbound)
%
% - We assume that -
%       x-dim is left to right
%       y-dim is anterior to posterior (rostral to caudal)
%       z-dim is superior to inferior (dorsal to ventral)
% note: sometimes, y dim and z dim are inverted in an image file, such that 
% y-dim is superior to inferior and z-dim is anterior to posterior.
% Use flipper to permute the dimensions z <-> y
%
% Author: Nicolas Darmanthe
%%------------------------------------------------------------------------

function maskmake(infile, threshold, crop, xbound, ybound, zbound, outfile, voxelsize, origin, datatype)
% Inputs:
%   infile - is a string for input filename
%   
%   threshold - takes mask creation threshold value
%
%   crop - takes value 0 (no crop) or 1 (crop). This removes elements which are not
%          part of the brain out of the mask
%
%   xbound - is a 2-element vector which contains the brain boundaries within
%           the image volume in the x dimension ... [left right] (as seen in VINCI)
%   
%   ybound - is a 2-element vector which contains the brain boundaries within
%           the image volume in the y dimension [anterior posterior]
%
%   zbound - is a 2-element vector which contains the brain boundaries within
%           the image volume in the z dimension... [superior inferior]
%
%   outfile - is a string for output filename
%
%   voxelsize - is a vector containing the size of voxels in each dimension
%               x, y and z e.g. [0.2 0.2 0.2]
%
%   origin - is a vector containing the matrix coordinates of the origin (e.g. Bregma)
%   x, y and z e.g. [60 64 62]
%
%   datatype - is data format of outputfile. Takes values
%             2 - uint8,  4 - int16,  8 - int32,  16 - float32,
%             64 - float64,  128 - RGB24
%
% e.g. PETmaskmake('/Users/nicolas/Desktop/test.nii', 400000, 1, [47 70], [46 76], [32 49], '/Users/nicolas/Desktop/processedtest.nii', [0.949 0.949 0.796], [60 67 62], 4)
%%-------------------------------------------------------------------------

V = spm_vol(infile); %loads nifti as structure
Vm = spm_read_vols(V); %extracts matrix from structure such that left-right (X), 
                       %anterior-posterior (Y) and up-down (Z) are flipped
                       %relative to image coordinates/how the images are 
                       %displayed in VINCI, hence*
%
%      ?     /-------/|
%           /       / |
%  128 0   /-------/  |
%          |       |  |   ?
%          |       |  / 95 95
%  --Y--   |       | /
%          |       |/         --Z--
%   0 128   -------/ 0  0
%
%    ?  0      128
%      128      0
%
%           --X--
%
% ARROWS ? indicate how VINCI displays coordinates
% The associated row/column represents where the VINCI coordinates are
% mapped in the imported matrix/array
%
%              z                   z                                                    
%  --> x       ^                   ^  
% |            |                   |                VINCI orthoview pannels
% ?             --> x               --> y
% y
%
%
%      y                          z
%      ^       x <--              ^
%      |            |             |                 MATLAB array
% x <--             ?              -->y   
%                   z              
%------------------------------------------------------------------------                       
                       
mask = Vm >= threshold; %this is a logical
mask = double(mask); %convert to double

if crop == 1
    mask(1:size(mask,1)-xbound(1,2),:,:) =0;            %makes right side of mask outside boundary 0
                                                        %(by cropping left side of matrix *)
    mask(size(mask,1)-xbound(1,1):size(mask,1),:,:)=0;  %makes left side of mask outside boundary 0
    
    mask(:,1:size(mask,2)-ybound(1,2),:)=0;             %makes posterior side of mask outside boundary 0
    mask(:,size(mask,2)-ybound(1,1):size(mask,2),:)=0;  %makes anterior side of mask outside boundary 0
  
    mask(:,:,1:zbound(1,1))=0;                          %makes inferior side of mask outside boundary 0
    mask(:,:,zbound(1,2):size(mask,3))=0;               %%makes superior side of mask outside boundary 0
  
end

vincitomatorigin = [size(Vm,1) size(Vm,2)] - origin(1:2); %x and y are flipped but z stays the same. note: in elastix, z takes a negative sign.
vincitomatorigin = [vincitomatorigin origin(3)];
maskvolume = make_nii(mask, voxelsize, vincitomatorigin, datatype); %creates mask image volume
nii.hdr.dime.datatype = datatype; %sets datatype for save_nii function
nii.hdr.dime.bitpix = datatype; %sets datatype for save_nii function
save_nii(maskvolume, outfile)


