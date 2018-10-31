%%-------------------------------------------------------------------------
% This script crops images based on boundary coordinates. These coordinates
% are visually determined using VINCI. 
% - We assume that - (see maskmake.m notes for details)
%       x-dim is left to right
%       y-dim is anterior to posterior (rostral to caudal)
%       z-dim is superior to inferior (dorsal to ventral)
% note: sometimes, y dim and z dim are inverted in an image file, such that 
% y-dim is superior to inferior and z-dim is anterior to posterior.
% Use flipper to permute/reorient the dimensions z <-> y
%
% Author: Nicolas Darmanthe
%%-------------------------------------------------------------------------

function cropper(infile, xbound, ybound, zbound, outfile, voxelsize, origin, datatype)
% Inputs:
%   infile - is a string for input filename
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
%   origin - is a vector containing the VINCI coordinates of the origin (e.g. Bregma)
%   x, y and z e.g. [60 64 62]
%
%   datatype - is data format of outputfile. Takes values
%             2 - uint8,  4 - int16,  8 - int32,  16 - float32,
%             64 - float64,  128 - RGB24
%%-------------------------------------------------------------------------

V = spm_vol(infile);
Vm = spm_read_vols(V);

Vmc = Vm(size(Vm,1)-xbound(1,2):size(Vm,1)-xbound(1,1),:,:);
Vmc = Vmc(:,size(Vm,2)-ybound(1,2):size(Vm,2)-ybound(1,1),:);
Vmc = Vmc(:,:,zbound(1,1):zbound(1,2));

vincitomatorigin = [size(Vm,1) size(Vm,2)] - origin(1:2); %x and y are flipped but z stays the same. note: in elastix, z takes a negative sign.
vincitomatorigin = [vincitomatorigin origin(3)];
V = make_nii(Vmc, voxelsize, vincitomatorigin, datatype);

nii.hdr.dime.datatype = datatype; %sets datatype for save_nii function
nii.hdr.dime.bitpix = datatype; %sets datatype for save_nii function
save_nii(V, outfile); %saves image