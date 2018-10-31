%%-------------------------------------------------------------------------
% This script sets the origin of our images to Bregma 
% (or whatever other origin you would like to set).
% 
% e.g. makebregma('/Users/nicolas/Desktop/baseline_workdir/Rat6Scan1_mask.nii', '/Users/nicolas/Desktop/baseline_workdir/Rat6Scan1_mask.nii', [0.949 0.949 0.796], [48 35 71], 16)
% The bregma coordinates [48 35 71] are given in VINCI space... The software takes care of converting them to matlab
% space
% Author: Nicolas Darmanthe
%%-------------------------------------------------------------------------

function makebregma(infile, outfile, voxelsize, origin, datatype)
V = spm_vol(infile);
Vm = spm_read_vols(V);
vincitomatorigin = [size(Vm,1) size(Vm,2)] - origin(1:2); %x and y are flipped but z stays the same. note: in elastix, z takes a negative sign.
vincitomatorigin = [vincitomatorigin origin(3)];
BregmaVolume = make_nii(Vm, voxelsize, vincitomatorigin, datatype);
nii.hdr.dime.datatype = datatype; %sets datatype for save_nii function
nii.hdr.dime.bitpix = datatype; %sets datatype for save_nii function
save_nii(BregmaVolume, outfile)