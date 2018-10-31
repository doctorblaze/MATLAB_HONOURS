%------------------------------------------------------------------------
%   This script converts binary microPET images to nifti images           
%                                                                        
%   Author: Nicolas Darmanthe                                             
%   
%   Usage: microPETtoNifti(infile, outfile, dimensions, voxelsize, origin, datatype)
%           infile: a string containing output file path 
%           outfile: a string containing output file path 
%           dimensions: a vector containing image volume dimensions along [X Y Z]
%           voxelsize: a vector containing image voxel size along [X Y Z]
%           origin: a vector containing the coordinates of the origin [X Y Z]. 
%                   coordinates are input as seen in VINCI (Max Planck Institute)
%           datatype: defines data format of outputfile. Takes values
%           2 - uint8,  4 - int16,  8 - int32,  16 - float32,
%           64 - float64,  128 - RGB24
%           
%   e.g. microPETtoNifti('/Users/nicolas/Desktop/workdir/Rat6Scan1.img', '/Users/nicolas/Desktop/workdir/Rat6Scan1.nii', [128 128 95], [0.949 0.949 0.796], [48 35 71], 16)
%------------------------------------------------------------------------

%%
function microPETtoNifti(infile, outfile, dimensions, voxelsize, origin, datatype)

xdim = dimensions(1);
ydim = dimensions(2);
zdim = dimensions(3);

% Load input image binary file
fid = fopen(infile);

% Read binary and shape into matrix
I = reshape(fread(fid, Inf, 'float32'), xdim, ydim, zdim);

% Cleanup for supervoxels outside of brain.... we set voxels outside of brain area to 0                                            
xbound = [41 89]; %empirically determined, but brain always seems to fall within these boundaries
ybound = [6 80];
zbound = [12 71];
I(1:size(I,1)-xbound(1,2),:,:) = 0;
I(:,1:size(I,2)-ybound(1,2),:) = 0;
I(:,:,1:zbound(1,1)) = 0;
I(size(I,1)-xbound(1,1):end,:,:) = 0;
I(:,size(I,2)-ybound(1,1):end,:) = 0;
I(:,:,zbound(1,2):end) = 0;

% Close input image binary file
fclose(fid);

% Make image volume
nii = make_nii(I, voxelsize, origin, datatype);

%Save image volume
nii.hdr.dime.datatype = datatype; %sets datatype for save_nii function
nii.hdr.dime.bitpix = datatype; %sets datatype for save_nii function
save_nii(nii, outfile) 