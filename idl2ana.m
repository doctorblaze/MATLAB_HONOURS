function [Out] = idl2ana(infile, voxels, outfile, type, range, format)
%
% Inputs: 
%   infile  - is a string for input filename
%   voxels  - is a vector which contains the number of voxel in each
%             dimension, e.g. [128 128 85 1]
%   outfile - is a string for output filename
%   type    - defines the type of image output. It takes two values and it 
%             can either 'dynamic' or 'static'
%   range   - is 2-element vector that defines the range of frames saved 
%             for dynamic or summed for static
%   format  - is data format of outputfile. Takes values
%             2 - uint8,  4 - int16,  8 - int32,  16 - float32,
%             64 - float64,  128 - RGB24
%
% Example:
% idl2ana('/Users/nicolas/Desktop/out.img',[40 39 48 22],'/Users/nicolas/Desktop/testana.img', 'static', [10 18], 4);
%
% ----------------------------------------
% Author: Giorgos Angelis
%% Parse input arguments ------------------------------------
   
    datatype = 4; % this is 4-bytes or 4*8=32 bits float
    vx = voxels(1);
    vy = voxels(2);
    vz = voxels(3);
    vt = voxels(4);
    
%% Load input image 
            
    fileInfo = dir(infile);
    fileSize = fileInfo.bytes;
    if(fileSize~=vx*vy*vz*vt*datatype)
        fprintf(2, ' ERROR! Inconsistent number of input voxels. Please check your input data\n');
        Out = 0;
        return;
    end

    fprintf(1, ' Loading image file: %s\n', infile);
    fid = fopen(infile);
    I = reshape(fread(fid, Inf, 'float'), vx, vy, vz, vt); %reshapes binary into image matrix
    fclose(fid);
    fprintf(1, ' Done!\n');

%% View loaded image ---------------
    figure; 
    for i=1:vt
        imagesc(I(:,:,20,i)'); axis image off; 
        title(sprintf('Frame %d',i)); drawnow
        pause(0.05); 
    end
    
%% Convert to Analyze ----------------------

    if(strcmp(type,'dynamic'))
        fprintf(1, ' Saving dynamic frames %d to %d\n', range(1), range(2));
        Out = I(:,:,:,range(1):range(2));
    elseif(strcmp(type,'static'))
        fprintf(1, ' Saving summed static frames %d to %d\n', range(1), range(2));
        Out = sum(I(:,:,:,range(1):range(2)),4);
    else
        fprintf(2, ' WARNING! Incorrect type. Saving entire image as dynamic!\n');
        Out = I;
    end

    fprintf(1, ' Exporting output image file %s\n',outfile);
    ana = make_ana(Out, [0.5 0.5 1], [0.0 0.0 0.0], format); %voxel size [0.5 0.5 1] - must match image matrix size, otherwise stretched images..
    save_untouch_nii(ana, outfile);
    fprintf(1, ' Done!\n');
    
end