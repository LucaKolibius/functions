% load template
nii=load_nii('Y:\Archive\derivatives\sub-0001\anat\wsub-0001_ses-postop_T1w.nii');

% get blank image
% img = zeros(size(nii.img));
img = nii.img;

% get transformation matrix
trans = [nii.hdr.hist.qoffset_x,nii.hdr.hist.qoffset_y,nii.hdr.hist.qoffset_z];

% add co-ordinates
mni_coords = [-19	0	-35;-23	0	-33;-30	1	-34;-36	2	-31;-43	2	-29;-49	3	-29];

% transform co-ordinates into matlab space
mat_coords = mni_coords + abs(trans); % look at why abs

% cycle through each coordinate
for i = 1 : size(mat_coords,1)
    
    % change coordinate to white
    img(mat_coords(i,1),mat_coords(i,2),mat_coords(i,3)) = 255;    
end

% add back into nifti
nii.img = img;
nii.hdr.dime.bitpix = 64;
nii.fileprefix = 'D:\brain' ;
save_nii(nii,'D:\brain.nii');