%get all images

clear all;


%get all filenames in directory 
mydir = 'scans';
filelist = getAllFileNames(mydir);


%parameters for all files
doplot = false;
storeM = 8;%size of subimages stored before cropping
M = 16; %resized M x M subimages after cropping
c =0;%amount to crop from EACH edge
mytol = 252; %tolerance for mean value of image

%intialize holder for all subimages
bigX = cell(4,1);

%loop thruogh all file names
for ii = 1:numel(filelist)
    filename = filelist{ii};
    
    %extract subimages
    process_image;
    % pause(0.1);
    %see if there are x's

    
    if size(X,1)==232
    bigX{ii} = X;
    end

    %information about the ith image is in the ith row of X
end

%save to file (images_scans & images_stef_scans)

    