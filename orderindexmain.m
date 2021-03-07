%% Here to calculate the order index of fiber-like structures in a 3D context
% This is the main program
% Part of the code is modified from Refs.: Liu et al., Biomedical Optics
% Express,6, 2015; Liu et al., Biomaterials, 116, 2017; Liu et al.,
% Biomaterials, 179, 2018
% Here to load images and create a 3D stack
clear;
sz = 512;
sz2 = 512;
sz3 = 80; % modify 1, define sz, sz2, and sz3 according to the size of the 3D image

fiberstack = zeros(sz,sz2,sz3);

for i = 1:sz3
    fiberstack(:,:,i) = imread(['D:\Examples\',num2str(i),'.tif']); % modify 2 
    % 'fiberstack' is the 3D image used for order index analysis.
    % Define the path of the folder where the images are saved
end

%% Here to create the mask selecting fiber regions

threshint = 3;  % modify 3 
% 'threshint' is the raw background intensity acquired by averaging intensity of several regions identified as background
maska3d = fiberstack > threshint; % maska3d is a raw binary mask based on the raw background intensity

meanint = mean(fiberstack(maska3d)); % calculate the mean intensity within the raw binary mask
threshinta = 0.45*meanint; % modify 4 
% 'threshinta' is an inmproved intensity threshold, and the factor '0.45' can be adjusted according to different sample
maska3da = fiberstack > threshinta; % 'maska3da' is an improved binary mask

% Here to create the finalmask
finalmaska = maska3d.*maska3da; % 'finalmaska' is the binary mask acquired based on the above two masks

finalmask = zeros(sz,sz2,sz3);
for j = 1:sz3
    Maskf1 = finalmaska(:,:,j); 
    Maskf2 = tinyremoval(Maskf1,15); % 'tinyremoval.m' is a function which is able to remove very tiny structures
    finalmask(:,:,j) = Maskf2;
end

finalmask = logical(finalmask); % 'finalmask' is the final binary mask selecting the fiber-only regions

%% Here to calculate the voxel-wise 3D orientation
armdx = 4; % modify 5 
% the width of window in both 'x' and 'y' dimensions is both '2*armdx+1'.
% 'armdx' is chosen based on the diamater of fibers. Typically '2*armdx+1'
% is 2 to 3 times the diameter of fiber so as to provide optimal accuracy
armdz = 5; % modify 6
% '2*armdz+1' is the width of the window in 'z' dimension. We would expect
% the actual width of the window in 'x', 'y' dimensions is equal to the one
% in 'z' dimension. Therefore, 'armdz' is determined by the sampling 
% frequency in the 'xy' dimension and 'z' dimension  
filton = 1; % modify 7
% 'filton' is a binary variable that tells program whether to filter the
% data prior to analysis. 1: filtering; 0: no filtering
para = 0.65; % modify 8
% 'para' is the ratio of sampling frequency between 'xy' dimension and 'z'
% dimension
[aSDE,pSDER,bSDER,gSDER] = fiberangle3D(fiberstack,armdx,armdz,filton,para); 
% 'calcfibang3D.m' is the function which calculates the voxel-wise 3D
% orientation of a 3D image. 'aSDE' is the calculated theta stack, 'pSDER'
% is the calculated phi stack, 'bSDER' is the calculated beta stack, and
% 'gSDER' is the calculated gamma stack

%% Here to calculate the voxel-wise 3D order index based on a localized window
bwhole = sqrt(1./(tan(2*bSDER*pi/180).^2)+1./(tan(2*gSDER*pi/180).^2));
Cwhole = (bwhole./sqrt(1+bwhole.^2)).*cos(2*aSDE*pi/180);
Swhole = (bwhole./sqrt(1+bwhole.^2)).*sin(2*aSDE*pi/180);
Zwhole = zeros(sz,sz2,sz3);
% 'bwhole', 'Cwhole', 'Swhole' and 'Zwhole' are all assisting variables
% used to acquire the order index results.  
for m = 1:sz
   for n = 1:sz2
     for o = 1:sz3
        if bSDER(m,n,o) <= 90   
            Zwhole(m,n,o) = 1/sqrt(1+bwhole(m,n,o)^2);
        else
            Zwhole(m,n,o) = -1/sqrt(1+bwhole(m,n,o)^2);
        end
     end
   end
end

Cwholemean = zeros(sz,sz2,sz3);
Swholemean = zeros(sz,sz2,sz3);
Zwholemean = zeros(sz,sz2,sz3);
ll = 0;

armdxx = 4; % modify 9
armdzz = 5; % modify 10
% 'armdxx' and 'armdzz' are defining the size of the window used to
% acquired the 3D order index (xy: 2*armdxx+1; z: 2*armdzz+1).
% Typically, these two parameters are defined the same as 'armdx' and
% 'armdz', so that the window used to calculate the voxel-wise orientation 
% and the one used to acquire localized order index have the same 
% size. However, they can be set to other values according to the purpuse 
% of analysis 

h = waitbar(0,'Please Wait');
ijk = 0;
nijk = (2*armdxx+1)*(2*armdxx+1)*(2*armdzz+1);
tic

for  i = -armdxx:armdxx
    for j = -armdxx:armdxx
        for k = -armdzz:armdzz
            ijk = ijk+1;
            waitbar(ijk/nijk)
            Cim = circshift(Cwhole,[i j k]);
            Sim = circshift(Swhole,[i j k]);
            Zim = circshift(Zwhole,[i j k]);
                        
            Cwholemean = Cwholemean+Cim;
            Swholemean = Swholemean+Sim;
            Zwholemean = Zwholemean+Zim;
            ll = ll+1;
        end
    end
end
    
toc
close(h)
Cwholemean = Cwholemean/ll;
Swholemean = Swholemean/ll;
Zwholemean = Zwholemean/ll;
Omatr = sqrt(Cwholemean.^2+Swholemean.^2+Zwholemean.^2);
Omatr(isnan(Omatr)) = 0;
% 'Omatr' is the voxel-wise 3D order index stack, with every voxel
% filled with an order index value showing the fiber organization
% within the localized window region surrounding it
%% Here to do post-processing

% For post-processing, we prepare the 'pretty' images of orientation and
% order index. In these 'pretty' images, the raw intensity image 
% is used to provide the contrast of fiber features, and the orientation or
% order index maps are labeled by different colors to show the orientation 
% or order index information

% first, prepare the pretty theta orientation stack
fiberstackre = fiberstack;
fiberstackre = fiberstackre/max(max(max(fiberstackre)));
prettytheta = zeros(sz,sz2,3,sz3);

uplim = 180;
botlim = 0;
bright = 0.99;
dark = 0.01;

for mm = 1:sz3
    fiberima = fiberstackre(:,:,mm);
    thetaima = aSDE(:,:,mm);
    thetaprettyima = prettymap(thetaima,fiberima,'none',hsv(64),uplim,botlim,bright,dark);
    % 'prettymap.m' is the function used to create 'pretty' maps.
    % 'thetaima' is the theta orientation map, 'fiberima' is the raw fiber
    % image. 'hsv(64)' designates the color scheme. 'uplim' and 'botlim' 
    % are the upper and bottom limits of the orientation metric. The range 
    % of both theta and phi is from 0 to 180. 'bright' and 'dark' are used 
    % to enhance the contrast of the image
    prettytheta(:,:,:,mm) = thetaprettyima;
end

% second, prepare the pretty phi orientation stack
prettyphi = zeros(sz,sz2,3,sz3);
for mm = 1:sz3
    fiberima = fiberstackre(:,:,mm);
    phiima = pSDER(:,:,mm);
    phiprettyima = prettymap(phiima,fiberima,'none',hsv(64),uplim,botlim,bright,dark);
    % Here 'phiima' is the phi orientation map 
    prettyphi(:,:,:,mm) = phiprettyima;
end

% third, prepare the pretty order index stack, with order index acquired with localized window
uplim = 1;
botlim = 0;
prettyorder = zeros(sz,sz2,3,sz3);
for mm = 1:sz3
    fiberima = fiberstackre(:,:,mm);
    orderima = Omatr(:,:,mm);
    orderprettyima = prettymap(orderima,fiberima,'none',jet(64),uplim,botlim,bright,dark);
    % Here 'orderima' is the order index map. For order
    % index, the range is from 0 to 1. Therefore, 'uplim' and 'botlim'
    % are modified accordingly
    prettyorder(:,:,:,mm) = orderprettyima;
end





%% Here to calculate the order index of fiber-like structures in a 2D context
fiberima = imread(['D:\Examples\1.tif']); % modify 1
% Define the path of the folder where the images are saved

fiberima = double(fiberima);

sz = size(fiberima,1);
sz2 = size(fiberima,2);

% Here to create the mask selecting fiber regions

threshint = 3;  % modify 2 
% 'threshint' is the raw background intensity acquired by averaging intensity of several regions identified as background
maska = fiberima > threshint; % maska is a raw binary mask based on the raw background intensity
meanint = mean(fiberstack(maska)); % calculate the mean intensity within the raw binary mask
threshinta = 0.45*meanint; % modify 3 
% 'threshinta' is an inmproved intensity threshold, and the factor '0.45' can be adjusted according to different sample
maskaa = fiberima > threshinta; % 'maskaa' is an improved binary mask
% Here to create the finalmask
finalmaska = maska.*maskaa; % 'finalmaska' is the binary mask acquired based on the above two masks
finalmask = tinyremoval(finalmaska,15);
finalmask = logical(finalmask); % 'finalmask' is the final binary mask selecting the fiber-only regions

armd = 4; % modify 4
[aS] = fiberanglespeed(fiberima,armd,1);

Ctotal = cos(2*aS);
Stotal = sin(2*aS);

Omatr2d = zeros(sz,sz2);
    for i = (armd+1):(sz-armd)
        for j = (armd+1):(sz2-armd)
            localCw = Ctotal((i-armd):(i+armd),(j-armd):(j+armd));
            localSw = Stotal((i-armd):(i+armd),(j-armd):(j+armd));
            maskre = logical(finalmask((i-armd):(i+armd),(j-armd):(j+armd)));
            Cvalue = mean(localCw(maskre));
            Svalue = mean(localSw(maskre));
            Ovalue = sqrt(Cvalue^2+Svalue^2);
            Omatr2d(i,j) = Ovalue;
        end
    end
    
Omatr2d(isnan(Omatr2d)) = 0;

% Here to do the post-processing
uplim = 180;
botlim = 0;
bright = 0.99;
dark = 0.01; % modify 5, including uplim, botlim, bright, and dark
aSDE = aS*180/pi;
thetamappretty = prettymap(aSDE,fiberima,'none',hsv(64),uplim,botlim,bright,dark);

uplim = 1;
botlim = 0; % modify 6, including uplim, botlim
ordermappretty = prettymap(Omatr2d,fiberima,'none',jet(64),uplim,botlim,bright,dark);



















