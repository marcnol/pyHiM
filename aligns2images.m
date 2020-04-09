tic
%% sets example image files and folder names
rootFolder='/home/marcnol/Documents/Images/Embryo_debug_dataset/rawImages';
file1='scan_001_DAPI_017_ROI_converted_decon_ch01.tif';
%file1='scan_004_RT18_017_ROI_converted_decon_ch00.tif';
file2='scan_004_RT20_017_ROI_converted_decon_ch00.tif';
% rootFolder='/home/marcnol/Documents/Images/Experiment15_embryo001/rawImages';
% file1='scan_001_RT14_002_ROI_converted_decon_ch01.tif';
% file2='scan_001_RT10_002_ROI_converted_decon_ch01.tif';


fullfile1=[rootFolder,'/',file1];
fullfile2=[rootFolder,'/',file2];
%% loads images

tiff_info = imfinfo(fullfile1); % return tiff structure, one element per image
I1 = imread(fullfile1, 1) ; % read in first image
%concatenate each successive tiff to tiff_stack
for ii = 2 : size(tiff_info, 1)
    temp_tiff = imread(fullfile1, ii);
    I1 = cat(3 , I1, temp_tiff);
end


tiff_info = imfinfo(fullfile2); % return tiff structure, one element per image
I2= imread(fullfile2, 1) ; % read in first image
%concatenate each successive tiff to tiff_stack
for ii = 2 : size(tiff_info, 1)
    temp_tiff = imread(fullfile2, ii);
    I2 = cat(3 , I2, temp_tiff);
end

%% 2D projections

I1_2d=sum(I1,3);

I2_2d=sum(I2,3);

%% shows images

figure,imshowpair(I1_2d,I2_2d,'falsecolor')

%% contrast 
higher_threshold=.999
lower_threshold=0.99
I1_2d_scaled=mat2gray(max(I1_2d,[],3));
I1_2d_scaled=imadjust(I1_2d_scaled,stretchlim(I1_2d_scaled, [lower_threshold higher_threshold]),[]);
I1_2d_scaled(I1_2d_scaled<lower_threshold)=0;


I2_2d_scaled=mat2gray(max(I2_2d,[],3));
I2_2d_scaled=imadjust(I2_2d_scaled,stretchlim(I2_2d_scaled, [lower_threshold higher_threshold]),[]);
I2_2d_scaled(I2_2d_scaled<lower_threshold)=0;

figure, imshowpair(I1_2d_scaled,I2_2d_scaled,'falsecolor')
%%
I1_1d=reshape(I1_2d,[2048*2048,1]);
I1_scaled_1d=reshape(I1_2d_scaled,[2048*2048,1]);

I2_1d=reshape(I2_2d,[2048*2048,1]);
I2_scaled_1d=reshape(I2_2d_scaled,[2048*2048,1]);

figure,
subplot(2,2,1)
histogram(I1_1d)
subplot(2,2,3)
histogram(I1_scaled_1d)
subplot(2,2,2)
histogram(I2_1d)
subplot(2,2,4)
histogram(I2_scaled_1d)
%% aligment

% image crosscorrelation
settings.searchWindow=[];
settings.fittingWindow=10;
[Xcorr,~]=merfish_finds2DShiftusingXcorrelation(I1_2d_scaled, I2_2d_scaled, settings.searchWindow, settings.fittingWindow,'convn')


toc