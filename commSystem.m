% run the system 
clear all, close all
%% Image pre-processing
filename = 'image_example.png';
qbits = 8;
[Ztres,r,c,m,n,minval,maxval]=ImagePreProcess_gray(filename,qbits); % returns array of 8x8 blocks
%% Convert to bit stream

%% Modulate bitstream

%% Pass modulated bitstream through channel

%% Add noise

%% Pass signal through matched filter

%% Pass signal through equalizer

%% Sample and recover bitstream

%% Convert back to image

%% Image post processing