% Example 1: Synthesis of a "text" texture image, using
% Portilla-Simoncelli texture analysis/synthesis code, based on
% alternate projections onto statistical constraints in a complex
% overcomplete wavelet representation.
%
% See Readme.txt, and headers of textureAnalysis.m and
% textureSynthesis.m for more details.
%
% Javier Portilla (javier@decsai.ugr.es).  March, 2001

close all
clear all

% select image directory

ImDir = uigetdir;
cd(ImDir);

listing = dir('*.pgm');

AllImages = {listing.name};

for nIm = 1:numel(AllImages)
    %im0 = pgmRead(AllImages{nIm});	% im0 is a double float matrix!
    im0 = double(imread(AllImages{nIm}));
    
    Nsc = 4; % Number of scales
    Nor = 4; % Number of orientations
    Na = 9;  % Spatial neighborhood is Na x Na coefficients
	 % It must be an odd number!

    params = textureAnalysis(im0, Nsc, Nor, Na);

    Niter = 25;	% Number of iterations of synthesis loop
    Nsx = 512;	% Size of synthetic image is Nsy x Nsx
    Nsy = 512;	% WARNING: Both dimensions must be multiple of 2^(Nsc+2)

    res = textureSynthesis(params, [Nsy Nsx], Niter);
    
    close all
    figure(1)
    showIm(im0, 'auto', 1, 'Original texture');
    figure(2)
    showIm(res, 'auto', 1, 'Synthesized texture');
    
    imwrite(uint8(res), ['4', AllImages{nIm}], 'pgm');
end
