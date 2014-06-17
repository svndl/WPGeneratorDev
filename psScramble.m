function [scrambled, image] = psScramble( inputImg )
% Portilla-Simoncelli scrambling tool 
% Original image has to be returned, because it might have been altered

    Nsc = 4; % Number of scales
    Nor = 4; % Number of orientations
    Na = 11;  % Spatial neighborhood is Na x Na coefficients It must be an odd number!

    Niter = 15;	% Number of iterations of synthesis loop
    
    %copy input image, because we might have to modify it 
    image = inputImg;
    
    Nx = size(image, 1);
    Ny = size(image, 2);
    
    %for scrambling algorithm, image dimentions have to be  multiple of 2^(Nsc+2)
    req_mult = 2^(Nsc + 2);
    req_Nx = req_mult*floor(Nx/req_mult);
    req_Ny = req_mult*floor(Ny/req_mult);
    % if current wallpaper's size is less than required 
    
    if (req_Nx > Nx || req_Ny > Ny)
        added_X = cat(2, image, image(:, 1:req_mult));
        added_XY = cat(1, added_X, added_X(1:req_mult, :));
        image = added_XY;
    
    end
    Nsx = req_Nx;
    Nsy = req_Ny;
    mask = zeros(Nsx, Nsy);
    imKeep = zeros(Nsx*Nsy, 1);
    im0 = image(1:Nsx, 1:Nsy);
    try 
        params = textureAnalysis(im0, Nsc, Nor, Na);
    
        imKeep(:, 1) = reshape(mask, [Nsx*Nsy, 1]);
        imKeep(:, 2) = reshape(im0, [Nsx*Nsy, 1]);
        scrambled = textureSynthesis(params, [Nsx Nsy], Niter, [], imKeep);
    catch err
        disp('psScramble:Error  ');
        disp(err.message);
        disp(err.stack(1));
        disp(err.stack(2));
        scrambled = 0;
    end
end
