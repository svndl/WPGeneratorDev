function [scrambled, freqScrambled] = psScramble( inputSet )
% Portilla-Simoncelli scrambling tool
% Separate implementation

    Nsc = 4; % Number of scales
    Nor = 4; % Number of orientations
    Na = 11;  % Spatial neighborhood is Na x Na coefficients It must be an odd number!

    Niter = 15;	% Number of iterations of synthesis loop
    

    x = numel(inputSet);
    scrambled = cell(x, 1);
    freqScrambled = cell(x, 1);
    %first, go through the original image set and  
    
    for i = 1:x
        y = numel(inputSet{i});
        % two cell arrays to store freq-scrambled and raw-scrambled
        scrambled{i} = cell(y, 1);
        freqScrambled{i} = cell(y, 1);
        
        for j = 1:y
            image = inputSet{i}{j};
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
            
            params = textureAnalysis(im0, Nsc, Nor, Na);

            imKeep(:, 1) = reshape(mask, [Nsx*Nsy, 1]);
            imKeep(:, 2) = reshape(im0, [Nsx*Nsy, 1]);

            res = textureSynthesis(params, [Nsx Nsy], Niter, [], imKeep);
            scrambled{i}{j} = res;
            freqScrambled{i}{j} = fft2(res);
        end
    end
end
