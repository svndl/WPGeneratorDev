# WPGeneratorDev
Matlab code for wallpaper generation

Instructions for running:
To generate the wallpers, run *generateWPTImages_multi.m* or *generateWPTImages_multi(groupNames,nGroups,tileSize,outDir)*

Default values for the arguments:

* *groupNames = {'P1', 'P2', 'PM' ,'PG', 'CM', 'PMM', 'PMG', 'PGG', 'CMM', 'P4', 'P4M', 'P4G', 'P3', 'P3M1', 'P31M', 'P6', 'P6M'}*;    
* *nGroups = 10* -- number of exemplars per group;    
* *tileSize = 100* -- the area of a smallest repeating tile in pixels;    
* *outDir = '~/Desktop/WPSet'* -- directory where images are saved.     

Default image size *wpSize* is 600 pixels.

# Wallpaper generation pipeline

Function *image = new_SymmetricNoise(type, N, n, optTexture)* generates a single wallpaper exemplar. 
Arguments: 
*type* -- wallpaper type;
*N* -- the size of the output image (for complex groups returned image will have size at least NxN);
*n* -- the size of repeating pattern for all groups;
*optTexture* -- optional base texture.

For most of the group types, generation is done within the switch block. Groups 'P3', 'P3M1', 'P31M', 'P6', 'P6M' are generated by separate scripts.
Repeating tile area is kept constant for all groups, we achieve it by cutting out the templates of a different sizes from fundamental region.  
__________________________ ____________________ _______________________ ____________________ _____________  

Tile sizes by groups      | fundamental region |  repeating tile       |  square ratio      | width ratio  
__________________________ ____________________ _______________________ ____________________ _____________  
P1:                       |       (n, n)       |      (n, n)           |   1                |  1  
__________________________ ____________________ _______________________ ____________________ _____________  

P2:                       |       (n, n)       |      (n, 2n)          |   2                |  0.5  
__________________________ ____________________ _______________________ ____________________ _____________  

PM, PG:                   |       (n, n)       |      (2n, n)          |   2                |  1  
__________________________ ____________________ _______________________ ____________________ _____________  

PMG, PMM, P4, P4M:        |       (n, n)       |      (2n, 2n)         |   4                |  0.5  
__________________________ ____________________ _______________________ ____________________ _____________  

PGG:                      |       (n, 2n)      |      (2n, 2n)         |   4                |  0.5  
__________________________ ____________________ _______________________ ____________________ _____________  

CM:                       |       (n, n)       |      (n, 2n)          |   2                |  0.25  
__________________________ ____________________ _______________________ ____________________ _____________  

P4G:                      |       (n, 2n)      |      (4n, 4n)         |   16               |  0.25  
__________________________ ____________________ _______________________ ____________________ _____________  

CMM:                      |       (n, n)       |      (4n, 4n)         |   16               |  0.25  
__________________________ ____________________ _______________________ ____________________ _____________  

P3:                       |       (n, n)       |      (3n, n sqrt(3))  |   3sqrt(3)         |  1/sqrt(3)  
__________________________ ____________________ _______________________ ____________________ _____________  

P31M: s = round(n*2.632); |    (s, s)          |      (3s, s sqrt(3))  |   3s^2 sqrt(3)/n^2 |  1/(2.632*sqrt(3))  
__________________________ ____________________ _______________________ ____________________ _____________  

P3M1, P6, P6M:            |    (s, s)          |      (s, s sqrt(3))   |   s^2 sqrt(3)/n^2  |  1/(2.632*sqrt(3))   
__________________________ ____________________ _______________________ ____________________ _____________  

# Post-processing

Raw wallpapers generated by *new_SymmetricNoise.m* are averaged within the group, each is lowpass-filtered and masked (circle).

