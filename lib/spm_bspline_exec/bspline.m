%{
===============================================================================
dMRIharmonization (2018) pipeline is written by-

TASHRIF BILLAH
Brigham and Women's Hospital/Harvard Medical School
tbillah@bwh.harvard.edu, tashrifbillah@gmail.com

===============================================================================
See details at https://github.com/pnlbwh/dMRIharmonization
Submit issues at https://github.com/pnlbwh/dMRIharmonization/issues
View LICENSE at https://github.com/pnlbwh/dMRIharmonization/blob/master/LICENSE
===============================================================================
%}

function []= bspline(inPrefix)

    % load lowResImg
    dataFile= [inPrefix '_data.mat'];
    load(dataFile);

    % load sp_low, sp_high, s0rder, size
    gridFile= [inPrefix '_sp.mat'];
    load(gridFile);

    step = double(sp_high ./ sp_low);
    imgDim = double(imgDim);
    sOrder = double(sOrder);

    [x,y,z]=ndgrid(1:step(1):(imgDim(1)+step(1)+0.01), ...
                   1:step(2):(imgDim(2)+step(2)+0.01), ...
                   1:step(3):(imgDim(3)+step(3)+0.01));

    v1=spm_bsplinc(double(lowResImg),[sOrder sOrder sOrder 0 0 0]);
    highResImg = single(spm_bsplins(v1,x,y,z,[sOrder sOrder sOrder 0 0 0]));

    % save highResImg
    save([inPrefix '_resampled.mat'], 'highResImg');


end
