function []= bspline(inPrefix)

    % load lowResImg
    dataFile= [inPrefix '_data.mat'];
    load(dataFile);

    % load x,y,z, s0rder
    gridFile= [inPrefix '_grid.mat'];
    load(gridFile);

    v1=spm_bsplinc(double(lowResImg),[sOrder sOrder sOrder 0 0 0]);
    highResImg = single(spm_bsplins(v1,x,y,z,[sOrder sOrder sOrder 0 0 0]));

    % save highResImg
    save([inPrefix '_resampled2.mat'], 'highResImg');


end
