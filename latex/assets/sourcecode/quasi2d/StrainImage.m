function [alpha, tau] = StrainImage( I1, I2, L, deltaAx, maxX, numPoints, columnRadius, maxAlpha )

    columns = size(I1, 2);

    alpha = zeros(maxX / deltaAx, columns);
    tau = zeros(maxX / deltaAx, columns);
    h = 0;
    
    tStart = tic;

    parfor_progress(columns);
    parfor i = 1:columns
        [alpha(:, i), tau(:, i)] = Strain(I1(:, i), I2, L, deltaAx, maxX, numPoints, h, i, columns, columnRadius, maxAlpha);
        parfor_progress;
    end
    parfor_progress(0);
    
    processingTime = toc(tStart);
    fprintf('Elastogram took %s!\n', datestr(processingTime / 3600 / 24, 'HH:MM:SS'));

end
