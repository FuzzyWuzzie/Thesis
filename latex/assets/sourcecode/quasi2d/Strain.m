function [alpha, tau] = Strain( I1, I2, L, deltaAx, maxX, numPoints, h, column, columns, columnRadius, maxAlpha )

    M = maxX / deltaAx;
    alpha = [];
    tau = [];
    
    for m = 1:M
        pr1 = Pr1(m, deltaAx);
        pr2 = Pr2(alpha, deltaAx);
        
        r1 = R1( I1, pr1, L, numPoints );
        alphas = maxAlpha:-0.001:1;
        taus = zeros(size(alphas));
        correlations = zeros(size(alphas));
        
        for i = 1:length(alphas)
            columnCorrelations = [];
            for c = (column-columnRadius):(column+columnRadius)
                % make sure we have a valid column
                if (c < 1) || (c > columns)
                    continue;
                end
                
                % get the part of the image we need
                r2 = R2( I2(:, c), pr2, L / alphas(i), numPoints );
                
                % and correlate it
                columnCorrelations = [columnCorrelations, Correlate(r1, r2)];
            end
            
            % now pick the best column that had the lowest correlation for
            % this alpha
            c = -columnRadius:column+columnRadius;
            [correlations(i), mindex] = min(columnCorrelations);
            taus(i) = c(mindex);
        end

        % only use the estimated alpha if the slope isn't too steep
        [~, mindex] = min(correlations);
        if m > (M/3)
            if abs(alphas(mindex) - alpha(length(alpha))) > 0.02
                alpha = [alpha; alpha(length(alpha))];
                tau = [tau; tau(length(tau))];
            else
                alpha = [alpha; alphas(mindex)];
                tau = [tau; taus(mindex)];
            end
        else
            alpha = [alpha; alphas(mindex)];
            tau = [tau; taus(mindex)];
        end

    end

end
