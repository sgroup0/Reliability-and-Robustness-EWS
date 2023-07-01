function [H, H_error] = DFA_func(TimeSeries, scale_list, m_val)
    
    %==================================================
    % ASSIGN TIMESERIES AS A ROW VECTOR
    
    X = TimeSeries;
    X = reshape(X, 1, length(X));

    
    %==================================================
    % CONVERT NOISE LIKE TIME SERIES TO RANDOM WALK

    RW_X = cumsum(X - mean(X));
    % RW_X = X;
    
    
    %==================================================
    % MONOFRACTAL DETRENDED FLUCTUATION ANALYSIS

    scale = scale_list;
    m = m_val;
    scale = flip(scale);

    for ns = 1:length(scale)
        segments(ns) = floor(length(X) / scale(ns));

        for v = 1:segments(ns)
            Idx_start = ((v-1) * scale(ns)) + 1;
            Idx_stop = v * scale(ns);
            Index{v, ns} = Idx_start:Idx_stop;
            X_Idx = RW_X(Index{v, ns});
            C = polyfit(Index{v, ns}, RW_X(Index{v, ns}), m);
            trendfit{v, ns} = polyval(C, Index{v,ns});
            RMS{ns}(v) = sqrt(mean((X_Idx - trendfit{v, ns}).^2));

            RMS_display{ns}(Index{v, ns}) = RMS{ns}(v) * ones(1, length(Index{v, ns}));
        end
        
        % Repeat for flipped vector
        Y = flip(RW_X);
        for v = 1:segments(ns)
            Idx_start = ((v-1) * scale(ns)) + 1;
            Idx_stop = v * scale(ns);
            Y_Index{v, ns} = Idx_start:Idx_stop;
            Y_Idx = Y(Y_Index{v, ns});
            C_Y = polyfit(Y_Index{v, ns}, Y(Y_Index{v, ns}), m);
            Y_trendfit{v, ns} = polyval(C_Y, Y_Index{v,ns});
            RMS{ns}(v + segments(ns)) = sqrt(mean((Y_Idx - Y_trendfit{v, ns}).^2));

            RMS_display{ns}(Index{v, ns}) = RMS{ns}(v) * ones(1, length(Index{v, ns}));
        end

        F(ns) = sqrt(mean(RMS{ns}.^2));

    end

    %==================================================
    % PLOT DFA LINE

    C = polyfit(log2(scale), log2(F), 1);
    H = C(1);
    RegLine = polyval(C, log2(scale));

%     fig = figure(100);
%     clf(fig);
%     hold on
%     plot(log2(scale), log2(F), 'LineStyle', 'none', 'LineWidth', 1.5, 'Marker', 'o', 'MarkerFaceColor', 'blue', 'MarkerEdgeColor', 'black');
%     plot(log2(scale), RegLine, 'LineWidth', 1.5, 'Color', 'blue');
% 
%     xlim([log2(min(scale)), max(log2(scale))]);
%     ylim([min([min(log2(F)), min(RegLine)]), max([max(log2(F)), max(RegLine)])]);
    
    % Get Error in H
    fitobject = fit(log2(scale)', log2(F)', 'poly1');
    ci = confint(fitobject, 0.90);
    H_avg = (ci(1, 1) + ci(2, 1)) / 2;
    H_error = ci(2, 1) - H_avg;
    
    
end