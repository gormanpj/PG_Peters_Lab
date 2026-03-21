
col = [0.8902 0.301 0.301; 0.66 0.266 0.612; 0.5804 0.427 0.204];

dayBins = {
    struct('mask', learnDayIdx_allQui <= -3,                         'title', '<=-3')
    struct('mask', learnDayIdx_allQui >= -2 & learnDayIdx_allQui <= -1, 'title', '-2:-1')
    struct('mask', learnDayIdx_allQui >= 0 & learnDayIdx_allQui <= 1, 'title', '0:1')
    struct('mask', learnDayIdx_allQui >= 2,                          'title', '>=2')
};

figure;
tiledlayout(2,2,'TileSpacing','compact','Padding','compact')

for g = 1:numel(dayBins)
    nexttile; hold on
    
    for orientation = 1:3
        subset = psthData_allQui( ...
            orientationIdx_allQui == orientation & dayBins{g}.mask, :);

        % Drop sparse trials
        trialCoverage = mean(~isnan(subset), 2);
        subset = subset(trialCoverage >= 0.99, :);

        avg = mean(subset, 1, 'omitnan');
        nValid = sum(~isnan(subset), 1);
        sem = std(subset, 0, 1, 'omitnan') ./ sqrt(nValid);

        bline = avg(21);
        avgsub = avg - bline;

        x = 1:numel(avgsub);
        upper = avgsub + sem;
        lower = avgsub - sem;

        fill([x fliplr(x)], [upper fliplr(lower)], col(orientation,:), ...
            'FaceAlpha', 0.2, 'EdgeColor', 'none');
        plot(x, avgsub, 'Color', col(orientation,:), 'LineWidth', 4)
    end

    ylim([-0.3, 0.3]);
    xline(21, 'k--')
    title(dayBins{g}.title)
    hold off
end






%deriv

col = ['b', 'k', 'r'];
figure; tiledlayout
for day = -2:2
    nexttile;hold on
    
    for orientation = 1:3
        subset = psthData_all(orientationIdx_all == orientation & learnDayIdx_all == day, :);
        
        % Drop sparse trials
        trialCoverage = mean(~isnan(subset), 2);
        subset = subset(trialCoverage >= 0.99, :);
        
        % Derivative per trial
        subsetDeriv = diff(subset, 1, 2);
        
        % Mean and SEM of derivative traces
        avg = mean(subsetDeriv, 1, 'omitnan');
        nValid = sum(~isnan(subsetDeriv), 1);
        sem = std(subsetDeriv, 0, 1, 'omitnan') ./ sqrt(nValid);
        
        bline = avg(21);
        avgsub = avg-bline;

        x = 1:numel(avgsub);
        
        % Shaded SEM band
        fill([x fliplr(x)], ...
             [avgsub + sem, fliplr(avgsub - sem)], ...
             col(orientation), ...
             'FaceAlpha', 0.2, 'EdgeColor', 'none');
        
        % Mean derivative line
        plot(x, avgsub, col(orientation), 'LineWidth', 4)
    end
    ylim([-0.03, 0.03]);

    xline(20.5, 'k--')   % diff shifts by one sample
    title(['Derivative, day ' num2str(day)])
    hold off
end