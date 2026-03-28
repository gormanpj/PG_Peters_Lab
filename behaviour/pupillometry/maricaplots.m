
col = [0.066666666666667	0.266666666666667	0.611764705882353;  0.580392156862745	0.427450980392157	0.203921568627451; 0.890196078431372	0.309803921568627	0.309803921568627];

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

figure;
tiledlayout(2,2,'TileSpacing','compact','Padding','compact')

for g = 1:numel(dayBins)
    nexttile; hold on
    
    for orientation = 1:3
        subset = psthData_allQui( ...
            orientationIdx_allQui == orientation & dayBins{g}.mask, :);
        
        % Drop sparse trials
        trialCoverage = mean(~isnan(subset), 2);
        subset = subset(trialCoverage >= 0.00, :);

        % Derivative per trial
        subsetDeriv = diff(subset, 1, 2);

        % Mean and SEM of derivative traces
        avg = mean(subsetDeriv, 1, 'omitnan');
        nValid = sum(~isnan(subsetDeriv), 1);
        sem = std(subsetDeriv, 0, 1, 'omitnan') ./ sqrt(nValid);

        % Baseline derivative using pre-stim bins
        bline = mean(avg(20:21), 'omitnan');
        avgsub = avg - bline;

        x = 1:numel(avgsub);
        upper = avgsub + sem*1.96;
        lower = avgsub - sem*1.96;

        c = col(orientation,:);

        fill([x fliplr(x)], [upper fliplr(lower)], c, ...
            'FaceAlpha', 0.2, 'EdgeColor', 'none');
        plot(x, avgsub, 'Color', c, 'LineWidth', 4)
    end

    ylim([-0.03, 0.03]);
    xline(20.5, 'k--')
    title(dayBins{g}.title)
    hold off
end

% auc line plot 

aucMeans = [];
aucSems = [];
peakwindow = 30:40;

figure;
tiledlayout(2,2,'TileSpacing','compact','Padding','compact')
for g = 1:numel(dayBins)
    nexttile; hold on
    
    for orientation = 1:3
        subset = psthData_allQui( ...
            orientationIdx_allQui == orientation & dayBins{g}.mask, :);
        mouseIdxSubset = mouseIdx_allQui( ...
            orientationIdx_allQui == orientation & dayBins{g}.mask, :);

        % Drop sparse trials
        trialCoverage = mean(~isnan(subset), 2);
        subset = subset(trialCoverage >= 0.99, :);
        mouseIdxSubset = mouseIdxSubset(trialCoverage >= 0.99, :);
        mouseSubset = unique(mouseIdxSubset);

        % Derivative per trial
        subsetDeriv = diff(subset, 1, 2);

        for auc = 1:numel(mouseSubset)
            mouse = diff(subset(mouseIdxSubset == mouseSubset(auc),:), 1, 2);
            aucs = trapz(mouse(:,peakwindow),2);
            aucMeans(g,orientation) = mean(aucs,1,'omitnan');
            aucSems(g, orientation) = std(aucs, 0, 1, 'omitnan') ./ sqrt(numel(aucs));
        end    
            
       
        plot()
    end

    ylim([-0.03, 0.03]);
    xline(20.5, 'k--')
    title(dayBins{g}.title)
    hold off
end



% AUC quantification from derivative traces, averaged across mice


basewindow = 20:21;   % pre-stimulus derivative bins
peakwindow = 30:50;  % AUC window on derivative trace

mouseIDs = unique(mouseIdx_allQui);

% one value per mouse x day-bin x orientation
aucMouse = nan(numel(mouseIDs), numel(dayBins), 3);

% group summaries across mice
aucMeans = nan(numel(dayBins), 3);
aucSems  = nan(numel(dayBins), 3);

for g = 1:numel(dayBins)
    for orientation = 1:3
        for m = 1:numel(mouseIDs)
            trialMask = orientationIdx_allQui == orientation & ...
                        dayBins{g}.mask & ...
                        mouseIdx_allQui == mouseIDs(m);

            subset = psthData_allQui(trialMask, :);

            % Drop sparse trials
            trialCoverage = mean(~isnan(subset), 2);
            subset = subset(trialCoverage >= 0.99, :);

            % Derivative per trial
            subsetDeriv = diff(subset, 1, 2);

            % Baseline subtract each trial in derivative space
            % (pre-stimulus baseline; cleaner than using a single noisy point)
            trialBase = mean(subsetDeriv(:, basewindow), 2, 'omitnan');
            subsetDeriv = subsetDeriv - trialBase;

            % Trial-wise AUC in the window
            trialAUC = trapz(peakwindow, subsetDeriv(:, peakwindow), 2);

            % Mouse-level mean AUC
            aucMouse(m, g, orientation) = mean(trialAUC, 'omitnan');
        end

        vals = aucMouse(:, g, orientation);
        aucMeans(g, orientation) = mean(vals, 'omitnan');
        nMouse = sum(~isnan(vals));
        aucSems(g, orientation) = std(vals, 0, 'omitnan') ./ sqrt(nMouse);
    end
end

% Plot across day bins 
figure; hold on
x = 1:numel(dayBins);

for orientation = 1:3
    y = aucMeans(:, orientation);
    e = aucSems(:, orientation);
    c = col(orientation,:);

    errorbar(x, y, e, ...
        'Color', c, ...
        'LineWidth', 3, ...
        'CapSize', 0);
end

xlim([1 numel(dayBins)])
ylim([-0.3 0.3])   % adjust if needed
set(gca, 'XTick', x, ...
         'XTickLabel', cellfun(@(s) s.title, dayBins, 'UniformOutput', false))

xlabel('Day bin')
ylabel('Mean derivative AUC')
title('Derivative AUC by day bin and orientation')

hold off









for g = 1:numel(dayBins)
    vals = squeeze(peakMouse(:, g, :));  % [mouse x orientation]

    % remove mice with any NaNs
    valid = all(~isnan(vals), 2);
    vals = vals(valid, :);

    % pairwise comparisons
    [~, p12] = ttest(vals(:,1), vals(:,2));
    [~, p13] = ttest(vals(:,1), vals(:,3));
    [~, p23] = ttest(vals(:,2), vals(:,3));

    fprintf('Day %d: p12=%.3g, p13=%.3g, p23=%.3g\n', g, p12, p13, p23);
end
% 
% hex = {'#E34F4F', '#11449C', '#946D34'};
% 
% col = cell2mat(cellfun(@(h) sscanf(h(2:end),'%2x%2x%2x',[1 3])/255, ...
%     hex, 'UniformOutput', false));

%figure; histogram(cov)



