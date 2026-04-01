%% Script for videographic plots 

col = [0.066666666666667	0.266666666666667	0.611764705882353;  0.580392156862745	0.427450980392157	0.203921568627451; 0.890196078431372	0.309803921568627	0.309803921568627];

dayBins = {
    struct('mask', learnDayIdx_allQui <= -3,                         'title', '<=-3')
    struct('mask', learnDayIdx_allQui >= -2 & learnDayIdx_allQui <= -1, 'title', '-2:-1')
    struct('mask', learnDayIdx_allQui >= 0 & learnDayIdx_allQui <= 1, 'title', '0:1')
    struct('mask', learnDayIdx_allQui >= 2,                          'title', '>=2')
};


%% deriv with avg and sem across mice, no baseline subtraction

figure;
tiledlayout(1,4,'TileSpacing','compact','Padding','compact')

mouseIDs = unique(mouseIdx_allQui);

for g = 1:numel(dayBins)
    nexttile; hold on
    axis square

    for orientation = 1:3
        mouseTraces = [];   

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

            % Baseline subtract each trial
            % trialBase = mean(subsetDeriv(:, 20:21), 2, 'omitnan');
            % subsetDeriv = subsetDeriv - trialBase;

            % Mouse-level mean derivative trace
            mouseTrace = mean(subsetDeriv, 1, 'omitnan');
            mouseTraces = [mouseTraces; mouseTrace]; %#ok<AGROW>
        end

        % Mean and SEM across mice
        avg = mean(mouseTraces, 1, 'omitnan');
        nValid = sum(~isnan(mouseTraces), 1);
        sem = std(mouseTraces, 0, 1, 'omitnan') ./ sqrt(nValid);
        
        c = col(orientation,:); 
        % plot individual traces
        % for m = 1:size(mouseTraces,1)
        %     plot(mouseTraces(m,:), ...
        %         'Color', [c 0.4], ...   % transparency
        %         'LineWidth', 1);
        % end

        x = 1:numel(avg);
        upper = avg + sem;
        lower = avg - sem;

        fill([x fliplr(x)], [upper fliplr(lower)], c, ...
            'FaceAlpha', 0.2, 'EdgeColor', 'none');
        plot(x, avg, 'Color', c, 'LineWidth', 4)
    end

    ylim([-0.03, 0.035]);
    xline(15.5, 'k--')
    title(dayBins{g}.title)
    hold off
end


%% AUC quantification from derivative traces, averaged across mice, no baseline subtraction


basewindow = 15:16;   % pre-stimulus derivative bins
peakwindow = 25:45;  % AUC window on derivative trace

mouseIDs = unique(mouseIdx_allQui);

aucMouse = nan(numel(mouseIDs), numel(dayBins), 3);
aucMeans = nan(numel(dayBins), 3);
aucSems  = nan(numel(dayBins), 3);

for g = 1:numel(dayBins)
    for orientation = 1:3
        for m = 1:numel(mouseIDs)
            trialMask = orientationIdx_allQui == orientation & dayBins{g}.mask & mouseIdx_allQui == mouseIDs(m);

            subset = psthData_allQui(trialMask, :);

            % Drop sparse trials
            trialCoverage = mean(~isnan(subset), 2);
            subset = subset(trialCoverage >= 0.99, :);

            % Derivative per trial
            subsetDeriv = diff(subset, 1, 2);

            % Baseline subtract each trial 
            % trialBase = mean(subsetDeriv(:, basewindow), 2, 'omitnan');
            % subsetDeriv = subsetDeriv - trialBase;

            % Trial-wise AUC in the window
            trialAUC = trapz(peakwindow, subsetDeriv(:, peakwindow), 2);

            % Mouse-level mean AUC
            aucMouse(m, g, orientation) = mean(trialAUC, 'omitnan');
            %figure;hold on;for n = 1:size(subsetDeriv,1) plot(subsetDeriv(n,:));end
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

    % plot individual mice
    % for g = 1:numel(dayBins)
    %     vals = aucMouse(:, g, orientation);
    %     xj = g + (rand(size(vals))-0.5)*0.15; % jitter
    % 
    %     scatter(xj, vals, 30, c, ...
    %         'filled', 'MarkerFaceAlpha', 0.5);
    % end

    errorbar(x, y, e, ...
        'Color', c, ...
        'LineWidth', 3, ...
        'CapSize', 0);
end

xlim([1 numel(dayBins)])
ylim([-0.3 0.3])   % adjust if needed
set(gca, 'XTick', x, ...
         'XTickLabel', cellfun(@(s) s.title, dayBins, 'UniformOutput', false))
axis square
axis padded
xlabel('Day bin')
ylabel('Mean derivative AUC')

hold off



%% pairwise comparisons

for g = 1:numel(dayBins)
    vals = squeeze(aucMouse(:, g, :));  % [mouse x orientation]

    % remove mice with any NaNs
    valid = all(~isnan(vals), 2);
    vals = vals(valid, :);

    % pairwise comparisons
    [~, p12] = ttest(vals(:,1), vals(:,2));
    [~, p13] = ttest(vals(:,1), vals(:,3));
    [~, p23] = ttest(vals(:,2), vals(:,3));

    fprintf('Day %d: p12=%.3g, p13=%.3g, p23=%.3g\n', g, p12, p13, p23);
end