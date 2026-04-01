
col = [0.066666666666667	0.266666666666667	0.611764705882353;  0.580392156862745	0.427450980392157	0.203921568627451; 0.890196078431372	0.309803921568627	0.309803921568627];

dayBins = {
    struct('mask', learnDayIdx_allQui <= -3,                         'title', '<=-3')
    struct('mask', learnDayIdx_allQui >= -2 & learnDayIdx_allQui <= -1, 'title', '-2:-1')
    struct('mask', learnDayIdx_allQui >= 0 & learnDayIdx_allQui <= 1, 'title', '0:1')
    struct('mask', learnDayIdx_allQui >= 2,                          'title', '>=2')
};

%% diameter plot
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


%% deriv (across all trials)

figure;
tiledlayout(1,4,'TileSpacing','compact','Padding','compact')

for g = 1:numel(dayBins)
    nexttile; hold on
    axis square
    for orientation = 1:3
        subset = psthData_allQui( ...
            orientationIdx_allQui == orientation & dayBins{g}.mask, :);
        
        % Drop sparse trials
        trialCoverage = mean(~isnan(subset), 2);
        subset = subset(trialCoverage >= 0.99, :);

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
        upper = avgsub + sem;
        lower = avgsub - sem;

        c = col(orientation,:);

        fill([x fliplr(x)], [upper fliplr(lower)], c, ...
            'FaceAlpha', 0.2, 'EdgeColor', 'none');
        plot(x, avgsub, 'Color', c, 'LineWidth', 4)
    end

    ylim([-0.03, 0.035]);
    xline(20.5, 'k--')
    title(dayBins{g}.title)
    hold off
end


%% AUC quantification from derivative traces, averaged across mice


basewindow = 20:21;   % pre-stimulus derivative bins
peakwindow = 30:50;  % AUC window on derivative trace

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
            trialBase = mean(subsetDeriv(:, basewindow), 2, 'omitnan');
            subsetDeriv = subsetDeriv - trialBase;

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

    errorbar(x, y, e*1.96, ...
        'Color', c, ...
        'LineWidth', 3, ...
        'CapSize', 0);
end

xlim([1 numel(dayBins)])
ylim([-0.3 0.3])   % adjust if needed
set(gca, 'XTick', x, ...
         'XTickLabel', cellfun(@(s) s.title, dayBins, 'UniformOutput', false))
axis square
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
% 
% hex = {'#E34F4F', '#11449C', '#946D34'};
% 
% col = cell2mat(cellfun(@(h) sscanf(h(2:end),'%2x%2x%2x',[1 3])/255, ...
%     hex, 'UniformOutput', false));

%figure; histogram(cov)



% %% to compare peaks instead of AUC
% 
% basewindow = 20:21;   % pre-stimulus derivative bins
% peakwindow = 30:50;   % peak window on derivative trace
% 
% mouseIDs = unique(mouseIdx_allQui);
% 
% % one value per mouse x day-bin x orientation
% peakMouse = nan(numel(mouseIDs), numel(dayBins), 3);
% 
% % group summaries across mice
% peakMeans = nan(numel(dayBins), 3);
% peakSems  = nan(numel(dayBins), 3);
% 
% for g = 1:numel(dayBins)
%     for orientation = 1:3
%         for m = 1:numel(mouseIDs)
%             trialMask = orientationIdx_allQui == orientation & ...
%                         dayBins{g}.mask & ...
%                         mouseIdx_allQui == mouseIDs(m);
% 
%             subset = psthData_allQui(trialMask, :);
% 
%             % Drop sparse trials
%             trialCoverage = mean(~isnan(subset), 2);
%             subset = subset(trialCoverage >= 0.99, :);
% 
%             % Derivative per trial
%             subsetDeriv = diff(subset, 1, 2);
% 
%             % Baseline subtract each trial in derivative space
%             trialBase = mean(subsetDeriv(:, basewindow), 2, 'omitnan');
%             subsetDeriv = subsetDeriv - trialBase;
% 
%             % Trial-wise peak in the window
%             trialPeak = max(subsetDeriv(:, peakwindow), [], 2);
% 
%             % Mouse-level mean peak
%             peakMouse(m, g, orientation) = mean(trialPeak, 'omitnan');
%         end
% 
%         vals = peakMouse(:, g, orientation);
%         peakMeans(g, orientation) = mean(vals, 'omitnan');
%         nMouse = sum(~isnan(vals));
%         peakSems(g, orientation) = std(vals, 0, 'omitnan') ./ sqrt(nMouse);
%     end
% end
% 
% % Plot across day bins
% figure; hold on
% x = 1:numel(dayBins);
% 
% for orientation = 1:3
%     y = peakMeans(:, orientation);
%     e = peakSems(:, orientation);
%     c = col(orientation,:);
% 
%     errorbar(x, y, e, ...
%         'Color', c, ...
%         'LineWidth', 3, ...
%         'CapSize', 0);
% end
% 
% xlim([1 numel(dayBins)])
% set(gca, 'XTick', x, ...
%          'XTickLabel', cellfun(@(s) s.title, dayBins, 'UniformOutput', false))
% axis square
% xlabel('Day bin')
% ylabel('Mean derivative peak')
% title('Derivative peak by day bin and orientation')
% 
% hold off


%% deriv with avg and sem across mice

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
            trialBase = mean(subsetDeriv(:, 20:21), 2, 'omitnan');
            subsetDeriv = subsetDeriv - trialBase;

            % Mouse-level mean derivative trace
            mouseTrace = mean(subsetDeriv, 1, 'omitnan');
            mouseTraces = [mouseTraces; mouseTrace]; %#ok<AGROW>
        end

        % Mean and SEM across mice
        avg = mean(mouseTraces, 1, 'omitnan');
        nValid = sum(~isnan(mouseTraces), 1);
        sem = std(mouseTraces, 0, 1, 'omitnan') ./ sqrt(nValid);

        % Baseline derivative using pre-stim bins
        bline = mean(avg(20:21), 'omitnan');
        avgsub = avg - bline;

        x = 1:numel(avgsub);
        upper = avgsub + sem;
        lower = avgsub - sem;

        c = col(orientation,:);

        fill([x fliplr(x)], [upper fliplr(lower)], c, ...
            'FaceAlpha', 0.2, 'EdgeColor', 'none');
        plot(x, avgsub, 'Color', c, 'LineWidth', 4)
    end

    ylim([-0.03, 0.035]);
    xline(20.5, 'k--')
    title(dayBins{g}.title)
    hold off
end

%% AUC taken across all trials

% % AUC quantification from derivative traces, averaged across all trials
% 
% basewindow = 20:21;   % pre-stimulus derivative bins
% peakwindow = 30:50;   % AUC window on derivative trace
% 
% aucTrialMeans = nan(numel(dayBins), 3);
% aucTrialSems  = nan(numel(dayBins), 3);
% 
% for g = 1:numel(dayBins)
%     for orientation = 1:3
%         trialMask = orientationIdx_allQui == orientation & dayBins{g}.mask;
% 
%         subset = psthData_allQui(trialMask, :);
% 
%         % Drop sparse trials
%         trialCoverage = mean(~isnan(subset), 2);
%         subset = subset(trialCoverage >= 0.99, :);
% 
%         % Derivative per trial
%         subsetDeriv = diff(subset, 1, 2);
% 
%         % Baseline subtract each trial
%         trialBase = mean(subsetDeriv(:, basewindow), 2, 'omitnan');
%         subsetDeriv = subsetDeriv - trialBase;
% 
%         % Trial-wise AUC in the window
%         trialAUC = trapz(peakwindow, subsetDeriv(:, peakwindow), 2);
% 
%         % Mean and SEM across all trials
%         aucTrialMeans(g, orientation) = mean(trialAUC, 'omitnan');
%         nTrial = sum(~isnan(trialAUC));
%         aucTrialSems(g, orientation) = std(trialAUC, 0, 'omitnan') ./ sqrt(nTrial);
%     end
% end
% 
% % Plot across day bins
% figure; hold on
% x = 1:numel(dayBins);
% 
% for orientation = 1:3
%     y = aucTrialMeans(:, orientation);
%     e = aucTrialSems(:, orientation);
%     c = col(orientation,:);
% 
%     errorbar(x, y, e*1.96, ...
%         'Color', c, ...
%         'LineWidth', 3, ...
%         'CapSize', 0);
% end
% 
% xlim([1 numel(dayBins)])
% ylim([-0.3 0.3])
% set(gca, 'XTick', x, ...
%          'XTickLabel', cellfun(@(s) s.title, dayBins, 'UniformOutput', false))
% axis square
% xlabel('Day bin')
% ylabel('Mean derivative AUC')
% 
% hold off
% 


%% diameter plot across animals

figure;
tiledlayout(1,4,'TileSpacing','compact','Padding','compact')

mouseIDs = unique(mouseIdx_allQui);

for g = 1:numel(dayBins)
    nexttile; hold on
    axis square
    axis padded
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

            % Mouse-level mean diameter trace
            mouseTrace = mean(subset, 1, 'omitnan');

            % Baseline subtract this mouse trace
            mouseTrace = mouseTrace - mouseTrace(16);

            mouseTraces = [mouseTraces; mouseTrace]; 
        end

        % Mean and SEM across mice
        avg = mean(mouseTraces, 1, 'omitnan');
        nValid = sum(~isnan(mouseTraces), 1);
        sem = std(mouseTraces, 0, 1, 'omitnan') ./ sqrt(nValid);

        x = 1:numel(avg);
        upper = avg + sem;
        lower = avg - sem;
        
        c = col(orientation,:);
        % plot individual traces
        % for m = 1:size(mouseTraces,1)
        %     plot(mouseTraces(m,:), ...
        %         'Color', [c 0.4], ...   % transparency
        %         'LineWidth', 1);
        % end

        fill([x fliplr(x)], [upper fliplr(lower)], c, ...
            'FaceAlpha', 0.2, 'EdgeColor', 'none');
        plot(x, avg, 'Color', c, 'LineWidth', 4)
    end

    ylim([-0.3, 0.3]);
    xlim([0, 61])
    xline(16, 'k--')
    title(dayBins{g}.title)
    hold off
end


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
        for m = 1:size(mouseTraces,1)
            plot(mouseTraces(m,:), ...
                'Color', [c 0.4], ...   % transparency
                'LineWidth', 1);
        end

        x = 1:numel(avg);
        upper = avg + sem;
        lower = avg - sem;

        fill([x fliplr(x)], [upper fliplr(lower)], c, ...
            'FaceAlpha', 0.2, 'EdgeColor', 'none');
        plot(x, avg, 'Color', c, 'LineWidth', 4)
    end

    ylim([-0.03, 0.035]);
    xline(20.5, 'k--')
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
    for g = 1:numel(dayBins)
        vals = aucMouse(:, g, orientation);
        xj = g + (rand(size(vals))-0.5)*0.15; % jitter

        scatter(xj, vals, 30, c, ...
            'filled', 'MarkerFaceAlpha', 0.5);
    end

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




%% for completeness: AUC for diameter plot too

%% diameter AUC plot across animals

baseIdx   = 21;      % baseline frame
aucWindow = 30:50;   % window to quantify

mouseIDs = unique(mouseIdx_allQui);

diamAUCmouse = nan(numel(mouseIDs), numel(dayBins), 3);
diamAUCmean  = nan(numel(dayBins), 3);
diamAUCsem   = nan(numel(dayBins), 3);

for g = 1:numel(dayBins)
    for orientation = 1:3
        mouseAUCs = [];

        for m = 1:numel(mouseIDs)
            trialMask = orientationIdx_allQui == orientation & ...
                        dayBins{g}.mask & ...
                        mouseIdx_allQui == mouseIDs(m);

            subset = psthData_allQui(trialMask, :);

            % Drop sparse trials
            trialCoverage = mean(~isnan(subset), 2);
            subset = subset(trialCoverage >= 0.99, :);

            % Mouse-level mean diameter trace
            mouseTrace = mean(subset, 1, 'omitnan');

            % Baseline subtract this mouse trace
            mouseTrace = mouseTrace - mouseTrace(baseIdx);

            % AUC in the chosen window
            mouseAUC = trapz(aucWindow, mouseTrace(aucWindow));

            mouseAUCs = [mouseAUCs; mouseAUC]; %#ok<AGROW>
            diamAUCmouse(m, g, orientation) = mouseAUC;
        end

        % Mean and SEM across mice
        diamAUCmean(g, orientation) = mean(mouseAUCs, 'omitnan');
        nValid = sum(~isnan(mouseAUCs));
        diamAUCsem(g, orientation) = std(mouseAUCs, 0, 'omitnan') ./ sqrt(nValid);
    end
end

% Plot across day bins
figure; hold on
x = 1:numel(dayBins);

for orientation = 1:3
    y = diamAUCmean(:, orientation);
    e = diamAUCsem(:, orientation);
    c = col(orientation,:);

    % % plot individual mice
    % for g = 1:numel(dayBins)
    %     vals = diamAUCmouse(:, g, orientation);
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
set(gca, 'XTick', x, ...
         'XTickLabel', cellfun(@(s) s.title, dayBins, 'UniformOutput', false))
axis square
axis padded
xlabel('Day bin')
ylabel('Diameter AUC')

hold off