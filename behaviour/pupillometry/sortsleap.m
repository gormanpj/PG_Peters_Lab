% sort_sleap_outputs_by_workflow.m
% Copy SLEAP outputs (.pred.slp and .analysis.h5) into workflow folders
% discovered under the original data tree (bonsai folder).
%
% For each output file named like:
%   AM026_2024-07-27_Recording_1229_mousecam.pred.slp
% the script infers ID=AM026, DATE=2024-07-27, RECORDING=Recording_1229,
% then checks Z:\Data\<ID>\<DATE>\<RECORDING>\bonsai\* for workflow folders.
% It then copies files to:
%   SORTED_ROOT\<workflow>\<ID>\...
%
% Configure the USER CONFIG below before running.

%% -------------------- USER CONFIG --------------------
OUTPUT_ROOT = 'C:\Users\pgorman\Documents\SLEAP\Outputs_AV';   % original outputs root (per-ID folders)
DATA_ROOT   = 'Z:\Data';                                    % original recordings root containing bonsai folders
SORTED_ROOT = 'C:\Users\pgorman\Documents\SLEAP\Outputs_AV_Sorted'; % where to copy sorted files
WORKFLOW_WHITELIST = { 'lcr_passive_size60', 'lcr_passive', 'wheel_aversive_wn_stage1_longer_noincorrect_size60', 'wheel_aversive_wn_stage1_longer_noincorrect', 'wheel_aversive_wn_stage1_noincorrect'}; % optional; empty = accept any workflow found
DRY_RUN = false;   % true = print actions but don't copy
VERBOSE = true;
LOG_CSV = fullfile(SORTED_ROOT, 'sort_summary.csv');
% ------------------------------------------------------

if ~exist(OUTPUT_ROOT,'dir')
    error('OUTPUT_ROOT does not exist: %s', OUTPUT_ROOT);
end
if ~exist(DATA_ROOT,'dir')
    error('DATA_ROOT does not exist: %s', DATA_ROOT);
end
if ~exist(SORTED_ROOT,'dir')
    mkdir(SORTED_ROOT);
end

% find output files recursively (only the analysis.h5 files)
fileList = dir(fullfile(OUTPUT_ROOT, '**', '*.analysis.h5'));

fprintf('Found %d output files to consider.\n', numel(fileList));

% Prepare results container
results = cell(numel(fileList), 6); % src, filename, ID, DATE, RECORDING, workflow(s), dest(s), status
row = 1;

for k = 1:numel(fileList)
    f = fileList(k);
    srcPath = fullfile(f.folder, f.name);
    try
        % Parse filename to get ID/DATE/RECORDING
        [idToken, dateToken, recToken] = parse_output_filename(f.name);
        
        if isempty(idToken) || isempty(dateToken) || isempty(recToken)
            warning('Could not parse tokens from filename: %s', f.name);
            workflowNames = {'unknown'};
            destDirs = { fullfile(SORTED_ROOT,'unknown','unknownID') };
            status = 'parse_failed';
        else
            % locate bonsai folder for that recording
            recordingPath = fullfile(DATA_ROOT, idToken, dateToken, recToken);
            bonsaiDir = fullfile(recordingPath, 'bonsai');
        
            if ~exist(bonsaiDir,'dir')
                workflowNames = {'no_bonsai'};
                destDirs = { fullfile(SORTED_ROOT,'no_bonsai', idToken) };
                status = 'no_bonsai';
            else
                wdirs = dir(bonsaiDir);
                wdirs = wdirs([wdirs.isdir] & ~ismember({wdirs.name},{'.','..'}));
                wnames = {wdirs.name};
        
                if ~isempty(WORKFLOW_WHITELIST)
                    selected = intersect(wnames, WORKFLOW_WHITELIST, 'stable');
                else
                    selected = wnames;
                end
        
                if isempty(selected)
                    workflowNames = {'no_workflow_found'};
                    destDirs = { fullfile(SORTED_ROOT,'no_workflow_found', idToken) };
                    status = 'no_workflow_found';
                else
                    workflowNames = selected;
                    destDirs = cellfun(@(w) fullfile(SORTED_ROOT, w, idToken), ...
                                       workflowNames, 'UniformOutput', false);
                    status = 'ok';
                end
            end
        end

        % Ensure destination dirs exist and perform copy
        destsMade = {};
        for di = 1:numel(destDirs)
            destDir = destDirs{di};
            if ~exist(destDir,'dir') && ~DRY_RUN
                mkdir(destDir);
            end
            destPath = fullfile(destDir, f.name);
            if DRY_RUN
                if VERBOSE, fprintf('[DRY] Would copy\n  %s\n-> %s\n', srcPath, destPath); end
            else
                copyfile(srcPath, destPath);
                if VERBOSE, fprintf('Copied:\n  %s\n-> %s\n', srcPath, destPath); end
            end
            destsMade{di} = destPath; %#ok<SAGROW>
        end

        % Record result
        results{row,1} = srcPath;
        results{row,2} = f.name;
        results{row,3} = idToken;
        results{row,4} = dateToken;
        results{row,5} = recToken;
        results{row,6} = strjoin(workflowNames, ';');
        results{row,7} = strjoin(destsMade, ';');
        results{row,8} = status;
    catch ME
        results{row,1} = srcPath;
        results{row,2} = f.name;
        results{row,3} = '';
        results{row,4} = '';
        results{row,5} = '';
        results{row,6} = '';
        results{row,7} = '';
        results{row,8} = ['ERROR: ' ME.message];
        warning('Error processing %s: %s', srcPath, ME.message);
    end
    row = row + 1;
end

% Write summary CSV (trim empty trailing rows)
results = results(1:row-1, :);
T = cell2table(results, 'VariableNames', {'src','filename','ID','DATE','RECORDING','workflows','dests','status'});
writetable(T, LOG_CSV);
fprintf('Wrote summary to %s\n', LOG_CSV);

%% -------- helper: parse output filename --------
function [idToken, dateToken, recToken] = parse_output_filename(fname)
    % Expected format:
    % ID_DATE_RECORDING_mousecam.analysis.h5

    idToken = ''; dateToken = ''; recToken = '';

    if ~endsWith(fname, '.analysis.h5')
        return;
    end

    base = extractBefore(fname, '.analysis.h5');
    tokens = split(base, '_');

    if numel(tokens) < 4
        return;
    end

    idToken   = tokens{1};
    dateToken = tokens{2};
    recToken  = strjoin(tokens(3:end-1), '_');
end
