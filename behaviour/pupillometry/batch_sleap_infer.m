% batch_sleap_infer.m
% Batch-run SLEAP-NN inference on a list of videos and export analysis HDF5s.
% - Calls the Python in venv to ensure the correct environment.
% - Converts .slp -> .analysis.h5 using sleap-convert
%
% USAGE:
%    Edit USER CONFIG block and run this script in MATLAB.
%
% NOTE: This script runs jobs serially (safe for single-GPU). 

%% ------------------------ USER CONFIG ------------------------
% Path to the Python executable in venv 
PYTHON_EXE = 'C:\Users\pgorman\sleap-env-py312\Scripts\python.exe';

% Model directory (the folder that contains training_config.json & best.ckpt)
MODEL_DIR = 'C:\Users\pgorman\Documents\SLEAP\Projects\models\251117_190801.centroid.n=625';

% Second model 
MODEL_DIRS = { MODEL_DIR, 'C:\Users\pgorman\Documents\SLEAP\Projects\models\251117_212458.centered_instance.n=625' };

% GPU device to use
GPU_DEVICE = 'cuda:0';

% Output root (script will create per-ID subfolders)
OUTPUT_ROOT = 'C:\Users\pgorman\Documents\SLEAP\Outputs_AV';

% Temporary work/log folders
LOG_DIR = fullfile(OUTPUT_ROOT, 'logs');
if ~exist(LOG_DIR,'dir'), mkdir(LOG_DIR); end

% Batch size and tracking options for inference
BATCH_SIZE = 10; % Default in SLEAP is 4
DO_TRACKING = true;

% Path list mode or build-from-IDs mode?
% MODE = 'paths'    -> specify a cell array of full video paths (VIDEO_LIST)
% MODE = 'build'    -> supply a table/struct of ID, date, time and use PATH_TEMPLATE
MODE = 'build';

% If MODE == 'paths', list video files here:
%VIDEO_LIST = {
    %'Z:\Data\AM022\2024-04-06\Recording_1428\mousecam\mousecam.mj2'
    % 'Z:\Data\AM022\2024-04-07\Recording_1429\mousecam\mousecam.mj2'
%};

% If MODE == 'build', supply a table 'VID_TABLE' with columns ID, dateStr, timeStr
% and set PATH_TEMPLATE. Example template where you can use placeholders:
%   {ID}, {DATE}, {TIME}, {VIDEO_NAME}
% Example:
% PATH_TEMPLATE = 'Z:\Data\{ID}\{DATE}\{RECORDING}\mousecam\mousecam.mj2';

%table( {'AM012'}', {'2023-12-08'}', {'Recording_1403'}', 'VariableNames', {'ID','DATE','RECORDING'} );

VID_TABLE = vid_table;

PATH_TEMPLATE = 'Z:\Data\{ID}\{DATE}\{RECORDING}\mousecam\mousecam.mj2';

% -----------------------------------------------------------------

% Make sure output root exists
if ~exist(OUTPUT_ROOT,'dir'), mkdir(OUTPUT_ROOT); end

% Build list of videos
videos = {};
if strcmpi(MODE,'paths')
    videos = VIDEO_LIST(:);
elseif strcmpi(MODE,'build')
    for r = 1:height(VID_TABLE)
        p = PATH_TEMPLATE;
        p = strrep(p, '{ID}', VID_TABLE.ID{r});
        p = strrep(p, '{DATE}', VID_TABLE.DATE{r});
        p = strrep(p, '{RECORDING}', VID_TABLE.RECORDING{r});
        videos{end+1,1} = p; %#ok<SAGROW>
    end
else
    error('Unknown MODE: %s', MODE);
end

% Basic checks
nJobs = numel(videos);
fprintf('Found %d videos to process.\n', nJobs);

% Test whether sleap-convert exists in PATH (system 'where' -- Windows)
venvScripts = 'C:\Users\pgorman\sleap-env-py312\Scripts';
sleapConvertExe = fullfile(venvScripts,'sleap-convert.exe');

if exist(sleapConvertExe,'file') == 2
    have_sleap_convert = true;
    fprintf('Found sleap-convert at: %s\n', sleapConvertExe);
else
    have_sleap_convert = false;
    fprintf('sleap-convert not found in venv Scripts (%s). \n', venvScripts);
end

% Prepare results log
results = cell(nJobs,6); % columns: status(0/1), video, out_slp, out_h5, elapsed_sec, message
resultsCSV = fullfile(LOG_DIR, ['batch_results_' datestr(now,'yyyymmdd_HHMMSS') '.csv']);
fidlog = fopen(fullfile(LOG_DIR,'batch_run.log'),'a');
fprintf(fidlog, '=== Batch run start %s ===\n', datestr(now));
fclose(fidlog);

% Loop over videos
for i=1:nJobs
    vidPath = videos{i};
    fprintf('\n[%d/%d] Processing: %s\n', i, nJobs, vidPath);
    tStart = tic;
    try
        if ~isfile(vidPath)
            error('Video file does not exist: %s', vidPath);
        end

        % Determine per-ID output folder: try to extract an ID from the path or fallback
        % attempt to extract tokens from path
        [idToken, dateToken, recToken] = extract_tokens_from_path(vidPath); % helper function below
        
        % Get explicit ID from VID_TABLE (always correct)
        if ismember('ID', VID_TABLE.Properties.VariableNames)
            idToken = VID_TABLE.ID{i};
        else
            idToken = get_id_from_path(vidPath);   % use regex-based path parser
        end
        
        % Final safety fallback
        if isempty(idToken)
            idToken = 'unknownID';
        end
    
        % fallback tokens if not found
        if isempty(dateToken), dateToken = 'unknownDate'; end
        if isempty(recToken), recToken = 'unknownRecording'; end
        
        % sanitize tokens for filenames (remove spaces, special chars)
        idSafe  = sanitize_filename(idToken);
        dateSafe = sanitize_filename(dateToken);
        recSafe  = sanitize_filename(recToken);
        
        % basename from file (e.g., mousecam)
        [~, baseName, ~] = fileparts(vidPath);
        baseSafe = sanitize_filename(baseName);
        
        % build a composite filename 
        outBase = sprintf('%s_%s_%s_%s', idSafe, dateSafe, recSafe, baseSafe);
        
        % limit filename length to avoid Windows issues
        maxLen = 200;
        if numel(outBase) > maxLen
            outBase = outBase(1:maxLen);
        end
        
        % per-ID directory
        outDirID = fullfile(OUTPUT_ROOT, idSafe);
        if ~exist(outDirID, 'dir'), mkdir(outDirID); end
        
        % final output file paths
        outSLP = fullfile(outDirID, [outBase, '.pred.slp']);
        outH5  = fullfile(outDirID, [outBase, '.analysis.h5']);

        % Build sleap-nn command using venv python -m sleap_nn.cli track
        % This guarantees the venv python is the interpreter.
        modelArgs = '';
        for mi = 1:numel(MODEL_DIRS)
            modelArgs = [modelArgs ' --model_paths "' MODEL_DIRS{mi} '"']; %#ok<AGROW>
        end

        trackFlag = '';
        if DO_TRACKING
            trackFlag = ' --tracking';
        end

        cmdTrack = sprintf('"%s" -m sleap_nn.cli track --data_path "%s" %s --output_path "%s" --device %s --batch_size %d %s --max_instances 1 --ensure_grayscale --robust_best_instance 1 --tracking_window_size 3 --max_tracks 1 --track_matching_method hungarian --post_connect_single_breaks' , ...
            PYTHON_EXE, vidPath, modelArgs, outSLP, GPU_DEVICE, BATCH_SIZE, trackFlag);

        % Clear environment variables that might force CPU
        % (This affects only the spawned process)
        % Windows: use 'set VAR=' before the command in one line; easier: call via system after clearing in MATLAB environment
        setenv('CUDA_VISIBLE_DEVICES','');
        setenv('SLEAP_DEVICE','');

        fprintf('Running inference command:\n%s\n', cmdTrack);
        [st, cmdOut] = system([cmdTrack ' 2>&1']); % capture stderr and stdout
        elapsed = toc(tStart);

        if st ~= 0
            % Failure
            fprintf('Inference failed for %s\nOutput:\n%s\n', vidPath, cmdOut);
            results{i,1} = 0;
            results{i,6} = sprintf('Inference failed: exit %d', st);
            results{i,5} = elapsed;
            results{i,2} = vidPath;
            results{i,3} = outSLP;
            results{i,4} = outH5;
            continue;
        else
            fprintf('Inference completed. Time: %.1f s\n', elapsed);
        end

        % Convert SLP -> analysis.h5
        convStart = tic;
        if have_sleap_convert
            convCmd = sprintf('"%s" --format analysis --output "%s" "%s"', sleapConvertExe, outH5, outSLP);
            [cst, cout] = system([convCmd ' 2>&1']);
        else
            fprintf('sleap-convert not found in venv. Unable to convert to HDF5 file. \n');
        end
        convElapsed = toc(convStart);

        if cst ~= 0
            fprintf('Conversion failed for %s\nOutput:\n%s\n', outSLP, cout);
            results{i,1} = 0;
            results{i,6} = sprintf('Conversion failed: exit %d', cst);
            results{i,5} = elapsed;
            results{i,2} = vidPath;
            results{i,3} = outSLP;
            results{i,4} = outH5;
            continue;
        else
            fprintf('Conversion succeeded (%.1fs). Output: %s\n', convElapsed, outH5);
            results{i,1} = 1;
            results{i,6} = 'OK';
            results{i,5} = elapsed + convElapsed;
            results{i,2} = vidPath;
            results{i,3} = outSLP;
            results{i,4} = outH5;
        end

    catch ME
        fprintf('EXCEPTION processing %s: %s\n', vidPath, ME.message);
        results{i,1} = 0;
        results{i,6} = ['EXCEPTION: ' ME.message];
        results{i,5} = toc(tStart);
        results{i,2} = vidPath;
        results{i,3} = '';
        results{i,4} = '';
    end

    % write an incremental CSV row so long runs can be inspected
    append_results_csv(resultsCSV, results(1:i,:));
    % small pause to let system settle (optional)
    pause(0.5);
end

% Finalize
fprintf('\nBatch finished. Writing final CSV to %s\n', resultsCSV);
append_results_csv(resultsCSV, results);

%% -------------------- helper functions -----------------------
function id = get_id_from_path(p, idPattern)
% get_id_from_path - robustly infer an ID token from a path
% Usage:
%   id = get_id_from_path(path)
%   id = get_id_from_path(path, idPattern)
%
% idPattern (optional) - a regular expression string to prefer (e.g. '^AM\d+')
% If idPattern is provided, tokens matching that pattern are preferred.
% Otherwise function searches for the first token matching generic ID pattern.
    if nargin < 2
        idPattern = '';
    end

    id = '';

    try
        % split on both slash/backslash
        parts = regexp(p, '[\\/]+', 'split');
        parts = parts(~cellfun(@isempty, parts));

        % prefer user-specified pattern if given
        if ~isempty(idPattern)
            for k = 1:numel(parts)
                if ~isempty(regexp(parts{k}, idPattern, 'once'))
                    id = parts{k};
                    return;
                end
            end
        end

        % Generic ID heuristic: token that is letters followed by digits, e.g. AM023, AP123, S1, M12
        genericIDpat = '^[A-Za-z]{1,4}\d{1,5}$';
        for k = 1:numel(parts)
            tok = parts{k};
            if ~isempty(regexp(tok, genericIDpat, 'once'))
                id = tok;
                return;
            end
        end

        % If nothing matched, consider tokens *after* the drive/root segment,
        % but avoid picking common folder names like 'Data', 'raw', 'videos', 'archive'
        commonNames = {'data','raw','videos','archive','recordings','datasets'};
        for k = 1:numel(parts)
            tok = parts{k};
            % skip drive-like tokens (e.g., 'Z:' or 'C:')
            if length(tok) <= 2 && endsWith(tok, ':')
                continue;
            end
            if ~any(strcmpi(tok, commonNames))
                % choose the first non-common token as fallback
                id = tok;
                return;
            end
        end

    catch
        id = '';
    end
end

function append_results_csv(csvPath, resultsCell)
    % Append or create CSV summarizing resultsCell rows
    % resultsCell: Nx6 cell with columns as described earlier
    headers = {'ok','video','out_slp','out_h5','seconds','message'};
    % convert to table
    n = size(resultsCell,1);
    T = cell2table(resultsCell, 'VariableNames', headers);
    % write
    if exist(csvPath,'file')
        % overwrite with full contents
        writetable(T, csvPath);
    else
        writetable(T, csvPath);
    end
end

function s = sanitize_filename(s)
    % Replace characters that are unsafe in filenames
    if isempty(s), s = 'NA'; return; end
    % replace spaces and slashes
    s = regexprep(s, '[\\/:"*?<>|]+', '_');
    s = regexprep(s, '\s+', '_');            % spaces -> _
    s = regexprep(s, '[^\w\-\._]', '_');     % anything not alnum, _, -, . -> _
    % collapse multiple underscores
    s = regexprep(s, '_+', '_');
    % trim underscores
    s = regexprep(s, '^_+|_+$', '');
end

function [idToken, dateToken, recToken] = extract_tokens_from_path(p, idPattern)
% Robust extraction of ID, DATE, and RECORDING from a path
% usage: [id,date,rec] = extract_tokens_from_path(path) or with preferred idPattern
    if nargin < 2, idPattern = ''; end
    idToken = ''; dateToken = ''; recToken = '';

    try
        parts = regexp(p, '[\\/]+', 'split');
        parts = parts(~cellfun(@isempty, parts));

        % 1) ID - try preferred pattern, then generic pattern
        if ~isempty(idPattern)
            for k=1:numel(parts)
                if ~isempty(regexp(parts{k}, idPattern, 'once'))
                    idToken = parts{k}; break;
                end
            end
        end
        if isempty(idToken)
            genericIDpat = '^[A-Za-z]{1,4}\d{1,5}$';
            for k=1:numel(parts)
                if ~isempty(regexp(parts{k}, genericIDpat, 'once'))
                    idToken = parts{k}; break;
                end
            end
        end

        % 2) DATE - look for YYYY-MM-DD
        datepat = '^\d{4}-\d{2}-\d{2}$';
        for k=1:numel(parts)
            if ~isempty(regexp(parts{k}, datepat, 'once'))
                dateToken = parts{k}; break;
            end
        end

        % 3) RECORDING - look for folder names starting with 'Recording' or 'rec'
        for k=1:numel(parts)
            tok = parts{k};
            if startsWith(lower(tok), 'record', 'IgnoreCase', true) || startsWith(lower(tok), 'rec', 'IgnoreCase', true)
                recToken = tok; break;
            end
        end

        % If rec not found, but date found, maybe next segment after date is recording
        if isempty(recToken) && ~isempty(dateToken)
            idx = find(strcmp(parts, dateToken),1);
            if ~isempty(idx) && idx < numel(parts)
                candidate = parts{idx+1};
                if startsWith(lower(candidate), 'record', 'IgnoreCase', true) || startsWith(lower(candidate), 'rec', 'IgnoreCase', true)
                    recToken = candidate;
                end
            end
        end

        % Final fallback: leave idToken empty rather than return 'Data'
        % (do NOT set idToken = parts{2} as before)

    catch
        idToken=''; dateToken=''; recToken='';
    end
end
