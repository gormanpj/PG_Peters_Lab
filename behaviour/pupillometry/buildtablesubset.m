ids = {'AP033','AP034'};
vid_table = build_vid_table_for_ids('Z:\Data', ids);

function VID_TABLE = build_vid_table_for_ids(rootDir, idsOrCsv, dateStart, dateEnd, recordingRegexp, outputRoot, skipAlreadyProcessed)
% VID_TABLE = build_vid_table_for_ids(rootDir, idsOrCsv, dateStart, dateEnd, recordingRegexp, outputRoot, skipAlreadyProcessed)
%
% rootDir: top-level data folder (e.g. 'Z:\Data')
% idsOrCsv: cell array of ID strings OR path to a CSV/text file listing IDs (one per line)
% dateStart, dateEnd: optional date strings 'yyyy-mm-dd' filter (inclusive). Use [] for no filter.
% recordingRegexp: optional regexp to filter recording folder names (use [] to disable)
% outputRoot: optional path to outputs root (used to skip already processed). If empty, skip check disabled.
% skipAlreadyProcessed: logical (true/false). If true and output exists, the video is omitted.
%
% Returns a table with columns: ID, DATE, RECORDING, video_path

if nargin < 2 || isempty(idsOrCsv)
    error('You must provide a list of IDs or path to an ID CSV file.');
end
if nargin < 3, dateStart = []; end
if nargin < 4, dateEnd = []; end
if nargin < 5, recordingRegexp = []; end
if nargin < 6, outputRoot = ''; end
if nargin < 7, skipAlreadyProcessed = false; end

% Load IDs if CSV path provided
if ischar(idsOrCsv) || isstring(idsOrCsv)
    csvPath = char(idsOrCsv);
    if exist(csvPath,'file')==2
        % simple one-ID-per-line reading
        fid = fopen(csvPath,'r');
        ids = {};
        tline = fgetl(fid);
        while ischar(tline)
            tline = strtrim(tline);
            if ~isempty(tline)
                ids{end+1,1} = tline; %#ok<AGROW>
            end
            tline = fgetl(fid);
        end
        fclose(fid);
    else
        error('IDs file not found: %s', csvPath);
    end
else
    ids = idsOrCsv(:)';
end

% convert date strings to serial if provided
if ~isempty(dateStart), ds = datenum(dateStart,'yyyy-mm-dd'); else ds = []; end
if ~isempty(dateEnd), de = datenum(dateEnd,'yyyy-mm-dd'); else de = []; end

rows = {};
row = 1;
for i = 1:numel(ids)
    ID = ids{i};
    idDir = fullfile(rootDir, ID);
    if ~exist(idDir,'dir')
        warning('ID folder not found: %s (skipping)', idDir);
        continue;
    end
    dateDirs = dir(idDir);
    dateDirs = dateDirs([dateDirs.isdir] & ~ismember({dateDirs.name},{'.','..'}));
    for j = 1:numel(dateDirs)
        dateName = dateDirs(j).name;
        % filter date if needed (expecting YYYY-MM-DD)
        if ~isempty(ds) || ~isempty(de)
            try
                dnum = datenum(dateName,'yyyy-mm-dd');
            catch
                % if not parseable, skip if filtering active
                continue;
            end
            if ~isempty(ds) && dnum < ds, continue; end
            if ~isempty(de) && dnum > de, continue; end
        end
        % list recording subfolders
        recDirs = dir(fullfile(idDir, dateName));
        recDirs = recDirs([recDirs.isdir] & ~ismember({recDirs.name},{'.','..'}));
        for k = 1:numel(recDirs)
            recName = recDirs(k).name;
            % optional: ensure recording name matches regexp
            if ~isempty(recordingRegexp) && isempty(regexp(recName, recordingRegexp, 'once'))
                continue;
            end
            % try known video path patterns inside the recording
            % (extend this list if you have multiple camera names)
            camCandidates = { fullfile(idDir, dateName, recName, 'mousecam', 'mousecam.mj2'), ...
                              fullfile(idDir, dateName, recName, 'mousecam', 'mousecam.mp4'), ...
                              fullfile(idDir, dateName, recName, 'mousecam', 'mousecam.avi') };
            found = false;
            for c = 1:numel(camCandidates)
                vp = camCandidates{c};
                if exist(vp,'file')
                    % skip if output already exists and skipAlreadyProcessed true
                    if skipAlreadyProcessed && ~isempty(outputRoot)
                        outDir = fullfile(outputRoot, ID);
                        baseName = sanitize_for_filename([ID '_' dateName '_' recName '_' get_base_name(vp)]);
                        outH5 = fullfile(outDir, [baseName '.analysis.h5']);
                        if exist(outH5,'file')
                            fprintf('Skipping (already processed): %s\n', vp);
                            found = true; % treat as handled and skip
                            break;
                        end
                    end
                    rows(row,:) = {ID, dateName, recName, vp}; %#ok<AGROW>
                    row = row + 1;
                    found = true;
                    break;
                end
            end
            % if not found, optionally log
            if ~found
                % video not found in expected locations; skip silently
            end
        end
    end
end

if isempty(rows)
    VID_TABLE = table([],[],[],[],'VariableNames',{'ID','DATE','RECORDING','video_path'});
else
    VID_TABLE = cell2table(rows, 'VariableNames', {'ID','DATE','RECORDING','video_path'});
end

end

%% Helper small functions embedded (or put in separate files)
function s = sanitize_for_filename(s)
    s = regexprep(s, '[\\/:"*?<>|]+', '_');
    s = regexprep(s, '\s+', '_');
    s = regexprep(s, '[^\w\-\._]', '_');
    s = regexprep(s, '_+', '_');
    s = regexprep(s, '^_+|_+$', '');
end

function bn = get_base_name(p)
    [~,bn,~] = fileparts(p);
end
