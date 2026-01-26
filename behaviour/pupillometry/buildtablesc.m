rootDir = 'Z:\Data';

IDs = {}; DATEs = {}; RECORDINGs = {}; VIDEOs = {};

dIDs = dir(rootDir);
dIDs = dIDs([dIDs.isdir] & ~startsWith({dIDs.name}, '.'));

row = 1;

for i = 1:numel(dIDs)
    ID = dIDs(i).name;
    dateDirs = dir(fullfile(rootDir, ID));

    dateDirs = dateDirs([dateDirs.isdir] & ~startsWith({dateDirs.name}, '.'));

    for j = 1:numel(dateDirs)
        DATE = dateDirs(j).name;

        recDirs = dir(fullfile(rootDir, ID, DATE));
        recDirs = recDirs([recDirs.isdir] & startsWith({recDirs.name}, 'Recording'));

        for k = 1:numel(recDirs)
            RECORDING = recDirs(k).name;

            % look for a single mousecam video
            vidPath = fullfile(rootDir, ID, DATE, RECORDING, 'mousecam', 'mousecam.mj2');
            if exist(vidPath, 'file')
                IDs{row,1} = ID;
                DATEs{row,1} = DATE;
                RECORDINGs{row,1} = RECORDING;
                VIDEOs{row,1} = vidPath;
                row = row + 1;
            end
        end
    end
end

VID_TABLE = table(IDs, DATEs, RECORDINGs, VIDEOs, ...
    'VariableNames', {'ID','DATE','RECORDING','video_path'});