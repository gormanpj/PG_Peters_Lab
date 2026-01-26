%% Load behavior, add 'learned_days' and 'days_from_learning' fields
data_path = fullfile(plab.locations.server_path,'Lab','Papers','Marica_2025','data');
% Set stat and p-value to define learning day
use_stat = 'firstmove_mean';
learn_stat_p = 0.05;
% Load behavior
load(fullfile(data_path,'bhv'));
% Set "learned_days" and "days_from_learning"
bhv.learned_days = cellfun(@(x) x < learn_stat_p,bhv.(['stimwheel_pval_',use_stat]));
for curr_animal = unique(bhv.animal)'
    curr_rows = strcmp(bhv.animal,curr_animal);
    bhv.days_from_learning(curr_rows) = ...
        (1:sum(curr_rows))' - ...
        max([NaN,find(bhv.learned_days(curr_rows),1)]);
end