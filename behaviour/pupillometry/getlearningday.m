
% See unique animals in bhv
aniUniq = unique(bhv.animal);

% initialize cell to store learning days
aniDates = cell(numel(aniUniq),2);

for n = 1:numel(aniUniq); 
    id = aniUniq{n};

    aniDates(n,1) = {id};
    
    % Get only recordings for this animal
    daysID = bhv(strcmp(bhv.animal, id),:);

    if any(daysID.days_from_learning == 0) % check if there was a learning day
        ld = daysID(daysID.days_from_learning == 0, :);
        aniDates(n,2) = ld.rec_day; % set learning date
    else
        aniDates(n,2) = {NaN}; % if no learning date, NaN
    end
end

