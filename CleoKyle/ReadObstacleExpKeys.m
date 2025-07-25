% Read the file into a cell array of lines
filename = 'ExptKeys.txt';
fileLines = readlines(filename);

% Initialize empty arrays
start_frames = [];
end_frames = [];

% Loop through each line to extract frame numbers
for i = 1:length(fileLines)
    line = strtrim(fileLines(i));
    if startsWith(line, 'trial ')
        % Extract numbers between square brackets
        tokens = regexp(line, '\[(\d+),\s*(\d+)\]', 'tokens');
        if ~isempty(tokens)
            nums = str2double(tokens{1});
            start_frames(end+1) = nums(1);
            end_frames(end+1) = nums(2);
        end
    elseif startsWith(line, 'trial_types')
       data_str = extractBetween(line, '[', ']');
       trial_types = str2num(data_str{1});
    end
end