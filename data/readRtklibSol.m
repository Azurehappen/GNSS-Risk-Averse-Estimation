function [pos_ecef, sow, Q, dtime] = readRtklibSol(filename)

% Open the file
fid = fopen(filename, 'r');

% Check if file opened successfully
if fid == -1
    error('Error opening the file. Please check the filename and path.');
end

% Initialize empty arrays for the data
sow = [];
pos_ecef = [];
Q = [];
dtime = [];

% Read the file line by line
while ~feof(fid)
    line = fgetl(fid);
    
    % Skip header lines and empty lines
    if isempty(line) || line(1) == '%'
        continue;
    end
    
    % Extract data from the line
    tokens = strsplit(line);
    time_str = [tokens{1} ' ' tokens{2}];
    date_time = datetime(time_str, 'InputFormat', 'yyyy/MM/dd HH:mm:ss.SSS');
    date_vector = datevec(date_time);
    [~, ~, sow_i] = date2gnsst(date_vector);
    sow = [sow,sow_i];
    dtime = [dtime,date_time];
    pos_ecef = [pos_ecef, [str2double(tokens(3));str2double(tokens(4));str2double(tokens(5))]];
    Q = [Q,str2double(tokens(6))];
end

% Close the file
fclose(fid);

end