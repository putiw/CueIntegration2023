function filename = get_info
BASE  = [pwd '/results/'];
keepasking = true;
while keepasking
    prompt = {'Subject ID: ', 'Session: ', 'Run: '};
    numlines = 1;
    title = 'Experiment Parameters';
    defAns = {'0248', '01', '1'};
    info = inputdlg(prompt, title, numlines, defAns);
    
    subID = info{1};
    session = str2double(info{2});
    run = str2double(info{3});
    
    % First, create the data file which will store the results and
    % find out if the file is already present, so you don't overwrite
    % and existing file
    filename = [BASE,'sub-' subID,'_ses-',sprintf('%02d',session),'_task-cue_run-', info{3},'_' datestr(now,30) '.mat'];
    % check response and ask again
    if ~ischar(subID) || length(subID)~=4
        disp('error! please enter 4 digit numerical for ID')
        keepasking = true;
    elseif isnan(session) || length(info{2})<1 || length(info{2})>3
        disp('error! please enter 2 digit numerical for session')
        keepasking = true;
    elseif isnan(run)
        disp('error! run should be a number')
        keepasking = true;
    elseif exist(filename, 'file')== 2
        % ask if experimenter wants to overwrite or not
        overwrite = input('File already exists, Do you want to overwrite (yes or no): ', 's');
        if strcmp(overwrite, 'yes')
            keepasking = false;
        elseif strcmp(overwrite, 'no')
            keepasking = true;
        end
    else
        keepasking = false;
    end
end
