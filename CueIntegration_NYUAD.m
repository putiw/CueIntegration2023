function [VP, pa] = CueIntegration_NYUAD(display);
addpath(genpath('HelperToolbox'));
%change
% subject = input('Enter subject [test]: ','s');
% if isempty(subject)
%     % error('Enter subject name to start experiment');
%     subject =
%     'test';
% end
filename = get_info;
%% Setup parameters and viewing geometry
data = [];
global GL; % GL data structure needed for all OpenGL demos
backGroundColor = [0.5 0.5 0.5].*255; % Gray-scale - calibrate for display so white and black dots have same contrast with background
skipSync = 1; % skip Sync to deal with sync issues (should be for debugging only)
VP = SetupDisplay_NYUAD(skipSync, backGroundColor, display);

if VP.stereoMode == 8 && display ~=2
    Datapixx('SetPropixxDlpSequenceProgram',1); % 1 is for RB3D mode, 3 for setting up to 480Hz, 5 for 1440Hz
    Datapixx('RegWr');
     
    Datapixx('SetPropixx3DCrosstalkLR', 0); % minimize the crosstalk
    Datapixx('SetPropixx3DCrosstalkRL', 0); % minimize the crosstalk
end

VP.backGroundColor = backGroundColor;
priorityLevel=MaxPriority(VP.window);
Priority(priorityLevel);
pa = SetupParameters_NYUAD(VP);
pa.response = zeros(pa.numberOfTrials,2);
kb = SetupKeyboard();
pa.trialNumber = 0;
fn = 1; %Frame 1
dontClear = 0; % Don't clear,on flip, we manually clear screen

kb = SetupKeyboard();
VP = MakeTextures(pa,VP);

%% Generate new dot matrices for quick drawing rather than doing the calculations between frames
Screen('SelectStereoDrawbuffer', VP.window, 0);
Screen('DrawText', VP.window, 'Preparing Experiment...L',VP.Rect(3)/2-130,VP.Rect(4)/2);
Screen('SelectStereoDrawbuffer', VP.window, 1);
Screen('DrawText', VP.window, 'Preparing Experiment...R',VP.Rect(3)/2-130,VP.Rect(4)/2);
VP.vbl = Screen('Flip', VP.window, [], dontClear);
create_stim_NYUAD(VP,pa)
s_size = strrep(num2str(pa.stimulusSizeDeg),'.','_'); % for indexing into the dot matrix structure created above.
load('DotBank.mat')
StateID = 0;
OnGoing = 1;
skip = 0;
% Preload mex files
GetSecs; KbCheck;
kbIdx = GetKeyboardIndices;

count = 1;

%% Experiment Starts
while ~kb.keyCode(kb.escKey) && OnGoing
    
    [kb.keydown, ~, kb.keyCode] = KbCheck(-1);
    
    if kb.keyCode(kb.escKey) % Quit on escape
        % error('Escape key pressed');
        OnGoing = 0;
    end
    
    %% States control the experimental flow (e.g., inter trial interval, stimulus, response periods)
    switch StateID
        case 0
            % Draw blank window until button pressed
            Screen('SelectStereoDrawbuffer', VP.window, 0);
            Screen('DrawTexture', VP.window, VP.bg(VP.curBg));
            Screen('DrawText', VP.window, 'Press a button to begin...L',VP.Rect(3)/2-120,VP.Rect(4)/2);
            
            Screen('SelectStereoDrawbuffer', VP.window, 1);
            Screen('DrawTexture', VP.window, VP.bg(VP.curBg));
            Screen('DrawText', VP.window, 'Press a button to begin...R',VP.Rect(3)/2-120,VP.Rect(4)/2);
            VP.vbl = Screen('Flip', VP.window, [], dontClear);
            
            kb.keyIsDown = 0;
            while kb.keyIsDown == 0
                [kb,~] = CheckKeyboard(kb); % if response with keyboard
                                   %[kb,~] = CheckResponseButton_MRI(kb); % if response with response button MRI
            end
            pause(0.3)
            % Draw blank window until MRI triggers
            Screen('SelectStereoDrawbuffer', VP.window, 0);
            Screen('DrawTexture', VP.window, VP.bg(VP.curBg));
            Screen('DrawText', VP.window, 'Waiting for trigger...',VP.Rect(3)/2-95,VP.Rect(4)/2);
            
            Screen('SelectStereoDrawbuffer', VP.window, 1);
            Screen('DrawTexture', VP.window, VP.bg(VP.curBg));
            Screen('DrawText', VP.window, 'Waiting for trigger...',VP.Rect(3)/2-95,VP.Rect(4)/2);
            
            VP.vbl = Screen('Flip', VP.window, [], dontClear);
            
            %waiting for trigger
            switch display
                case 1 %nyuad
                                    
                    Datapixx('RegWrRd');
                    triggerStart = dec2bin(Datapixx('GetDinValues'));
                    kb.keyIsDown = 0;  
                    while ~kb.keyIsDown
                        [kb,~] = CheckTrigger_MRI(kb,triggerStart); % if response with response button MRI
                        %[kb,~] = CheckKeyboard(kb); % if response with keyboard                       
                    end
                case 2 %puti laptop
                    kb.keyIsDown = 0;
                    pause(0.3)
                    while kb.keyIsDown == 0;
                        [kb,~] = CheckKeyboard(kb); % if response with keyboard
                    end
                case 3 %cbi
                    kb.keyIsDown = 0;
                    pause(0.5)
                    while ~kb.keyIsDown
                        [kb,~] = CheckTrigger_MRI_CBI(kb); % if response with response button MRI
                        [kb,~] = CheckKeyboard(kb); % if response with keyboard
                        fprintf('>>>>>>>>>>> waiting for the trigger from the scanner.... \n')
                    end                    
                    fprintf('>>>>>>>>>>> trigger detected \n')  
                   
            end
            
             %pause(4)
             StateID = 1; % Inter trial interval
            
        case 1 %% inter trial interval
            
            pa.trialNumber = pa.trialNumber + 1;
            
            % Get this trial's parameters
            if pa.trialNumber > pa.numberOfTrials
                pause(pa.pause)
                OnGoing = 0; % End experiment
                break;
            else
                pa.trial = pa.design(pa.trialNumber,:); % Get this trial's parameters
            end
            
            % Get the stimulus location
            if size(pa.allPositions,1)>1
                [x, y] = pol2cart(d2r(pa.allPositions(pa.trial(1),1)), tand(pa.allPositions(pa.trial(1),2))*VP.screenDistance);
            else
                [x, y] = pol2cart(d2r(pa.allPositions(1)), tand(pa.allPositions(2))*VP.screenDistance);
            end
            VP.dstCenter = [x, y];
            
            % Get condition and load a stimulus
            pa.current_condition = char(pa.conditionNames(pa.trial(5)));
            pa.whichCondition = char(pa.current_condition);
            % We have a large stimulus "bank" with many repeats of the same stimulus type: pick a random one
            dotMatrix.blank = dotMatrix.monoR;
            rand_stim = randi([1,size(dotMatrix.comb,6)]);
            pa.current_stimulus = squeeze(dotMatrix.(char(pa.current_condition))(pa.trial(1),pa.trial(4),:,:,:,rand_stim));
            
            % The stimuli are identical for towards and away motions,
            % the only thing that changes is the order in which the
            % frames are presented.
            if pa.directions(pa.trial(3)) == 1
                pa.current_stimulus = flip(pa.current_stimulus,3);
            end
            
            
            StateID = 3; % send to fixation point
            fn = 1;
            
            
        case 3 % Begin drawing stimulus
            
           
            colors = pa.current_stimulus(:,5:7,fn);
            
            for view = 0:1 %VP.stereoViews
                Screen('SelectStereoDrawbuffer', VP.window, view);
                
                %% Draw dots based on condition and viewing parameters
                if view == 0
                    pa.dotPosition = [pa.current_stimulus(:,1,fn), pa.current_stimulus(:,3,fn)].*VP.pixelsPerMm;
                    
                else
                    pa.dotPosition = [pa.current_stimulus(:,2,fn), pa.current_stimulus(:,3,fn)].*VP.pixelsPerMm;
                    
                end
                
                switch pa.whichCondition
                    case{'comb'}
                        Screen('DrawDots',VP.window, pa.dotPosition', pa.current_stimulus(:,4,fn), colors', [VP.Rect(3)/2, VP.Rect(4)/2], 2);
                        
                        if view == 0 && pa.photo_align
                            Screen('DrawDots',VP.window, [VP.Rect(1)+25 VP.Rect(4)-25], 50, [255 255 255],[0 0],2);
                        elseif view == 1 && pa.photo_align
                            Screen('DrawDots',VP.window, [VP.Rect(3)-25 VP.Rect(4)-25], 50, [255 255 255],[0 0],2);
                        end
                    case {'monoL'}
                        if view == 0 % only draw on this trial's eye
                            Screen('DrawDots',VP.window, pa.dotPosition', pa.current_stimulus(:,4,fn), colors', [VP.Rect(3)/2, VP.Rect(4)/2], 2);
                            if pa.photo_align
                                Screen('DrawDots',VP.window, [VP.Rect(1)+25 VP.Rect(4)-25], 50, [255 255 255],[0 0],2);
                            end
                        end
                        
                    case {'monoR'}
                        if view == 1 % only draw on this trial's eye
                            Screen('DrawDots',VP.window, pa.dotPosition', pa.current_stimulus(:,4,fn), colors', [VP.Rect(3)/2, VP.Rect(4)/2], 2);
                            if pa.photo_align
                                Screen('DrawDots',VP.window, [VP.Rect(3)-25 VP.Rect(4)-25], 50, [255 255 255],[0 0],2);
                            end
                        end
                        
                    case {'bino'}
                        %                             pa.dotPosition(:,1:2) = (pa.dotPosition(:,1:2) + [VP.dstCenter(1), VP.dstCenter(2)]).*VP.pixelsPerMm;
                        Screen('DrawDots',VP.window, pa.dotPosition', pa.current_stimulus(:,4,fn), colors', [VP.Rect(3)/2, VP.Rect(4)/2], 2);
                        if view == 0 && pa.photo_align
                            Screen('DrawDots',VP.window, [VP.Rect(1)+25 VP.Rect(4)-25], 50, [255 255 255],[0 0],2);
                        elseif view == 1 && pa.photo_align
                            Screen('DrawDots',VP.window, [VP.Rect(3)-25 VP.Rect(4)-25], 50, [255 255 255],[0 0],2);
                        end
                        
                    case {'blank'}
                        colors = repmat(backGroundColor,size(colors,1),1);
                        if view == 1 % only draw on this trial's eye
                            Screen('DrawDots',VP.window, pa.dotPosition', pa.current_stimulus(:,4,fn), colors', [VP.Rect(3)/2, VP.Rect(4)/2], 2);
                            if pa.photo_align
                                Screen('DrawDots',VP.window, [VP.Rect(3)-25 VP.Rect(4)-25], 50, [255 255 255],[0 0],2);
                            end
                        end
                        
                end                


                Screen('DrawTexture', VP.window, VP.bg(VP.curBg));
                 Screen('DrawDots', VP.window, [0 0],pa.fixationDotSize,pa.fixationDotColor,[VP.Rect(3)./2 VP.Rect(4)./2],2);
                
            end
                        
            VP.vbl = Screen('Flip', VP.window, [], dontClear); % Draw frame
            
            if fn == 1 && pa.trialNumber ==1 && skip == 0
                pa.firstFrame = VP.vbl;
                skip = 1;
            else
%                 
                timeSoFar = GetSecs - pa.firstFrame;
                pa.fn(count,1) = timeSoFar;
                pa.fn(count,3) = round(timeSoFar/(1/VP.frameRate));
                pa.fn(count,2) = pa.trialNumber;
%                 
                fn = pa.timeStamps(floor(timeSoFar/(1/VP.frameRate)),1);
             end
%                         
             pa.fn(count,4) = fn;
             count = count +1;


            if fn> pa.numFlips    % If we have exceeded the frames then the stimulus is complete
                StateID = 4;              
            end
            
            
        case 4  %% Get your response
            
            StateTime = GetSecs - 0.025;
            
            % Wait for response....
            
            Screen('SelectStereoDrawbuffer', VP.window, 0);
            Screen('DrawTexture', VP.window, VP.bg(VP.curBg));
            Screen('DrawDots', VP.window, [0 0],pa.fixationDotSize,[255/5 255/5 255/5],[VP.Rect(3)./2 VP.Rect(4)./2],2);
            
            Screen('SelectStereoDrawbuffer', VP.window, 1);
            Screen('DrawTexture', VP.window, VP.bg(VP.curBg));
            Screen('DrawDots', VP.window, [0 0],pa.fixationDotSize,[255/5 255/5 255/5],[VP.Rect(3)./2 VP.Rect(4)./2],2);
            
            VP.vbl = Screen('Flip', VP.window, [], dontClear);
            
            EndT = 0;
            kb.keyIsDown = 0;
            stop = 0;
            while ~kb.keyIsDown %&& stop ==0
                KbCheck;
               % [kb,stop] = CheckKeyboard(kb); % if response with keyboard
                switch display
                    case 1
                        [kb,stop] = CheckResponseButton_MRI(kb); % if response with response button MRI
                        pa.response(pa.trialNumber,1) = kb.resp;
                        pa.response(pa.trialNumber,2) = kb.secs;
                      
                    case 3                        
                        [kb,stop] = CheckResponseButton_MRI_CBI(kb); % if response with response button MRI
                        pa.response(pa.trialNumber,:) = kb.resp;
                        
                end

                if  EndT >  VP.frameRate*(pa.ITI+pa.trialDuration)*pa.trialNumber 
                    break
                end
                EndT = (GetSecs - pa.firstFrame)./(1/VP.frameRate);

            end

            while 2>1
                if  EndT > VP.frameRate*(pa.ITI+pa.trialDuration)*pa.trialNumber+1
                    break
                end
                EndT = (GetSecs - pa.firstFrame)./(1/VP.frameRate);
            end

            StateID = 1;                      
    end
end
%% Save your data

accuracy = pa.design(:,3)==pa.response(:,1);
for condition = 1:pa.exp_mat(5)
    acc=sum(accuracy(pa.design(:,5)==condition))/sum(pa.design(:,5)==condition)*100;
    disp(['Accuracy - ' char(pa.conditionNames{condition}) ': ' num2str(acc) '%'])
end
%filename = fullfile([pwd, '/results/' subject '-' datestr(now,30), '.mat']);
save(filename,'pa','VP');

%% Clean up
RestrictKeysForKbCheck([]); % Reenable all keys for KbCheck:
ListenChar; % Start listening to GUI keystrokes again
ShowCursor;
clear moglmorpher;
Screen('CloseAll');%sca;
clear moglmorpher;
Priority(0);
if VP.stereoMode == 8 && display ~=2
    Datapixx('SetPropixxDlpSequenceProgram',0);
    Datapixx('RegWrRd');
end
Datapixx('Close');

end