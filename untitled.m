display = 1;
global GL; % GL data structure needed for all OpenGL demos
backGroundColor = [0.5 0.5 0.5].*255; % Gray-scale - calibrate for display so white and black dots have same contrast with background
skipSync = 1; % skip Sync to deal with sync issues (should be for debugging only)
VP = SetupDisplay_NYUAD(skipSync, backGroundColor, display);
if VP.stereoMode == 8
    Datapixx('SetPropixxDlpSequenceProgram',1); % 1 is for RB3D mode, 3 for setting up to 480Hz, 5 for 1440Hz
    Datapixx('RegWr');
    
    Datapixx('SetPropixx3DCrosstalkLR', 0); % minimize the crosstalk
    Datapixx('SetPropixx3DCrosstalkRL', 0); % minimize the crosstalk
end

VP.backGroundColor = backGroundColor;
priorityLevel=MaxPriority(VP.window);
Priority(priorityLevel);
pa = SetupParameters_NYUAD(VP);
pa.response = zeros(pa.numberOfTrials,1);
kb = SetupKeyboard();
pa.trialNumber = 0;
fn = 1; %Frame 1
dontClear = 0; % Don't clear,on flip, we manually clear screen

kb = SetupKeyboard();
VP = MakeTextures(pa,VP);

    [ monitorFlipInterval nrValidSamples stddev ] =Screen('GetFlipInterval', VP.window);

%%

kb.keyIsDown = 0;
count = 1;
aa = [];
curren = GetSecs;
while kb.keyIsDown == 0;
    for view = 0:1
        Screen('SelectStereoDrawbuffer', VP.window, view);
        %     DrawBackground(VP);
        %
        Screen('DrawText', VP.window, [num2str(view)],VP.Rect(3)/2,VP.Rect(4)/2); %pw
        Screen('DrawDots', VP.window, [0 0],pa.fixationDotSize,pa.fixationDotColor,[VP.Rect(3)./2 VP.Rect(4)./2],2);
    end
   

    VP.vbl = Screen('Flip', VP.window, [], dontClear);
    aa(count,1) = VP.vbl -curren;
    
    [kb,~] = CheckKeyboard(kb); % if response with keyboard
    [kb,~] = CheckResponseButton_MRI_CBI(kb); % if response with response button MRI
    count = count+1;
end

aa(end)/numel(aa)
%%
    
    
    
             
            RestrictKeysForKbCheck([]); % Reenable all keys for KbCheck:
ListenChar; % Start listening to GUI keystrokes again
ShowCursor;
clear moglmorpher;
Screen('CloseAll');%sca;
clear moglmorpher;
Priority(0);
if VP.stereoMode == 8
Datapixx('SetPropixxDlpSequenceProgram',0);
Datapixx('RegWrRd');
end
Datapixx('Close');