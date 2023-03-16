function pa = SetupParameters_NYUAD(VP)

% 10/25/2022 less dots, no fixation during stimulus, check disparity

rng('shuffle'); % shuffle the random number generator seeds so you don't repeat!
%rng('default');
%% Stimulus Parameters
% These are the locations of the stimuli
thetas = [0]+45;    %,180,270                                                       % Polar angle(s) of stimulus
pa.thetaDirs = thetas;
pa.rDirs = 0;                                                              % Eccentricity of circle
radius_stim = 0;                                                            % Eccentricity of stimulus
pa.stimX_deg = radius_stim*cosd(thetas);                                 
pa.stimY_deg = radius_stim*sind(thetas);
pa.stimulusSizeDeg = 6.6667;                                                  % Radius
pa.apertureLipConst = 1.2;                                                   % Can shrink the stimulus to have small border: 1 means there is no border
pa.screenAperture = pa.apertureLipConst*pa.stimulusSizeDeg;                % aperture after considering border
pa.numberOfDots = 80;     %22                                                 % number of dots

conditionNow = 1;

switch conditionNow
    
    case 1
pa.numberOfRepeats = 5;                                                % number of blocks to complete
pa.trialDuration = 1;                                                      % duration of stimulus
pa.ITI = 8;                                                            % duration between stimuli
pa.numberOfBlanks = 0; %
pa.conditionNames   = {'monoL','monoR','bino','comb','blank'};          % Stimulus conditions
pa.pause = 0; % seconds of blank screen at the end

    case 2 %change pause time in cueintegration
pa.numberOfRepeats = 30;                                             % number of blocks to complete
pa.trialDuration = 1.5;                                                      % duration of stimulus
pa.ITI = 3;                                                                % duration between stimuli
pa.numberOfBlanks = 30; %
pa.conditionNames   = {'comb','blank'};          % Stimulus conditions
pa.pause = 15;

     case 3 %change pause time in cueintegration
pa.numberOfRepeats = 40;                                             % number of blocks to complete
pa.trialDuration = 1.5;                                                      % duration of stimulus
pa.ITI = 3;                                                                % duration between stimuli
pa.numberOfBlanks = 10; %
pa.conditionNames   = {'comb','blank'};          % Stimulus conditions
pa.pause = 15;
 
    case 4 %change pause time in cueintegration
pa.numberOfRepeats = 10;                                             % number of blocks to complete
pa.trialDuration = 10.5;                                                      % duration of stimulus
pa.ITI = 10.5;                                                                % duration between stimuli
pa.numberOfBlanks = 0; %
pa.conditionNames   = {'comb','blank'};          % Stimulus conditions
pa.pause = 0;
 
    case 5 %change pause time in cueintegration
pa.numberOfRepeats = 22;                                             % number of blocks to complete
pa.trialDuration = 3;                                                      % duration of stimulus
pa.ITI = 6;                                                                % duration between stimuli
pa.numberOfBlanks = 0; %
pa.conditionNames   = {'comb','blank'};          % Stimulus conditions
pa.pause = 24;

end
      
pa.fixationAcqDura = 0;                                                    % duration of fixation prior to stimulus
pa.disparityLimit = 0.3;  %1                                               % using the same front and back disparity, what is the limit?
pa.loops = 2;   % 1                                                           % # of times dots travel the volume (determines speed)
pa.reversePhi = 1;                                                         % dots change color on wrapping to reduce apparent motion
pa.directions = [-1 1];                                                    % experiment directions (+1:towards, -1:away)
pa.coherence = 1;                                                          % Motion coherence levels
pa.photo_align = 0;
if VP.stereoMode == 1
    pa.numFlips = floor(pa.trialDuration*VP.frameRate/2);                  % every other frame for each eye when in interleaved stereo mode
else
    pa.numFlips = floor(pa.trialDuration*VP.frameRate);                    % each frame to both eyes
end


[pTH,pR] = cart2pol(pa.stimX_deg, pa.stimY_deg);
pTH = -rad2deg(pTH);
if pTH < 0
    pTH = 360 + pTH;
end
pa.allPositions = [pTH; pR]';

%% Dot Parameterss
pa.fixationRadius   = 0.6; %0.7;    (in ???)
pa.fixationDotSize  = 7; %4;                    % fixation dot size (in pixels)
pa.saccadeDotSize = 12; %in pixels
pa.fixationDotColor  = [0 0 0]; % red
pa.dotDiameterinDeg = pa.fixationDotSize/VP.pixelsPerDegree;
pa.dotDiameter = pa.dotDiameterinDeg * VP.pixelsPerDegree;
pa.dotColor = [255, 255, 255, 255]; % white
pa.dotSpacing = (pa.dotDiameterinDeg*1.5)/VP.degreesPerMm;  % in mm since dots are in mm

%% Initial dot motion calculations
% What is the rate at which the disparity must change?
pa.deltaDisp = pa.disparityLimit*2*pa.loops/(VP.frameRate/2);  %pa.disparityLimit*2*pa.loops/(pa.trialDuration*VP.frameRate/2);
% Maximum depth extent (to calc the rear radius and max dot size)
xLocFixationPoint = 0;
ConvergenceAngleFixation = atand(((VP.IOD/2)+xLocFixationPoint)/VP.screenDistance) + atand(((VP.IOD/2)-xLocFixationPoint)/VP.screenDistance); %left eye + right eye
f = @(xDepth) abs(((atand((VP.IOD/2)/(xDepth)) + atand((VP.IOD/2)/(xDepth))) - ConvergenceAngleFixation)-pa.disparityLimit);
[X,FVAL] = fminsearch(@(xDepth) f(xDepth), VP.screenDistance);%, ops); % Finds a minimum of the anonymous function evaluated at xDepth and near the screen
pa.RetrievedDepth(1) = X; % Relative to the screen, which is Z = 0;
f = @(xDepth) abs(((atand((VP.IOD/2)/(xDepth)) + atand((VP.IOD/2)/(xDepth))) - ConvergenceAngleFixation)-(-pa.disparityLimit));
[X,FVAL] = fminsearch(@(xDepth) f(xDepth), VP.screenDistance);%, ops); % Finds a minimum of the anonymous function evaluated at xDepth and near the screen
pa.RetrievedDepth(2) = X; % Relative to the sscreen, which is Z = 0;
pa.depthExtent = abs(pa.RetrievedDepth(2) - pa.RetrievedDepth(1));
pa.dz = pa.depthExtent*pa.loops/(pa.trialDuration*VP.frameRate); %.*pa.modeConstant); % (in mm)

%% Design Structure
pa.exp_mat = [size(pa.allPositions,1) 1 length(pa.directions) length(pa.coherence) length(pa.conditionNames)]; % last index is for stimulation. 1 = no. 2 = yes
pa.repeat_design = fullfact(pa.exp_mat); 
% Remove the extra 0's if we have 0% coherence
if ismember(0,pa.coherence)
    extraZeroEntries = find(~(pa.repeat_design(:,4) == 1 & pa.repeat_design(:,3) == 2));
    pa.repeat_design = pa.repeat_design(extraZeroEntries,:);
end

% Repeat the trial structure with random permutations in each block
pa.design = [];
for r = 1:pa.numberOfRepeats
    pa.temp_design = pa.repeat_design(randperm(size(pa.repeat_design,1)),:);
    pa.design = [pa.design; pa.temp_design];
end


 whichBlank = find(pa.design(:,5)==max(pa.design(:,5)));
 pa.design(whichBlank(randperm(numel(whichBlank),numel(whichBlank)-pa.numberOfBlanks)),:) = [];


pa.repeats_completed = 0;
pa.numberOfTrials = size(pa.design,1);

pa.totalFrames = pa.numberOfTrials * (pa.ITI+pa.trialDuration)*VP.frameRate; %pw
pa.timeStamps = [repmat(1:(pa.ITI+pa.trialDuration)*VP.frameRate,1, ...
    pa.numberOfTrials); repelem(1:pa.numberOfTrials,...
    1,(pa.ITI+pa.trialDuration)*VP.frameRate); zeros(1,pa.totalFrames)]';
pa.fn = [];


%condition wise design matrix
pa.tr = 1;
TRcon = 1:((pa.ITI+pa.trialDuration)/pa.tr):(pa.numberOfTrials * (pa.ITI+pa.trialDuration))./pa.tr;

    pa.dsCon = zeros((pa.numberOfTrials * (pa.ITI+pa.trialDuration)+pa.pause)./1.5,2*(numel(pa.conditionNames)-1));
    for iCue = 1:numel(pa.conditionNames)-1
        pa.dsCon(TRcon,iCue*2-1) = pa.design(:,3)==1&pa.design(:,5)==iCue; %away
        pa.dsCon(TRcon,iCue*2) = pa.design(:,3)==2&pa.design(:,5)==iCue; %toward
    end

 
%trial wise design matrix
    onset = sum(pa.dsCon,2);
    pa.dsTrial = zeros(size(onset,1),numel(find(onset)));     
    ind = sub2ind(size(pa.dsTrial),find(onset),(1:size(pa.dsTrial,2))');
    pa.dsTrial(ind) = 1;



%% Savefile parameters
pa.subjectName = 'Jim';
pa.movie = 0;
pa.screenDump = 0; % Get Image?
pa.baseDir = pwd;



