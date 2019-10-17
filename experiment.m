function experiment

%-%-%-%-%-
%- INIT %-
%-%-%-%-%-

% ask for condition information (stimulus type, N, etc)
settings = getExperimentSettings();

% set some condition-independent variables
settings.screen_width    = 40;   % in cm (Dell@T115A: ~48cm; Dell@T101C: ~40 cm)
settings.barwidth        = .3;   % width of stimulus bar (deg)
settings.barheight       = .8;   % height of stimulus bar (deg)
settings.ellipseArea     = settings.barwidth*settings.barheight; % ellipse size (deg^2)
settings.jitter          = .6;   % amount of x/y-jitter (deg)
settings.fgdac           = 200;  % foreground grayvalue (RGB)
settings.stimecc         = 7;    % stimulus eccentricity (deg)
settings.ITT             = 500;  % inter stimulus time (ms)
settings.targetort       = -45;  % target orientation
settings.cueCont         = 50;   % contrast of bar target cue
settings.cueEcc          = .9;   % eccentricity of ellipse target cue
settings.stimtime        = 60;   % stimulus presentation time (ms)
settings.barContrasts    = [20 40]; % low and high contrast
settings.cueCont         = 15;   % contrast of bar target cue
settings.ellipseEcc      = [.7 .9]; % low and high eccentricity

% Set background luminance
settings.bgdac = 128;

% Set possible distractor orientations
if strcmp(settings.disttype,'homo')
    settings.distorts = -30;
else
    settings.distorts = -90:90;
end

% screen info
screenNumber=max(Screen('Screens'));       % use external screen if exists
[w h]=Screen('WindowSize', screenNumber);  % screen resolution
screen_resolution = [w h];                 % screen resolution
screen_center = screen_resolution/2;       % screen center
screen_distance = 60;                      % distance between observer and screen (in cm)
screen_angle = 2*(180/pi)*(atan((settings.screen_width/2) / screen_distance)) ; % total visual angle of screen
screen_ppd = screen_resolution(1) / screen_angle;  % pixels per degree
screen_fixposxy = screen_resolution .* [.5 .5]; % fixation position

settings.ellipseArea = settings.ellipseArea * screen_ppd^2;

% open screen
gray=GrayIndex(screenNumber);
windowPtr = screen('OpenWindow',screenNumber,gray,[],32,2);

% create cue texture
if strcmp(settings.stimtype,'bar')
    w = round(settings.barwidth * screen_ppd);
    h = round(settings.barheight * screen_ppd);
    im = ones(w,h) * settings.bgdac * (1 + settings.cueCont/100);
    cuetex = screen('MakeTexture',windowPtr,im);
    cuerot = settings.targetort;
else
    b = sqrt(settings.ellipseArea * sqrt(1 - settings.cueEcc^2) / pi);
    a = settings.ellipseArea / (pi*b);
    im = drawEllipse(2*b,2*a,settings.targetort,settings.fgdac,settings.bgdac);
    cuetex = screen('MakeTexture',windowPtr,im);
    cuerot = 0;
end
cuesrcrect = [0 0 size(im)];

% determine breakpoints (i.e., trialnrs after which to show progress + target reminder)
breakpoints = round([1 2 3] .* (settings.nTrials/4));

% show start screen
screen('TextSize',windowPtr,20);
textx = 100;
texty = screen_center(2) - 50;
dy = 30;
screen('DrawText',windowPtr,'Your task is to detect whether the target (45 deg)',textx,texty,[255 255 255]); texty = texty + dy;
screen('DrawText',windowPtr,'is present among distractors',textx,texty,[255 255 255]); texty = texty + 2*dy;
screen('DrawText',windowPtr,'Next, enter your confidence on a scale of 1-3 (1=not confident,',textx,texty,[255 255 255]); texty = texty + dy;
screen('DrawText',windowPtr,'3=very confident)',textx,texty,[255 255 255]); texty = texty + 2*dy;
screen('DrawText',windowPtr,'Spread your responses across all ratings',textx,texty,[255 255 255]); texty = texty + 2*dy;
[newx newy] = screen('DrawText',windowPtr,'The ORIENTATION of the target looks like this:',textx,texty,[255 255 255]); texty = texty + 3*dy;
destrect = centerRectOnPoint(cuesrcrect,newx+20,newy+20);
screen('drawtexture',windowPtr,cuetex,cuesrcrect,destrect,cuerot);

screen('Flip', windowPtr);
waitForKey;
screen('FillRect', windowPtr, gray);
screen('Flip', windowPtr);
screen('TextSize', windowPtr, 15);
HideCursor;

yesKey = kbName('y');
noKey = kbName('u');
escKey = 27;

datafilename = [settings.subjid '_' settings.stimtype '_' settings.disttype '_' num2str(settings.nStim) '_' datestr(now,'yyyymmddTHHMMSS') '.mat'];

%-%-%-%-%-%-%-%-%-%-%-%-%
%- LOOP THROUGH TRIALS %-
%-%-%-%-%-%-%-%-%-%-%-%-%
fixsize = 4;
fixcol = 0;
nextFlipTime = 0; % just to initialize...
aborted = 0;
trialnr = 0;

while (trialnr < settings.nTrials) && ~aborted
    trialnr = trialnr + 1;
    
    target_present = rand > .5;
    if (target_present)
        targetidx = ceil(rand*settings.nStim);
    else
        targetidx = -1;
    end
    
    % create stimulus patches
    clear stimtex;
    posTheta = (225/360)*2*pi;
    for ii=1:settings.nStim
        [x y] = pol2cart(posTheta,screen_ppd * settings.stimecc);
        stimPos(ii,:) = [x y] + screen_center;
        posTheta = posTheta + (2*pi)/settings.nStim;
        
        if strcmp(settings.stimtype,'bar')
            w = round(settings.barwidth * screen_ppd);
            h = round(settings.barheight * screen_ppd);
            barConts(ii) = settings.barContrasts(ceil(rand*length(settings.barContrasts)));
            im = ones(w,h) * settings.bgdac * (1 + barConts(ii)/100);
            if (ii==targetidx)
                rotation(ii) = settings.targetort;
                orientations(ii) = settings.targetort;
            else
                rotation(ii) = settings.distorts(ceil(rand*length(settings.distorts)));
                orientations(ii) = rotation(ii);
            end
        else
            rotation(ii) = 0;
            ellipseEccs(ii) = settings.ellipseEcc(ceil(rand*length(settings.ellipseEcc)));
            b = sqrt(settings.ellipseArea * sqrt(1 - ellipseEccs(ii)^2) / pi);
            a = settings.ellipseArea / (pi*b);
            if (ii==targetidx)
                im = drawEllipse(2*b,2*a,settings.targetort,settings.fgdac,settings.bgdac);
                orientations(ii) = settings.targetort;
            else
                randort = settings.distorts(ceil(rand*length(settings.distorts)));
                im = drawEllipse(2*b,2*a,randort,settings.fgdac,settings.bgdac);
                orientations(ii) = randort;
            end
        end
        stimtex(ii) = screen('MakeTexture',windowPtr,im);
        patchsize(ii,:) = size(im);
        patchsize(ii,:) = patchsize(ii,[2 1]);
    end
    
    % Jitter stimulus positions
    stimPos = stimPos + round(rand(settings.nStim,2)*settings.jitter*screen_ppd/2 - settings.jitter*screen_ppd/2);
    
    % SCREEN 1a: FIXATION
    screen('fillRect',windowPtr,settings.bgdac);
    drawfixation(windowPtr,screen_fixposxy(1),screen_fixposxy(2),fixcol,fixsize);
    screen('flip',windowPtr,nextFlipTime);
    nextFlipTime = getSecs + .5;
        
    % SCREEN 2: STIMULUS
    screen('fillRect',windowPtr,settings.bgdac);
    drawfixation(windowPtr,screen_fixposxy(1),screen_fixposxy(2),fixcol,fixsize);
    for ii=1:settings.nStim
        srcrect = [0 0 patchsize(ii,:)];
        destrect = centerRectOnPoint(srcrect,stimPos(ii,1),stimPos(ii,2));
        screen('drawtexture',windowPtr,stimtex(ii),srcrect,destrect,-rotation(ii));
    end
    screen('flip',windowPtr,nextFlipTime);
    stimStartTime = getSecs;
    nextFlipTime = getSecs + settings.stimtime/1000;

    screen('fillRect',windowPtr,settings.bgdac);
    drawfixation(windowPtr,screen_fixposxy(1),screen_fixposxy(2),fixcol,fixsize);
    screen('flip',windowPtr,nextFlipTime);   
    nextFlipTime = getSecs + 0.3;
    
    % SCREEN 3a: YES/NO RESPONSE
    screen('fillRect',windowPtr,settings.bgdac);
    drawfixation(windowPtr,screen_fixposxy(1),screen_fixposxy(2),fixcol,fixsize);
    screen('DrawText',windowPtr,'Was the target present? ("y"=yes, "u"=no)',50,screen_center(2)+50,settings.bgdac-50);
    screen('flip',windowPtr,nextFlipTime);
    REALSTIMTIME = round(1000*(getSecs - stimStartTime));
    
    responseStartTime = getSecs;
    done=0;
    while ~done
        keyCode = waitForKey;
        if (keyCode==yesKey)
            YESNORESP = 1;
            done=1;
        elseif (keyCode==noKey)
            YESNORESP = 0;
            done=1;
        elseif (keyCode == escKey)
            aborted=1;
            break;
        end
        
    end
    
    RT = round(1000*(getSecs - responseStartTime));
    
    % SCREEN 3b: CONFIDENCE RATING
    screen('fillRect',windowPtr,settings.bgdac);
    drawfixation(windowPtr,screen_fixposxy(1),screen_fixposxy(2),fixcol,fixsize);
    screen('DrawText',windowPtr,'Rate your confidence on a scale of 1-3 (1=not confident, 3=very confident)',50,screen_center(2)+50,settings.bgdac-50);
    screen('flip',windowPtr,nextFlipTime);
    
    done=0;
    while ~done
        keyCode = waitForKey;
        if (keyCode==49 || keyCode==35)
            done=1;
            CONFIDENCE = 1;
        elseif (keyCode==50 || keyCode==40)
            done=1;
            CONFIDENCE = 2;
        elseif (keyCode==51 || keyCode==34)
            done=1;
            CONFIDENCE = 3;
        elseif (keyCode == escKey)
            aborted=1;
            break;
        end
    end
    
    if strcmp(settings.stimtype,'bar')
        datamatrix(trialnr,:) = [target_present YESNORESP CONFIDENCE REALSTIMTIME RT orientations barConts];
    else
        datamatrix(trialnr,:) = [target_present YESNORESP CONFIDENCE REALSTIMTIME RT orientations ellipseEccs];
    end
    save(datafilename,'settings','datamatrix');
    
    % show progress info + target reminder
    if ~isempty(intersect(breakpoints,trialnr))
        screen('fillRect',windowPtr,settings.bgdac);
        screen('DrawText',windowPtr,['Your have finished ' num2str(round(100*trialnr/settings.nTrials)) '% of the trials'],100,screen_center(2)-80,[255 255 255]);
        screen('DrawText',windowPtr,'Just as a reminder, the ORIENTATION of the target is like this:',100,screen_center(2)-50,[255 255 255]);
        destrect = centerRectOnPoint(cuesrcrect,screen_center(1),screen_center(2));
        screen('drawtexture',windowPtr,cuetex,cuesrcrect,destrect,cuerot);
        screen('Flip', windowPtr);
        waitForKey;
    end
    
    % SCREEN 4: INTER TRIAL DISPLAY
    screen('fillRect',windowPtr,settings.bgdac);
    drawfixation(windowPtr,screen_fixposxy(1),screen_fixposxy(2),fixcol,fixsize);
    screen('flip',windowPtr);
    nextFlipTime = getSecs + (settings.ITT/1000);
    
end

% show start screen
screen('fillRect',windowPtr,settings.bgdac);
screen('DrawText',windowPtr,'End of this block. Press <ENTER> to continue',250,screen_center(2) - 50,[255 255 255]);
screen('Flip', windowPtr);
key = 0;
while (key ~= 13)
    key = waitForKey;
end

%-%-%-%-%-%-%-
%- FINALIZE %-
%-%-%-%-%-%-%-
if (aborted)
    delete(datafilename);
end
ShowCursor;
screen('closeall');



%-%-%-%-%-%-%-%-%-%-%-%-%- HELPER FUNCTIONS %-%-%-%-%-%-%-%-%-%-%-%-%-%-%-

function keyCode = waitForKey
keyCode = ones(1,256);
while sum(keyCode(1:254))>0
    [keyIsDown,secs,keyCode] = KbCheck;
end
while sum(keyCode(1:254))==0
    [keyIsDown,secs,keyCode] = KbCheck;
end
keyCode = min(find(keyCode==1));

function drawfixation(windowPtr,x,y,color,size)
Screen('DrawLine',windowPtr,color,x-size,y,x+size,y,2);
Screen('DrawLine',windowPtr,color,x,y-size,x,y+size,2);

function settings = getExperimentSettings()
stimtype='';
disttype='';
clc;
subjid   = input('Subject initials...........................: ','s');
while ~(strcmp(stimtype,'bar') || strcmp(stimtype,'ellipse'))
    stimtype = input('Stimulus type (bar/ellipse)................: ','s');
end
while ~(strcmp(disttype,'homo') || strcmp(disttype,'hetero'))
    disttype = input('Distractor type (homo/hetero)..............: ','s');
end
nStim    = input('Number of stimuli..........................: ');
nTrials  = input('Number of trials...........................: ');
settings.subjid = upper(subjid);
settings.stimtype = stimtype;
settings.disttype = disttype;
settings.nStim = nStim;
settings.nTrials = nTrials;

function im = drawEllipse(d1,d2,rot,fgcol,bgcol)
rot = -rot-90;  
d1 = round(d1);
d2 = round(d2);
% make sure that d1 is the minor axis
if (d1>d2)
    d3=d1;
    d1=d2;
    d2=d3;
end
% draw ellipse
rot = -rot/180*pi;
im = ones(2*d2,2*d2)*bgcol;
minX = -d2;
maxX = minX + 2*d2 - 1; 
[X Y] = meshgrid(minX:maxX,minX:maxX);
X_new = X * cos(rot) - Y * sin(rot);
Y = X * sin(rot) + Y * cos(rot);
X = X_new;
idx = (X.^2/(d1/2)^2 + Y.^2/(d2/2)^2)<1;
idx_low = (X.^2/((1.05*d1)/2)^2 + Y.^2/((1.025*d2)/2)^2)<1;
idx_super_low = (X.^2/((1.075*d1)/2)^2 + Y.^2/((1.05*d2)/2)^2)<1;
im(idx_super_low) = mean([fgcol bgcol bgcol]);
im(idx_low) = mean([fgcol bgcol]);
im(idx) = fgcol;
% crop
while im(:,1)==bgcol
    im = im(:,2:end);
end
while im(1,:)==bgcol
    im = im(2:end,:);
end
while im(end,:)==bgcol
    im = im(1:end-1,:);
end
while im(:,end)==bgcol
    im = im(:,1:end-1);
end
