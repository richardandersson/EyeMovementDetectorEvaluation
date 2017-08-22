clear all, close all, clc

%--------------------------------------------------------------------------
% Init parameters
%--------------------------------------------------------------------------
% start path should be the Algorithms folder
addpath(genpath(cd))
% set your R path in gtwrapper.m also, if you use Windows.
% load([cd, filesep 'testData' filesep '1250Hz_3_Participants.mat']);
% load([cd, filesep 'testData' filesep '05UL27_trial17_labelled_RA.mat']);
% inputFolder = 'testData';
inputFolder = 'allData';
coderNames = {'MN', 'RA'}; % data file coder initials, only two used.
nameOfAlgorithm = {'coderMN', 'coderRA', 'CDT','EK','IDT','IDTk',...
    'IKF', 'IMST', 'IHMM', 'IVT', 'NH', 'BIT', 'LNS'};
% nameOfAlgorithm = {'IDT', 'IDTk'};
% nameOfAlgorithm = {'coderMN', 'IDT'};%, 'gazetoolsV', 'gazetoolsVA', 'gazetoolsVI'};
% "algorithm" name of coders must have this format of 'coder' + two initals
params = getParams(); % set general parameters
params.ploton = 0;
eventTypes = {'Fixation', 'Saccade', 'PSO', 'Pursuit', 'Blink', 'Undefined'};
% events - classification result (1- fixation, 2 - saccade,
% 3 - PSO,4 - pursuit,5 - blink,6 - undefined)
stimuliTypes = {'_img_', '_trial\d+_', '_video_'}; % a regexp pattern of filenames to batch by

OPTIMIZING = 1;

% addpath(genpath('/home/richard/Dropbox/Matlab'));

%--------------------------------------------------------------------------
% Set up input data
%--------------------------------------------------------------------------

inputFiles = getEtFileNames(inputFolder, coderNames);

% return % TEST: stop here if you only want a list of filenames.

%--------------------------------------------------------------------------
% Begin detection
%--------------------------------------------------------------------------
                 
% events - classification result (1- fixation, 2 - saccade,
% 3 - glissade,4 - pursuit,5 - blink,6 - undefined)

topResult = struct(); % Will store all results based stimuliTypes

for s = 1:length(stimuliTypes)+1 % +1 for storing 'all' as last category
    events = struct();
    fnames = {}; % empty files to evalute
    for i = 1:length(inputFiles) % for all input file candidates
        if s == length(stimuliTypes)+1 % all files as last category
            fnames = {fnames{:}, inputFiles{i}}; % add this file
        elseif ~isempty(regexp(inputFiles{i}, stimuliTypes{s})) % if file from correct stimulus
            fnames = {fnames{:}, inputFiles{i}}; % add file to evalute
        end
    end
    for i = 1:length(fnames) % for every data file (events(i) will match list in inputFiles
        fname = fnames{i}
%         if isempty(stimuliTypes) % if no types specified, test all
%             fname = inputFiles{i}
%         else
%             if ~isempty(regexp(inputFiles{i}, stimuliTypes{s})) % if file is of particular stimulus type.
%                 fname = inputFiles{i}
%             else
%                 continue
%             end            
%         end
%         [x,y] = ETload([inputFolder filesep fname]);
        ET = load( [inputFolder filesep fname] );
        ET.ETdata.pos = POScleaner(ET.ETdata.pos); % amends tagging-GUI bug
        x = ET.ETdata.pos(:,4);
        y = ET.ETdata.pos(:,5);
        
        for j = 1:length(nameOfAlgorithm)
            algoName = nameOfAlgorithm{j};
            fprintf('%s',['Detecting events: ',algoName])
            tic
            switch algoName
                 case 'coderMN'
                     events(i).coderMN = ET.ETdata.pos(:,6);
                 case 'coderRA'
                     temp = [fname(1:end-6) 'RA.mat'];
                     temp = load(temp);
                     temp.ETdata.pos = POScleaner(temp.ETdata.pos);
                     events(i).coderRA = temp.ETdata.pos(:,6);
                case 'CDT' % FDT
                    events(i).CDT = veneri(x',y',params);             
                case 'EK' % Engbert & Kliegl
                    events(i).EK = engbert(x',y',params)';            
                case 'IDT'
                    events(i).IDT = IDT(x',y',params);
                case 'IDTk'
                    events(i).IDTk = IDTk(x',y',params);                
                case 'IKF'
                    events(i).IKF = IKF(x',y',params);         
                case 'IMST'
                    events(i).IMST = IMST(x',y',params); 
                case 'IHMM'
                    events(i).IHMM = IHMM(x',y',params); 
                case 'IVT'
                    events(i).IVT = IVT(x',y',params);                 
                case 'NH' % Nyström & Holmqvist
                    events(i).NH = nyst(x',y',params);
                case 'BIT' % BIT
                    events(i).BIT = vdLans(x',y',params);   
                case 'gazetoolsV'
                    events(i).gazetoolsV = gtwrapper('V',x,y,params);
                case 'gazetoolsVA'
                    events(i).gazetoolsVA = gtwrapper('VA',x,y,params);
                case 'gazetoolsVI'
                    events(i).gazetoolsVI = gtwrapper('VI',x,y,params);
                case 'LNS'
                    events(i).LNS = Larsson(fname);
                otherwise
                    error('Wrong algo name.')
            end
            t = toc;
            fprintf(', %0.1f s\n',t);

        end
        
        % Optional plotting of events after each data file
        if params.ploton == 1;
            scarfPlots(events, inputFolder, inputFiles, i)
            pause
        end
        
    end
    topResult.(['stimtype_' num2str(s)]) = events;
end
fprintf('%s\n','Done!')

% Plotting function
% scarfPlots(events, inputFolder, inputFiles, 12)

% Tidy things up
% clear ET x y algoName fname i f j

% NO FILTERING OF EVENTS
% STANDARD OUTPUT SHOULD BE A VECTOR OF EVENT NUMBERS (0 and 1:6)
% SUBSELECT WITHIN BLOCKS/FUNCTIONS ONLY

% temp here
[eventSummary, eventObject] = eventStats(topResult.stimtype_1, eventTypes, nameOfAlgorithm, params);
[e_img, ~] = eventStats(topResult.stimtype_1, eventTypes, nameOfAlgorithm, params);
[e_dot, ~] = eventStats(topResult.stimtype_2, eventTypes, nameOfAlgorithm, params);
[e_vid, ~] = eventStats(topResult.stimtype_3, eventTypes, nameOfAlgorithm, params);

[e_all, ~] = eventStats(topResult.stimtype_4, eventTypes, nameOfAlgorithm, params);
    
return % stops main detection loop, leaving the rest optional.

%--------------------------------------------------------------------------
% EXTRACTING META INFORMATION
%--------------------------------------------------------------------------
noSamples = countSamples(inputFolder, inputFiles);
noStimSamp = [countStimuliSamples(inputFolder, inputFiles, stimuliTypes{1})
    countStimuliSamples(inputFolder, inputFiles, stimuliTypes{2})
    countStimuliSamples(inputFolder, inputFiles, stimuliTypes{3}) ]'% use this to get noSample for each stimuli type.
participantList = unique(arrayfun(@(a)a{1}(3:4), inputFiles, 'UniformOutput', false))
empiricalEvents(events) % What events to they actually detect? Return codes.
rmsd = empRMSD('allData', inputFiles, events, params); % precision of the data

%--------------------------------------------------------------------------
% PLOTTING EVENT SIMILARITIES
%--------------------------------------------------------------------------
bubblePlot2(e_img, e_dot, e_vid, nameOfAlgorithm)
% export_fig('bubbleplot2.pdf', '-transparent')

%--------------------------------------------------------------------------
% GENERALIZED COHENS KAPPA FOR ALL EVENTS
%--------------------------------------------------------------------------
kappaObject = callKappa1(events, eventTypes, nameOfAlgorithm);
% should have another Kappa that aggregates data from alla data files
% N.B. this may be a massive memory hog
%kappaSummary = callKappa2(events, eventTypes, nameOfAlgorithm);
% Kappa separate for algorithms and stimulus type, but aggregated over data files:
% Images:
kappaMat = callKappa3(events, eventTypes, inputFiles, nameOfAlgorithm, stimuliTypes, 'i')
% Dots:
kappaMat = callKappa3(events, eventTypes, inputFiles, nameOfAlgorithm, stimuliTypes, 'd')
% Videos
kappaMat = callKappa3(events, eventTypes, inputFiles, nameOfAlgorithm, stimuliTypes, 'v')
% Full table (sort of)
kappaMat = callKappa6(events, eventTypes, inputFiles, nameOfAlgorithm, stimuliTypes, {'coderMN', 'coderRA'})
printFieldTables(kappaMat, nameOfAlgorithm, stimuliTypes)

% You can look at Kappa against individual coders:
kappaMat = callKappa6(events, eventTypes, inputFiles, nameOfAlgorithm, stimuliTypes, {'coderMN'})
printFieldTables(kappaMat, nameOfAlgorithm, stimuliTypes)


%--------------------------------------------------------------------------
% EVENT DISTRIBUTIONS PRODUCED BY ALGORITHMS
%--------------------------------------------------------------------------
% [eventSummary, eventObject] = eventStats(topResult.stimtype_2, eventTypes, nameOfAlgorithm, params);
% eventVisual1(eventSummary, eventObject); % Print summary table of all events
% eventVisual2(eventSummary, eventObject); % Ugly summary histogram of all events.
% % rows are algorithms; columns are number types
% % Fixations
% [e_img, ~] = eventStats(topResult.stimtype_1, eventTypes, nameOfAlgorithm, params);
% [e_dot, ~] = eventStats(topResult.stimtype_2, eventTypes, nameOfAlgorithm, params);
% [e_vid, ~] = eventStats(topResult.stimtype_3, eventTypes, nameOfAlgorithm, params);
% mat2tex(e_img.Fixation, nameOfAlgorithm, {'Algorithm', 'Mean', 'SD', 'Count'}, 0)
% mat2tex(e_dot.Fixation, nameOfAlgorithm, {'Algorithm', 'Mean', 'SD', 'Count'}, 0)
% mat2tex(e_vid.Fixation, nameOfAlgorithm, {'Algorithm', 'Mean', 'SD', 'Count'}, 0)
% mat2tex(e_img.Saccade, nameOfAlgorithm, {'Algorithm', 'Mean', 'SD', 'Count'}, 0)
% mat2tex(e_dot.Saccade, nameOfAlgorithm, {'Algorithm', 'Mean', 'SD', 'Count'}, 0)
% mat2tex(e_vid.Saccade, nameOfAlgorithm, {'Algorithm', 'Mean', 'SD', 'Count'}, 0)
% mat2tex(e_img.PSO, nameOfAlgorithm, {'Algorithm', 'Mean', 'SD', 'Count'}, 0)
% mat2tex(e_dot.PSO, nameOfAlgorithm, {'Algorithm', 'Mean', 'SD', 'Count'}, 0)
% mat2tex(e_vid.PSO, nameOfAlgorithm, {'Algorithm', 'Mean', 'SD', 'Count'}, 0)
% mat2tex(e_img.Pursuit, nameOfAlgorithm, {'Algorithm', 'Mean', 'SD', 'Count'}, 0)
% mat2tex(e_dot.Pursuit, nameOfAlgorithm, {'Algorithm', 'Mean', 'SD', 'Count'}, 0)
% mat2tex(e_vid.Pursuit, nameOfAlgorithm, {'Algorithm', 'Mean', 'SD', 'Count'}, 0)
% mat2tex(e_img.Blink, nameOfAlgorithm, {'Algorithm', 'Mean', 'SD', 'Count'}, 0)
% mat2tex(e_dot.Blink, nameOfAlgorithm, {'Algorithm', 'Mean', 'SD', 'Count'}, 0)
% mat2tex(e_vid.Blink, nameOfAlgorithm, {'Algorithm', 'Mean', 'SD', 'Count'}, 0)

%--------------------------------------------------------------------------
% SIMPLE AGREEMENT SCORE (PROPORTION OF IDENTICAL CODES BY TWO PARTIES)
%--------------------------------------------------------------------------
simpAgreeMat = simpleAgreement(events, nameOfAlgorithm);
num2str(simpAgreeMat, '  %0.2f')

% per stimuli type and event type against humans, in proportions
resultsMatrix = simpleAgreement2(topResult, nameOfAlgorithm, eventTypes)

% per stimuli, algorithm and event type against humans, as Cohen's Kappa
resultsMatrix2 = simpleAgreement3(topResult, nameOfAlgorithm, eventTypes)
resultsMatrix2 = cell2mat(resultsMatrix2)
num2str(resultsMatrix2(:,1:9), '   %0.2f')
%--------------------------------------------------------------------------
% VISUALIZATION OF COHENS KAPPA RESULTS
%--------------------------------------------------------------------------

%--------------------------------------------------------------------------
% ALL-EVENT CONFUSION MATRIX COMPARISONS
%--------------------------------------------------------------------------
confMatrix = confusion(events, eventTypes, 'coderMN', 'coderRA')
% first algo is row algo in matrix, the second is column algo.

% "Source accumulation" of confusion -- areas to improve algorithms?
confMatrix2 = confusion3(topResult, eventTypes, nameOfAlgorithm)

% (Yes, this is a pain in the ass to manually transfer to a table - but no
% function for it yet). Below is a partial solution: remember to set both
% indices
conftemp = ([confMatrix2{:,6}])'; % index number is event type
conftemp2 = conftemp(3:3:end,:); % first index: 1 for images, 2 for dots, 3 for vid
conftemp2 = conftemp2';
conftemp2(:) % displays one column with each algorithms over- and underclassification
% the column should be 2*length(nameOfAlgorithm) in length

tempname = repmat(nameOfAlgorithm,2,1);
tempname(:)

noTemp = cell2mat(confMatrix2(:,end)); % get columns of total disagreeing samples
bsxfun(@rdivide, noTemp, noStimSamp)*100 % Percentage disagreeing samples; noStimSamp defined above


%--------------------------------------------------------------------------
% RECONSTRUCTED ALGORITHM PARAMETERS FROM HUMANS
%--------------------------------------------------------------------------
recParams = reconstructParams(inputFolder, inputFiles, params);

%--------------------------------------------------------------------------
% ALGORITHM RANKINGS
%--------------------------------------------------------------------------
% rankMatrix = simpleMinSS(e_dot.Fixation, e_dot.Fixation(1:2,:))
% horzcat(nameOfAlgorithm', num2cell(rankMatrix))
simpleMinSS(e_all.Fixation, e_all.Fixation(1:2,:), nameOfAlgorithm)
simpleMinSS(e_all.Saccade, e_all.Saccade(1:2,:), nameOfAlgorithm)


% Rankings assumes the baseline "the average" human.
simpleMinSS(e_img.Fixation, e_img.Fixation(1:2,:), nameOfAlgorithm)
simpleMinSS(e_dot.Fixation, e_dot.Fixation(1:2,:), nameOfAlgorithm)
simpleMinSS(e_vid.Fixation, e_vid.Fixation(1:2,:), nameOfAlgorithm)

simpleMinSS(e_img.Saccade, e_img.Saccade(1:2,:), nameOfAlgorithm)
simpleMinSS(e_dot.Saccade, e_dot.Saccade(1:2,:), nameOfAlgorithm)
simpleMinSS(e_vid.Saccade, e_vid.Saccade(1:2,:), nameOfAlgorithm)

simpleMinSS(e_img.PSO, e_img.PSO(1:2,:), nameOfAlgorithm)
simpleMinSS(e_dot.PSO, e_dot.PSO(1:2,:), nameOfAlgorithm)
simpleMinSS(e_vid.PSO, e_vid.PSO(1:2,:), nameOfAlgorithm)
    % When comparing against a hypothetical "average human", then both
    % humans are (unfairly) more similar to this person than any algorithm.

% --- misc events
simpleMinSS(e_img.PSO, e_img.PSO(1:2,:), nameOfAlgorithm)
simpleMinSS(e_dot.PSO, e_dot.PSO(1:2,:), nameOfAlgorithm)
simpleMinSS(e_vid.PSO, e_vid.PSO(1:2,:), nameOfAlgorithm)

simpleMinSS(e_img.Pursuit, e_img.Pursuit(1:2,:), nameOfAlgorithm)
simpleMinSS(e_dot.Pursuit, e_dot.Pursuit(1:2,:), nameOfAlgorithm)
simpleMinSS(e_vid.Pursuit, e_vid.Pursuit(1:2,:), nameOfAlgorithm)

simpleMinSS(e_img.Blink, e_img.Blink(1:2,:), nameOfAlgorithm)
simpleMinSS(e_dot.Blink, e_dot.Blink(1:2,:), nameOfAlgorithm)
simpleMinSS(e_vid.Blink, e_vid.Blink(1:2,:), nameOfAlgorithm)

% ---------------
    
    
% Rankings against Human MN (1)
simpleMinSS(e_img.Fixation, e_img.Fixation(1,:), nameOfAlgorithm)
simpleMinSS(e_dot.Fixation, e_dot.Fixation(1,:), nameOfAlgorithm)
simpleMinSS(e_vid.Fixation, e_vid.Fixation(1,:), nameOfAlgorithm)
    % Note that for fix in vid, MN (rank 1) and NH (rank 2) are better
    % than RA.

simpleMinSS(e_img.Saccade, e_img.Saccade(1,:), nameOfAlgorithm)
simpleMinSS(e_dot.Saccade, e_dot.Saccade(1,:), nameOfAlgorithm)
simpleMinSS(e_vid.Saccade, e_vid.Saccade(1,:), nameOfAlgorithm)
    % LNS more similar to MN than RA for vid and img, but not dot.

simpleMinSS(e_img.PSO, e_img.PSO(1,:), nameOfAlgorithm)
simpleMinSS(e_dot.PSO, e_dot.PSO(1,:), nameOfAlgorithm)
simpleMinSS(e_vid.PSO, e_vid.PSO(1,:), nameOfAlgorithm)
    % Note that LNS (2) and NH (3) are closer (than RA) to MN for img and vid, but not dot.

% Rankings against Human RA
simpleMinSS(e_img.Fixation, e_img.Fixation(2,:), nameOfAlgorithm)
simpleMinSS(e_dot.Fixation, e_dot.Fixation(2,:), nameOfAlgorithm)
simpleMinSS(e_vid.Fixation, e_vid.Fixation(2,:), nameOfAlgorithm)
    % MN closest all the time.

simpleMinSS(e_img.Saccade, e_img.Saccade(2,:), nameOfAlgorithm)
simpleMinSS(e_dot.Saccade, e_dot.Saccade(2,:), nameOfAlgorithm)
simpleMinSS(e_vid.Saccade, e_vid.Saccade(2,:), nameOfAlgorithm)
    % MN closest all the time.

simpleMinSS(e_img.PSO, e_img.PSO(2,:), nameOfAlgorithm)
simpleMinSS(e_dot.PSO, e_dot.PSO(2,:), nameOfAlgorithm)
simpleMinSS(e_vid.PSO, e_vid.PSO(2,:), nameOfAlgorithm)
    % LNS closer than MN all the time. NH always worst.
    
%--------------------------------------------------------------------------
% INTRACODER RELIABILITY
%--------------------------------------------------------------------------
MN_irr_events = intraCoderReliability('MN', eventTypes);