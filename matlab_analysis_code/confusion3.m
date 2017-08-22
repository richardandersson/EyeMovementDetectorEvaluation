function [confMatrix] = confusion2(topResult, eventTypes, nameOfAlgorithm)
% Creates a confusion matrix of two different algorithms/coders and all
% events. Assumes two first algorithms are the human coders! Only uses
% first three stimuli types in topResult.
% Same as confusion2.m except outputs and additional column with total
% proportion of disagreeing samples.

% Create a confusion matrix visavi human coders, and quantify for every
% algorithm the source of confusion for every event.
% Distinguish between overclassification and underclassification.

% Rows: algorithms
% Columns: events
% Cells: over_stim1, over_stim2, over_stim3
%        under_stim1, under_stim2, under_stim3
%
% XXX   Fixation       Saccade
% Alg1  0.1/0.3/0.2    0.1/0.4/0.2
%       0.1/0.4/0.2    0.1/0.3/0.2
% Alg2
%
% Example outcome:
% If human coder scores a sample as belonging to a fixation, but Alg1
% scores it as belonging to a saccade, then for Alg1 the saccade column
% gets +1 for overclassification and the fixation column gets +1 for
% underclassification. These scores are later transformed to ratios.
%
% Interpretation: Alg1 is more conservative than humans, with regard to
% fixations, but slightly more liberal concerning saccades. Algorithm 2 is
% slightly liberal with fixations, but show no bias for saccades.
% Furthermore, more errors are related to fixations for algorithm 2, so
% efforts should be concentrated there.
% Alg1 has more misclassified samples (compared to human experts) than
% Alg2, so overall Alg2 is better. Alg1's primary disadvantage is in
% distguishing between fixations and saccades, and an improvement 0.5 %
% correct classified samples dinguishing between fixation and saccade would
% be enough to make this algorithm the winner.
%

% Events are slightly recoded to:
% Fixation = 1
% Saccade = 2
% PSO = 3
% Pursuit = 4
% Blink = 5
% Other = 6 (including undefined and not coded)

stims = fieldnames(topResult); % name of stimuli types
coders = {nameOfAlgorithm{1:2}}; % names of human coders.

confMatrix = cell(length(nameOfAlgorithm), length(eventTypes)); % +1 for null classifications
confMatrix(:) = {[0 0 0; 0 0 0]}; % over_stim1--3; under_stim1--3


% Step 1: count samples
for s = 1:length(stims)-1 % for every stimuli type except 'all'.
    for i = 1:length(nameOfAlgorithm) % for every algorithm (econfMatrix2 = confusion2(topResult, eventTypes, nameOfAlgorithm)
        a = nameOfAlgorithm{i};
        for j = 1:length(topResult.(stims{s})) % for every data file
            for k = 1:length(topResult.(stims{s})(j).(a)) % for every data sample
                for c = 1:length(coders) % for every human coder
                    if strcmp(nameOfAlgorithm{i}, coders{c})
                        continue
                    end
                    algEvent = topResult.(stims{s})(j).(a)(k);
                    humEvent = topResult.(stims{s})(j).(coders{c})(k);
                    if algEvent ~= humEvent
                        if algEvent == 0
                            algEvent = 6;
                        end
                        % store overclassfication
                        confMatrix{i,algEvent}(1,s) = confMatrix{i,algEvent}(1,s) + 1; % minus 2 for no of human coders
                        % store underclassification
                        confMatrix{i,humEvent}(2,s) = confMatrix{i,humEvent}(2,s) + 1;
                    end
                end
            end
        end
    end
end


% Step 2: Create a ratio of the numbers
cellColumn = cell(size(confMatrix,1),1);
for i = 1:size(confMatrix,1) % for every algorithm
    sumMatrix = [0,0,0;0,0,0];
    for j = 1:size(confMatrix,2) % for every event class
        sumMatrix = sumMatrix + confMatrix{i,j};
    end
    for j = 1:size(confMatrix,2) % for every event class (again)
        confMatrix{i,j} = confMatrix{i,j}./sumMatrix;
    end
    cellColumn{i} = sumMatrix; % store the number of all disagreeing samples
end

% Note, Double counting because each sample is both over- and
% under-classified. Each algorithm is also counted once for each human
% coder (two of them).
for i = 1:size(cellColumn)
    if i <= length(coders)
        cellColumn{i} = cellColumn{i}(1,:);
    elseif i > length(coders)
        cellColumn{i} = cellColumn{i}(1,:)/2;
    end
end
confMatrix = [confMatrix cellColumn];



