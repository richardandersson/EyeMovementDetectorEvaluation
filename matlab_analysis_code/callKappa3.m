function [kappaMat] = callKappa3(events, eventTypes, inputFiles, nameOfAlgorithm, stimuliTypes, stimchar)
% Provides a Cohen's Kappa that is aggregated across individual data files.
% Note that this still provide a Kappa value for every algorithm against
% every other algorithm, and separately for data from different stimuli
% types.
% stimchar is a character code ('i','d','v') specifying what stimuli data
% to include (image, dot, video).


% SELECT RELEVANT SUBSET OF DATA
tc = file2stimtype(inputFiles, stimuliTypes);
events = events(tc==stimchar); % subset



genKappa = struct(); % For storing Cohens Kappa by algorithm and event.
% Example: genKappa.eventtype(datafileindexnumber) = [matrix of comparisons]

% find out the total number of events and prepare to store all the data
f = fields(events(1));
numTot = 0;
for j = 1:length(events) % for every data file
    numTot = numTot + numel(events(j).(f{1}));
end
event_vector = zeros(numTot,length(nameOfAlgorithm), 'uint8'); % store data vector of all data (for two c


% extract all the data
idx_start = 1;
idx_stop = 0;
for j = 1:length(events) % for every data file
    idx_stop = idx_start - 1 + numel(events(j).(nameOfAlgorithm{1}));
    for i = 1:length(nameOfAlgorithm)
        event_vector(idx_start:idx_stop, i) = cast(events(j).(nameOfAlgorithm{i}), 'uint8');
    end
    idx_start = idx_start + numel(events(j).(nameOfAlgorithm{1})); % update index
end


% Initialize comparison matrix
for e = 1:length(eventTypes)
    kappaMat.(eventTypes{e}) = zeros( length(nameOfAlgorithm), length(nameOfAlgorithm) );
end
    
X = zeros(2,2);

% This variant repeats the coding vectors for each event type
for e = 1:length(eventTypes)
    evMat = event_vector == e; % has current event or not
    for i = 1:length(nameOfAlgorithm) % row algorithm
        for j = 1:length(nameOfAlgorithm) % column algorithm
            X(1,1) = sum( evMat(:,i) &  evMat(:,j) ); % 11 (true algo1, true algo2)
            X(1,2) = sum( evMat(:,i) & ~evMat(:,j) ); % 10
            X(2,1) = sum(~evMat(:,i) &  evMat(:,j) ); % 01
            X(2,2) = sum(~evMat(:,i) & ~evMat(:,j) ); % 00
            kappaMat.(eventTypes{e})(i,j) = kappa(X);
        end
    end
end


% -------------------------------------------------------------------------

end % end main function