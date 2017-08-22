function [kappaMat] = callKappa6(events, eventTypes, inputFiles, nameOfAlgorithm, stimuliTypes, standard)
% Provides a Cohen's Kappa that is aggregated across individual data files.
% Note that this still provide a Kappa value for every algorithm against
% every other algorithm, and separately for data from different stimuli
% types.
%
% Version 6 should aid in providing a table in the form of:
% EVENT1             EVENT2           EVENT3
% image dot vid      image dot vid    image dot vid
% 
% The function will actually produce a struct with a field for each event,
% but that field contains the different stimuli as columns.

% each human against each human
% selectable whether algorithms are compared again "average human" or
% a particular human.

% Argument "standard" is a cell array of algorithm/human names that are the
% comparison standard. If array size > 1, then the standard is the average
% of the two (or more) algorithms/humans in the array.

%standard = {'coderMN', 'coderRA'};
stimchars = ['i', 'd', 'v'];
kappaMat = struct();

tc = file2stimtype(inputFiles, stimuliTypes); % what file is what stimuli type, give one-letter vector

for i = 1:length(eventTypes) % for every event (number)
    kappaMat.(eventTypes{i}) = zeros( length(nameOfAlgorithm), length(stimuliTypes) );
    for j = 1:length(stimuliTypes) % for every stimuli type (string)
        % preallocate memory
        events2 = events(tc==stimchars(j)); % select by logical vector
        f = fields(events2(1));
        numTot = 0;
        for q = 1:length(events2) % for every data file
            numTot = numTot + numel(events2(q).(f{1})); % number of elements for files of this stimuli type
        end
        event_vector = zeros(numTot,length(nameOfAlgorithm), 'uint8'); % store data vector of all data
        
        % extract all the data
        idx_start = 1;
        idx_stop = 0;
        for q = 1:length(events2) % for every data file (of relevant stimuli type)
            idx_stop = idx_start - 1 + numel(events2(q).(nameOfAlgorithm{1}));
            for r = 1:length(nameOfAlgorithm)
                event_vector(idx_start:idx_stop, r) = cast(events2(q).(nameOfAlgorithm{r}), 'uint8');
            end
            idx_start = idx_start + numel(events2(q).(nameOfAlgorithm{1})); % update index
        end
        
        % multiply if several coders constitute the baseline
        % Up until here, the event vector is a matrix of rows of samples x
        % algorithms. Cells contain the classification numbers.
%         event_vector = repmat(event_vector, [], length(standard));
        event_vector = repmat(event_vector, [length(standard),1]);
        comp_vector = zeros(size(event_vector,1), size(event_vector,2) ,1); % create comparison vector
        for z = 1:length(standard) % populate comparison vector
            [~, c_id] = ismember(standard(z), nameOfAlgorithm); % gets index of coders constituting standard
            idx_start = ((z-1) * length(event_vector)/length(standard)) + 1;
            idx_stop = z * ( length(event_vector)/length(standard) );
            comp_vector(idx_start:idx_stop,1) = event_vector(idx_start:idx_stop,c_id);
        end
       
        
        for k = 1:length(events2) % for every relevant data file
            evMat = event_vector == i; % has current event or not
            cpMat = comp_vector == i; % same for comparison vector
            
            X = zeros(2,2);
            for l = 1:length(nameOfAlgorithm) % row algorithm
                %for m = 1:length(nameOfAlgorithm) % column algorithm
                X(1,1) = sum( evMat(:,l) &  cpMat(:,1) ); % 11 (true algo1, true algo2)
                X(1,2) = sum( evMat(:,l) & ~cpMat(:,1) ); % 10
                X(2,1) = sum(~evMat(:,l) &  cpMat(:,1) ); % 01
                X(2,2) = sum(~evMat(:,l) & ~cpMat(:,1) ); % 00
                kappaMat.(eventTypes{i})(l,j) = kappa(X); % <- XXX
                %end
            end
        end
    end
end




% % SELECT RELEVANT SUBSET OF DATA
% 
% events = events(tc==stimchar); % subset
% 
% 
% 
% genKappa = struct(); % For storing Cohens Kappa by algorithm and event.
% % Example: genKappa.eventtype(datafileindexnumber) = [matrix of comparisons]
% 
% % find out the total number of events and prepare to store all the data
% f = fields(events(1));
% numTot = 0;
% for j = 1:length(events) % for every data file
%     numTot = numTot + numel(events(j).(f{1}));
% end
% event_vector = zeros(numTot,length(nameOfAlgorithm), 'uint8'); % store data vector of all data (for two c
% 
% 
% % extract all the data
% idx_start = 1;
% idx_stop = 0;
% for j = 1:length(events) % for every data file
%     idx_stop = idx_start - 1 + numel(events(j).(nameOfAlgorithm{1}));
%     for i = 1:length(nameOfAlgorithm)
%         event_vector(idx_start:idx_stop, i) = cast(events(j).(nameOfAlgorithm{i}), 'uint8');
%     end
%     idx_start = idx_start + numel(events(j).(nameOfAlgorithm{1})); % update index
% end
% 
% 
% % Initialize comparison matrix
% for e = 1:length(eventTypes)
%     kappaMat.(eventTypes{e}) = zeros( length(nameOfAlgorithm), length(nameOfAlgorithm) );
% end
%     
% X = zeros(2,2);
% 
% for e = 1:length(eventTypes)
%     evMat = event_vector == e; % has current event or not
%     for i = 1:length(nameOfAlgorithm) % row algorithm
%         for j = 1:length(nameOfAlgorithm) % column algorithm
%             X(1,1) = sum( evMat(:,i) &  evMat(:,j) ); % 11 (true algo1, true algo2)
%             X(1,2) = sum( evMat(:,i) & ~evMat(:,j) ); % 10
%             X(2,1) = sum(~evMat(:,i) &  evMat(:,j) ); % 01
%             X(2,2) = sum(~evMat(:,i) & ~evMat(:,j) ); % 00
%             kappaMat.(eventTypes{e})(i,j) = kappa(X);
%         end
%     end
% end
% 


% -------------------------------------------------------------------------

end % end main function