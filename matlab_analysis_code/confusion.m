function [confMatrix] = confusion(events, eventTypes, algo1, algo2)
% Creates a confusion matrix of two different algorithms/coders and all
% events

confMatrix = zeros(length(eventTypes), length(eventTypes));

for i = 1:length(events) % for every data file
    for j = 1:length(events(i).(algo1)) % for every sample
        confMatrix( events(i).(algo1)(j), events(i).(algo2)(j) ) = ...
            confMatrix( events(i).(algo1)(j), events(i).(algo2)(j) ) + 1;
    end
end