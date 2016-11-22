%[dates, depths, temp] = list2table(data, tspan)
%Function converting dates, depths and temperature data from a three-column
%list "data" to a table "temp" with entry vectors "dates" and "depths".
%
% date1  depth1  temp11
% date1  depth2  temp12                    | date1  date2  ...
% date1  ...     ...                -------|------------------
% date2  depth1  temp21     -->     depth1 | temp11 temp21
% date2  depth2  temp22             depth2 | temp21 temp22
% date2  ...     ...                ...    | ...    ...
% ...
%
%The measurements are truncated according to the two elements in "tspan";
%if no period is specified, all measurements are considered.
%
%In case the depths vary from date to date, temperatures are linearly
%interpolated to obtain values at all measurement depths.

function [dates, depths, temp] = list2table(data, tspan)
    
    data = sortrows(data,1);
    dates = datenum(data(:,1));
    tspan = datenum(tspan);
    if nargin == 1 || any(isnan(tspan))
        tspan = [dates(1) dates(end)];
    elseif length(tspan)~=2
        error('The timespan vector must contain two values (start and end).')
    end
    data = data(dates>=tspan(1) & dates<=tspan(2),:); %Selecting data
    if isempty(data)
        warning('Timespan vector not consistent with data time series.')
    end
    
    dates = unique(datenum(data(:,1))); %Dates vector
    depths = sort(unique(data(:,2)),'descend'); %Depths vector
    
    temp = zeros(length(depths),length(dates));
    for j = 1:length(dates)
        D1 = data(datenum(data(:,1))==dates(j),2:3); %Depths and temperatures at that date
        D1 = sortrows(D1,-1); %Correcting order
        if sum(D1(:,1)==0)>1
            D1 = D1(2:end,:); %Removing double measurements at 0m
        end
        if length(D1(:,1)) < length(depths) %fI measurements are missing
            d = length(depths)-length(D1(:,1));
            D1(end+1:end+d,:) = [setdiff(depths,D1(:,1)) nan(d,1)]; %Filling missing values with NaNs
            D1 = sortrows(D1,-1); %Correcting order
            val_ok = ~isnan(D1(:,2));
            if sum(val_ok)>1
                %Interpolating linearly to replace the NaNs
                D1(:,2) = interp1(depths(val_ok),D1(val_ok,2),depths);
            end
        end
        temp(:,j) = D1(:,2);
    end
    
end