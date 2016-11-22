%This script reads the output of the k-epsilon model, extracts the data
%that can be compared to actual measurements, rewrites this data into a new
%text file and writes the PEST instruction file that allows PEST to
%understand it. This script is run by PEST after each kepsilon simulation.
if ~exist('wdir','var'), wdir='';
elseif wdir(end)~='\', wdir = [wdir '\']; end

%Path to kepsilon model directory (for function GetResults)
addpath('..');
%Types, times, depths and NaNs of available measurements from script Control.m
load('obsprop.mat');

%Extract zobs x tobs matrices
fid=fopen([wdir 'kepsilon_PEST.par']);
for i=1:7; line=fgetl(fid); end
fclose(fid);
for i=1:length(obs)
    [time{i},z{i},data(i),label(i)]=GetResults(line,obs(i),[],zobs{i});
end
for i=1:length(obs)
    if (min(tobs{i})<min(time{i}) || max(tobs{i})>max(time{i}))
        warning(['Some given %s measurements fall out of simulated time period. '...
                 'NaN values are to be expected in the model observations file, which PEST won''t accept. '...
                 'This should be taken care of in the script Control.m.'],obs{i})
    end
    if (min(zobs{i})<min(z{i}) || max(zobs{i})>max(z{i}))
        warning(['Some given %s measurements fall out of simulated depth range. '...
                 'NaN values are to be expected in the model observations file, which PEST won''t accept. '...
                 'This should be taken care of in the script Control.m.'],obs{i})
    end
end

for i=1:length(data)
    data{i} = interp1(time{i},data{i}',tobs{i})';
end

%Rewrite model observations in one file
fid = fopen([wdir 'ModelObs.dat'],'w');
for i=1:length(data)
    for it=1:length(tobs{i})
        for iz=1:length(zobs{i})
            if isempty(nanobs{i}) || ~any(it==nanobs{i}(:,1) & iz==nanobs{i}(:,2))
                fprintf(fid,'%12.4e',data{i}(iz,it));
            end
        end
        fprintf(fid,'\r\n');
    end
end
fclose(fid);

%Write corresponding instructions file for PEST
fid = fopen([wdir 'keps_obs.ins'],'w');
fprintf(fid,'pif @\r\n');
for i=1:length(data)
    for it=1:length(tobs{i})
        strf = 'l1';
        for iz=1:length(zobs{i})
            if isempty(nanobs{i}) || ~any(it==nanobs{i}(:,1) & iz==nanobs{i}(:,2))
                strf = [strf ' [' label{i}(1:4) sprintf('_%d_%d',it,iz) ']' sprintf('%d:%d',12*iz-11,12*iz)];
            end
        end
        fprintf(fid,[strf '\r\n']);
    end
end
fclose(fid);
fclose('all');