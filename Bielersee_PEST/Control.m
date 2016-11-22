%This script writes the PEST control file, given matrices of measured
%data and given a set of parameters with their related properties. This
%script should be run once, to prepare an optimization with PEST.
clear all, close all

%Measurements, in time x depth matrices
fid = fopen('../LacLeman.txt');
zobs{1} = str2double(strsplit(fgetl(fid),','));
zobs{1} = zobs{1}(~isnan(zobs{1}));
obsval = textscan(fid,['%s' repmat('%f',1,25)],'Delimiter',',');
tobs{1} = datenum(obsval{1}(:,1),'yyyy-mm-dd HH:MM');
obsl = {'Temperature'}; obs = {'T'};
obsdat{1} = cell2mat(obsval(2:26));
sim_T = datenum({'1-Jan-1981','1-Jan-2014'}); %Remove measurements out of simulation period
obsdat{1} = obsdat{1}(tobs{1}>=sim_T(1) & tobs{1}<=sim_T(2),:);
tobs{1} = tobs{1}(tobs{1}>=sim_T(1) & tobs{1}<=sim_T(2));
return
%Define parameter names, group, values, min, max
parname = {'lat','pair','alpha','qNN','fwind','C10','CD','fgeo','kmin','p1','p2','beta','albsw'};
partype = {'fixed','fixed','none','none','none','none','none','none','fixed','none','none','none','none'};
parlim = {'factor','factor','factor','factor','factor','factor','factor','relative','factor','factor','factor','relative','relative'};
parval = [47.20  990 0.030 1.25 1.0 0.0017 0.002 0.5 1e-15 1.1 0.9 0.3 0.2];
parmin = [40.00  900 0.001 0.7 1.0 0.001 0.001 0.0 1e-30 0.5 0.5 0.0 0.0];
parmax = [50.00 1000 0.100 1.3 2.0 0.003 0.005 1.0 1e-09 1.5 1.5 0.4 0.3];
pargrp = {'none','none','fit','fit','fit','fit','fit','fit','fit','fit','fit','fit','fit'};
npar = length(parname);

%Total number of NaNs
nnan = 0;
for i=1:length(obsdat), nnan=nnan+sum(isnan(obsdat{i}(:))); end

%Write corresponding control file for PEST
fid = fopen('keps_calib.pst','w');
fprintf(fid,'pcf\r\n');
fprintf(fid,'* control data\r\n');
fprintf(fid,'restart estimation\r\n');
fprintf(fid,'%d %d 1 0 1\r\n',npar,sum(cellfun(@numel,obsdat))-nnan);
fprintf(fid,'1 1 single nopoint 1 0 0\r\n');
fprintf(fid,'5.0 2.0 0.3 0.01 10\r\n');
fprintf(fid,'5.0 5.0 0.001\r\n');
fprintf(fid,'0.1\r\n');
fprintf(fid,'30 0.005 4 3 0.01 3\r\n');
fprintf(fid,'1 1 1\r\n');
fprintf(fid,'* parameter groups\r\n');
fprintf(fid,' fit\trelative\t0.01\t0.00001\tswitch\t2.0\tparabolic\r\n');
fprintf(fid,'* parameter data\r\n');
for i=1:npar
    fprintf(fid,'%6s\t%s\t%s\t%10.4e\t%10.4e\t%10.4e\t%4s\t1.0\t0.0\t1\r\n',parname{i},partype{i},parlim{i},parval(i),parmin(i),parmax(i),pargrp{i});
end
fprintf(fid,'* observation groups\r\n');
fprintf(fid,'obs\r\n');
fprintf(fid,'* observation data\r\n');
for i=1:length(obsdat)
    k = 1;
    for it=1:length(tobs{i})
        for iz=1:length(zobs{i})
            if ~isnan(obsdat{i}(it,iz))
                fprintf(fid,'%s_%d_%d\t%12.4e\t1.0\tobs\r\n',obsl{i}(1:4),it,iz,obsdat{i}(it,iz));
            else %Store location of missing measurements so as not to write corresponding model observations
                nanobs{i}(k,:) = [it iz];
                k = k+1;
            end
        end
    end
end
fprintf(fid,'* model command line\r\n');
fprintf(fid,'model.bat\r\n');
fprintf(fid,'* model input/output\r\n');
fprintf(fid,'keps_par.tpl\tkepsilon_PEST.par\r\n');
fprintf(fid,'keps_obs.ins\tModelObs.dat\r\n');
fprintf(fid,'* prior information\r\n');
fclose(fid);

%Save measurement properties for script OutputInstructions.m
save('obsprop.mat','obs','tobs','zobs','nanobs');
fclose('all');