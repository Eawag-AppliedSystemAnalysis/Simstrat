%Run the k-epsilon model in a sequenced way i.e. changing initial
%conditions 
clear all; close all; clc;
cd D:\kepsilon\LakeBiel\Love\Run1\Bin
addpath(genpath('D:\kepsilon\LakeBiel\Love\Run1\Files_FinalRun'))

%% Activate certain model parts
% Inflow or not
% Run model with no inflow (=0)
% Run model with inflow (=1)
zero=1;

% Sane Redirection
% Run model with redirection (=1)
% Run model without redirection (=0)
removal=0;


 %  Muhlberg heating or not
%  %Mühleberg heat pressent (=1)
%  %without Mühleberg (=0)
    Muhl=1;
    
%     Run with climat change forcing (=1) 
Climat=1;
if Climat==1
JJJJ=1;KKKK=1;
end




%% Time step, Time output and sequence length
dt=600;% model time step [sec]
% OBS change Output_time_Biel.dat
% dt=600 should be 144 time steps for dayly output
% dt=60 should be 1440 time steps for dayly output

time_res=6/24;%time resolution for river input [days] also for model out put


%How long the model should run (in days) befor new initial conditions are
%made (min req. 15 days see batch run report)
%30 days maximum to make winter and summer alfa work

 SequenceVect=20;%[Days]


%% parameters
% % RMS best fit (rivers applied)
% %     P1        P2         q        AlfaS      AlfaW    CD        C10       K
% %     1.3000    1.2000    1.3000    0.0060    0.0040    0.0050    0.0016    0.7000
% % 
% %     RMS      Corr.     Max. Diff.
% %     0.7333    0.9900    4.7507
% % %          best fit values new (RMS) Use this since it gives an overal
% % %          better model and RMS is more used than RMSD
     p1  = 1.3;
     p2  = 1.2;
     q   = 1.3;
     aS  = 0.006; 
     aW  = 0.004;
     CD  = 0.005;
     C10 = 0.0016;
     K   = 0.7;   
    
% % RMS best fit Calibration without river input
% %     P1        P2         q        AlfaS      AlfaW    CD        C10       K
% %     1.3000    1.5000    1.3000    0.0090    0.0010    0.0030    0.0014    0.8000
% % 
% %     RMS      Corr.     Max. Diff.
% %     0.9616    0.9780    5.2123
% %          best fit values new (RMS) Use this since it gives an overal
% %          better model and RMS is more used than RMSD
%      p1  = 1.3;
%      p2  = 1.5;
%      q   = 1.3;
%      aS  = 0.009; 
%      aW  = 0.001;
%      CD  = 0.003;
%      C10 = 0.0014;
%      K   = 0.8;

% % RMS best fit Calibration with river input and COSMO2 forcing data 
% RMS best fit
%     P1        P2         q        AlfaS      AlfaW    CD        C10       K
%     1.2000    0.9000    1.3000    0.0060    0.0010    0.0020    0.0018    0.6000
% 
%     RMS      Corr.     Max. Diff.
%     0.8758    0.9811    4.5904     



     
%          %best fit values new (RMSD) Use this since it gives an overal
%          %better model
%      p1  = 1.3;
%      p2  = 2.1;
%      q   = 1.3;
%      aS  = 0.006; 
%      aW  = 0.004;
%      CD  = 0.002;
%      C10 = 0.0010;
%      K   = 0.6;
%      

    %best fit values new (correlation)
%      p1  = 1.3;
%      p2  = 1.8;
%      q   = 1.6;
%      aS  = 0.006; 
%      aW  = 0.002;
%      CD  = 0.01;
%      C10 = 0.0012;
%      K   = 0.7;
 

% Calibration of model (5 values/variable to fit design)
% P1= 0.9:0.1:1.3;
% P2=0.9:0.3:2.1;
% q=1:0.3:2.2;%q=1.8:0.2:2.6;
% AlfaS=0.006:0.003:0.018;
% AlfaW=0.001:0.001:0.005;
% CD=0.001:0.001:0.005;
% C10=0.001:0.0002:0.0018;
% K= [0.6:0.1:1.1 1.3:0.2:1.7];
% 
% %Fractional factorial design (7 variables, 196 runs, 5 levels (values) in each variable)
% x=sortrows(rowexch(7,100,'linear','cat',1:7,'levels',5,'tries',10));
% CalMatrix=zeros(length(x),7);
% for I=1:length(x)
%     CalMatrix(I,1)=P1(x(I,1));
%     CalMatrix(I,2)=P2(x(I,2));
%     CalMatrix(I,3)=q(x(I,3));
%     CalMatrix(I,4)=AlfaS(x(I,4));
%     CalMatrix(I,5)=AlfaW(x(I,5));
%     CalMatrix(I,6)=CD(x(I,6));
%     CalMatrix(I,7)=C10(x(I,7));
% end
% clear P1 P2 q AlfaS AlfaW Cd C10
%         file=['D:\kepsilon\LakeBiel\Love\Run1\Files_FinalRun\Output\Result\Sequence\Run_Files\CalibrationMatrix.mat'];
%         save(file,'CalMatrix','K');
%%Restart run with the calibration matrix, 
% load(file)

%Factorial design of model calibration
% CalMatrix=zeros(length(P1)*length(P2)*length(q)*length(AlfaS)*length(AlfaW)*length(CD)*length(C10),7);
% GG=0;
% for I=1:length(P1)
%     for II=1:length(P2)
%         for III=1:length(q)
%             for IIII=1:length(AlfaS)
%                  for IIIII=1:length(AlfaW)
%                      for IIIIII=1:length(CD)
%                          for IIIIIII=1:length(C10)
% 
%                    CalMatrix(GG*length(C10)+IIIIIII,1)=P1(I);
%                    CalMatrix(GG*length(C10)+IIIIIII,2)=P2(II);
%                    CalMatrix(GG*length(C10)+IIIIIII,3)=q(III);
%                    CalMatrix(GG*length(C10)+IIIIIII,4)=AlfaS(IIII);
%                    CalMatrix(GG*length(C10)+IIIIIII,5)=AlfaW(IIIII);
%                    CalMatrix(GG*length(C10)+IIIIIII,6)=CD(IIIIII);
%                    CalMatrix(GG*length(C10)+IIIIIII,7)=C10(IIIIIII);
% 
%                    if IIIIIII==length(C10)
%                        GG=GG+1;
%                    end
%                          end
%                      end
%                  end
%             end
%         end
%     end
% end        
% 
% 
% save('D:\kepsilon\LakeBiel\Love\Run1\Files_FinalRun\Output\Result\Sequence\Run_Files\CalibrationMatrix.mat',...
%     'CalMatrix','K');
%  load('D:\kepsilon\LakeBiel\Love\Run1\Files_FinalRun\Output\Result\Sequence\Run_Files\CalibrationMatrix.mat');


% % Best fit 2 for model with old abs file
% p1 = 1.2;
% p2 = 1.08;
% %Winter Alpha
% aW = 0.00302;
% % summer Alpha
% % aSS = [0.004 0.005 0.006 0.007 0.008 0.009 0.010 0.011 0.012];
% %Best Fit
% aS=0.006; %Old abs file
% % aS=0.011;% New abs file

% % try different values for q as in Goudsmith et al. 2002 Table 1 (q_NN in model)
% % q_val_best=0.99;
% % q_val = [0.77 0.88 0.99 1.02];%other lake values
% q_val = [0.8 0.9 1 1.1 1.2 1.3 1.4 1.6 1.7];
% % to se if we get more surcafe turbulance compaired to fitted value q=0.99

                



%% Startdate
%%OBS use the right initial condition file

% months considerd to be winter season for seson Alfa
winter_months=[01 02 03 10 11 12];


% %Calibration
%Abs file ABS_use_Daily
% dateini= datenum('28-Feb-1994')-datenum('01-Jan-1904');
% % stop date
% dateend = datenum('01-Jan-2005')-datenum('01-Jan-1904');

%    Calibration and Validation run
%  INITIAL_COND='INITCON_19940228.dat';
%  dateini = datenum('28-Feb-1994')-datenum('01-Jan-1904');
%  dateend = datenum('01-Jun-2014')-datenum('01-Jan-1904');% river data ends

%Cosmo2 forcing and added Absorption throughe Kompensation depth to compair
%INITCON_2008_02_20.dat SimForceCOSMO2_2008_2013.dat3 ABS_KompDeph_2008_2013.dat 
% dateini = datenum('20-Feb-2008')-datenum('01-Jan-1904');
% dateend = datenum('31-Dec-2013')-datenum('01-Jan-1904');

%       Muhleberg heat analys and 
%        INITIAL_COND='INITCON_2008_12_10.dat';                                                   %INITCON_2008_02_20.dat  cosmo 2 forcing
% dateini= datenum('10-Dec-2008')-datenum('01-Jan-1904');%For when we have data for al rivers (Suze limited)
% dateend = datenum('01-Apr-2014')-datenum('01-Jan-1904');%For when we have data for al rivers (Zhil limited)


%       Maxumal river data cover over climate change period
%        INITIAL_COND=%'INITCON_2004_02_18.dat';
%remember in data analyis remove the first year as spin up
dateini= datenum('01-Mar-2004')-datenum('01-Jan-1904');%For when we have data for al rivers (Suze limited) 
dateend = datenum('31-Dec-2009')-datenum('01-Jan-1904');%For when we have data for al rivers (Zhil limited)




%% Forcing files used  

%        
% COSMO2 data impact chek (for publication)
%        INITIAL_COND= 'INITCON_2008_02_20.dat';
%        Absorption_File='ABS_use_Daily.dat';AbsDepth=-20;%Max depth given in ABS file    
% %       COSMO2 forcing
%         Forcing_File='SimForceCOSMO2_2008_2014.dat3'; 
% %      %  Ordenary forcing file
% %         Forcing_File='SimForceStationComb_1994_2014.dat3';%Met stations combined as in Forcing_data_construct.mat

        
      %Files for climate run  
     INITIAL_COND= 'INITCON_2004_02_18.dat';
     Forcing_File='SimForceStationComb_1994_2014.dat3';%Met stations combined as in Forcing_data_construct.mat
       Absorption_File='ABS_use_Daily.dat';AbsDepth=-20;%Max depth given in ABS file    

%       INITIAL_COND='INITCON_19940228.dat';
%        INITIAL_COND='INITCON_2008_12_10.dat';  %'INITCON_1994_02_28.dat';%                                                     %INITCON_2008_02_20.dat  cosmo 2 forcing
%         Forcing_File='SimForceStationComb_1994_2014.dat3';%Met stations combined as in Forcing_data_construct.mat
%        Forcing_File='Cosmo2Force_1994_2014.dat3';                                                %Forcing_New.dat3           
%        Absorption_File='ABS_KompDeph_2008_2013.dat';AbsDepth=-20;%Max depth given in ABS file                
%        Absorption_File='ABS_use_Daily.dat';AbsDepth=-20;%Max depth given in ABS file    
%        Absorption_File='ABSSecci.dat';AbsDepth=-70;%Max depth given in ABS file ABS_Secci_Construction.mat   
       load('D:\kepsilon\LakeBiel\BielData\River\River_Comb_Data_1980_2014');%River Data
       lake=load('D:\kepsilon\LakeBiel\BielData\Lake\Bieler_See_volume_area');
       LakeArea='D:\kepsilon\LakeBiel\Love\Run1\Files_FinalRun\Inputfiles\MORPHBiel.dat';
       


  
 %% the model part  
   run = 1;%parameter to give each run a name (increasing with each run)
  for KK=1:19% forloop to run only once (ordinary condition), or with changed forcing

      
   kk=1;
%    for kk=1:length(K)% forloop to calibrate K (Inactivate kk above)
  
     % Absobtion file construction 
        file_out='D:\kepsilon\LakeBiel\Love\Run1\Files_FinalRun\Inputfiles\ABSSecci.dat';
        ABS_Secci_Construction(K(kk),file_out,dateini,dateend)
   
%        for II=1:length(CalMatrix)%forloop to repeat the model Ie run with calibration matrix (except K)
        fprintf('RUN NUMBER %d\n',run)


        % diferent run parameter setups

        %Calibration of model with new ALPHA file and combined forcing file
%      p1 = CalMatrix(II,1);
%      p2 = CalMatrix(II,2);
%      q =  CalMatrix(II,3);
%      aS =  CalMatrix(II,4); 
%      aW = CalMatrix(II,5);
%      CD = CalMatrix(II,6);
%      C10 = CalMatrix(II,7);


        %% Parameters
        
        %used for the matlabe iterations, no impact on actual model
        k = 1;
        cont=1;
        Stop_Next=0;

        %How long the model should run (in days) befor new initial conditions are made
        Sequence=SequenceVect;%To get a compleat day 
        t1=dateini;
        t2=Sequence+dateini;

        %sumer and winter Alpa (OBS model must be run with batch smaler than
        %180 days (if start date is given in at Oct-01 or Apr-01)
        %winter
        if  ismember(str2num(datestr(t1+datenum('01-Jan-1904'),'mm')),winter_months)==1
            alpha=aW;
        else %summer
            alpha=aS;
        end

        %Parameters to store compleat model results
        Result={};
        Time=[];
        Temp=[];
        U=[];
        V=[];
        
        Intrusion_depth=[];
        flow_hagneck=[];
        flow_zhil=[];
        flow_suze=[];
        flow_brugg=[];
        
        temp_hag=[];
        temp_zhil=[];
        temp_suze=[];
        temp_brugg=[];
        River_time=[];
        
        
        TKE_k=[];
        diss=[];

%% make the M file
%% adjust for climat change 
if Climat==1
% load river Temp, Flow and Air temp
    load D:\kepsilon\LakeBiel\Love\Run1\Files_FinalRun\Inputfiles\ClimatChange\RiverTemp_Air2_With_CCHydro.mat
    
% Change of Air temp in Flow_Forcing file (ClimatAirTempAdjust.m)
    

% Change River Temperature and flow to match CCHydro climat forcing
  %KK Number of runs
  %Filt to make 24 H samples
  LL=length(RiverTempClimate.Suze.Y2021.Date);
  Filt=zeros(24,LL)*NaN;
  Filt(1,:)=RiverTempClimate.Suze.Y2021.Date';Filt=naninterp(reshape(Filt,1,24*LL));
  River.Suze.Temp_time=Filt;
  
  LL=length(RiverTempClimate.Hagneck.Y2021.Date);
  Filt=zeros(24,LL)*NaN;
  Filt(1,:)=RiverTempClimate.Hagneck.Y2021.Date';Filt=naninterp(reshape(Filt,1,24*LL));
  River.Hagneck.Temp_time=Filt;
  
  LL=length(RiverTempClimate.Zihlkanal.Y2021.Date);
  Filt=zeros(24,LL)*NaN;
  Filt(1,:)=RiverTempClimate.Zihlkanal.Y2021.Date';Filt=naninterp(reshape(Filt,1,24*LL));
  River.Zhil.Temp_time=Filt;
  
  LL=length(RiverTempClimate.Aegerten.Y2021.Date); 
  Filt=zeros(24,LL)*NaN;
  Filt(1,:)=RiverTempClimate.Aegerten.Y2021.Date';Filt=naninterp(reshape(Filt,1,24*LL));
  River.Brugg.Temp_time=Filt;
 
  River.Suze.Flow_time=River.Suze.Temp_time;
  River.Hagneck.Flow_time=River.Hagneck.Temp_time;
  River.Zhil.Flow_time=River.Zhil.Temp_time;
  River.Brugg.Flow_time=River.Brugg.Temp_time;
  
  
  if KK==1 % controll run as forcing
  LL=length(RiverTempClimate.Suze.Control.Temp');
  Filt=zeros(24,LL)*NaN;
  Filt(1,:)=RiverTempClimate.Suze.Control.Temp';Filt=naninterp(reshape(Filt,1,24*LL));
  River.Suze.Temp=Filt;
  Filt=zeros(24,LL)*NaN;
  Filt(1,:)=RiverTempClimate.Suze.Control.Flow';Filt=naninterp(reshape(Filt,1,24*LL));
  River.Suze.Flow=Filt;
  
  LL=length(RiverTempClimate.Hagneck.Control.Temp');
  Filt=zeros(24,LL)*NaN;
  Filt(1,:)=RiverTempClimate.Hagneck.Control.Temp';Filt=naninterp(reshape(Filt,1,24*LL));
  River.Hagneck.Temp=Filt;
  Filt=zeros(24,LL)*NaN;
  Filt(1,:)=RiverTempClimate.Hagneck.Control.Flow';Filt=naninterp(reshape(Filt,1,24*LL));
  River.Hagneck.Flow=Filt;
  
  LL=length(RiverTempClimate.Zihlkanal.Control.Temp');
  Filt=zeros(24,LL)*NaN;
  Filt(1,:)=RiverTempClimate.Zihlkanal.Control.Temp';Filt=naninterp(reshape(Filt,1,24*LL));
  River.Zhil.Temp=Filt;
  Filt=zeros(24,LL)*NaN;
  Filt(1,:)=RiverTempClimate.Zihlkanal.Control.Flow';Filt=naninterp(reshape(Filt,1,24*LL));
  River.Zhil.Flow=Filt;

  LL=length(RiverTempClimate.Aegerten.Control.Temp');
  Filt=zeros(24,LL)*NaN;
  Filt(1,:)=RiverTempClimate.Aegerten.Control.Temp';Filt=naninterp(reshape(Filt,1,24*LL));
  River.Brugg.Temp=Filt;
  Filt=zeros(24,LL)*NaN;
  Filt(1,:)=RiverTempClimate.Aegerten.Control.Flow';Filt=naninterp(reshape(Filt,1,24*LL));
  River.Brugg.Flow=Filt;           
     
  elseif KK>1 & KK<=10%2021
      
  LL=length(RiverTempClimate.Suze.Y2021.RiverTemp{KKKK,JJJJ});
  Filt=zeros(24,LL)*NaN;
  Filt(1,:)=RiverTempClimate.Suze.Y2021.RiverTemp{KKKK,JJJJ};Filt=naninterp(reshape(Filt,1,24*LL));
  River.Suze.Temp=Filt;
  Filt=zeros(24,LL)*NaN;
  Filt(1,:)=RiverTempClimate.Suze.Y2021.Flow{KKKK,JJJJ};Filt=naninterp(reshape(Filt,1,24*LL));
  River.Suze.Flow=Filt;
      
  LL=length(RiverTempClimate.Hagneck.Y2021.RiverTemp{KKKK,JJJJ});
  Filt=zeros(24,LL)*NaN;
  Filt(1,:)=RiverTempClimate.Hagneck.Y2021.RiverTemp{KKKK,JJJJ};Filt=naninterp(reshape(Filt,1,24*LL));
  River.Hagneck.Temp=Filt;
  Filt=zeros(24,LL)*NaN;
  Filt(1,:)=RiverTempClimate.Hagneck.Y2021.Flow{KKKK,JJJJ};Filt=naninterp(reshape(Filt,1,24*LL));
  River.Hagneck.Flow=Filt;

  LL=length(RiverTempClimate.Zihlkanal.Y2021.RiverTemp{KKKK,JJJJ});
  Filt=zeros(24,LL)*NaN;
  Filt(1,:)=RiverTempClimate.Zihlkanal.Y2021.RiverTemp{KKKK,JJJJ};Filt=naninterp(reshape(Filt,1,24*LL));
  River.Zhil.Temp=Filt;
    % No CCHydro predictions for Zihlkanal
%   Filt=zeros(24,LL)*NaN;
%   Filt(1,:)=RiverTempClimate.Zihlkanal.Y2021.Flow{KKKK,JJJJ};Filt=naninterp(reshape(Filt,1,24*LL));
%   River.Zhil.Flow=Filt;  
     
  
  LL=length(RiverTempClimate.Aegerten.Y2021.RiverTemp{KKKK,JJJJ});
  Filt=zeros(24,LL)*NaN;
  Filt(1,:)=RiverTempClimate.Aegerten.Y2021.RiverTemp{KKKK,JJJJ};Filt=naninterp(reshape(Filt,1,24*LL));
  River.Brugg.Temp=Filt;
  Filt=zeros(24,LL)*NaN;
  Filt(1,:)=RiverTempClimate.Aegerten.Y2021.Flow{KKKK,JJJJ};Filt=naninterp(reshape(Filt,1,24*LL));
  River.Brugg.Flow=Filt;  
     
  
               KKKK=KKKK+1;
        if KKKK==4
            KKKK=1;
            JJJJ=JJJJ+1;
        end
        if JJJJ==4
            KKKK=1;
            JJJJ=1;
        end
  
  
  elseif KK>10  %2070
  LL=length(RiverTempClimate.Suze.Y2070.RiverTemp{KKKK,JJJJ});
  Filt=zeros(24,LL)*NaN;
  Filt(1,:)=RiverTempClimate.Suze.Y2070.RiverTemp{KKKK,JJJJ};Filt=naninterp(reshape(Filt,1,24*LL));
  River.Suze.Temp=Filt;
  Filt=zeros(24,LL)*NaN;
  Filt(1,:)=RiverTempClimate.Suze.Y2070.Flow{KKKK,JJJJ};Filt=naninterp(reshape(Filt,1,24*LL));
  River.Suze.Flow=Filt;
      
  LL=length(RiverTempClimate.Hagneck.Y2070.RiverTemp{KKKK,JJJJ});
  Filt=zeros(24,LL)*NaN;
  Filt(1,:)=RiverTempClimate.Hagneck.Y2070.RiverTemp{KKKK,JJJJ};Filt=naninterp(reshape(Filt,1,24*LL));
  River.Hagneck.Temp=Filt;
  Filt=zeros(24,LL)*NaN;
  Filt(1,:)=RiverTempClimate.Hagneck.Y2070.Flow{KKKK,JJJJ};Filt=naninterp(reshape(Filt,1,24*LL));
  River.Hagneck.Flow=Filt;

  LL=length(RiverTempClimate.Zihlkanal.Y2070.RiverTemp{KKKK,JJJJ});
  Filt=zeros(24,LL)*NaN;
  Filt(1,:)=RiverTempClimate.Zihlkanal.Y2070.RiverTemp{KKKK,JJJJ};Filt=naninterp(reshape(Filt,1,24*LL));
  River.Zhil.Temp=Filt;
    % No CCHydro predictions for Zihlkanal
%   Filt=zeros(24,LL)*NaN;
%   Filt(1,:)=RiverTempClimate.Zihlkanal.Y2070.Flow{KKKK,JJJJ};Filt=naninterp(reshape(Filt,1,24*LL));
%   River.Zhil.Flow=Filt;  
     
  
  LL=length(RiverTempClimate.Aegerten.Y2070.RiverTemp{KKKK,JJJJ});
  Filt=zeros(24,LL)*NaN;
  Filt(1,:)=RiverTempClimate.Aegerten.Y2070.RiverTemp{KKKK,JJJJ};Filt=naninterp(reshape(Filt,1,24*LL));
  River.Brugg.Temp=Filt;
  Filt=zeros(24,LL)*NaN;
  Filt(1,:)=RiverTempClimate.Aegerten.Y2070.Flow{KKKK,JJJJ};Filt=naninterp(reshape(Filt,1,24*LL));
  River.Brugg.Flow=Filt;      
 
              KKKK=KKKK+1;
        if KKKK==4
            KKKK=1;
            JJJJ=JJJJ+1;
        end
        if JJJJ==4
            KKKK=1;
            JJJJ=1;
        end
  
  
  end
         

  %If removal of Mühleberg is activ during climate change
  if Muhl==0 % folowing method in Muhleberg_heat_removed_from_Aare_Filemaking.mat
    load('D:\kepsilon\LakeBiel\Love\Run1\Files_FinalRun\Inputfiles\ClimatChange\MuhlbergRemovalParameters.mat')
    addpath(genpath('D:\kepsilon\LakeBiel\Love\Mfiles_Love\seawater'))
    
    
    % Water mass passing a certain point during 1 sec          
    WaterMass=sw_dens_lake(River.Hagneck.Temp).*River.Hagneck.Flow;
    % Energy content using temp in [K]
    river_energy=Cp*(River.Hagneck.Temp+273.15).*WaterMass;%[J]

    
    % find places
    [~,places]=histc(River.Hagneck.Temp_time,Muhleberg_energy_date);
          River_place=find(places>0);
          Heat_place=places(River_place);
    
    % New temp
    River.Hagneck.Temp=((river_energy(River_place)-Muhleberg_energy(Heat_place)')./(WaterMass(River_place)*Cp))-273.15;%[°C]
    River.Hagneck.Temp_time=River.Hagneck.Temp_time(River_place);
  end
end




%% river dates
%create cell structure with all dates for river Hagneck
date.Hagneck.Temp = River.Hagneck.Temp_time-datenum('01-Jan-1904');
date.Hagneck.Flow = River.Hagneck.Flow_time-datenum('01-Jan-1904');
date.Hagneck.Temp_NoMuhl = River.Hagneck.Temp_NoMuhl_time-datenum('01-Jan-1904');

date.Suze.Temp = River.Suze.Temp_time-datenum('01-Jan-1904');
date.Suze.Flow = River.Suze.Flow_time-datenum('01-Jan-1904');

date.Zhil.Temp = River.Zhil.Temp_time-datenum('01-Jan-1904');
date.Zhil.Flow = River.Zhil.Flow_time-datenum('01-Jan-1904');



            %% make the river intrusion file
            file=['D:\kepsilon\LakeBiel\Love\Run1\Files_FinalRun\Inputfiles\',INITIAL_COND];
            % Creates dayly mean for the river intrusion
            if zero==1
            [Mean,DateVect]=Mean_calc(file,River,date,t1,t2,time_res,removal,Muhl); %[m3/s]
            %         removing extra date
            Mean.date=Mean.date(1:end-1);
            Mean.Hagneck.mean_temp=Mean.Hagneck.mean_temp(1:end-1);
            Mean.Hagneck.mean_flow=Mean.Hagneck.mean_flow(1:end-1);
            Mean.Zhil.mean_temp=Mean.Zhil.mean_temp(1:end-1);
            Mean.Zhil.mean_flow=Mean.Zhil.mean_flow(1:end-1);
            Mean.Suze.mean_temp=Mean.Suze.mean_temp(1:end-1);
            Mean.Suze.mean_flow=Mean.Suze.mean_flow(1:end-1);
            DateVect=DateVect(1:end-1);
            else
                Mean=[];DateVect=[];
            end
            
            
            % Calculates intrusion depth and new ABS and Forcing files
            [Inflow]=intrusion(file,t1,t2,time_res,zero,Mean,lake,DateVect,Forcing_File,Absorption_File,AbsDepth,LakeArea,dt,Climat,KK); 

            First=1;

    while cont==1 %run model with batch setup untill end date is reached
            %% Generate kepsilon.par file

            fid = fopen('D:\kepsilon\LakeBiel\Love\Run1\bin\kepsilon.par','w');
            fprintf(fid,'%s\n','*** Input Files ************************************');
            if k==1
                fprintf(fid,'%s\n',['D:\kepsilon\LakeBiel\Love\Run1\Files_FinalRun\Inputfiles\',INITIAL_COND]);
            else
                fprintf(fid,'%s\n','D:\kepsilon\LakeBiel\Love\Run1\Files_FinalRun\Inputfiles\Changin_Initial_Cond.dat');
            end
            fprintf(fid,'%s\n','D:\kepsilon\LakeBiel\Love\Run1\Files_FinalRun\Inputfiles\Grid_Biel.dat');
            fprintf(fid,'%s\n', LakeArea);
            fprintf(fid,'%s\n','D:\kepsilon\LakeBiel\Love\Run1\Files_FinalRun\Inputfiles\Forcing_use.dat3');
            fprintf(fid,'%s\n','D:\kepsilon\LakeBiel\Love\Run1\Files_FinalRun\Inputfiles\ABS_use.dat');
            fprintf(fid,'%s\n','D:\kepsilon\LakeBiel\Love\Run1\Files_FinalRun\Output\Result\Sequence\Results_Batch.dat');
            fprintf(fid,'%s\n','D:\kepsilon\LakeBiel\Love\Run1\Files_FinalRun\Output\Output_depth_Biel.dat');
            fprintf(fid,'%s\n','D:\kepsilon\LakeBiel\Love\Run1\Files_FinalRun\Output\Output_time_Biel.dat');
            if zero ==0
            fprintf(fid,'%s\n','D:\kepsilon\LakeBiel\Love\Run1\Files_FinalRun\Inputfiles\Qinp_Biel_Zero.dat');
            fprintf(fid,'%s\n','D:\kepsilon\LakeBiel\Love\Run1\Files_FinalRun\Inputfiles\Qout_Biel_Zero.dat');
            fprintf(fid,'%s\n','D:\kepsilon\LakeBiel\Love\Run1\Files_FinalRun\Inputfiles\Tinp_Biel_Zero.dat');
            fprintf(fid,'%s\n','D:\kepsilon\LakeBiel\Love\Run1\Files_FinalRun\Inputfiles\Sinp_Biel_Zero.dat');
            else
            fprintf(fid,'%s\n','D:\kepsilon\LakeBiel\Love\Run1\Files_FinalRun\Inputfiles\Qinp_Biel.dat');
            fprintf(fid,'%s\n','D:\kepsilon\LakeBiel\Love\Run1\Files_FinalRun\Inputfiles\Qout_Biel.dat');
            fprintf(fid,'%s\n','D:\kepsilon\LakeBiel\Love\Run1\Files_FinalRun\Inputfiles\Tinp_Biel.dat');
            fprintf(fid,'%s\n','D:\kepsilon\LakeBiel\Love\Run1\Files_FinalRun\Inputfiles\Sinp_Biel.dat');    
            end
            fprintf(fid,'%s\n','*** Model set up: time steps and grid resolution ***');
            fprintf(fid,'%s\n',sprintf('%.0f       dt (Time step [s])',dt));
            fprintf(fid,'%s\n',sprintf('%.0f       t_start [d]',t1));
            fprintf(fid,'%s\n',sprintf('%.0f       t_end [d]',t2));
            fprintf(fid,'%s\n','*** Select Model and Boundary conditions ***********');
            fprintf(fid,'%s\n','1           Model (1: kepsilon, 2: MY)');
            fprintf(fid,'%s\n','1           Stability function (1: const, 2: quasi-equilibrium)');
            fprintf(fid,'%s\n','1           Flux condition (1: no-flux, 0: Dirichlet condition)');
            fprintf(fid,'%s\n','3           Forcing (1: SST+WIND+SolRad, 2: Meteo, 3: Meteo+Cloud, 4: Heat fluxes)');
            fprintf(fid,'%s\n','0           Pressure gradients (0: no pressure gradients)');
            fprintf(fid,'%s\n','0           Display simulation (1: if data are saved)');
            fprintf(fid,'%s\n','0           Display diagnose (1: all settings, 2: all setting and forcing each iteration)');
            fprintf(fid,'%s\n','*** Seiche Parameters ******************************');
            fprintf(fid,'%s\n',sprintf('%.5f         alpha_Seiche',alpha));
            fprintf(fid,'%s\n',sprintf('%.2f            q_NN',q));%'0.99            q_NN');
            fprintf(fid,'%s\n',sprintf('%.3f           CD',CD));
            fprintf(fid,'%s\n','*** Other Parameters *******************************');
            fprintf(fid,'%s\n','47.2       Lat (latitude to calculate Coriolis parameter)');
            fprintf(fid,'%s\n',sprintf('%.4f     C10',C10));
            fprintf(fid,'%s\n','0.080      fgeo [W/m2] Geothermal heat flux');
            fprintf(fid,'%s\n','1e-9       k_min');
            fprintf(fid,'%s\n','951        air_press [mbar]');
            fprintf(fid,'%s\n','2          norm_ctr Seiche normalization (1: max N^2; 2: integral)');
            fprintf(fid,'%s\n',sprintf('%f   p1',p1));
            fprintf(fid,'%s\n',sprintf('%f   p2',p2));
            fprintf(fid,'%s\n','10         igoal');
            fprintf(fid,'%s\n','1          ModSal');
            fclose(fid);
            
                % Geothermal heat flux = 80 mW/m2 
                % 1979
                % 978-3-642-95359-0
                % Terrestrial Heat Flow in Europe
                % ?ermák, Vladimír
                % Rybach, Ladislaus
                % 10.1007/978-3-642-95357-6_1
                % Heat Flow Map of Europe
                % http://dx.doi.org/10.1007/978-3-642-95357-6_1
                % Springer Berlin Heidelberg
                % 3-40
         
                
            % to stop model to chek river intrusion 
%             if t1>datenum('15-Nov-2005')-datenum('01-Jan-1904');
%                 disp(1)
%             end       
%             
            %% Run the model
            fprintf(' Running model from %s to %s...: ',datestr(t1+datenum('01-Jan-1904')),datestr(t2+datenum('01-Jan-1904')));tic;
            !start /I /min /low /wait kepsmodel.exe
            toc;fprintf('Done.\n');

            %%
            %Write initial conditions for following period
            [z,ini] =  get_final_Biel('D:\kepsilon\LakeBiel\Love\Run1\Files_FinalRun\Output\Result\Sequence\Results_Batch.dat');
            fini = fopen('D:\kepsilon\LakeBiel\Love\Run1\Files_FinalRun\Inputfiles\Changin_Initial_Cond.dat','w');
            fprintf(fini,'z [m]\tV [m/s]\tU [m/s]\tT [°C]\tS [‰]\tk [J/kg]\teps [W/kg]');
            for i=2:length(z)%do not use cell above normal water level
                fprintf(fini,'\n%.0f',z(i));
                fprintf(fini,'\t%.3e',ini(i,:));                     
            end
            fprintf(fini,'\n'); 
            fclose(fini);

            %
                    %% Store compleat model results 
       if First==1 %for the first batch sequence
            [z,~,Model] = get_model_result('D:\kepsilon\LakeBiel\Love\Run1\Files_FinalRun\Output\Result\Sequence\Results_Batch.dat');
            Time= [Time (Model.Results.Time+round(datenum('01-Jan-1904')))];%datenumber
            Temp= [Temp Model.Results.T(2:end,:)];%[C°]
            U= [U Model.Results.U(2:end,:)];%[m/s]
            V= [V Model.Results.V(2:end,:)];%[m/s]
            TKE_k= [TKE_k Model.Results.k(2:end,:)];%[C°]
            diss= [diss Model.Results.diss(2:end,:)];%[C°]

            Result{k} = Model;%all results

           if zero==1
           ha=Inflow.Hagneck_Depth(1:end-1)';
           su=Inflow.Suze_Depth(1:end-1)';
           zh=Inflow.Zhil_Depth(1:end-1)';
           Intrusion_depth= [Intrusion_depth; [ha su zh]];
           flow_hagneck=[flow_hagneck Mean.Hagneck.mean_flow(1:end-1)];
           flow_zhil=[flow_zhil Mean.Zhil.mean_flow(1:end-1)];
           flow_suze=[flow_suze Mean.Suze.mean_flow(1:end-1)];
           flow_brugg=[flow_brugg (Mean.Hagneck.mean_flow(1:end-1)+Mean.Zhil.mean_flow(1:end-1)+Mean.Suze.mean_flow(1:end-1))];
           temp_hag=[temp_hag Mean.Hagneck.mean_temp(1:end-1)];
           temp_zhil=[temp_zhil Mean.Zhil.mean_temp(1:end-1)];
           temp_suze=[temp_suze Mean.Suze.mean_temp(1:end-1)];
           River_time=[River_time Mean.date(1:end-1)+round(datenum('01-Jan-1904'))];
           temp_brugg=[temp_brugg Model.Results.T(2,:)];

           else
           Intrusion_depth=[];
           flow_hagneck=[];
           flow_zhil=[];
           flow_suze=[];
           flow_brugg=[];
           temp_hag=[];
           temp_zhil=[];
           temp_suze=[];
           temp_brugg=[];
           River_time=[];
           end       
       else %For all other batch runs (to not store the same date twice due to model having to rerun first date)
            [z,~,Model] = get_model_result('D:\kepsilon\LakeBiel\Love\Run1\Files_FinalRun\Output\Result\Sequence\Results_Batch.dat');
            BatchTime=((Model.Results.Time(1:end-1))+round(datenum('01-Jan-1904')));
            Place=ismember(round(BatchTime*24),round(Time*24));
            PlaceUse=find(Place==0);
            Time= [Time BatchTime(PlaceUse)];%datenumber
            Temp= [Temp Model.Results.T(2:end,PlaceUse)];%[C°]
            U= [U Model.Results.U(2:end,PlaceUse)];%[m/s]
            V= [V Model.Results.V(2:end,PlaceUse)];%[m/s]
            TKE_k= [TKE_k Model.Results.k(2:end,PlaceUse)];%[C°]
            diss= [diss Model.Results.diss(2:end,PlaceUse)];%[C°]

            Result{k} = Model;%all results

           if zero==1 %inflows are recalculated for first day and extends one day exrtra, the model however stopps on t2 while the intrusion for this day have to be calculated.
           
           BatchRiverTime=Mean.date(1:end-2)+round(datenum('01-Jan-1904'));
           Place=ismember(round(BatchRiverTime*24),round(River_time*24));
           PlaceUse=find(Place==0);
%          do not go outside model time span
           if BatchRiverTime(end)>dateend+round(datenum('01-Jan-1904'))
               
              PPP1=find(BatchRiverTime(PlaceUse)<dateend+round(datenum('01-Jan-1904')));
              PlaceUse=PlaceUse(PPP1); 
           end
           
           River_time=[River_time BatchRiverTime(PlaceUse)];

           ha=Inflow.Hagneck_Depth(PlaceUse)';
           su=Inflow.Suze_Depth(PlaceUse)';
           zh=Inflow.Zhil_Depth(PlaceUse)';
           Intrusion_depth= [Intrusion_depth; [ha su zh]];
           flow_hagneck=[flow_hagneck Mean.Hagneck.mean_flow(PlaceUse)];
           flow_zhil=[flow_zhil Mean.Zhil.mean_flow(PlaceUse)];
           flow_suze=[flow_suze Mean.Suze.mean_flow(PlaceUse)];
           flow_brugg=[flow_brugg (Mean.Hagneck.mean_flow(PlaceUse)+Mean.Zhil.mean_flow(PlaceUse)+Mean.Suze.mean_flow(PlaceUse))];
           temp_hag=[temp_hag Mean.Hagneck.mean_temp(PlaceUse)];
           temp_zhil=[temp_zhil Mean.Zhil.mean_temp(PlaceUse)];
           temp_suze=[temp_suze Mean.Suze.mean_temp(PlaceUse)];
           temp_brugg=[temp_brugg Model.Results.T(2,PlaceUse)];
           else
           Intrusion_depth=[];
           flow_hagneck=[];
           flow_zhil=[];
           flow_suze=[];
           flow_brugg=[];
           temp_hag=[];
           temp_zhil=[];
           temp_suze=[];
           temp_brugg=[];
           River_time=[];           
           end
       end
       
%             % to stop model to chek river intrusion 
%             if t1>datenum('15-Nov-2005')-datenum('01-Jan-1904');
%                 disp(1)
%             end       

    %% for next sequence
          k = k+1;    
          t1=round(t2);
          t2=t1+Sequence;
           t1=t1-1;%to start closer to the initial condition
           
           
           
          %If next batch run is outside stop date    
          %run one more time then stop
          if Stop_Next==1
              k=k-1;
              break
          end
          if t2>dateend
              t2=dateend;
              Stop_Next=1;
          end          


           %to switch from first batch
           First=0;


            %make the river intrusion file
            file='D:\kepsilon\LakeBiel\Love\Run1\Files_FinalRun\Inputfiles\Changin_Initial_Cond.dat';
            
            % Creates dayly mean for the river intrusion
            if zero==1
            [Mean,DateVect]=Mean_calc(file,River,date,t1,t2,time_res,removal,Muhl);
            else
                Mean=[];DateVect=[];
            end            
            % Calculates intrusion depth       
            [Inflow]=intrusion(file,t1,t2,time_res,zero,Mean,lake,DateVect,Forcing_File,Absorption_File,AbsDepth,LakeArea,dt,Climat,KK); 


        if ismember(str2num(datestr(t1+datenum('01-Jan-1904'),'mm')),winter_months)==1
            alpha=aW;
        else %summer
            alpha=aS;
        end




    end

        file=['D:\kepsilon\LakeBiel\Love\Run1\Files_FinalRun\Output\Result\Sequence\Run_Files\' sprintf('Run_%d.mat',run)];
        save(file,'Time','Temp','U','V','TKE_k','diss','Intrusion_depth','flow_hagneck','flow_zhil','flow_suze','flow_brugg',...
            'temp_hag','temp_zhil','temp_suze','temp_brugg','River_time','Result');
        run = run+1;


%        end
 end



