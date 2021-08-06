% POST-PROCESS QUALITY CONTROL BGC ARGO RADIOMETRY 
% CALCULATE SENSOR-TEMPERATURE DEPENDENT CORRECTION
% NIGHT METHOD
% UPDATED 8.2.2021
% author: Terence O'Brien
% terence.obrien@maine.edu

% 5 SECTIONS
% CAN RUN WHOLE SCRIPT OR SECTION-BY-SECTION (TO CHECK & RECORD STATS)
% Output is "G(wavelength)cors" array of x0 & x1 per float (see section 5)

%% 1.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%%% 1. LOAD DATA, ASSESS NIGHT PROFILES FOLLOWING QUALITY CONTROL PROCEDURES FROM ORGANELLI 2016
%%%% SECTION 1 ADAPTED FROM EMMANUEL BOSS

clear all
close all
%Irradiance data in W m^-2 nm^-1
% PAR in micromol quanta m^-2 s^-1

load data490.mat %- use for 490 & 412
%load data380.mat
%data490 = data380; % for ease- if using data 380 don't
%                           need to change variable name throughout script

% ^need to change variables throughout script when changing from 380 to 490
% check before replacing
%% 

A=data490.Properties.VariableNames;
for i=1:length(A)
    I(1,i)=isempty(find(string(A{1,i})=='dt'));
    I(2,i)=isempty(find(string(A{1,i})=='lat'));
    I(3,i)=isempty(find(string(A{1,i})=='lon'));
    I(4,i)=isempty(find(string(A{1,i})=='PRES'));
    I(5,i)=isempty(find(string(A{1,i})=='TEMP'));
    I(6,i)=isempty(find(string(A{1,i})=='PSAL'));
    I(7,i)=isempty(find(string(A{1,i})=='FLUORESCENCE_CHLA'));
    I(8,i)=isempty(find(string(A{1,i})=='CHLA'));
    I(9,i)=isempty(find(string(A{1,i})=='FLUORESCENCE_CDOM'));
    I(10,i)=isempty(find(string(A{1,i})=='CDOM'));
    I(11,i)=isempty(find(string(A{1,i})=='CP'));
    I(12,i)=isempty(find(string(A{1,i})=='CP660'));
    I(13,i)=isempty(find(string(A{1,i})=='BBP'));
    I(14,i)=isempty(find(string(A{1,i})=='BBP470'));
    I(15,i)=isempty(find(string(A{1,i})=='BBP532'));
    I(16,i)=isempty(find(string(A{1,i})=='BBP700'));
    I(17,i)=isempty(find(string(A{1,i})=='DOWNWELLING_PAR'));
    I(18,i)=isempty(find(string(A{1,i})=='DOWN_IRRADIANCE380')); %**
    I(19,i)=isempty(find(string(A{1,i})=='DOWN_IRRADIANCE412')); %**
    I(20,i)=isempty(find(string(A{1,i})=='DOWN_IRRADIANCE443')); 
    I(21,i)=isempty(find(string(A{1,i})=='DOWN_IRRADIANCE490')); %**
    I(22,i)=isempty(find(string(A{1,i})=='WMO')); %float ID
end

%find the position of the variables of interest.
ind=mod(find(I'==0),47);
ind(find(ind==0))=47;


outNight490=zeros(1,8);
counter=0;
y=1;
Flag1=0; %number of profiles.
K=0; %Ed profile is Nans
KK=0; %profile not monotonic in pressure
KKK=0; %lat or long not in real bounds
KKKK=0; %solar azumith angle test
M = 0; % <5 data points in profile

for i=1:length(data490{:,1}) %analysis is performed a profile at a time
    i;
    Flag1;
    % getting the profiles and test that we want to keep them
    Ed490=cell2mat(data490{i,ind(19)})*100;  %miltiply by 100 for comparison with Organelli's units
                            % ^^^^^ !!! index is wavelength specific 
                            
    p=cell2mat(data490{i,ind(4)}); T=cell2mat(data490{i,ind(5)});S=cell2mat(data490{i,ind(6)});ID=data490{i,ind(22)};
    lat=data490{i,ind(2)};lon=data490{i,ind(3)};time=data490{i,ind(1)};
    %plot(Ed490,-p,'*')
    %pause
    
    %if profile is all NaN, remove it
    if isempty(find(isfinite(Ed490)))==1 K=K+1; continue; end 
    %Remove NaN from profile in not all NaN
    J=find(isfinite(Ed490));Ed490=Ed490(J);p=p(J); T=T(J); S=S(J);
    %if profile is not monotonic in pressure, remove it
    if sum(sort(p)-p) KK=KK+1; continue; end 
    %testing that reasonable LAT and LON, if not skip
    II=find(lat<=90 & lat>=-90 & lon<=180 & lon>=-180);
    if isempty(II)  KKK=KKK+1; continue; end  
    %testing that solar elevation is above 15degrees
    DateTimeUTC=time;
    [Az,El] = SolarAzEl(DateTimeUTC,lat,lon,0); %(Koblick 2021)
    %detect only reasonable values
    if El >0   KKKK=KKKK+1; continue; end    % 
   
    % if there are less than 5 data point in a profile, skip it.
    if length(Ed490) < 5 M=M+1; continue; end
     Ed=Ed490;pr=p;
  
         location.latitude=lat;location.longitude=lon;location.altitude = 0;
    t.year=str2num(datestr(time,'yyyy'));t.month=str2num(datestr(time,'mm'));t.day=str2num(datestr(time,'dd'));
    t.hour=str2num(datestr(time,'hh'));t.min=str2num(datestr(time,'MM'));t.sec=0; t.UTC=0;

       L=length(Ed);
       outNight490(counter+1:counter+L,1)=location.latitude;
         outNight490(counter+1:counter+L,2)=location.longitude;
         outNight490(counter+1:counter+L,3)=DateTimeUTC;
         outNight490(counter+1:counter+L,4)=p; 
         outNight490(counter+1:counter+L,5)=T;       
         outNight490(counter+1:counter+L,6)=S;
         outNight490(counter+1:counter+L,7)=Ed;    
         outNight490(counter+1:counter+L,8)=ID;
         counter=counter+L;
       clear p pp pr Ed Ed490 II location time JJ
       if mod(i,5000)==0
           disp(i)
       end
end

%% Check for NaNs in temperature profiles:

a = find(isnan(outNight490(:,5)));
if ~isempty(a)
outNight490(a,:)=NaN;
for ii=1:8
    cc = outNight490(:,ii);
    c = ~isnan(cc);
    outNight490(:,ii) = NaN;
    cc = cc(c);
    for j=1:length(cc)
        outNight490(j,ii) = cc(j);
    end
    
end

outNight490 = outNight490(1:find(isnan(outNight490(:,1)),1),:);

else 
    outNight490(length(outNight490)+1,1:8) = NaN;  % Intentionally setting last line as NaNs
end
%% clear vars 
clearvars -except outNight490 data490

%% 2
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%%% 2. COMPUTE SENSOR TEMPERATURE
%%%% CALLS ON FUNCTION(SENSORTEMP)
%%%% TENV TO TINS


cns=outNight490; 

Tsen3=NaN(length(cns),9);   %Final Output w corrected temps, and everything else.
x=1; %for Tenv index. resets to one with new lt
Tenv=NaN(10,2);  %to extract Tenv from cns per profile
Ts3 = NaN(2,1); %for use inside loops
lt=cns(1,1); %current latitute
plt=cns(1,1); %previous latitude
cnt=1; %do not reset. row counter for Tsen
n=1; %counter for Tsen(x,temp)=Ts(n,1) 

 K=200; %Time constant for function [s]

for i=1:length(cns)
    lt=cns(i,1);
    
    if lt==plt && ~isnan(lt) %same latitude
       Tenv(x,1)=cns(i,5); %temperature
       Tenv(x,2)=cns(i,4); %pressure for time calculation
        x=x+1;
       
         
        plt = lt;  
        
    end
       if ~(lt == plt) || i==length(cns)       % lt=/=plt
                              % we are at a new profile. Convert Tenv to Ts
                            % make sure to save the first iteration of the new profile. 
        temps = Tenv(:,1); 
       
        c = ~isnan(temps);
        temps = temps(c);
        pres = Tenv(:,2); 
        pres2 = pres(c);
        if min(pres2)>250  %Tconversion only on p<=250
               j=i-length(pres2);
            k=i-1;
            for iii=j:k                
                Tsen3(cnt,1)=cns(iii,1);
                Tsen3(cnt,2)=cns(iii,2);
                Tsen3(cnt,3)=cns(iii,3);
                Tsen3(cnt,4)=cns(iii,4);
                Tsen3(cnt,5)=cns(iii,5);     
                Tsen3(cnt,6)=cns(iii,5);
                Tsen3(cnt,7)=cns(iii,6);
                Tsen3(cnt,8)=cns(iii,7);
                Tsen3(cnt,9)=cns(iii,8);
                cnt=cnt+1;
           
                
            end
        else
        cc = pres2<= 250; % consistent Ed readings begin @250db 
        temps = temps(cc);
        pres3 = pres2(cc);
        temps = flipud(temps);  % so T(deepest) is index 1
        time = (max(pres3)-min(pres3))/0.1;
        t = 0:(time/(length(temps)-1)):time; %array of times, same length as temps
       % To=temps(1);         %T(0)
        ct=zeros(1);          %"current temp"
        Ts3=zeros(1,length(t));
        pres3 = flipud(pres3); %to fit w temps
         %now convert Tenv to Ts
         %if max(temps) - min(temps) > 3
             % calculate av dTdz for first 20db
             on = pres3(1);
             w = find((pres2 < on - 20),1);
             if ~isempty(w)
                 dT = temps(1:w);
                 dz = pres3(1:w);
                 fit = polyfit(dz,dT,1);
                 dTdt = fit(1);   
                To = temps(1)+dTdt*20;
             else
                 dT = temps(1:length(temps));
                 dz = pres3(1:length(pres3));
                 fit = polyfit(dz,dT,1);
                 dTdt = fit(1);
                 To = temps(1)+dTdt*(max(pres3)-min(pres3));
             end    
            for ii=1:length(t)      %to calculate Tins 
             ct(ii)=temps(ii);
    
              intg=SensorTemp(ct,t(ii),K);           %function to return integral over current t range
    
              % Ts(t)=T(0)e^(-t/K)+(e^(-t/K)int(0,t')Tenv(t')e^(t'/K)dt)/K
              Ts3(ii)=To*exp(-t(ii)/K)+exp(-t(ii)/K)*intg; %creates Ts array of Sensor temps
            end
            
            %Now fill Tsen with Ts and all other data.            
            Ts3 = fliplr(Ts3); %inverted compared to rest of data
            j=i-length(pres2);
            k=i-1;
            for iii=j:k
                if cns(iii,4)<= 250
                Tsen3(cnt,1)=cns(iii,1);
                Tsen3(cnt,2)=cns(iii,2);
                Tsen3(cnt,3)=cns(iii,3);
                Tsen3(cnt,4)=cns(iii,4);
                Tsen3(cnt,5)=cns(iii,5);     %Ts is a row array
                Tsen3(cnt,6)=Ts3(1,n);
                Tsen3(cnt,7)=cns(iii,6);
                Tsen3(cnt,8)=cns(iii,7);
                Tsen3(cnt,9)=cns(iii,8);
                cnt=cnt+1;
                n=n+1;
                else 
                Tsen3(cnt,1)=cns(iii,1);
                Tsen3(cnt,2)=cns(iii,2);
                Tsen3(cnt,3)=cns(iii,3);
                Tsen3(cnt,4)=cns(iii,4);
                Tsen3(cnt,5)=cns(iii,5);     %Ts not calculated >250db
                Tsen3(cnt,6)=cns(iii,5);
                Tsen3(cnt,7)=cns(iii,6);
                Tsen3(cnt,8)=cns(iii,7);
                Tsen3(cnt,9)=cns(iii,8);
                cnt=cnt+1;
                end
            end
        end
          n=1; %n is reset 
         plt = lt;  
        x=1; %x is reset;
        Tenv=NaN(10,2);
        Ts3=NaN(2,1);
        
        %now take out new Ts since we have lt=/=plt
       Tenv(x,1)=cns(i,5); %temperature
       Tenv(x,2)=cns(i,4);
       x=x+1;
       end  
end

Tsen3 = Tsen3(1:find(isnan(Tsen3(:,1)),1),:);

clearvars -except outNight490 data490 Tsen3



%% 3
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%%% 3. DOWNWELLING IRRADIANCE TEST 
%%%% DETECT AND REMOVE SECTIONS OF PROFILES DETECTING DOWNWELLING IRRADIANCE
%%%% TESTS 150m, 100m, 50m TO SURFACE

ns=Tsen3;         

lt=1;   %latitude
plt=0;  %must start different than lt for first profile
Id = 0;
fail=1;
fail2=0; %second l-l test
fail3=0; %third l-l test
maxFail=0;
ctL150=0; %light detected 150db-0
ctL100=0; %light detected 100db-0
ctL50=0; %light detected 50db-0
ctT=0; %total profiles

ctl=1;
for i=1:length(ns)
    lt=ns(i,1);  %current latitude
    
    if not(lt==plt) %when @new lat/profile
        ctT = ctT+1;
        ctl = 1; %reset ctl for counter at end
        pr=NaN(2,1); %reset
        e=NaN(2,1);  %reset
        
        for l=0:5000 %determine if new profile is good or bad
                    % all profiles < 5000 readings
                    
           if ns(l+i,1)==lt && l+i < length(ns)-1 %if latitude is the same
             %extract pr and e for current profile
             pr(l+1)=ns(l+i,4); %pressure
             e(l+1)=ns(l+i,8); % Ed
             Id = ns(l+i,9);
             
             else  %=@new latitude/profile. run calculation on current profile. 
                %1. Filter downwelling light
             
                minEd = min(e);
                 for r=1:length(e)             %shift data so none are negative for log test
                     e(r)=e(r)-minEd + 0.001;
                 end
               
                    cc= pr<150; %150 m to surface for test
                    pr1=pr(cc);
                    e1=e(cc);
                    p=polyfit(pr1,log(e1),1);
                    %begin test
                    
                    if p(1)<-0.01   %first log-lin test. 150m-surface
                          fail=1;  %Fail - trend detected
                          ns(i,:)=NaN; %remove surfacemost measurement
                           y=pr<150; %
                           ctL150=ctL150+1;
                            for x=1:length(pr)
                              if y(x)==1
                                 pr(x)=NaN;  %remove upper 150 m of data
                                 e(x)=NaN; 
                              end
                            end
                            
                    else %test 100m - surface
                          y = 0; %reset
                          y=pr<100; % to test 100m to surface
                          %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%  
                          pr2=pr(y);
                            e2=e(y);
                             p=polyfit(pr2,log(e2),1);
                                if p(1)<-0.01   %SECOND log-lin test. 100m-0
                                     fail2=1;
                                     ns(i,:)=NaN; %remove surfacemost measurement
                                 ctL100=ctL100+1;
                                for x=1:length(pr)
                                     if y(x)==1
                                     pr(x)=NaN;  %remove upper 100 m of data
                                     e(x)=NaN; 
                                     end
                                end
                                else
                                  fail2=0;
                                  y3=pr<50; %get pressures less than than 50
                               
                                  pr3=pr(y3);
                                  e3=e(y3);
                                   p=polyfit(pr3,log(e3),1);
                                   if p(1)<-0.01   %THIRD log-lin test. 
                                      fail3=1; %light detected
                                     ctL50=ctL50+1;
                                        for xxx=1:length(pr)
                                         if y3(xxx)==1
                                          pr(xxx)=NaN;  %remove upper 50 m of data
                                          e(xxx)=NaN; 
                                         end
                                        end
                                   else
                                               fail3=0;
                                   end
                                end
                       end  %end of light filter poly test
                     
             break 
           end %if ns(l+i,1) == lt
        end % for l=0:1000
        
   else 
       % output amended ns
     %  ns(i,:)=NaN;  
      if isnan(pr(ctl)) %if pressure is NaN, make row NaN  
         ns(i,:)=NaN;
      end
      ctl=ctl+1; %ctl to count pressure. resets when at new profile. 
        
       
      end     %if not(lt==plt)
    
    plt=lt;     
end




%% nans

for i=1:9
    cc = ns(:,i);
    c = ~isnan(cc);
    ns(:,i) = NaN;
    cc = cc(c);
    for j=1:length(cc)
        ns(j,i) = cc(j);
    end
    
end

ns = ns(1:find(isnan(ns(:,1)),1),:);

clearvars -except outNight490 data490 Tsen3 ns



%% 4
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%%% 4. MAX IRRADIANCE FILTER
%%%% 
%%%% REMOVE VALUES OUTSIDE [-0.3 0.3] W m-2 nm-1


cns = ns;

ctH = 0;
ctL = 0;
ctT = 0;
fns = cns;
for i = 1:length(cns)
    if cns(i,8)> 0.03
        fns(i,:) = NaN;
        ctH = ctH +1;
        ctT = ctT +1;
    end
    if cns(i,8) < -0.03
        fns(i,:) = NaN;
        ctL = ctL +1;
        ctT = ctT +1;
    end
end


%% Condense NaNs

for i=1:9
    cc = fns(:,i);
    c = ~isnan(cc);
    fns(:,i) = NaN;
    cc = cc(c);
    for j=1:length(cc)
        fns(j,i) = cc(j);
    end
    
end

fns = fns(1:find(isnan(fns(:,1)),1),:);
clearvars -except outNight490 data490 Tsen3 ns fns


%% 5
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%%% 5. ROBUST FIT AND dE/dT CONSTRAINT
%%%% ROBUST LINEAR FIT ON COMPILED PROFILES PER FLOAT PER WAVELENGTH
%%%% 1.5*IQR CONSTRAINT ACROSS ALL FLOATS dE/dT AT SAME WAVELENGTH

%fns: 1. lat 2. long 3. time 4. pres 5. Tm 6. Ts 7. Sal 8. Ed 9. ID

cns = fns;  
N490_cors=NaN(1,8);
short = 0; %counter 
tt = 1;
i = 1;
Trangefail = 0;
rhofail = 0;
while tt < length(cns)
    Id=cns(tt,9);
    
    % extract temp, ed, pres for Id
    a = find(cns(:,9)==Id);
    t = cns(a,6);
    Ed = cns(a,8);
    p = cns(a,4); %to sort 
    
    % robust fit on T vs Ed, pres 50-250 db
    ccc = (p>=50);   %Calculate for Pressure>=50 db; 
    t = t(ccc); %use t for median
    Ed = Ed(ccc);
    p = p(ccc);
    
    cc = (p <= 250);  %calculate for pressure <= 250
    t2 = t(cc);         %use t2, Ed2 for dE/dT
    Ed2 = Ed(cc);
    p = p(cc);
    
    maxT = max(t2);
    minT = min(t2);
    if mod(length(Ed2),2)==0 %Even # - median isn't real
        Ed3 = sort(Ed2);
        c = find(Ed3 > median(Ed2),1);
        c = find(Ed2 == Ed3(c),1); %index of median+ 
    else                    %odd number; median index exists
        c = find(Ed2 == median(Ed2),1);
    end
    TmedEd = t2(c);
   % t2=t2-minT; %sets y intercept @ Tmin
    
    if length(t)>5 %else too few to characterize blank  
     if maxT - minT > 2.5   %if Trange < 2.5, no dE/dT; use median
     [p2,stats]=robustfit(t2,Ed2);
     q=stats.se;
     s=stats.s;
     rho = corr(t2,Ed2,'Type','Spearman');
         if abs(rho) < 0.3 % No temp correlation found
           p2(1) = median(Ed2);
           p2(2) = 0;
           rhofail = rhofail+1;
         end
     else
         Trangefail = Trangefail +1; 
         p2(1) = median(Ed2);
         p2(2) = 0;
         q(1) = iqr(Ed2);
         q(2) = iqr(Ed2);
         rho = corr(t2,Ed2,'Type','Spearman');
     end
    
    N490_cors(i,1)=Id;
    N490_cors(i,2)=p2(1); %yintercept x0
    N490_cors(i,3)=p2(2); %slope  x1
    N490_cors(i,4)=q(1); %yintercept se
    N490_cors(i,5)=q(2); %slope se
    N490_cors(i,6)=maxT;
    N490_cors(i,7)=minT;
    N490_cors(i,8)=rho;   %Spearman's Rank Coefficient
    N490_cors(i,9)=median(Ed2);
    N490_cors(i,10)=prctile(Ed,2);
    N490_cors(i,11)=prctile(Ed,98);
    N490_cors(i,12) = TmedEd;
    N490_cors(i,13)=0;
    i = i+1;
    
    else
        short = short+1;
    end

tt = tt+length(a);
clear a
end




%% Filter for 1.5IQR 
N = N490_cors;
tt = find(~(N490_cors(:,3)==0));
N = N(tt,:);            %IQR test for those w/ dE/dT =/= 0 

  x1iqr = iqr(N(:,3));
  x1med = median(N(:,3));

  x1lb = x1med - 1.5*x1iqr;
  x1ub = x1med + 1.5*x1iqr;
  
  ct = 0;
  ct2 = 0;
  ctt = 0;
  G490cors = zeros(length(N490_cors),13); %the goods
  % G490cors is:
   % 1. Id 2. x0 3. x1 4. x0se 5. x1se 6. maxT 7. minT 8. rho 9. med(Ed2)
    % 10. 2% 11. 98% 12. T(medEd) 13. bool
    
for i=1:length(N490_cors)  %check IQR, save only those that pass
    x1 = N490_cors(i,3);
    
    if x1 > x1ub
        N490_cors(i,3) = x1ub;
        N490_cors(i,2) = N490_cors(i,9) - (x1ub*(N490_cors(i,12))); %x0
        N490_cors(i,13) = 1;
        ct = ct+1;
        ctt = ctt+1;
   
    else
         
        if x1 < x1lb
        N490_cors(i,3) = x1lb;
       N490_cors(i,2) = N490_cors(i,9) - (x1lb*(N490_cors(i,12))); %x0
        N490_cors(i,13) = 1;
        ct2 = ct2+1;
         ctt = ctt+1;
         
        end
    end

 G490cors(i,1:13) = N490_cors(i,1:13);

end



%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% FUNCTION
function [Tins] = SensorTemp(temps,time,K)   %
% Return integral int(0,t)Tenv(t')e^(t'/K)dt
 
t = 0:(time/(length(temps)-1)):time; %divide time into #units = #units temps
y = temps;           %array of Tenv
%y = fliplr(y);          %so temp array ascends from depth     

ex=exp(t/K);
 
F = y .* ex/K;        %Tenv(t')e^(-t'/K)dt; F   

if length(t)>1
    Tins = trapz(t,F);      %trapezoidal integration of F over t
else
    Tins = 0;         % <there is no integration if length(t) =1
end
end


