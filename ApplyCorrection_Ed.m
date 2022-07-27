% POST-PROCESS QUALITY CONTROL BGC ARGO RADIOMETRY 
% APPLY SENSOR-TEMPERATURE DEPENDENT CORRECTION
% PRODUCES QUALITY-CONTROLLED RADIOMETRY PROFILES "outc" (see section 3)
% 
% 
% FOLLOWS DAYMETHOD.m AND NIGHTMETHOD.m
% UPDATED 07.27.2022
% author: Terence O'Brien
% terence.obrien@maine.edu

% apply Ed correction
% E('actual') = E(meas) - [x0 + x1 (Tins)]

%steps: 
% - clean profiles per organelli method
% - don't remove the deep portion
% - compute Ts
% - apply correction
% - night if available
% - day if available
% - save both corrections

%% 1.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%%% 1.LOAD DATA, ASSESS PROFILES FOLLOWING QUALITY CONTROL PROCEDURES FROM ORGANELLI 2016
%%%% SCRIPT 1 ADAPTED FROM E.BOSS

clear all
close all
%Irradiance data is in W m^-2 nm^-1
%PAR in micromol quanta m^-2 s^-1
load Ed490_RAW_20220328.mat    
%load data490.mat %use for 490 and 412
%data380 = data490; %Set this to avoid changing variable all throughout script
%% 
A=data380.Properties.VariableNames;
for i=1:length(A)           %creates a logical array(isempty is yes/no)
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
    I(18,i)=isempty(find(string(A{1,i})=='DOWN_IRRADIANCE380'));
    I(19,i)=isempty(find(string(A{1,i})=='DOWN_IRRADIANCE412'));
    I(20,i)=isempty(find(string(A{1,i})=='DOWN_IRRADIANCE443'));
    I(21,i)=isempty(find(string(A{1,i})=='DOWN_IRRADIANCE490'));
    I(22,i)=isempty(find(string(A{1,i})=='WMO')); %float ID
end
%find the position of the variables of interest.
%find the position of the variables of interest.
ind=mod(find(I'==0),55);
ind(find(ind==0))=55;

 %% or load workspace unpacked 
load('UnpackedWorkspaceMay2022.mat')

%% 
data380 = data490;
out380=zeros(1,9);
counter=0;
Flag1=0; %number of profiles.
K=0;        %Nans
KK=0;       %bad pressure 
KKK=0;      %bad lat/long
KKKK=0;     % Solar Angle < 15
M=0;        % length profile < 10
MM=0;       % no blank profile
MMM=0;      % profile above dark < 5 (after first polytest)
L=0;        %'type 2' r2< 0.998  (after 2nd polytest)
LL=0;       % filter after 2nd polytest--removed profiles if isfinite(R(1,2))==0


for i=1:length(data380{:,1}) %analysis is performed a profile at a time
    i;
    Flag1;    
    % getting the profiles and test that we want to keep them
    Ed380=cell2mat(data380{i,ind(17)})*100;  %miltiply by 100 for comparison with Organelli's units
      % ^^^^^ !!! index is wavelength specific 14=380
                            % 15 = 412 17 = 490  
    Ed3802 = Ed380; 
    p=cell2mat(data490{i,ind(4)}); T=cell2mat(data490{i,ind(5)});S=cell2mat(data490{i,ind(6)});ID=data490{i,ind(18)};
    lat=data490{i,ind(2)};lon=data490{i,ind(3)};time=data490{i,ind(1)};
    %plot(Ed380,-p,'*')
    %pause
    
    %if profile is all NaN, remove it
    if isempty(find(isfinite(Ed380)))==1 K=K+1; continue; end 
    %Remove NaN from profile if not all NaN
    J=find(isfinite(Ed380));Ed380=Ed380(J);p=p(J); T=T(J); S=S(J);
    %if profile is not monotonic in pressure, remove it
    if sum(sort(p)-p) KK=KK+1; continue; end %how does this work %sort puts elements in ascending order
    %testing that reasonable LAT and LON, if not skip
    II=find(lat<=90 & lat>=-90 & lon<=180 & lon>=-180);
    if isempty(II)  KKK=KKK+1; continue; end  
    %testing that solar elevation is above 15degrees
    DateTimeUTC=time;
    [Az,El] = SolarAzEl(DateTimeUTC,lat,lon,0); 
    %detect only reasonable values
    if El < 15 KKKK=KKKK+1; continue; end
    % if there are less than 10 data point in a profile, skip it.
    if length(p) < 10 M=M+1; continue; end
    %compute blank based on Lilliefors test with alpha=0.01;
     blankEd380=0; %default 
     Ed=Ed380;pr=p;                         %Organelli a) dark signal
     for j=1:length(Ed)
       B=lilliefors(Ed'); 
       while B(2)<0.01 
           J=find(pr>min(pr));
           if J<10 break; end  %if too few deep values to characterize the blank
           pr=pr(J);
           minp=min(pr); %needed for future fits
           Ed=Ed(J);
           blankEd380=median(Ed);          % also storing dark value mean to compare
           B=lilliefors(Ed');               % ^ consider using median instead 
       end
       break
     end
     p_dark=min(pr); %keep for stats of depth below which we cannot compute K_L
    if B(2)<0.01 blankEd380=0; MM=MM+1; end % did not find a blank
    
    %now only work with data above noise signal
    JJ=find(p<minp);
    JK = find(p>=minp); pd = p(JK); p=p(JJ);
    Ed380d = Ed380(JK);
    Ed380=Ed380(JJ);
    pp = polyfit(p,log(Ed380),4);           %1st polytest
    res=log(Ed380)-polyval(pp,p);           
    for k=1:length(res)
        if abs(res(k)-mean(res))>2*std(res)     %flag3 values removed/NaN'd
            p(k)=NaN;
            Ed380(k)=NaN;
        end
    end
    III=isfinite(Ed380);
    if length(III) < 5 MMM=MMM+1; continue; end  %less than 5 values remain after 1st polytest--skip
    pp = polyfit(p(III),log(Ed380(III)),4);     %2nd polytest
    R=corrcoef(log(Ed380(III)),polyval(pp,p(III)));
    if R(1,2)^2<0.998  L=L+1; continue; end  %filter profiles that are not well behaved -'type 2' r2< 0.998
    if isfinite(R(1,2))==0  LL=LL+1; continue; end  %filter profiles that are not well behaved
    %compare to clear sky irradiance    %^skipping any type 2 profiles
    
   % kd380=pp(1)*4*pp(1).^3+pp(2)*3*p.^2+pp(3)*2*p+pp(4);    %diffuse attenuation coefficient
    location.latitude=lat;location.longitude=lon;location.altitude = 0;
    t.year=str2num(datestr(time,'yyyy'));t.month=str2num(datestr(time,'mm'));t.day=str2num(datestr(time,'dd'));
    t.hour=str2num(datestr(time,'hh'));t.min=str2num(datestr(time,'MM'));t.sec=0; t.UTC=0;
    [es,wl] = fullspec_es(t,location,-1,-1,-1);


    
   % figure(1)
   % title('Downwelling irradiance(Wm^-2 nm^-1) vs depth(pressure)')
   % plot(Ed380,-p,'*',exp(polyval(pp,p(III))),-p(III),'r'); 
   % hold on
   % plot(es380*10,0,'*r')
   % hold off
   % figure(2)
   % subplot(2,1,1)
   % plot(p,log(Ed380),'*',p(III),polyval(pp,p(III)),'r');
   % subplot(2,1,2)
   % plot(p(III),log(Ed380(III))-polyval(pp,p(III)),'*b');

  % keep whole profile, above 300 db
    Ed380 = union(Ed380,Ed380d,'stable');
    p = union(p,pd,'stable');
        aok = isfinite(Ed380);
        Ed380 = Ed380(aok);
        p = p(aok);


     JL = p<300;       %only taking p<300. 
     Ed380 = Ed380(JL); %oligotrophic waters may have light past 250db
     p = p(JL);


    IV = isfinite(Ed380);
    Ed380 = Ed380(IV);
    

    Flag1=Flag1+1;
    %Lat,Lon,YYYY,MM,DD,HH,MM,Pressure,ED_380,BBP700(m^-1),TEMP,SALINITY,CHL(ug/L),CDOM(ppb),ID
    
    L=length(Ed380(IV));
         out380(counter+1:counter+L,1)=location.latitude;
         out380(counter+1:counter+L,2)=location.longitude;
         out380(counter+1:counter+L,3)=DateTimeUTC;
         out380(counter+1:counter+L,4)=p(IV);
         out380(counter+1:counter+L,5)=T(IV);    % to make sure T & S align   
         out380(counter+1:counter+L,6)=S(IV);          %w/ saved Ed380
         out380(counter+1:counter+L,7)=Ed380;      
         out380(counter+1:counter+L,8)=ID; 
         %out380(counter+1:counter+L,9) = blankEd380;
         counter=counter+L;
% 
     

     blank(Flag1)=blankEd380;
      lat(Flag1)=location.latitude;
      lon(Flag1)=location.longitude;
      dark_pres(Flag1)=p_dark;
    clear p pp pr pd Ed Ed380 Ed380d JK JL JJ I II III IV aok location time blankEd380 p_dark
    if mod(i,10000)==0
    disp(i)
    end
end

%% remove missing temperatures if any:
a = find(isnan(out380(:,5)));
if ~isempty(a)
out380(a,:)=NaN;
for ii=1:9
    cc = out380(:,ii);
    c = ~isnan(cc);
    out380(:,ii) = NaN;
    cc = cc(c);
    for j=1:length(cc)
        out380(j,ii) = cc(j);
    end
    
end
out380 = out380(1:find(isnan(out380(:,1)),1),:);
end
%remove any w Ed values <-100 (a few exist..._)
b = find(out380(:,7) < -100);
if ~isempty(b) 
out380(b,:)=NaN;
for ii=1:9
    cc = out380(:,ii);
    c = ~isnan(cc);
    out380(:,ii) = NaN;
    cc = cc(c);
    for j=1:length(cc)
        out380(j,ii) = cc(j);
    end
    
end

out380 = out380(1:find(isnan(out380(:,1)),1),:);
end
 if isempty(a) && isempty(b)
    out380(length(out380)+1,1:9) = NaN;  % Intentionally set the last line as NaNs
 end
%% 
clearvars -except out380 data380 data490


 
%% 2 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%%% 2. COMPUTE SENSOR TEMPERATURE
%%%% CALLS ON FUNCTION(SENSORTEMP)
%%%% TENV TO TINS


%1.lat 2.long 3.date 4.pres 5.Ed 6.temp 7.salinity 8.ID 9.Eddarkmedian

cns=out380; 


Tsen3=NaN(length(cns),9);   %Final Output w corrected temps, and everything else.
x=1; %for Tenv index. resets to one with new lt
Tenv=NaN(10,2);  %to extract Tenv from cns per profile
Ts3 = NaN(2,1); %for use inside loops
lt=cns(1,3); %current time
plt=cns(1,3); %previous time
cnt=1; %do not reset. row counter for Tsen
n=1; %counter for Tsen(x,temp)=Ts(n,1) 
qrs = 0;
K = 200;
i = 1;
while i<length(cns)
    a = find(cns(:,1)==cns(i,1) & cns(:,2) == cns(i,2));
    cp = cns(a,:);
    
    temps = cp(:,5);
    pres = cp(:,4);
        
       if ~isempty(temps)
           
        c = ~isnan(temps);
        temps = temps(c);
        pres3 = pres(c);

        temps = flipud(temps);  % so T(deepest) is index 1
        t = zeros(length(pres3),1);
        for tt = 1:length(pres3)
            t(tt) = (max(pres3)-pres3(tt))/0.1;
        end
       t = flipud(t);
        ct=zeros(1);          %"current temp"
        Ts3=zeros(length(t),1);
        pres3 = flipud(pres3); %to fit w temps
   
       To = temps(1)-0.1;
       
       %for profiles with deep values
       if max(pres3)>260 %&& length(find(pres3>250)) >5
              pr = find(pres3>248);
              tempr = temps(pr);
              timer = t(pr);
              pr = pres3(pr);
              filledpres = min(pr):1:max(pr);
              if (max(tempr)-min(tempr))>0
              filledtemp = min(tempr):(max(tempr)-min(tempr))/(length(filledpres)-1):max(tempr);
              else 
                  for ppp = 1:length(filledpres)     %some short deep profiles show const temp.
                      filledtemp(ppp) = tempr(1);
                  end
              end
              filledtime = min(timer):(max(timer)-min(timer))/(length(filledpres)-1):max(timer);
              filledtime(1) = 1;
              preslong = sort([filledpres'; pres3],'descend');
              templong = sort([filledtemp';temps]);
              timelong = sort([filledtime';t]);
%               if length(pr)<5
%                   pause()
%               end

               q = zeros(length(pres3),1);
              for iv = 1:length(pres3)
                   q(iv) = find(preslong == pres3(iv),1);
               end

            for ii=1:length(timelong)      %to calculate Tins 
             ct(ii)=templong(ii); %current temp
    
              intg=SensorTemp(ct,timelong(ii),K);           %function to return integral over current t range
    
              % Ts(t)=T(0)e^(-t/K)+(e^(-t/K)int(0,t')Tenv(t')e^(t'/K)dt)/K
              Ts3(ii)=To*exp(-timelong(ii)/K)+exp(-timelong(ii)/K)*intg; %creates Ts array of Sensor temps
              
            end
            Ts3 = Ts3(q);
      
       else 
           
           on = pres3(1);
             w = find((pres3 < on - 20),1);
             
                 dT = temps(1:w);
                 dz = pres3(1:w);
                 fit = polyfit(dz,dT,1);
                 dTdt = fit(1);   
                 To = temps(1)+dTdt*20; 
                
          if To > temps(1)
              To = temps(1)-0.1; %if To is warmer than T(1), revert to T(1)
              qrs = qrs+1;
          end
                
            for ii=1:length(t)      %to calculate Tins 
             ct(ii)=temps(ii);
    
              intg=SensorTemp(ct,t(ii),K);           %function to return integral over current t range
    
              % Ts(t)=T(0)e^(-t/K)+(e^(-t/K)int(0,t')Tenv(t')e^(t'/K)dt)/K
              Ts3(ii,1)=To*exp(-t(ii)/K)+exp(-t(ii)/K)*intg; %creates Ts array of Sensor temps
            end
            
           
        
            %Now fill Tsen with Ts and all other data.   
        end
            Ts3 = flipud(Ts3); %inverted compared to rest of data

            for iii=i:i+length(a)-1
                %if cns(iii,4)>= 50
                Tsen3(cnt,1)=cns(iii,1);
                Tsen3(cnt,2)=cns(iii,2);
                Tsen3(cnt,3)=cns(iii,3);
                Tsen3(cnt,4)=cns(iii,4);
                Tsen3(cnt,5)=cns(iii,5);     %Ts is a row array
                if Tsen3(cnt,4)>250 && Ts3(n,1)> cns(iii,5)
                   Tsen3(cnt,6) = cns(iii,5); %if a deep Tsen value is>Tm revert to Tm
                else 
                    
                     Tsen3(cnt,6)=Ts3(n,1);
                end
%                 if Tsen3(cnt,6)-Tsen3(cnt,5)>2
%                     pause()
%                 end
                   
                Tsen3(cnt,7)=cns(iii,6);
                Tsen3(cnt,8)=cns(iii,7);
                Tsen3(cnt,9)=cns(iii,8);
                cnt=cnt+1;
                n=n+1;

            end

          n=1; %n is reset 
         plt = lt;  
        x=1; %x is reset;
        Tenv=NaN(10,2);
        %Ts3=NaN(2,1);
        
        %now take out new Ts since we have lt=/=plt
       Tenv(x,1)=cns(i,5); %temperature
       Tenv(x,2)=cns(i,4);
       x=x+1;
       end 
       i = i+length(a);
end

out380Ts = Tsen3;

%% 
dif = out380Ts(:,6)-out380Ts(:,5);
hist(dif)
%% Nans out
for i=1:9
    cc = out380Ts(:,i);
    c = ~isnan(cc);
    out380Ts(:,i) = NaN;
    cc = cc(c);
    for j=1:length(cc)
        out380Ts(j,i) = cc(j);
    end
    
end

clearvars -except out380 data380 out380Ts

%% 3 change load

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%%% 3. APPLY CORRECTION TO GOOD PROFILES
%%%% LOAD NIGHT & DAY CORRECTIONS 
%%%% TENV TO TINS
% keep both corrections or if no correction make nan 


load G490corsDay.mat % WHATEVER DAY CORRECTIONS ARE SAVED AS 412,490...
Dcors = G490cors;
load G490corsNight.mat % WHATEVER NIGHT CORRECTIONS ARE SAVED AS
Ncors = G490cors;
%% 
out = out380Ts;
%clear out380Ts G380cors
pID = 1; %arbitrary

x = 1;
outc = zeros(10000,13);
% out is: 1.lat 2.long 3.date 4.pres 5.tm 6.ts 7.Edm 8.sal 9.Id 10.med_dark 11.cor_value_M1 12.Ed_cor_M1 13.cor_val_M2 14. Ed_cor_M2 

 ii = 1;
 while ii <= length(out)
     
    ID = out(ii,9);
    a = find(Ncors(:,1)==ID,1); %find Ncors 
      bb = find(out(:,9)==ID); %
    if isfinite(a) %Night correction found. apply night correction

           M1x0 = Ncors(a,2);
           M1x1 = Ncors(a,3);
           
       outc(x:x+length(bb),1:9) = out(ii:ii+length(bb),1:9);    
       outc(x:x+length(bb),10) = (M1x0 + M1x1*(out(ii:ii+length(bb),6)));%value of correction M1 Night
       outc(x:x+length(bb),11) = out(ii:ii+length(bb),8) - (M1x0 + M1x1*(out(ii:ii+length(bb),6))); %corrected Ed

    else 
                 outc(x:x+length(bb),1:9) = out(ii:ii+length(bb),1:9);    
              outc(x:x+length(bb),10) = NaN;%value of correction M1 Night
       outc(x:x+length(bb),11) = NaN; %corrected Ed
    end
    
     b = find(Dcors(:,1)==ID,1); %no night correction: use day
            if isfinite(b)
              M2x0 = Dcors(b,2);
              M2x1 = Dcors(b,3);
           
             outc(x:x+length(bb),12) = (M2x0 + M2x1*(out(ii:ii+length(bb),6)));%value of correction M2 day
              outc(x:x+length(bb),13) = out(ii:ii+length(bb),8) - (M2x0 + M2x1*(out(ii:ii+length(bb),6))); %corrected Ed
           
                 
            else %no correction at all
                   outc(x:x+length(bb),12) = NaN;%value of correction M2 day
              outc(x:x+length(bb),13) = NaN; %corrected Ed
            end
            
                   if isempty(bb)       % to prevent long searches if possible
                       bb = find(~(out(ii:length(out),9) == ID),1);
                    ii = ii+ bb;
                    disp(ii)
                   end
        
                 x = x+length(bb);
                 ii = ii+ length(bb);
                 disp(ii)
   
    end
 



%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% 
subplot(211)
hist(outc(:,10))
title('M1 correction')
subplot(212)
hist(outc(:,12))
title('M2 correction')
%%
outc490 = outc;
%% FUNCTIONS

function [Tins] = SensorTemp(temps,time,K)   %
% Return integral int(0,t)Tenv(t')e^(t'/K)/Kdt
 if time > 0
t = 0:(time/(length(temps)-1)):time; %divide time into #units = #units temps
y = temps;           %array of Tenv
%y = fliplr(y);          %so temp array ascends from depth     

ex=exp(t/K);
 
F = y .* ex/K;        %Tenv(t')e^(-t'/K)dt; F   

if length(t)>1      %CHANGED THIS from 1 to 4
    Tins = trapz(t,F);      %trapezoidal integration of F over t
else
    Tins = 0;         % <there is no integration if length(t) < 5
end
 else
    Tins = 0;  % <there is no integration if time = 0
end
end

function [Results]=lilliefors(x)
% KOLMOGOROV-SMIRNOV TEST - LILLIEFORS MODIFICATION
alpha=0.01;
n=length(x);
i=1:n;
y=sort(x);
fx=normcdf(zscore(y));
dplus=max(abs(fx-i/n));
dminus=max(abs(fx-(i-1)/n));
Dn=max(dplus,dminus);
% P = [n D20 D15]
P=[5 0.289 0.303;
   6 0.269 0.281;
   7 0.252 0.264;
   8 0.239 0.250;
   9 0.227 0.238;
   10 0.217 0.228;
   11 0.208 0.218;
   12 0.200 0.210;
   13 0.193 0.202;
   14 0.187 0.196;
   15 0.181 0.190;
   16 0.176 0.184;
   17 0.171 0.179;
   18 0.167 0.175;
   19 0.163 0.170;
   20 0.159 0.166;
   25 0.143 0.150;
   30 0.131 0.138;
   40 0.115 0.120;
   100 0.074 0.077;
   400 0.037 0.039;
   900 0.025 0.026];
 
aaa=P(:,1)';
subind=max(find(aaa<n));
upind=subind+1;
if (subind==22) upind=subind; subind=subind-1; end
xxx=P(subind:upind,:);
 
if aaa(upind)==n
   D20=xxx(2,2);
   D15=xxx(2,3);
else
    D20=xxx(1,2)+(n-aaa(subind))*((xxx(2,2)-xxx(1,2))/(xxx(2,1)-xxx(1,1)));
    D15=xxx(1,3)+(n-aaa(subind))*((xxx(2,3)-xxx(1,3))/(xxx(2,1)-xxx(1,1)));
end
 
a1=-7.01256*(n+2.78019);
b1=2.99587*sqrt(n+2.78019);
c1=2.1804661+0.974598/sqrt(n)+1.67997/n;
 
a2=-7.90289126054*(n^0.98);
b2=3.180370175721*(n^0.49);
c2=2.2947256;
 
if n>100
   D10=(-b2-sqrt(b2^2-4*a2*c2))/(2*a2);
   a=a2;
   b=b2;
   c=c2;
else
   D10=(-b1-sqrt(b1^2-4*a1*c1))/(2*a1);
   a=a1;
   b=b1;
   c=c1;
end
 
if Dn==D10
        pvalue=0.10;
        
    elseif Dn>D10
        pvalue=exp(a*Dn^2+b*Dn+c-2.3025851);
    elseif Dn>=D15
        pvalue=((0.10-0.15)/(D10-D15))*(Dn-D15)+0.15;
    elseif Dn>=D20
       pvalue=((0.15-0.20)/(D15-D20))*(Dn-D20)+0.20;
    else
       pvalue=0.20;
end
Results(1)=Dn;
Results(2)=pvalue;
end
