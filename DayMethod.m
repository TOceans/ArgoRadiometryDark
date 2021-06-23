% POST-PROCESS QUALITY CONTROL BGC ARGO RADIOMETRY 
% CALCULATE SENSOR-TEMPERATURE DEPENDENT CORRECTION
% DAY METHOD
% UPDATED 6.23.2021
% author: Terence O'Brien
% terence.obrien@maine.edu

% 5 SECTIONS
% CAN RUN WHOLE SCRIPT OR SECTION-BY-SECTION (TO CHECK & RECORD STATS)
% Output is "G(wavelength)cors" array of x0 & x1 per float (see section 5)

%% 1.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%%% 1. LOAD DATA, ASSESS PROFILES FOLLOWING QUALITY CONTROL PROCEDURES FROM ORGANELLI 2016
%%%% SCRIPT 1 ADAPTED FROM EMMANUEL BOSS
%%%% section 1 takes many minutes to run


clear all
close all
%Irradiance data in W m^-2 nm^-1
% PAR in micromol quanta m^-2 s^-1


%load data490.mat  %use for 490 & 412
load data380.mat
%% 
A=data380.Properties.VariableNames;
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
    I(18,i)=isempty(find(string(A{1,i})=='DOWN_IRRADIANCE380'));
    I(19,i)=isempty(find(string(A{1,i})=='DOWN_IRRADIANCE412'));
    I(20,i)=isempty(find(string(A{1,i})=='DOWN_IRRADIANCE443'));
    I(21,i)=isempty(find(string(A{1,i})=='DOWN_IRRADIANCE490'));
    I(22,i)=isempty(find(string(A{1,i})=='WMO')); %float ID
end

%find the position of the variables of interest.
ind=mod(find(I'==0),47);
ind(find(ind==0))=47;


darkout=zeros(1,8);
counter=0;
y=1;
Flag1=0; %number of profiles.
K=0; %Ed profile is Nans
KK=0; %pressure not monotonic
KKK=0; % lat or long not real
KKKK=0; %solar az angle test
M=0; % short profile
minp = 50;

for i=1:length(data380{:,1}) %analysis is performed a profile at a time
    i;
    Flag1;
    % getting the profiles and test that we want to keep them
    Ed380=cell2mat(data380{i,ind(18)})*100;  %miltiply by 100 for comparison with Organelli's units
                           % ^^^^^ !!! index is wavelength specific
    p=cell2mat(data380{i,ind(4)}); T=cell2mat(data380{i,ind(5)});S=cell2mat(data380{i,ind(6)});ID=data380{i,ind(22)};
    lat=data380{i,ind(2)};lon=data380{i,ind(3)};time=data380{i,ind(1)};
  %  plot(Ed380,-p,'*')
    %pause
    
    %if profile is all NaN, remove it
    if isempty(find(isfinite(Ed380)))==1 K=K+1; continue; end 
    %Remove NaN from profile in not all NaN
    J=find(isfinite(Ed380));Ed380=Ed380(J);p=p(J); T=T(J); S=S(J);
    %if profile is not monotonic in pressure, remove it
    if sum(sort(p)-p) KK=KK+1; continue; end 
    %testing that reasonable LAT and LON, if not skip
    II=find(lat<=90 & lat>=-90 & lon<=180 & lon>=-180);
    if isempty(II)  KKK=KKK+1; continue; end  
    %testing that solar elevation is above 15 degrees
    DateTimeUTC=time;
    [Az,El] = SolarAzEl(DateTimeUTC,lat,lon,0);
    %detect only reasonable values
    if El < 15  KKKK=KKKK+1; continue; end    % % %SAVE THESE, FOR EL<0
   
    % if there are less than 10 data point in a profile, skip it.
    if length(Ed380) < 10 M=M+1; continue; end
    %compute blank based on Lilliefors test with alpha=0.01;
     blankEd380=0; %default 
     Ed=Ed380;pr=p;
     for j=1:length(Ed)
       B=lilliefors(Ed');
       while B(2)<0.01 
           J=find(pr>min(pr));
          
           if J<10 break; end%if too few deep values to characterize the blank
           pr=pr(J);
           minp=min(pr); %needed for future fits
          % Ed=Ed(J);
          % blankEd380=mean(Ed);
           B=lilliefors(Ed');
       end
         location.latitude=lat;location.longitude=lon;location.altitude = 0;
    t.year=str2num(datestr(time,'yyyy'));t.month=str2num(datestr(time,'mm'));t.day=str2num(datestr(time,'dd'));
    t.hour=str2num(datestr(time,'hh'));t.min=str2num(datestr(time,'MM'));t.sec=0; t.UTC=0;
    JJ=find(p>minp);
    p = p(JJ); T = T(JJ); S = S(JJ); Ed = Ed(JJ);
    JJJ = find(p>50);
       L=length(JJJ);
   
         darkout(counter+1:counter+L,1)=location.latitude;
         darkout(counter+1:counter+L,2)=location.longitude;
         darkout(counter+1:counter+L,3)=DateTimeUTC;
         darkout(counter+1:counter+L,4)=p(JJJ); 
         darkout(counter+1:counter+L,5)=T(JJJ);       
         darkout(counter+1:counter+L,6)=S(JJJ);
         darkout(counter+1:counter+L,7)=Ed(JJJ);    
         darkout(counter+1:counter+L,8)=ID;
         counter=counter+L;
         clear p pp pr Ed Ed380 II location time JJ JJJ
     
        
       break
     end
     if (mod(i,5000)==0)
         disp(i)
     end
    
end

%%  remove any lines with missing temperatures
 a = find(isnan(darkout(:,5)));
 if ~isempty(a)
     darkout(a,:)=NaN;
     
for ii=1:8
    cc = darkout(:,ii);
    c = ~isnan(cc);
    darkout(:,ii) = NaN;
    cc = cc(c);
    for j=1:length(cc)
        darkout(j,ii) = cc(j);
    end
    
end

darkout = darkout(1:find(isnan(darkout(:,1)),1),:);
else 
    darkout(length(darkout)+1,1:8) = NaN;  % Intentionally setting the last line as NaNs
end
 
 %% clear vars 
clearvars -except darkout data380


%% 2
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%%% 2. COMPUTE SENSOR TEMPERATURE
%%%% CALLS ON FUNCTION(SENSORTEMP)
%%%% TENV TO TINS


cns=darkout; 
Tsen3=NaN(length(cns),9);   %Final Output w corrected temps, and everything else. 
x=1; %for Tenv index. resets to 1 with new lt. 
Tenv=NaN(10,2);  %to extract Tenv from cns per profile
Ts3 = NaN(2,1); %for use inside loops
lt=cns(1,1); %current latitute
plt=cns(1,1); %previous latitude
cnt=1; %do not reset. row counter for Tsen
n=1; %counter for Tsen(x,temp)=Ts(n,1); resets to 1
 K=200; %Time constant for function [s]

for i=1:length(cns)
    lt=cns(i,1);
    
    if lt==plt && ~isnan(lt) %same latitude. extract Tenv 
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
        if min(pres2)>250   %Tconverion only on p<=250
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
              % n=n+1;
              
            end
        
        else
        cc = pres2<= 250; % consistent Ed readings begin @250db  
        temps = temps(cc);
        pres3 = pres(cc);
        temps = flipud(temps);  % so T(deepest) is index 1
       % time = pres3(length(pres3))/.1;  %total time in seconds for float to rise, assuming 0.1m/s
        time = (max(pres3)-min(pres3))/0.1;
        t = 0:(time/(length(temps)-1)):time; %array of times, same length as temps
      %  To=temps(1);         %T(0)
        ct=zeros(1);          %"current temp"
        Ts3=zeros(1,length(t));
        pres3 = flipud(pres3); %to fit w temps
            
        if length(temps)>4 %max(temps)-min(temps) > 2
             % calculate av dTdz for first 20db
             on = pres3(1);
             w = find((pres3 < on - 20),1);
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
                 dTdt = abs(fit(1));
                 To = temps(1)-dTdt*(max(pres3)-min(pres3));
             end    
         %now convert Tenv to Ts
            for ii=1:length(t)      %to calculate Tins 
             ct(ii)=temps(ii);
    
              intg=SensorTemp(ct,t(ii),K);           %function to return integral over current t range
    
              %Ts(t)=T(0)e^(-t/K)+(e^(-t/K)int(0,t')Tenv(t')e^(t'/K)dt)/K
              Ts3(ii)=To*exp(-t(ii)/K)+exp(-t(ii)/K)*intg; %creates Ts array of Sensor temps
        
            end
        else 
            Ts3 = temps';
        end    
            %Now fill Tsen with Ts and all other data.            
            Ts3 = fliplr(Ts3); %inverted compared to rest of data
            j=i-length(pres2);
            k=i-1;
            for iii=j:k
                if pres2(n)<250
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
                n=n+1;
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

clearvars -except darkout data380 Tsen3


%% 3
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%%% 3. DOWNWELLING IRRADIANCE TEST 
%%%% DETECT AND REMOVE PROFILES DETECTING DOWNWELLING IRRADIANCE
%%%% 

ns=Tsen3;   


lt=1;   %latitude
plt=0;  %must start different than lt for first profile
fail=1;
ctf=0;  % counts failed
ctT = 0; %total profiles


 
for i=1:length(ns)
    pr=NaN(2,1); %reset
    e=NaN(2,1);  %reset
    lt=ns(i,1);  %current latitude
    
    if not(lt==plt) %when @new lat/profile, determine FAIL as 1 or zero. 
        ctT = ctT +1;
        for l=0:1000 %determine if new profile is good or bad
                    % also assuming no profile has more than 1000 values
                    
           if ns(l+i,1)==lt && l+i <length(Tsen3(:,1))%if latitude is the same
             %extract pr and e for current profile
             pr(l+1)=ns(l+i,4); %pressure
             e(l+1)=ns(l+i,8); % Ed
           
             
             else  %=@new latitude/profile. run calculation on current profile. 
           
           minEd = min(e);
              for r=1:length(e)             %shift data so none are negative for log test
                 e(r)=e(r)-minEd + 0.001;
              end
           
             p=polyfit(pr,log(e),1);
             
             if p(1) < -0.01
                  fail=1;
                ns(i,:)=NaN;
                ctf=ctf+1;
              
        %for visualizing if wanted
%              pp=polyval(p,pr); 
%              val=num2str(p);
%              plot(pr,log(e),'*',pr,pp,'-') 
%              set(gca,'LineWidth',1.5)
%              xlabel('Pressure(db)')
%              ylabel('log(E)')
%             % xlabel(val)
%              pause()
               
             else
               fail=0;
             end
                
             break
            
           end      % if ns(l+i,1)==lt
        end         %for l=1:1000
       
    else 
       if fail==1 
       ns(i,:)=NaN;  
       end
       
    end     %if not(lt==plt)
    
    plt=lt;     
end


%% Remove NaNs/condense


for i=1:9
    cc = ns(:,i);
    c = ~isnan(cc);
    ns(:,i) = NaN;
    cc = cc(c);
    for j=1:length(cc)
        ns(j,i) = cc(j);
    end
    
end

ns = ns(1:find(isnan(ns(:,2)),1),:); 


clearvars -except darkout data380 Tsen3 ns


%% 4
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%%% 4. MAX IRRADIANCE FILTER
%%%% 
%%%% REMOVE VALUES OUTSIDE [-0.3 0.3] W m-2 nm-1

fns = ns;

ctH = 0;
ctL = 0;

for i = 1:length(fns)   
    if fns(i,8)> 0.03
        fns(i,:) = NaN;
        ctH = ctH +1;
    end
    if fns(i,8) < -0.03
        fns(i,:) = NaN;
        ctL = ctL +1;
    end
end


%% Remove Nans

ns2 = fns;
for i=1:9
    cc = fns(:,i);
    c = ~isnan(cc);
    ns2(:,i) = NaN;
    cc = cc(c);
    for j=1:length(cc)
        ns2(j,i) = cc(j);
    end
    
end


ns2 = ns2(1:(find(isnan(ns2(:,1)),1)),:);

clearvars -except darkout data380 Tsen3 ns ns2

%% 5
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%%% 5. ROBUST FIT AND dE/dT CONSTRAINT
%%%% ROBUST LINEAR FIT ON COMPILED PROFILES PER FLOAT PER WAVELENGTH
%%%% 1.5*IQR CONSTRAINT ACROSS ALL FLOATS dE/dT AT SAME WAVELENGTH


%ns2: 1.lat 2.long 3.time 4.pres 5.Tm 6.Ts 7.Sal 8.Ed 9.ID

cns = ns2; 
D380_cors=NaN(1,8);
short = 0; %counter 
tt = 1;
i = 1;
while tt < length(cns)
    Id=cns(tt,9);
    
    % extract temp, ed, pres for Id
    a = find(cns(:,9)==Id);
    t = cns(a,6);
    Ed = cns(a,8);
    p = cns(a,4); %to sort 
    
    % robust fit on T vs Ed, pres > 50 db
    
     ccc = (p>=50);   %Calculate for Pressure>=50 db; 
     t = t(ccc);
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
    else                    %odd number; median index exists
        c = find(Ed2 == median(Ed2),1);
    end
    TmedEd = t2(c);
 %   t2=t2-minT; %sets y intercept @ Tmin
    
    if length(t2)>5 %else too few to characterize blank
      if (maxT-minT) > 2.5 %if Trange < 2.5 use median
        [p2,stats]=robustfit(t2,Ed2);
         q=stats.se;
         s=stats.s;
         rho = corr(t2,Ed2,'Type','Spearman');
         
         if abs(rho) < 0.3  %No temp correlation found
           p2(1) = median(Ed2);
           p2(2) = 0;
         end
      else
         p2(1) = median(Ed2);
         p2(2) = 0;
         q(1) = iqr(Ed2);
         q(2) = iqr(Ed2);
         rho = corr(t2,Ed2,'Type','Spearman');
     end
    D380_cors(i,1)=Id;
    D380_cors(i,2)=p2(1); %yintercept x0
    D380_cors(i,3)=p2(2); %slope  x1
    D380_cors(i,4)=q(1); %yintercept se
    D380_cors(i,5)=q(2); %slope se
    D380_cors(i,6)=maxT;
    D380_cors(i,7)=minT;
    D380_cors(i,8)=rho;   %Spearman's Rank coefficient
    D380_cors(i,9)=median(Ed2);
    D380_cors(i,10)=prctile(Ed,2); 
    D380_cors(i,11)=prctile(Ed,98);
    D380_cors(i,12) = TmedEd;
    D380_cors(i,13) = 0;
    i = i+1;
    % 1. Id 2. x0 3. x1 4. x0se 5. x1se 6. maxT 7. minT 8. rho 9. med(Ed2)
    % 10. 2% 11. 98% 12. T(medEd) 13. bool
    else
        short = short+1;
    end

tt = tt+length(a);
clear a
end


%% 1.5*IQR filter
N = D380_cors;
NN = D380_cors;
tt = find(~(D380_cors(:,3)==0));
N = N(tt,:);            %IQR test for those w/ dE/dT =/= 0 
  x1iqr = iqr(N(:,3));
  x1med = median(N(:,3));
  x1lb = x1med - 1.5*x1iqr;
  x1ub = x1med + 1.5*x1iqr;
  
  ct = 1;
  ct2 = 1;
  G380cors = zeros(length(D380_cors),13); %the goods
  % G380cors is:
   % 1. Id 2. x0 3. x1 4. x0se 5. x1se 6. maxT 7. minT 8. rho 9. med(Ed2)
    % 10. 2% 11. 98% 12. T(medEd) 13. bool
for i=1:length(D380_cors)  %check IQR, save only those that pass
    ID = NN(i,1);
    x1 = NN(i,3);
    %pass if x0,x1 <= |1.5*iqr() +/- med|
    if x1 > x1ub
        D380_cors(i,3) = x1ub;
        D380_cors(i,2) = D380_cors(i,9) - (x1ub*(D380_cors(i,12))); %x0
       D380_cors(i,13) = 1;
        ct = ct+1;
    end
    if x1<x1lb
        D380_cors(i,3) = x1lb;
        D380_cors(i,2) = D380_cors(i,9) - (x1lb*(D380_cors(i,12))); %x0
         D380_cors(i,13) = 1;
        ct2 = ct2+1;
    end

 G380cors(i,1:13) = D380_cors(i,1:13);

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











 %%
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