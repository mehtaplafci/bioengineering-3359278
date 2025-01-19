clear all 
close all

global  k12 k11 k10 k5 cIL1M0 k6 cIL10M0 k4 c2 cIL10C cIL1C cIL10F dIL10 dIL1 k2 k3 c1 k7 k8 k9 
global  yMc mMc muM yF mF muMc muS k14  
global  k13 muT muRT cRT k1 dR k22 dIL17 k18 k19 cIL17 cIL17M0 muN muRN cRN k15 cIL17Mc k16 muRMc cRMc k17 k21 k20  muRM0 cRM0 muRM1 cRM1 muRM2 cRM2 muRMS cRS omega gamma delta c3


% final time "tfinal" Oxygen (in seconds)
%tfinal=1/12; %1/6 1/4 1/2 1 ORT İS OXYGEN RESTORE TIME 
for tfinal=[1/12 1/6 1/4 1/2 1]

%for tfinal=1/12

for steminit1=[0 2*10^7]

TT=[];
PP=[];




% final time "tfinal1" Stem Cells (in seconds)
for tfinal1=(tfinal+0.01):1:30
    
k2=9*10^(-2);   %estimate
k3=0.0004;      %estimate    
k4=0.0005;          %from Wang 2012
k5=0.0005;          %from Wang 2012
k6=4;       %estimate
k7=0.7;     %estimate
k8=0.3;             %from Wang 2012
k10=26*10^(5);      %from Wang 2012
k11=0.0003;         %from Jin 2011
k12=0.25;           %from Jin 2011
cIL10M0=5;          %from Wang 2012
cIL1M0=10;          %from Wang 2012
dIL10=5;             %from Wang 2011
dIL1=10.5;           %from Wang 2012
muM=0.2;              %from Wang 2012
c1=25;                %from Wang 2012
c2=100;                %from Wang 2012
cIL10C=5;         %estimate
cIL1C=10;         %estimate
cIL10F=2.5;       %estimate

%treatment
     
%health initial conditions
mcinit=4*10^7; 
moinit=2*10^3;  
m1init=0;
m2init=0;
il10init=0.01;
il1init=0.1;
cinit=8395*10^8;
finit=1*10^8;
mdinit=0;  
tinit=0;    %change 
il17init=0.1;   %change 
rinit=210;    %change 
ninit=0;    %change 



tfinal2=30;    %change

mMc=5*10^(0); 
mF=5*10^(0);  
yF=0;
yMc=0;
muMc=0.3;      %0.3 mild MI estimate       %first model  4.0677 severe MI
k9=0.075;
muS=2;            %estimate
muN=0.3;          %from Dunster 2014 
k14=0;     
k1=0;       
dR=0; 
muRN=0; 
muRMc=0;  
muRM0=0; 
muRM1=0;
muRM2=0;
muRMS=0;
muRT=0; 
omega=1;
gamma=1;
delta=1; 


%
k13=4;                %change  
muT=2;                %change  
k22=1*10^(-6);        %change   k5=5*10^(-4);
dIL17=10.5;           %change  
k18=4;                %change   
k19=1;                %change                        ??????????????????similar???????????
k15=1*10^(-3) ; %1*10^(-4);       %change  muMc=3*10^(-1)  
cIL17Mc=10;           %change
k16=1*10^(-7);   % similar k2=9*10^(-2) (1*10^(-4))/(12*10^(4))
cRMc=1*10^2;              %change
k17=1*10^(-6);        %change  % similar k2=9*10^(-2)  ?????????????????????????????????????
k21=1*10^(-6);        %change  % similar k5=0.0005             ????????????????????????????????????
k20=0.7;              %change   k7=0.7;
cIL17=10^2;             %change                       ??????????????????similar???????????
cIL17M0=10;           %change
cRM0=1*10^2;     %change
cRM1=1*10^2;     %change
cRM2=1*10^2;     %change
cRS=1*10^2;      %change
cRT=1*10^2;      %change
cRN=1*10^2;      %change 
c3=100;          %change



sinit=4000;    
[T1,X1] = ode23s(@treatwithTcells_model,[0 tfinal],[mcinit moinit m1init m2init il10init il1init cinit finit mdinit sinit tinit rinit il17init ninit]); 

k9=0.075;
yMc=0.9; %severe MI 2.52;  %9*10^(-1) mild MI
yF=9*10^(-1);  
muMc=0; %second model
muS=2;
k14=10^(-5);   %estimate
%change the following parameters
k1=0.3;                          %??????????????????similar???????????  
dR=10.5;                          %??????????????????similar???????????     
muRN=1.5*10^(-3);    
muRMc=1.5*10^(-3);  
muRM0=1.5*10^(-3); 
muRM1=1.5*10^(-3);
muRM2=1.5*10^(-3);
muRMS=1.5*10^(-3);
muRT=1.5*10^(-3);
omega=0;            %change 
gamma=0;            %change 
delta=0;            %change 

%%%%%%%%%%%%%

[T2, X2] = ode23s(@treatwithTcells_model, [tfinal tfinal1], [X1(length (X1(:,1)),1) X1(length (X1(:,2)), 2) X1(length (X1(:,3)), 3) X1(length (X1(:,4)), 4) X1(length (X1(:,5)), 5) X1(length (X1(:,6)), 6) X1(length (X1(:,7)), 7) X1(length (X1(:,8)), 8) X1(length (X1(:,9)), 9) X1(length (X1(:,10)),10) X1(length (X1(:,11)),11) X1(length (X1(:,12)),12) X1(length (X1(:,13)),13) X1(length (X1(:,14)),14)] ); 

k9=0.15; %after stem cell injection

[T3, X3] = ode23s(@treatwithTcells_model, [tfinal1 tfinal2], [X2(length (X2(:,1)),1) X2(length (X2(:,2)), 2) X2(length (X2(:,3)), 3) X2(length (X2(:,4)), 4) X2(length (X2(:,5)), 5) X2(length (X2(:,6)), 6) X2(length (X2(:,7)), 7) X2(length (X2(:,8)), 8) X2(length (X2(:,9)), 9) X2(length (X2(:,10)),10)+steminit1 X2(length (X2(:,11)), 11) X2(length (X2(:,12)), 12) X2(length (X2(:,13)), 13) X2(length (X2(:,14)), 14)] ); 
 
% figure(25);


T=[T1; T2; T3];
X=[X1; X2; X3];



mcinit
X(length(X(:,1)),1)
tfinal1

FinalPercentage=(X(length(X(:,1)),1)/mcinit)*100;
%figure(1)
%hold on
%plot(tfinal1,FinalPercentage)
%xlabel('Stem cells injection time (days)'), ylabel('Cardiomyocytes recovery level (%)')

TT=[TT tfinal1];
PP=[PP FinalPercentage];
end

figure(2) 
plot(TT,PP)
xlabel('Stem cells injection time (days)'), ylabel('Cardiomyocytes recovery level (%)');
hold on


end
end

legend('ORT 2 hrs', 'ORT 2 hrs with SCI', 'ORT 4 hrs', 'ORT 4 hrs with SCI', 'ORT 6 hrs','ORT 6 hrs with SCI', 'ORT 12 hrs', 'ORT 12 hrs with SCI', 'ORT 24 hrs', 'ORT 24 hrs with SCI')