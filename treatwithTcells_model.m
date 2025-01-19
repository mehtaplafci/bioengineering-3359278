function dx = treatwithTcells_model(t,x)
global  k12 k11 k10 k5 cIL1M0 k6 cIL10M0 k4 c2 cIL10C cIL1C cIL10F dIL10 dIL1 k2 k3 c1 k7 k8 k9 
global  yMc mMc muM yF mF muMc muS k14  
global  k13 muT muRT cRT k1 dR k22 dIL17 k18 k19  muN muRN cRN k15 cIL17Mc k16 muRMc cRMc k17 k21 k20 cIL17 cIL17M0 muRM0 cRM0 muRM1 cRM1 muRM2 cRM2 muRMS cRS omega gamma delta c3

 Rc=200;  

if x(12)<Rc 
    gamma=0;
end
%gamma
dx = zeros(14,1);
dx(1)=yMc*x(10)*x(5)/(mMc+x(5))-muMc*x(1)-k15*x(1)*x(13)/(x(13)+cIL17Mc)-k16*x(1)*x(14)-gamma*muRMc*x(1)*x(12)/(x(12)+cRMc);  %1 Mc Cardiomyocytes
dx(2)=k6*x(9)-k7*x(2)*x(6) /(x(6)+cIL1M0)-k8*x(2)*x(5)/(x(5)+cIL10M0)-k20*x(2)*x(13)/(x(13)+cIL17M0)-muM*x(2)-gamma*muRM0*x(2)*x(12)/(x(12)+cRM0); %2 M0 Monocytes
dx(3)=k7*x(2)*x(6)/(x(6)+cIL1M0)+k20*x(2)*x(13)/(x(13)+cIL17M0)-k9*x(3)-muM*x(3)-gamma*muRM1*x(3)*x(12)/(x(12)+cRM1); %3 M1 Classically Activated Macrophages (Cleaners)
dx(4)=k8*x(2)*x(5)/(x(5)+cIL10M0)+k9*x(3)-muM*x(4)-gamma*muRM2*x(4)*x(12)/(x(12)+cRM2); %4 M2 Alternatively Activated Macrophages (Builders)
dx(5)=k21*x(11)*c3/(c3+x(5))+k5*x(4)*c2/(c2+x(5))-dIL10*x(5); %5 IL10 (good - synthesis)
dx(6)=k3*x(9)+k4*x(3)*c1 /(c1+x(5))-dIL1*x(6); %6 IL1 (bad - degredation)
dx(7)=k10*x(8)*x(5)/(x(5)+cIL10C)-k11*x(7)*x(6)/(x(6)+cIL1C); %7 C Collagen
dx(8)=yF*x(10)*x(5)/(mF+x(5))+k12*x(8)*x(5)/(x(5)+cIL10F); %8 F Fibroblasts
dx(9)=k15*x(1)*x(13)/(x(13)+cIL17Mc)+muMc*x(1)+gamma*muRMc*x(1)*x(12)/(x(12)+cRMc)-k2*x(3)*x(9)-k17*x(14)*x(9); %9 Md Dead Cardiomyocytes
dx(10)=-yMc*x(10)*x(5)/(mMc+x(5))-yF*x(10)*x(5) /(mF+x(5))-muS*x(10)-k14*x(10)*x(3)-gamma*muRMS*x(10)*x(12)/(x(12)+cRS); %10 Stem Cells
dx(11)=k13*x(9)-muT*x(11)-gamma*muRT*x(11)*x(12)/(x(12)+cRT); %11 T Cells
dx(12)=omega*k1*x(14)-dR*x(12);   %12 ROS
dx(13)=k22*x(11)*c3/(c3+x(5))-dIL17*x(13);  %13 IL17
dx(14)=delta*k18*x(9)+k19*x(13)/(x(13)+cIL17)-muN*x(14)-gamma*muRN*x(14)*x(12)/(x(12)+cRN);  %14 Neutrophil