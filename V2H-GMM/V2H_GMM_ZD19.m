function [Sv2h,sigmat,f]=V2H_GMM_ZD19(Rtype,func,m,r,ord,pdc,sc,mech,k)
% Usage:  [Sv2h,sigmat]=V2H_GMM_ZD19(Rtype,functionType,m, r, ord, pdc, sc, mech, k)
% this script implements the ground-motion models for predicting 5%-damped V/H 
% spectral acceleration in Iran as functions of moment magnitude, 
% distance, site condition and style of faulting. 
% 
% The general form of the functional form is as follows:
% for Rjb and Repi : 
%   logSv2h(T)=c1+m1.M+m2.M.^2 + r1.*(log10(sqrt(R.^2+h.^2)))
%                  +sCI.s1+sCII.s2+sCIII.s3
%                         + sR.SoF1+sS.SoF2
% for Rrup and Rhyp : 
%   logSv2h(T)=c1+m1.M+m2.M.^2 + r1.*(log10(R+h))
%                  +sI.s1+sII.s2+sIII.s3
%                         + sR.SoF1+sS.SoF2
% Alternative functional forms: 
% fun1-> excluding magnitude-dependent geometric attenuation term (r2=0 in (r1 +r2M))
% fun3-> m2=0 & r2=0
%
% Input:    Rtype: define distance type of interest: 'Rjb','Repi','Rhyp','Rrup' 
%           func: define functional form: 'fun3' , 'fun1' 
%           m :     Vector of Moment Magnitude
%           r :     distance (km) 
%           ord : index of period vector: [-2 -1 Vector of period up to 10 (s)]
%           pdc: load a text file including the period-dependent coefficients, 
%             choose based on desired distance type and functional form 
%             e.g.,IrV2H-Rrup-fun3.txt corresponding to Rrup and fun3 functional form 
%           sc : Soil condition of site that categorized based on the 
%               shear wave velocity (see Darzi et al. 2019) 
%             'I':VS30>750 m/s, 'II':(375<VS30<750 m/s), 'III': Vs30<375 m/s
%           mech: 'R' for reverse and thrust,  'S' for strike-slip
%           k :     term that represent standard deviation in calculation
%                  
% Output:   Sv2h : 5%-damped vertical -to-horizontal (V2H) spectral acceleration- PGV-PGD 
%           sigmat : Standard deviation 
% Note: The proposed attenuation relationship provide the best result for
%      distance between 4-200 km Moment magnitude between 4.5-7.5
% 
% Author: Atefe Darzi, PhD 
% Last Edit: 2019

% Reference:
% Zolfaghari MR, Darzi A (2019). A Prediction Model for Vertical-to-Horizontal 
% Ratios of PGA, PGV and 5%-damped Response Spectra (0.01-10 s) for Iran. 
% Journal of Seismology. 10.1007/s10950-019-09836-z 
% *************************************************************************    

%% define site coefficients
s1=0;s2=0;
switch sc
    case 'I'
       s1=0; 
    case 'II'
       s1=1;
    case 'III'
       s2=1;
end

%% define SOF
SoF1=0; SoF2=0;
switch mech
    case 'R'
        SoF1=1;
    case 'S'
        SoF2=1;
end

Tnga = load('Periods-nga.txt'); % T=0.01-20
per = [-2 -1 0 Tnga(1:60)']; % PGD,PGV,PGA, SA(0.01-10 s)

%% define coefficients at each period 
for i=1:length(ord)
 if strcmp(func,'fun1')
 % Coef={'fRV';'fSS';'C1';'m1';'m2';'h';'r1';'S(II)';'S(III-IV)';'phi';'tau';'sigma'};
sR = pdc(1,ord(i));sS = pdc(2,ord(i)); 
c1 = pdc(3,ord(i)); m1 = pdc(4,ord(i)); m2 = pdc(5,ord(i)); 
h = pdc(6,ord(i)); r1 = pdc(7,ord(i));
sII = pdc(8,ord(i)); sIII = pdc(9,ord(i)); 
phi = pdc(10,ord(i)); tau = pdc(11,ord(i)); sigma = pdc(12,ord(i)); 
     
%% define V2H-fun1 
    if length(ord)==1
        for j=1:length(m)
            for w=1:length(r)
                if strcmp(Rtype,'Rjb') | strcmp(Rtype,'Repi')
               Sv2h(j,w)=10.^(c1+m1.*m(j)+m2.*m(j).^2 + r1.*(log10(sqrt(r(w).^2+h.^2)))...
                    +sII.*s1+sIII.*s2 ...
                       +sR.*SoF1+sS.*SoF2+k.*sigma);
                elseif strcmp(Rtype,'Rhyp') | strcmp(Rtype,'Rrup') 
               Sv2h(j,w)=10.^(c1+m1.*m(j)+m2.*m(j).^2 + r1.*(log10(r(w)+h))...
                    +sII.*s1+sIII.*s2 ...
                       +sR.*SoF1+sS.*SoF2+k.*sigma);
                end
               sigmat(j,w,i)=sigma;
            end
        end
    else 
      if length(r)==1
        for j=1:length(m) 
            if strcmp(Rtype,'Rjb') | strcmp(Rtype,'Repi')
            Sv2h(i,j)=10.^(c1+m1*m(j)+m2.*m(j).^2 + r1.*(log10(sqrt(r.^2+h.^2)))...
                 +sII.*s1+sIII.*s2 ......
                        + sR.*SoF1+sS.*SoF2+k.*sigma);
            elseif strcmp(Rtype,'Rhyp') | strcmp(Rtype,'Rrup')
            Sv2h(i,j)=10.^(c1+m1*m(j)+m2.*m(j).^2 + r1.*(log10(r+h))...
                 +sII.*s1+sIII.*s2 ......
                        + sR.*SoF1+sS.*SoF2+k.*sigma);                
            end
             period(i,j)=per(ord(i));
             sigmat(i,j)=sigma;
        end
     else
            for w=1:length(r) 
            if strcmp(Rtype,'Rjb') | strcmp(Rtype,'Repi')                
                Sv2h(i,w)=10.^(c1+m1*m+m2.*m.^2 + r1.*(log10(sqrt(r(w).^2+h.^2)))...
                     +sII.*s1+sIII.*s2 ......
                            + sR.*SoF1+sS.*SoF2+k.*sigma);
            elseif strcmp(Rtype,'Rhyp') | strcmp(Rtype,'Rrup')
                Sv2h(i,w)=10.^(c1+m1*m+m2.*m.^2 + r1.*(log10(r(w)+h))...
                     +sII.*s1+sIII.*s2 ......
                            + sR.*SoF1+sS.*SoF2+k.*sigma);
            end
                 period(i,w)=per(ord(i));
                 sigmat(i,w)=sigma;
            end
      end
    end
 end

if strcmp(func,'fun3')
 % Coef={'fRV';'fSS';'C1';'m1';'h';'r1';'S(II)';'S(III-IV)';'phi';'tau';'sigma'};
sR = pdc(1,ord(i));sS = pdc(2,ord(i)); 
c1 = pdc(3,ord(i)); m1 = pdc(4,ord(i)); 
h = pdc(5,ord(i)); r1 = pdc(6,ord(i));
sII = pdc(7,ord(i)); sIII = pdc(8,ord(i));
phi = pdc(9,ord(i)); tau = pdc(10,ord(i)); sigma = pdc(11,ord(i)); 
     
%% define V2H-fun3 
    if length(ord)==1
        for j=1:length(m)
            for w=1:length(r)
                if strcmp(Rtype,'Rjb') | strcmp(Rtype,'Repi')
               Sv2h(j,w)=10.^(c1+m1.*m(j)+ r1.*(log10(sqrt(r(w).^2+h.^2)))...
                   +sII.*s1+sIII.*s2 ......
                       +sR.*SoF1+sS.*SoF2+k.*sigma);
                elseif strcmp(Rtype,'Rhyp') | strcmp(Rtype,'Rrup') 
               Sv2h(j,w)=10.^(c1+m1.*m(j)+ r1.*(log10(r(w)+h))...
                    +sII.*s1+sIII.*s2 ......
                       +sR.*SoF1+sS.*SoF2+k.*sigma);
                end
               sigmat(j,w,i)=sigma;
            end
        end
    else 
      if length(r)==1
        for j=1:length(m) 
            if strcmp(Rtype,'Rjb') | strcmp(Rtype,'Repi')
            Sv2h(i,j)=10.^(c1+m1*m(j) + r1.*(log10(sqrt(r.^2+h.^2)))...
                 +sII.*s1+sIII.*s2 ......
                        + sR.*SoF1+sS.*SoF2+k.*sigma);
            elseif strcmp(Rtype,'Rhyp') | strcmp(Rtype,'Rrup')
            Sv2h(i,j)=10.^(c1+m1*m(j) + r1.*(log10(r+h))...
                 +sII.*s1+sIII.*s2 ......
                        + sR.*SoF1+sS.*SoF2+k.*sigma);                
            end
             period(i,j)=per(ord(i));
             sigmat(i,j)=sigma;
        end
     else
            for w=1:length(r) 
            if strcmp(Rtype,'Rjb') | strcmp(Rtype,'Repi')                
                Sv2h(i,w)=10.^(c1+m1*m + r1.*(log10(sqrt(r(w).^2+h.^2)))...
                     +sII.*s1+sIII.*s2 ......
                            + sR.*SoF1+sS.*SoF2+k.*sigma);
            elseif strcmp(Rtype,'Rhyp') | strcmp(Rtype,'Rrup')
                Sv2h(i,w)=10.^(c1+m1*m + r1.*(log10(r(w)+h))...
                     +sII.*s1+sIII.*s2 ......
                            + sR.*SoF1+sS.*SoF2+k.*sigma);
            end
                 period(i,w)=per(ord(i));
                 sigmat(i,w)=sigma;
            end
      end
    end
end
 

end


end     
        
        
        
        
        
        
        

