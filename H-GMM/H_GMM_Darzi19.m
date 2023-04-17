function [Sa,sigmat]=H_GMM_Darzi19(m,r,ord,pdc,sc,mech,k)
% Usage:  [Sa,sigmat]=H_GMM_Darzi19(m,r,ord,pdc,sc,mech,k)
% An Empirical Ground-Motion Model for Horizontal PGV, PGA, and
% 5% Damped Elastic Response Spectra (0.01–10 s) in Iran
% Bulletin of the Seismological Society of America-April, (2019)

% The general form of their relation is as follows:
%
%   logSa(T)=c1+m1.M+m2.M.^2 + r1.*(log10(sqrt(R.^2+h.^2)))
%                   +sCI.s1+sCII.s2+sCIII.s3
%                         + sR.SoF1+sS.SoF2
% 
% applicable to 'Rjb','Repi','Rhyp','Rrup' distance metrics
% Input:   
%           m :     Vector of Moment Magnitude
%           r :     distance (km) 
%           ord : index of period vector: [-2 -1 Vector of period up to 10 (s)]
%           pdc: load a text file including the period-dependent coefficients, 
%             choose amongst IrH-jb.txt, IrH-rup, IrH-epi, IrH-hyp
%           sc :    Soil condition of site that categorized based on the 
%               shear wave velocity (see Darzi et al. 2019) 
%             'I':VS30>750 m/s, 'II':(375<VS30<750 m/s), 'III': Vs30<375 m/s
%           mech: 'R' for reverse and thrust,  'S' for strike-slip
%           k :     term that represent standard deviation in calculation
% Output:   sa :  5%-damped horizontal spectral acceleration in (cm/s^2)
%           sigmat : Standard deviation of the result
% Note: The proposed attenuation relationship provide the best result for
%      distance between 4-200 km Moment magnitude between 4.5-7.5
% 
% 
% Author: Atefe Darzi, KNTU, visiting PhD candidate at SED (ETH Zurich)  
% Last Edit: 2017

% Reference:
% Darzi A, Zolfaghari MR, Cauzzi C, Fäh D (2019). An Empirical Ground Motion Model 
% for Horizontal PGV, PGA and 5%-Damped Elastic Response Spectra (0.01-10 s) in Iran. 
% Bulletin of the Seismological Society of America. 10.1785/0120180196 

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
SoF1=0;SoF2=0;
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
% Coef={'fRV';'fSS';'C1';'m1';'m2';'h';'r1';'S(II)';'S(III-IV)';'phi';'tau';'sigma'};

sR = pdc(1,ord(i));sS = pdc(2,ord(i)); 
c1 = pdc(3,ord(i)); m1 = pdc(4,ord(i)); m2 = pdc(5,ord(i)); 
h = pdc(6,ord(i)); r1 = pdc(7,ord(i));
sCII = pdc(8,ord(i)); sCIII = pdc(9,ord(i));
phi = pdc(10,ord(i)); tau = pdc(11,ord(i)); sigma = pdc(12,ord(i)); 
     
%% define SA 
    if length(ord)==1
        for j=1:length(m)
            for w=1:length(r)
              Sa(j,w)=10.^(c1+m1.*m(j)+m2.*m(j).^2 + r1.*(log10(sqrt(r(w).^2+h.^2)))...
                    +sCII.*s1+sCIII.*s2...
                       +sR.*SoF1+sS.*SoF2+k.*sigma);
              sigmat(j,w,i)=sigma;
            end
        end
    else 
      if length(r)==1
        for j=1:length(m) 
              Sa(i,j)=10.^(c1+m1*m(j)+m2.*m(j).^2 + r1.*(log10(sqrt(r.^2+h.^2)))...
                   +sCII.*s1+sCIII.*s2...
                        + sR.*SoF1+sS.*SoF2+k.*sigma);
             period(i,j)=per(ord(i));
             sigmat(i,j)=sigma;
        end
     else
            for w=1:length(r) 
                Sa(i,w)=10.^(c1+m1*m+m2.*m.^2 + r1.*(log10(sqrt(r(w).^2+h.^2)))...
                    +sCII.*s1+sCIII.*s2...
                            + sR.*SoF1+sS.*SoF2+k.*sigma);
               period(i,w)=per(ord(i));
                 sigmat(i,w)=sigma;
            end
      end
    end
end
 
end     
        
        
        
        
        
        
        

