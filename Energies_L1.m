function [ESubNorm,ESubNormF,VarThreshold,std_flux]=Energies_L1(X, NObs, Lwind, OverLap, fcArt, Fs)
% Information used for detection (Energy)

%% Energy Computation
ESub = CalculEnergie(X, NObs, Lwind, OverLap, fcArt, Fs);

%% substraction of the Energy of negative frequency
ESubNorm(:,1)           = ESub(:,1)-ESub(:,2);
ESubNorm(ESubNorm<0,1)  = 0;
ESubNorm(:,2)           = ESub(:,2);

vSmooth = 25;
ESubNormF(:,1)=smooth(ESubNorm(:,1)',vSmooth,'mean');
ESubNormF(:,2)=smooth(ESubNorm(:,2)',vSmooth,'mean');

% [n1,bin1] = hist(ESubNormF(:,1),length(ESubNorm));
% [~,b1]    = max(n1(2:end));
% Mpv(1)    = bin1(b1+2);

[n2,bin2] = hist(ESubNormF(:,2),length(ESubNorm));
[~,b2]    = max(n2(10:end));
if b2+10>length(bin2)
    Mpv = bin2(b2);
else
    Mpv       = bin2(b2+10);
end

%%
Sat(1)  = mean(ESubNormF(:,1))+5*std(ESubNormF(:,1));
Sat(2)  = 10*Mpv;

ESubNormLim(:,1)                    = ESubNorm(:,1);
ESubNormLim(ESubNorm(:,1)>Sat(1),1) = Sat(1);
ESubNormF(:,1) = smooth(ESubNormLim(:,1)',vSmooth,'mean');

ESubNormLim(:,2)                    = ESubNorm(:,2);
ESubNormLim(ESubNorm(:,2)>Sat(2),2) = Sat(2);
ESubNormF(:,2) = smooth(ESubNormLim(:,2)',vSmooth,'mean');

%%
Fc = round((32/length(ESubNormF)*Fs)*4096/Fs);
Filtre        = [ones(1,Fc) zeros(1,4096-2*Fc) ones(1,Fc)];

Spectre1(1,:)  = fft(ESubNormF(:,1),4096);
temp           = real(ifft(Spectre1.*Filtre));
ESubNormF(:,1) = temp(1:length(ESubNormF));

Spectre2(1,:)  = fft(ESubNormF(:,2),4096);
temp           = real(ifft(Spectre2.*Filtre));
ESubNormF(:,2) = temp(1:length(ESubNormF));

STD(1)            = std(ESubNormF(:,1));
StdBruit(1)       = std(ESubNormLim(:,1)-ESubNormF(:,1));
VarThreshold(:,1) = 2*(2*STD(1)+StdBruit(1));

STD(2)            = std(ESubNormF(:,2));
StdBruit(2)       = std(ESubNormLim(:,2)-ESubNormF(:,2));
VarThreshold(:,2) = 2*(2*STD(2)+StdBruit(2));
%% new
D=ESubNorm-ESubNormF;

[a11,b11] = hist(D(D(:,1)<0,1),100);
c11       = find( cumsum(a11/sum(a11))/2 < 0.1/100);
if (isempty(c11) == 1)
    std_flux(1) = abs(b11(1));
else
    std_flux(1) = abs(b11(c11(end)));
end

[a22,b22] = hist(D(D(:,2)<0,2),100);
c22       = find(cumsum(a22/sum(a22))/2<0.1/100);
if (isempty(c22) == 1) 
    std_flux(2) = abs(b22(1));
else
    std_flux(2) = abs(b22(c22(end))); 
end






