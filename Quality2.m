function [DynamicSNR, Mpv_pos] = Quality2(ESubNorm, ESubNormF)

EEm = mean(ESubNormF(:,1));
EEc = ESubNormF(:,1)-EEm;

[indmin, indmax] = extr(EEc);

a = find(EEc(indmax)>0);

[n1,bin1]  = hist(EEc(indmax(a)),length(a));
[~,b1]     = max(n1(1:end));
Mpv_pos(1) = bin1(b1);

b = find(EEc(indmin)<0);
[n2,bin2]  = hist(EEc(indmin(b)),length(b));
[~,b2]     = max(n2(1:end));
Mpv_pos(2) = bin2(b2);

DynamicSNR = 10*log10( abs((Mpv_pos(1)-Mpv_pos(2))/mean( ESubNorm(:,1)-ESubNormF(:,1) )) );