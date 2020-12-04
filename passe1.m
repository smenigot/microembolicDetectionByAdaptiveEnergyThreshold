function [Mpv_apriori, Med_apriori] = passe1(X, Lwind, OverLap, fcArt, Fs)
% Pass 1 of the detection process
% X signal
% Lwind window size
% Overlap overlap
% fcArt cutoff frequency to remove artefact
% Fs sampling frequency

ESub = CalculEnergie(X, length(X), Lwind, OverLap, fcArt, Fs);

[n1,bin1]      = hist(ESub(:,1),10000);
[~,b1]         = max(n1);
Mpv_apriori(1) = bin1(b1);
Med_apriori(1) = median(ESub(:,1));

[n2,bin2]      = hist(ESub(:,2),10000);
[~,b2]         = max(n2);
Mpv_apriori(2) = bin2(b2);
Med_apriori(2) = median(ESub(:,2));
