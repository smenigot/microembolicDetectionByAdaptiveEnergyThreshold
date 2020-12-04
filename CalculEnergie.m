function ESub = CalculEnergie(X, Nmax, Lwind, OverLap, fcArt, Fs)
% Compute short term energy
% X signal
% Nmax sampling amount of the analysis
% Lwind window size
% Overlap overlap
% fcArt cutoff frequency to remove artefact
% Fs sampling frequency

Lf      = 256;
Delay   = round(Lwind-(Lwind*OverLap));

Ncp   = round(Fs/2/fcArt);  % 12 points
Ncn   = round(Fs/2/fcArt);
%%
l=1;
ESub = zeros(round((Nmax-Lwind)/Delay)+1,2);
while (l-1)*Delay+Lwind<Nmax
    x_wind = X((l-1)*Delay+1:(l-1)*Delay+Lwind,:);
    S=abs(fft(x_wind,Lf));
    
    %% Low-Pass filter to remove artefacts
    Subband11=Ncp:length(S)/2;
    ESub(l,1)=sum((S(Subband11).*hamming(length(Subband11))).^2)/length(Subband11);
    
    Subband12=length(S)/2+1 : length(S)-Ncn;
    ESub(l,2)=sum((S(Subband12).*hamming(length(Subband12))).^2)/length(Subband12);
    l=l+1;
end