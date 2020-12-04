
function [PositionArts,ThresholdArts,DurationArts,ZArts,EnergyArts]=DetectionArts(ESubNorm,ESubNormF,Mpv_apriori)

alphaA=0.8*2;%nouveau
ESubSmooth=smooth(ESubNorm(:,2),2,'mean');%nouveau

ThresholdArts=ESubNormF(:,2) + alphaA*Mpv_apriori*ones(length(ESubNormF),1);%nouveau

ZArts                              = zeros(length(ESubSmooth),1);
ZArts(ESubSmooth>ThresholdArts) = 1;

dZ = [0; diff(ZArts)];
PositionArts = find(dZ==1);
EnergyArts   = ESubSmooth(PositionArts) - ESubNormF(PositionArts,2);

Pp = find(dZ==1);
Pn = find(dZ==-1);
if (length(Pp) == length(Pn))
    DurationArts = abs(Pn-Pp);
else
    if (length(Pp)>length(Pn))
        Pn = [Pn; length(ZArts)];
        DurationArts = abs(Pn-Pp);
    else
        Pp = [Pp; length(ZArts)];
        DurationArts = abs(Pn-Pp);
    end
end



