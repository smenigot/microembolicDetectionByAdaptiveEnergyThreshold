
function [PositionHitsFinal, DurationHitsFinal, ZHits] = ...
    Detection3dB(ESubNorm,NewThresholdHits)

ESubNorm(:,1)=smooth(ESubNorm(:,1),10,'mean');%nouveau

Z = zeros(length(ESubNorm),1);
Z(ESubNorm(:,1)>NewThresholdHits) = 1;

dZ           = [0; diff(Z)];
c1           = find(dZ==1);
PositionHits = c1;
%EnergyHits   = ESubNorm(c1,1)-ESubNormF(c1,1);

Pp = find(dZ == 1);
Pn = find(dZ == -1);
if (length(Pp) == length(Pn))
    DurationHits = abs(Pn-Pp);
else
    if (length(Pp) > length(Pn))
        Pn = [Pn; length(Z)];
        DurationHits = abs(Pn-Pp);
    else
        Pp = [Pp; length(Z)];
        DurationHits = abs(Pn-Pp);
    end
end

a                 = find(DurationHits<15);
PositionHitsFinal = PositionHits(a);
% EnergyHitsFinal   = EnergyHits(a);
DurationHitsFinal = DurationHits(a);

ZHits = zeros(length(ESubNorm),1);
for k = 1:length(a)
    ZHits(PositionHitsFinal(k):PositionHitsFinal(k)+DurationHitsFinal(k),1) = 1;
end
