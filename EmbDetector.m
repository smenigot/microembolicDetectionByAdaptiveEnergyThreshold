
function [NEmb, EmbPos, EmbEner, EmbDur]=EmbDetector(PositionHits, PositionArts,DurationHits,DurationArts,ZHits,ZArts,ESubNorm,ESubNormF,Tolerance);
EmbPos=[];
EmbEner=[];
EmbDur=[];
Com = and(ZHits,ZArts);
% Tolerance=2*Tolerance;
% Com    = conv( double( and( ...
%     conv(ZHits,ones(1,Tolerance),'same')>0,conv(ZArts,ones(1,Tolerance),'same')>0)),...
%     ones(1,1),'same')>0;
IndCom = (find(Com==1))';

for k = 1:length(PositionHits);
    RHits = (repmat(PositionHits(k),1,length(PositionArts)))';
    AA = find(abs(PositionArts-RHits)<2*Tolerance,1);
    RHits2 = (repmat(PositionHits(k),length(IndCom),1))';
    BB = find(abs(IndCom-RHits2)<2*Tolerance,1);%find(IndCom-repmat(PositionHits(k),1,length(IndCom))==0,1);
    %CC = find(and((RHits<PositionArts-Tolerance),(RHits>PositionArts+DurationArts(1:length(PositionArts))+Tolerance))==1)
    if ( (isempty(AA)==1) && (isempty(BB)==1) )
        EmbPos  = [EmbPos, PositionHits(k)];
        EmbEner = [EmbEner, ESubNorm(PositionHits(k))/ESubNormF(PositionHits(k))];
    end
end
NEmb=length(EmbPos);

for k = 1:NEmb;
    VarInt = (repmat(EmbPos(k),1,length(PositionHits)))';
    [a] = find((PositionHits-VarInt)==0);
    EmbDur = [EmbDur DurationHits(a)];
end

EmbEner = EmbEner';
EmbPos = EmbPos';
EmbDur = EmbDur';





