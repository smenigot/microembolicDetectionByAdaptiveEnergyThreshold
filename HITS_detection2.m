function [ThresholdHits,SimilarityThresholdRatio,PositionHits,DurationHits,ZHits] =...
    HITS_detection2(ESubNorm, ESubNormF, VarThreshold, std_flux, Mpv_pos)

ThresholdHits(:,1) = ESubNormF(:,1) + std_flux;%nouveau
ThresholdHits(:,2) = mean(ESubNormF(:,1)) + 1.5*std_flux;

SimilarityThresholdRatio = 100-100*abs(mean(ThresholdHits(:,1))-ThresholdHits(1,2))/ThresholdHits(1,2);

if (SimilarityThresholdRatio<0)
    ThresholdHits(:,3) = ThresholdHits(1,2)*ones(1,length(ThresholdHits));
else
    ThresholdHits(:,3) = ThresholdHits(:,1);
end

[PositionHits,DurationHits,ZHits] = Detection3dB(ESubNorm,ThresholdHits(:,3));