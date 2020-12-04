#ifndef EMBOLEDETECTION_H_INCLUDED
#define EMBOLEDETECTION_H_INCLUDED

void pass1Detection(const double* inData, const int Nmax, const int fs, double fc, double* outMpvApriori);//, double* outMedApriori
void calculESub(const double* data, int Nmax, int lF, int Ncp, int Ncn, int delay, int lWind, double *ESub1, double *ESub2);
void pass2Detection(const double* data, const int Nmax, const int fs, const double Tobs, const double fc, const int tolerance, const double *mpvApriori, int *NEmbFinal, double *EmbPosFinal, double *EmbEnergFinal, int *EmbDuration, int *NArts);
void Energies_L1(const double* data, const int Nmax, const int delay, const int lWind, const int Ncp, const int Ncn, double* ESubNorm, double *ESubNormF, double *VarThreshold, double *std_flux);
void Quality(const double *ESubNormF, int lEsub, double *DynamicSNR, double *Mpv_pos);
void HITS_detection(const double *ESubNorm, const double *ESubNormF, const int lEsub, const double *VarThreshold, const double alphaH, const double *Mpv_pos, double *ThresholdHits, double *SimilarityThresholdRatio, int *PositionHits, int *DurationHits, int *ZHits, int *LHits);
void DetectionArts(const double *ESubNorm, const double *ESubNormF, const int lEsub, const double *Mpv_apriori, double *ThresholdArts, int *PositionArts, int *DurationArts, int *ZArts, double *EnergyArts, int *LArts);
void EmbDetector(const double *ESubNorm, const double *ESubNormF, const int lEsub, const int *PositionHits, const int *DurationHits, const int *ZHits, const int *lHits, const int *PositionArts, const int *ZArts, const int *lArts, const int Tolerance, int *NEmb, int* EmbPos, double *EmbEner, int *EmbDur);
void detectionEmbolus(const double* dataWave, const int Nmax, const int fs, const double Tobs, const double fc, const int tolerance, int *NEmbFinal, double *EmbPosFinal, double *EmbEnergFinal, int *EmbDuration, int *NArts);
#endif // EMBOLEDETECTION_H_INCLUDED
