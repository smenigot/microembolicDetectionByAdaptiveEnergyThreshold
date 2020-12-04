#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <string.h>
#include <omp.h>

#include "mathTools.h"
#include "emboleDetection.h"
#include "fftw3.h"

#define FREE(_p) do{ \
	 free( (_p) ); \
     (_p) = NULL; \
}while( (_p) != NULL)

#define FFTW_FREE(_p) do{ \
	 fftw_free( (_p) ); \
     (_p) = NULL; \
}while( (_p) != NULL)

#define FFTW_DESTROY_PLAN(_p) do{ \
	 fftw_destroy_plan( (_p) ); \
     (_p) = NULL; \
}while( (_p) != NULL)

// ---------------------------------------------------------------------------------------------------------------------------
// ---------------------------------------------------------------------------------------------------------------------------
// ---------------------------------------------------------------------------------------------------------------------------
// Passe 1 pour determiner les statistiques du signal complet
// ---------------------------------------------------------------------------------------------------------------------------
// ---------------------------------------------------------------------------------------------------------------------------
// ---------------------------------------------------------------------------------------------------------------------------
void pass1Detection(const double* data, const int Nmax, const int fs, double fc, double* mpvApriori)//, double* medApriori
{
    //Nmax: header->bytes_in_data
    //fs: header->frequency

    // declaration de variables
    int lWind=64;
    double overLap=20./100.;
    int delay=(int)round(lWind-(lWind*overLap));
    int lF=256;
    int Ncp=round((double)fs/2./fc), Ncn=Ncp;// selection de sousbandes // 206 Hz-> 12 samples
    int Nclasse=10000;
    int a1, b1, a2, b2;

    double *ESub1=NULL, *ESub2=NULL, *edges=NULL;
    int *histo=NULL;

    // affichage
    printf("Pass 1 : Extraction of a priori information\n");
    do
    {
        ESub1=(double*) malloc(sizeof(double) *(int)round((double)Nmax/(double)delay));
    }
    while(ESub1 == NULL);
    do
    {
        ESub2=(double*) malloc(sizeof(double) *(int)round((double)Nmax/(double)delay));
    }
    while(ESub2 == NULL);


    // calcul de l'energie a court terme
    calculESub(data, Nmax, lF, Ncp, Ncn, delay, lWind, ESub1, ESub2);

    // Histogramme pour la statistique
    do
    {
        histo = (int*) malloc(Nclasse *sizeof(int));
    }
    while(histo == NULL);
    do
    {
        edges = (double*) malloc(Nclasse *sizeof(double));
    }
    while(edges == NULL);

    fhist(ESub1, (int)round((double)Nmax/(double)delay), Nclasse, histo, edges);
    maximum(histo, Nclasse, &a1, &b1);
    mpvApriori[0]=edges[b1];
    // mediane(ESub1, (int)round((double)Nmax/(double)delay), medApriori);

    fhist(ESub2, (int)round((double)Nmax/(double)delay), Nclasse, histo, edges);
    maximum(histo, Nclasse, &a2, &b2);
    mpvApriori[1]=edges[b2];
    // mediane(ESub2, (int)round((double)Nmax/(double)delay), medApriori+1);

    FREE(ESub1);
    FREE(ESub2);
    FREE(histo);
    FREE(edges);
}

// ---------------------------------------------------------------------------------------------------------------------------
// ---------------------------------------------------------------------------------------------------------------------------
// ---------------------------------------------------------------------------------------------------------------------------
// Calcul de l'energie a court terme
// ---------------------------------------------------------------------------------------------------------------------------
// ---------------------------------------------------------------------------------------------------------------------------
// ---------------------------------------------------------------------------------------------------------------------------
void calculESub(const double* data, int Nmax, int lF, int Ncp, int Ncn, int delay, int lWind, double *ESub1, double *ESub2)
{
    int k=0, l=1;
    fftw_complex *in=NULL, *out=NULL;
    fftw_plan planFft=NULL;
    double *S=NULL;

    // definition de la fft
    do
    {
        in =(fftw_complex*) fftw_malloc(sizeof(fftw_complex) *lF);
    }
    while(in == NULL);
    do
    {
        out=(fftw_complex*) fftw_malloc(sizeof(fftw_complex) *lF);
    }
    while(out == NULL);
    do
    {
        planFft=fftw_plan_dft_1d(lF, in, out, FFTW_FORWARD, FFTW_ESTIMATE);
    }
    while(planFft == NULL);

    do
    {
        S=(double*) malloc(sizeof(double) *lF);
    }
    while(S == NULL);

    // passe 1
    while ((l-1)*delay+lWind<Nmax)
    {
        for (k=0; k<lF; k++)
        {
            //c + bytes_in_data*r
            if (k<lWind)
            {
                in[k][0]=data[(l-1)*delay+k]; // partie reelle
                in[k][1]=data[(l-1)*delay+k+Nmax*1]; // partie imaginaire
            }
            else
            {
                in[k][0]=0; // partie reelle
                in[k][1]=0; // partie imaginaire
            }
        }
        fftw_execute(planFft); //fft
        for (k=0; k<lF; k++)
        {
            S[k]=sqrt(pow(out[k][0],2) + pow(out[k][1],2)); // energie dans le domaine des frequences
        }

        // Energies de la partie positive
        ESub1[(l-1)]=0;
        for (k=0; k<(int)((double)lF/2.)-Ncp+1; k++)
        {
            ESub1[(l-1)]=ESub1[(l-1)] + pow(S[k+Ncp-1] * hamming(k, (int)((double)lF/2.)-Ncp),2);
        }
        ESub1[(l-1)] /= (double)lF/2.-(double)Ncp+1.;

        // Energies de la partie negative
        ESub2[(l-1)]=0;
        for (k=0; k<(int)((double)lF/2.)-Ncn; k++)
        {
            ESub2[(l-1)]=ESub2[(l-1)] + pow( S[k+(int)((double)lF/2.)-1] * hamming(k, (int)((double)lF/2.)-Ncn),2);
        }
        ESub2[(l-1)] /= (double)lF/2.-(double)Ncn;

        l++;
    }

    FFTW_FREE(in);
    FFTW_FREE(out);
    FFTW_DESTROY_PLAN(planFft);
    FREE(S);
}

// ---------------------------------------------------------------------------------------------------------------------------
// ---------------------------------------------------------------------------------------------------------------------------
// ---------------------------------------------------------------------------------------------------------------------------
// Passe 2 : Detection par morceaux de Tobs secondes
// ---------------------------------------------------------------------------------------------------------------------------
// ---------------------------------------------------------------------------------------------------------------------------
// ---------------------------------------------------------------------------------------------------------------------------
void pass2Detection(const double* data, const int Nmax, const int fs, const double Tobs, const double fc, const int tolerance, const double *mpvApriori, int *NEmbFinal, double *EmbPosFinal, double *EmbEnergFinal, int *EmbDuration, int *NArts)
{
    //Nmax: header->bytes_in_data
    //fs: header->frequency

    int k1=0,k2=0;
    int nObs=(int)round(Tobs)*fs;
    double *DynamicSNR=NULL, *SimilarityThresholdRatio=NULL;

    int Ncp=round((double)fs/2./fc), Ncn=Ncp;// selection de sousbandes // 206 Hz-> 12 samples
    int lWind=64;
    double overLap=80./100.;
    int delay=(int)round((double)lWind-((double)lWind*overLap));
    int lEsub=(int)round((double)(nObs-lWind)/(double)delay)+1;

    int *NEmbTab=NULL, **EmbPosTab=NULL, **EmbDurationTab=NULL;
    double **EmbEnergTab=NULL;

    do
    {
        DynamicSNR=(double*)malloc((int)((double)Nmax/(double)nObs)*sizeof(double));
    }
    while(DynamicSNR == NULL);
    do
    {
        SimilarityThresholdRatio=(double*)malloc((int)((double)Nmax/(double)nObs)*sizeof(double));
    }
    while(SimilarityThresholdRatio == NULL);

    do
    {
        NEmbTab=(int*)malloc((int)(Nmax/nObs)*sizeof(int));
    }
    while(NEmbTab == NULL);
    do
    {
        EmbPosTab=(int**)malloc((int)(Nmax/nObs)*sizeof(int*));
    }
    while(EmbPosTab == NULL);
    do
    {
        EmbDurationTab=(int**)malloc((int)(Nmax/nObs)*sizeof(int*));
    }
    while(EmbDurationTab == NULL);
    do
    {
        EmbEnergTab=(double**)malloc((double)(Nmax/nObs)*sizeof(double*));
    }
    while(EmbEnergTab == NULL);

    // affichage
    printf("Pass 2 : Emboli detection\nProgress\n");

    // init
    *NArts=0;// nombre total d'artefacts

    //boucle parallele
    #pragma omp parallel for  private(k1, k2) shared(DynamicSNR, SimilarityThresholdRatio, NArts, NEmbTab, EmbPosTab, EmbDurationTab, EmbEnergTab)
    for(k1=1; k1<=Nmax/nObs; k1++)//
    {
        double stdSignal;
        double *xWindow=NULL;

        // Selection de la fenetre d'analyse
        do
        {
            xWindow=(double*) malloc(2*nObs*sizeof(double));
        }
        while(xWindow == NULL);
        for(k2=0; k2<nObs; k2++)
        {
            xWindow[k2]=data[k2+(k1-1)*nObs];
            xWindow[k2+nObs]=data[k2+(k1-1)*nObs+(Nmax)*1];
            //if((k1==2)&&(k2<100)) printf("%d\t%f\t%f\n",k2, xWindow[k2],xWindow[k2+nObs]);
        }

        stdDevComplex(xWindow, nObs, &stdSignal);
        if (stdSignal>0.01)
        {
            // Variable locale
            double *ESubNorm=NULL, *ESubNormF=NULL, *VarThreshold=NULL;
            double *Mpv_pos=NULL;
            double *ThresholdHits=NULL;
            int *PositionHits=NULL, *DurationHits=NULL, *ZHits=NULL;
            int lHits=0, lArts=0;
            int *PositionArts=NULL, *DurationArts=NULL, *ZArts=NULL;
            double *ThresholdArts=NULL, *EnergyArts=NULL;
            int NEmb=0;
            int *EmbPos=NULL, *EmbDur;
            double *EmbEner=NULL;
            double *std_flux=NULL;

            do
            {
                ESubNorm=(double*)malloc(2*lEsub*sizeof(double));
            }
            while(ESubNorm == NULL);
            do
            {
                ESubNormF=(double*)malloc(2*lEsub*sizeof(double));
            }
            while(ESubNormF == NULL);
            do
            {
                VarThreshold=(double*)malloc(2*sizeof(double));
            }
            while(VarThreshold == NULL);
            do
            {
                Mpv_pos=(double*)malloc(2*sizeof(double));
            }
            while(Mpv_pos == NULL);
            do
            {
                ThresholdHits=(double*)malloc(3*lEsub*sizeof(double));
            }
            while(ThresholdHits == NULL);
            do
            {
                PositionHits=(int*)malloc(lEsub*sizeof(int));
            }
            while(PositionHits == NULL);
            do
            {
                DurationHits=(int*)malloc(lEsub*sizeof(int));
            }
            while(DurationHits == NULL);
            do
            {
                ZHits=(int*)malloc(lEsub*sizeof(int));
            }
            while(ZHits == NULL);
            do
            {
                ThresholdArts=(double*)malloc(lEsub*sizeof(double));
            }
            while(ThresholdArts == NULL);
            do
            {
                PositionArts=(int*)malloc(lEsub*sizeof(int));
            }
            while(PositionArts == NULL);
            do
            {
                DurationArts=(int*)malloc(lEsub*sizeof(int));
            }
            while(DurationArts == NULL);
            do
            {
                ZArts=(int*)malloc(lEsub*sizeof(int));
            }
            while(ZArts == NULL);
            do
            {
                EnergyArts=(double*)malloc(lEsub*sizeof(double));
            }
            while(EnergyArts == NULL);
            do
            {
                EmbPos=(int*)malloc(lEsub*sizeof(int));
            }
            while(EmbPos == NULL);
            do
            {
                EmbDur=(int*)malloc(lEsub*sizeof(int));
            }
            while(EmbDur == NULL);
            do
            {
                EmbEner=(double*)malloc(lEsub*sizeof(double));
            }
            while(EmbEner == NULL);
            do
            {
                std_flux=(double*)malloc(2*sizeof(double));
            }
            while(std_flux == NULL);

            // Analyse
            //printf("%d\t%d\n",k1,__LINE__);
            Energies_L1(xWindow, nObs, delay, lWind, Ncp, Ncn, ESubNorm, ESubNormF, VarThreshold, std_flux);
            Quality(ESubNormF, lEsub, &DynamicSNR[k1], Mpv_pos);
            HITS_detection(ESubNorm, ESubNormF, lEsub, VarThreshold, std_flux[0], Mpv_pos, ThresholdHits, &SimilarityThresholdRatio[k1], PositionHits, DurationHits, ZHits, &lHits);
            DetectionArts(ESubNorm, ESubNormF, lEsub, mpvApriori, ThresholdArts, PositionArts, DurationArts, ZArts, EnergyArts, &lArts);
            EmbDetector(ESubNorm, ESubNormF, lEsub, PositionHits, DurationHits, ZHits, &lHits, PositionArts, ZArts, &lArts, tolerance, &NEmb, EmbPos, EmbEner, EmbDur);
            //printf("%d\t%d\n",k1,__LINE__);
            *NArts+=lArts;

            // sauvegarde en tableau
            NEmbTab[k1-1]=NEmb;
            EmbPosTab[k1-1]=NULL;
            do
            {
                EmbPosTab[k1-1]=(int*)malloc(NEmb*sizeof(int));
            }
            while(EmbPosTab[k1-1] == NULL);
            EmbDurationTab[k1-1]=NULL;
            do
            {
                EmbDurationTab[k1-1]=(int*)malloc(NEmb*sizeof(int));
            }
            while(EmbDurationTab[k1-1] == NULL);
            EmbEnergTab[k1-1]=NULL;
            do
            {
                EmbEnergTab[k1-1]=(double*)malloc(NEmb*sizeof(double));
            }
            while(EmbEnergTab[k1-1] == NULL);

            for(k2=0; k2<NEmb; k2++)
            {
                EmbPosTab[k1-1][k2]=EmbPos[k2];
                EmbDurationTab[k1-1][k2]=EmbDur[k2];
                EmbEnergTab[k1-1][k2]=EmbEner[k2];
            }
            //printf("%d\t%d\n",k1,__LINE__);

            FREE(ESubNorm);
            FREE(ESubNormF);
            FREE(VarThreshold);
            FREE(Mpv_pos);
            FREE(ThresholdHits);
            FREE(PositionHits);
            FREE(DurationHits);
            FREE(ZHits);
            FREE(ThresholdArts);
            FREE(PositionArts);
            FREE(DurationArts);
            FREE(ZArts);
            FREE(EnergyArts);
            FREE(EmbPos);
            FREE(EmbEner);
            FREE(EmbDur);
            FREE(std_flux);
        }
        FREE(xWindow);
    }

    // Stockage finale
    *NEmbFinal=0;
    printf("Part\tEmbolus\tPos (s)\t\tEnergy\t\tDuration (ms)\n");
    for(k1=0; k1<Nmax/nObs; k1++)//
    {
        for(k2=0; k2<NEmbTab[k1]; k2++)
        {
            EmbPosFinal[k2+(int)(*NEmbFinal)]=(double)EmbPosTab[k1][k2]*(Tobs/((double)lEsub-1))+(double)k1*Tobs;
            EmbDuration[k2+(int)(*NEmbFinal)]=EmbDurationTab[k1][k2];
            EmbEnergFinal[k2+(int)(*NEmbFinal)]=EmbEnergTab[k1][k2];
            printf("%d\t%d\t%f\t%f\t%f\n",k1+1,k2+*NEmbFinal+1, EmbPosFinal[k2+(int)(*NEmbFinal)], EmbEnergFinal[k2+(int)(*NEmbFinal)], (double)EmbDuration[k2+(int)(*NEmbFinal)]/(double)fs*1000.);
        }
        *NEmbFinal+=NEmbTab[k1];
        FREE(EmbPosTab[k1]);
        FREE(EmbDurationTab[k1]);
        FREE(EmbEnergTab[k1]);

    }
    printf("Emboli detected : %d\n",*NEmbFinal);
    FREE(EmbPosTab);
    FREE(EmbDurationTab);
    FREE(EmbEnergTab);
    FREE(NEmbTab);

    FREE(DynamicSNR);
    FREE(SimilarityThresholdRatio);
}

// ---------------------------------------------------------------------------------------------------------------------------
// ---------------------------------------------------------------------------------------------------------------------------
// ---------------------------------------------------------------------------------------------------------------------------
// Calcul de l'energie a court terme pour chaque morceau et pre-calcul du seuil
void Energies_L1(const double* data, const int Nmax, const int delay, const int lWind, const int Ncp, const int Ncn, double* ESubNorm, double *ESubNormF, double *VarThreshold, double *std_flux)//,ESubNorm,ESub,ESubNormLim,ESubNormF,VarThreshold,Mpv0
{
    // declaration de variables
    int k;
    int lF=256;
    int lEsub=(int)round((double)(Nmax-lWind)/(double)delay)+1;
    int a2, b2;//a1, b1
    double Mpv2; //Mpv1
    int vSmooth=25;
    double sat1, sat2;
    int Fc;
    double std1, std2, stdBruit1, stdBruit2;
    int nDEsub1=0, nDEsub2=0, ic=0;

    double *edges=NULL;
    int *histo=NULL, *histoTemp=NULL;
    double *ESub1=NULL, *ESub2=NULL, *ESubNormN1=NULL, *ESubNormN2=NULL, *ESubNormF1=NULL, *ESubNormF2=NULL, *ESubNormLim1=NULL, *ESubNormLim2=NULL, *filtre=NULL;
    double *dEsub1=NULL, *dEsub2=NULL, *cumsum=NULL;
    int *c11=NULL;
    fftw_complex *in=NULL, *out=NULL;
    fftw_plan planFft=NULL;
    fftw_plan iplanFft=NULL;

    do
    {
        ESub1=(double*) malloc(sizeof(double) *lEsub);
    }
    while(ESub1 == NULL);
    do
    {
        ESub2=(double*) malloc(sizeof(double) *lEsub);
    }
    while(ESub2 == NULL);
    do
    {
        ESubNormN1=(double*) malloc(sizeof(double) *lEsub);
    }
    while(ESubNormN1 == NULL);
    do
    {
        ESubNormN2=(double*) malloc(sizeof(double) *lEsub);
    }
    while(ESubNormN2 == NULL);
    do
    {
        ESubNormF1=(double*) malloc(sizeof(double) *lEsub);
    }
    while(ESubNormF1 == NULL);
    do
    {
        ESubNormF2=(double*) malloc(sizeof(double) *lEsub);
    }
    while(ESubNormF2 == NULL);
    do
    {
        dEsub1=(double*) malloc(sizeof(double) *lEsub);
    }
    while(dEsub1 == NULL);
    do
    {
        dEsub2=(double*) malloc(sizeof(double) *lEsub);
    }
    while(dEsub2 == NULL);

    //---------------------------------------------//
    // Calcul du signal energetique
    //---------------------------------------------//
    calculESub(data, Nmax, lF, Ncp, Ncn, delay, lWind, ESub1, ESub2);

    // Lissage de l'energie
    //memcpy(ESubNormN1, ESub1, sizeof(double) *(int)round((double)Nmax/(double)delay));
    memcpy(ESubNormN2, ESub2, sizeof(double) *lEsub);
    for(k=0; k<lEsub; k++)
    {

        ESubNormN1[k]=0;
        ESubNormN1[k]=ESub1[k]-ESub2[k];
        if (ESubNormN1[k]<0)
            ESubNormN1[k]=0;
    }
    smoothMean(ESubNormN1, ESubNormF1, lEsub, vSmooth);
    smoothMean(ESubNormN2, ESubNormF2, lEsub, vSmooth);

    //---------------------------------------------//
    // Maximum de probabilite pour les frequences negatives
    //---------------------------------------------//
    do
    {
        histo = (int*) malloc(lEsub *sizeof(int));
    }
    while(histo == NULL);
    do
    {
        edges = (double*) malloc(lEsub *sizeof(double));
    }
    while(edges == NULL);

//    fhist(ESubNormF1, lEsub, lEsub, histo, edges);
//    histoTemp = (int*) malloc((lEsub-2) *sizeof(int));
//    memcpy(histoTemp, histo+1, sizeof(int) *(lEsub-2));
//    maximum(histoTemp, (lEsub-2), &a1, &b1);
//    FREE(histoTemp);
//    Mpv1=edges[b1+1];

    fhist(ESubNormF2, lEsub, lEsub, histo, edges);
    do
    {
        histoTemp = (int*) malloc((lEsub-10) *sizeof(int));
    }
    while(histoTemp == NULL);
    memcpy(histoTemp, histo+9, sizeof(int) *(lEsub-10));
    maximum(histoTemp, (lEsub-10), &a2, &b2);
    FREE(histoTemp);
    Mpv2=edges[b2+10];

    FREE(histo);
    FREE(edges);

    //---------------------------------------------//
    // Calcul des valeurs de saturation par les artefacts
    //---------------------------------------------//
    sat1 = fmean(ESubNormF1, lEsub) + 5.*stdDev(ESubNormF1, lEsub);
    sat2 = 10.*Mpv2;

    do
    {
        ESubNormLim1=(double*) malloc(sizeof(double) *lEsub);
    }
    while(ESubNormLim1 == NULL);
    do
    {
        ESubNormLim2=(double*) malloc(sizeof(double) *lEsub);
    }
    while(ESubNormLim2 == NULL);
    memcpy(ESubNormLim1, ESubNormN1, sizeof(double) *lEsub);
    memcpy(ESubNormLim2, ESubNormN2, sizeof(double) *lEsub);

    for (k=0; k<lEsub; k++)
    {
        if (ESubNormN1[k]>sat1)
        {
            ESubNormLim1[k]=sat1;
        }
        if (ESubNormN2[k]>sat2)
        {
            ESubNormLim2[k]=sat2;
        }
    }

    smoothMean(ESubNormLim1, ESubNormF1, lEsub, vSmooth);
    smoothMean(ESubNormLim2, ESubNormF2, lEsub, vSmooth);

    FREE(ESub1);
    FREE(ESub2);

    //---------------------------------------------//
    // Statistique pour le précalcul du seuil
    //---------------------------------------------//
    Fc=(int)round(32./(float)lEsub*4096.);
    do
    {
        filtre=(double*)malloc(sizeof(double) *(int)pow(2,round(log(lEsub)/log(2))));
    }
    while(filtre == NULL);
    for (k=0; k<(int)pow(2,round(log(lEsub)/log(2))); k++)
    {
        if (k<Fc)
        {
            filtre[k]=1.;
        }
        else if (k>=(int)pow(2,round(log(lEsub)/log(2)))-Fc)
        {
            filtre[k]=1.;
        }
        else
        {
            filtre[k]=0.;
        }
    }

    // Spectre 1
    do
    {
        in =(fftw_complex*) fftw_malloc(sizeof(fftw_complex) *(int)pow(2,round(log(lEsub)/log(2))));
    }
    while(in == NULL);
    do
    {
        out=(fftw_complex*) fftw_malloc(sizeof(fftw_complex) *(int)pow(2,round(log(lEsub)/log(2))));
    }
    while(out == NULL);
    for (k=0; k<(int)pow(2,round(log(lEsub)/log(2))); k++)
    {
        if (k<lEsub)
        {
            in[k][0]=ESubNormF1[k];
            in[k][1]=0.;
        }
        else
        {
            in[k][0]=0.;
            in[k][1]=0.;
        }
    }
    do
    {
        planFft=fftw_plan_dft_1d((int)pow(2,round(log(lEsub)/log(2))), in, out, FFTW_FORWARD, FFTW_ESTIMATE);
    }
    while(planFft == NULL);
    fftw_execute(planFft); //fft

    for (k=0; k<(int)pow(2,round(log(lEsub)/log(2))); k++)
    {
        out[k][0]*=filtre[k];
        out[k][1]*=filtre[k];
    }
    do
    {
        iplanFft=fftw_plan_dft_1d((int)pow(2,round(log(lEsub)/log(2))), out, in, FFTW_BACKWARD, FFTW_ESTIMATE);
    }
    while(iplanFft == NULL);
    fftw_execute(iplanFft); //fft
    for (k=0; k<lEsub; k++)
    {
        ESubNormF1[k]=in[k][0]/((double)pow(2,round(log(lEsub)/log(2))));
    }

    // Spectre 2
    for (k=0; k<(int)pow(2,round(log(lEsub)/log(2))); k++)
    {
        if (k<lEsub)
        {
            in[k][0]=ESubNormF2[k];
            in[k][1]=0;
        }
        else
        {
            in[k][0]=0;
            in[k][1]=0;
        }
    }
    fftw_execute(planFft); //fft
    for (k=0; k<(int)pow(2,round(log(lEsub)/log(2))); k++)
    {
        out[k][0]*=filtre[k];
        out[k][1]*=filtre[k];
    }
    fftw_execute(iplanFft); //fft
    for (k=0; k<lEsub; k++)
    {
        ESubNormF2[k]=in[k][0]/(double)pow(2,round(log(lEsub)/log(2)));
    }

    FREE(filtre);
    FFTW_FREE(in);
    FFTW_FREE(out);
    FFTW_DESTROY_PLAN(planFft);
    FFTW_DESTROY_PLAN(iplanFft);

    //---------------------------------------------//
    // Sortie de l'energie filtree
    //---------------------------------------------//
    for (k=0; k<lEsub; k++)
    {
        ESubNorm[k]=ESubNormN1[k];
        ESubNorm[k+lEsub]=ESubNormN2[k];

        ESubNormN1[k]=ESubNormLim1[k]-ESubNormF1[k];
        ESubNormN2[k]=ESubNormLim2[k]-ESubNormF2[k];


        ESubNormF[k]=ESubNormF1[k];
        ESubNormF[k+lEsub]=ESubNormF2[k];
    }
    std1=stdDev(ESubNormF1, lEsub);
    stdBruit1=stdDev(ESubNormN1, lEsub);
    VarThreshold[0]=2.*(2.*std1+stdBruit1);

    std2=stdDev(ESubNormF2, lEsub);
    stdBruit2=stdDev(ESubNormN2, lEsub);
    VarThreshold[1]=2.*(2.*std2+stdBruit2);

    //---------------------------------------------//
    // calcul Std Flux pour le seuil
    //---------------------------------------------//
    nDEsub1=0;
    nDEsub2=0;
    for (k=0; k<lEsub; k++)
    {
        if (ESubNorm[k]-ESubNormF[k]<0)
        {
            dEsub1[nDEsub1]=ESubNorm[k]-ESubNormF[k];
            nDEsub1++;
        }
        if (ESubNorm[k+lEsub]-ESubNormF[k+lEsub]<0)
        {
            dEsub2[nDEsub2]=ESubNorm[k+lEsub]-ESubNormF[k+lEsub];
            nDEsub2++;
        }
    }

    do
    {
        histoTemp = (int*) malloc(100 *sizeof(int));
    }
    while(histoTemp == NULL);
    do
    {
        edges = (double*) malloc(100 *sizeof(double));
    }
    while(edges == NULL);
    do
    {
        cumsum = (double*) malloc(100 *sizeof(double));
    }
    while(cumsum == NULL);
    do
    {
        c11 = (int*) malloc(100 *sizeof(int));
    }
    while(c11 == NULL);

    fhist(dEsub1, nDEsub1, 100, histoTemp, edges);
    cumsum[0]=0;
    ic=0;
    for (k=0; k<100; k++)
    {
        if (k==0)
            cumsum[k]=(double)histoTemp[k]/(double)nDEsub1/2.;
        else
            cumsum[k]=cumsum[k-1]+(double)histoTemp[k]/(double)nDEsub1/2.;

        if (cumsum[k]<0.1/100.)
        {
            c11[ic]=k;

            ic++;
        }
    }
    if (ic==0)
        std_flux[0]=fabs(edges[0]);
    else
    {
        ic--;
        std_flux[0]=fabs(edges[c11[ic]]);
    }

    fhist(dEsub2, nDEsub2, 100, histoTemp, edges);
    cumsum[0]=0;
    ic=0;
    for (k=0; k<100; k++)
    {
        if (k==0)
            cumsum[k]=(double)histoTemp[k]/(double)nDEsub2/2.;
        else
            cumsum[k]=cumsum[k-1]+(double)histoTemp[k]/(double)nDEsub2/2.;

        if (cumsum[k]<0.1/100.)
        {
            c11[ic]=k;
            ic++;
        }
    }
    if (ic==0)
        std_flux[1]=fabs(edges[0]);
    else
    {
        ic--;
        std_flux[1]=fabs(edges[c11[ic]]);
    }

    FREE(histoTemp);
    FREE(edges);
    FREE(cumsum);
    FREE(c11);

    FREE(ESubNormN1);
    FREE(ESubNormN2);
    FREE(ESubNormF1);
    FREE(ESubNormF2);
    FREE(ESubNormLim1);
    FREE(ESubNormLim2);
    FREE(dEsub1);
    FREE(dEsub2);
}

// ---------------------------------------------------------------------------------------------------------------------------
// ---------------------------------------------------------------------------------------------------------------------------
// ---------------------------------------------------------------------------------------------------------------------------
// Calcul de la qualite du signal dans le troncon
// ---------------------------------------------------------------------------------------------------------------------------
// ---------------------------------------------------------------------------------------------------------------------------
// ---------------------------------------------------------------------------------------------------------------------------
void Quality(const double *ESubNormF, int lEsub, double *DynamicSNR, double *Mpv_pos)
{
    int k;
    double EEm;
    double *EEc=NULL, *EEcA=NULL, *EEcB=NULL, *xHisto=NULL;
    int Lmin, Lmax, La, Lb=0;
    int *indmin=NULL, *indmax=NULL, *a=NULL, *b=NULL, *histo=NULL;
    int maxi;
    int b1, b2;

    do
    {
        EEc=(double*)malloc(lEsub*sizeof(double));
    }
    while(EEc == NULL);
    EEm=fmean(ESubNormF, lEsub);
    for (k=0; k<lEsub; k++)
    {
        EEc[k]=ESubNormF[k]-EEm;
    }

    do
    {
        indmin=(int*)malloc(lEsub*sizeof(int));
    }
    while(indmin == NULL);
    do
    {
        indmax=(int*)malloc(lEsub*sizeof(int));
    }
    while(indmax == NULL);
    extr(EEc, lEsub, indmin, &Lmin, indmax, &Lmax);

    do
    {
        a=(int*)malloc(Lmax*sizeof(int));
    }
    while(a == NULL);
    La=0;
    for(k=0; k<Lmax; k++)
    {
        if (EEc[indmax[k]]>0)
        {
            a[La]=k;
            La++;
        }
    }
    do
    {
        EEcA=(double*)malloc(La*sizeof(double));
    }
    while(EEcA == NULL);
    do
    {
        histo=(int*)malloc(La*sizeof(int));
    }
    while(histo == NULL);
    do
    {
        xHisto=(double*)malloc(La*sizeof(double));
    }
    while(xHisto == NULL);

    for(k=0; k<La; k++)
    {
        EEcA[k]=EEc[indmax[a[k]]];
    }
    fhist(EEcA, La, La, histo, xHisto);
    maximum(histo, La, &maxi, &b1);
    Mpv_pos[0]=xHisto[b1];
    FREE(a);
    FREE(EEcA);
    FREE(histo);
    FREE(xHisto);

    do
    {
        b=(int*)malloc(Lmin*sizeof(int));
    }
    while(b == NULL);
    for(k=0; k<Lmin; k++)
    {
        if (EEc[indmin[k]]<0)
        {
            b[Lb]=k;
            Lb++;
        }
    }
    do
    {
        EEcB=(double*)malloc(Lb*sizeof(double));
    }
    while(EEcB == NULL);
    do
    {
        histo=(int*)malloc(Lb*sizeof(int));
    }
    while(histo == NULL);
    do
    {
        xHisto=(double*)malloc(Lb*sizeof(double));
    }
    while(xHisto == NULL);
    for(k=0; k<Lb; k++)
    {
        EEcB[k]=EEc[indmin[b[k]]];
    }
    fhist(EEcB, Lb, Lb, histo, xHisto);
    maximum(histo, Lb, &maxi, &b2);
    Mpv_pos[1]=xHisto[b2];

    *DynamicSNR=10.*log10( abs( (Mpv_pos[0]-Mpv_pos[1])/Mpv_pos[1] ));

    FREE(b);
    FREE(EEcB);
    FREE(histo);
    FREE(xHisto);

    FREE(indmin);
    FREE(indmax);
}

// ---------------------------------------------------------------------------------------------------------------------------
// ---------------------------------------------------------------------------------------------------------------------------
// ---------------------------------------------------------------------------------------------------------------------------
// Detection des Hits (position, energie, duree)
// ---------------------------------------------------------------------------------------------------------------------------
// ---------------------------------------------------------------------------------------------------------------------------
// ---------------------------------------------------------------------------------------------------------------------------
void HITS_detection(const double *ESubNorm, const double *ESubNormF, const int lEsub, const double *VarThreshold, const double alphaH, const double *Mpv_pos, double *ThresholdHits, double *SimilarityThresholdRatio, int *PositionHitsFinal, int *DurationHitsFinal, int *ZHits, int *LHits)
{
    int k, l, nPp, nPn;
    int *Z=NULL, *dZ=NULL, *Pp=NULL, *Pn=NULL, *a=NULL;
    double moy=fmean(ESubNormF, lEsub);
    double *NewThresholdHits=NULL, *ESubSmooth=NULL;
    int *DurationHits=NULL;

    //---------------------------------------------//
    // Calcul du seuil
    //---------------------------------------------//
    do
    {
        NewThresholdHits=(double*)malloc(lEsub*sizeof(double));
    }
    while(NewThresholdHits == NULL);
    for (k=0; k<lEsub; k++)
    {
        ThresholdHits[k]=ESubNormF[k]+alphaH;
        ThresholdHits[k+lEsub]=moy+1.5*alphaH;
        ThresholdHits[k+2*lEsub]=ThresholdHits[k];
        NewThresholdHits[k]=ThresholdHits[k+2*lEsub];
    }

    *SimilarityThresholdRatio=100.-100.*fabs(fmean(ThresholdHits, lEsub)-ThresholdHits[lEsub])/ThresholdHits[lEsub];

    if (*SimilarityThresholdRatio<0.)
    {
        for (k=0; k<lEsub; k++)
        {
            ThresholdHits[k+2*lEsub]=ThresholdHits[lEsub];
            NewThresholdHits[k]=ThresholdHits[k+2*lEsub];
        }
    }

    //---------------------------------------------//
    //Detection de Hits
    //---------------------------------------------//

    //Smooth
    do
    {
        ESubSmooth=(double*)malloc(lEsub*sizeof(double));
    }
    while(ESubSmooth == NULL);
    smoothMean(ESubNorm, ESubSmooth, lEsub, 10);

    do
    {
        Z=(int*)malloc(lEsub*sizeof(int));
    }
    while(Z == NULL);
    do
    {
        dZ=(int*)malloc(lEsub*sizeof(int));
    }
    while(dZ == NULL);
    do
    {
        Pp=(int*)malloc(lEsub*sizeof(int));
    }
    while(Pp == NULL);
    do
    {
        Pn=(int*)malloc(lEsub*sizeof(int));
    }
    while(Pn == NULL);

    nPp=0;
    nPn=0;
    for (k=0; k<lEsub; k++)
    {
        if (ESubSmooth[k]>NewThresholdHits[k])
            Z[k]=1;
        else
            Z[k]=0;

        if (k==0)
            dZ[k]=0;
        else
            dZ[k]=Z[k]-Z[k-1];

        if (dZ[k]==1)
        {
            Pp[nPp]=k;
            nPp++;
        }
        else if (dZ[k]==-1)
        {
            Pn[nPn]=k;
            nPn++;
        }
    }

    // calcul de la duree des Hits
    l=0;
    if (nPp == nPn)
    {
        do
        {
            DurationHits=(int*)malloc(nPp*sizeof(int));
        }
        while(DurationHits == NULL);
        do
        {
            a=(int*)malloc(nPp*sizeof(int));
        }
        while(a == NULL);
        for (k=0; k<nPp; k++)
        {
            DurationHits[k]=abs(Pn[k]-Pp[k]);
            if (DurationHits[k]<15)
            {
                a[l]=k;
                l++;
            }
        }
    }
    else if (nPp > nPn)
    {
        do
        {
            DurationHits=(int*)malloc(nPp*sizeof(int));
        }
        while(DurationHits == NULL);
        do
        {
            a=(int*)malloc(nPp*sizeof(int));
        }
        while(a == NULL);

        for (k=0; k<nPp; k++)
        {
            if (k<nPn)
                DurationHits[k]=abs(Pn[k]-Pp[k]);
            else
                DurationHits[k]=abs(lEsub-Pp[k]);

            if (DurationHits[k]<15)
            {
                a[l]=k;
                l++;
            }
        }
    }
    else
    {
        do
        {
            DurationHits=(int*)malloc(nPn*sizeof(int));
        }
        while(DurationHits == NULL);
        do
        {
            a=(int*)malloc(nPn*sizeof(int));
        }
        while(a == NULL);

        for (k=0; k<nPn; k++)
        {
            if (k<nPp)
                DurationHits[k]=abs(Pn[k]-Pp[k]);
            else
                DurationHits[k]=abs(lEsub-Pp[k]);

            if (DurationHits[k]<15)
            {
                a[l]=k;
                l++;
            }
        }
    }

    *LHits=l;
    for (k=0; k<*LHits; k++)
    {
        PositionHitsFinal[k]=Pp[a[k]];
        DurationHitsFinal[k]=DurationHits[a[k]];
    }

    l=0;
    for (k=0; k<lEsub; k++)
    {

        if ((k >= PositionHitsFinal[l]) && (k < PositionHitsFinal[l]+DurationHitsFinal[l]))
        {
            //printf("%d\n", PositionHitsFinal[l]);
            ZHits[k]=1;
            if (k+1 == PositionHitsFinal[l]+DurationHitsFinal[l])
            {
                l++;
            }
        }
        else
            ZHits[k]=0;
    }

    FREE(Z);
    FREE(dZ);
    FREE(Pp);
    FREE(Pn);
    FREE(DurationHits);
    FREE(a);
    FREE(ESubSmooth);
    FREE(NewThresholdHits);
}

// ---------------------------------------------------------------------------------------------------------------------------
// ---------------------------------------------------------------------------------------------------------------------------
// ---------------------------------------------------------------------------------------------------------------------------
// Detection des artefacts
// ---------------------------------------------------------------------------------------------------------------------------
// ---------------------------------------------------------------------------------------------------------------------------
// ---------------------------------------------------------------------------------------------------------------------------
void DetectionArts(const double *ESubNorm, const double *ESubNormF, const int lEsub, const double *Mpv_apriori, double *ThresholdArts, int *PositionArts, int *DurationArts, int *ZArts, double *EnergyArts, int *LArts)
{
    int k, nPp=0, nPn=0;
    int *dZ=NULL, *Pn=NULL;
    double *ESubSmooth=NULL;
    double alphaA=0.8*2;

    do
    {
        ESubSmooth=(double*)malloc(lEsub*sizeof(double));
    }
    while(ESubSmooth == NULL);
    do
    {
        dZ=(int*)malloc(lEsub*sizeof(int));
    }
    while(dZ == NULL);
    do
    {
        Pn=(int*)malloc(lEsub*sizeof(int));
    }
    while(Pn == NULL);

    smoothMean(ESubNorm+lEsub, ESubSmooth, lEsub, 2);

    for(k=0; k<lEsub; k++)
    {
        ThresholdArts[k]=ESubNormF[k+lEsub] + alphaA*Mpv_apriori[1];

        if (ESubSmooth[k]>ThresholdArts[k])
        {
            ZArts[k]=1;
        }
        else
            ZArts[k]=0;

        if (k==0)
            dZ[k]=0;
        else
            dZ[k]=ZArts[k]-ZArts[k-1];

        if (dZ[k]==1)
        {
            PositionArts[nPp]=k;
            EnergyArts[nPp]=ESubSmooth[PositionArts[nPp]] - ESubNormF[PositionArts[nPp]+lEsub];
            nPp++;
        }
        else if (dZ[k]==-1)
        {
            Pn[nPn]=k;
            nPn++;
        }
    }

    if (nPp == nPn)
    {
        *LArts=nPp;
        for (k=0; k<nPp; k++)
            DurationArts[k]=abs(Pn[k]-PositionArts[k]);
    }
    else if (nPp > nPn)
    {
        *LArts=nPp;
        for (k=0; k<nPp; k++)
        {
            if (k<nPn)
                DurationArts[k]=abs(Pn[k]-PositionArts[k]);
            else
                DurationArts[k]=abs(lEsub-PositionArts[k]);
        }
    }
    else
    {
        *LArts=nPn;
        for (k=0; k<nPn; k++)
        {
            if (k<nPp)
                DurationArts[k]=abs(Pn[k]-PositionArts[k]);
            else
                DurationArts[k]=abs(lEsub-PositionArts[k]);
        }
    }

    FREE(dZ);
    FREE(Pn);
    FREE(ESubSmooth);
}

// ---------------------------------------------------------------------------------------------------------------------------
// ---------------------------------------------------------------------------------------------------------------------------
// ---------------------------------------------------------------------------------------------------------------------------
// Elimination des Hits qui coincident avec un artefact
// ---------------------------------------------------------------------------------------------------------------------------
// ---------------------------------------------------------------------------------------------------------------------------
// ---------------------------------------------------------------------------------------------------------------------------
void EmbDetector(const double *ESubNorm, const double *ESubNormF, const int lEsub, const int *PositionHits, const int *DurationHits, const int *ZHits, const int *lHits, const int *PositionArts, const int *ZArts, const int *lArts, const int Tolerance, int *NEmb, int* EmbPos, double *EmbEner, int *EmbDur)
{
    int k, k2, l=0, l2=0, l3=0, a;
    int *com=NULL, *indCom=NULL, *AA=NULL, *BB=NULL;

    do
    {
        com=(int*)malloc(lEsub*sizeof(int));
    }
    while(com == NULL);
    do
    {
        indCom=(int*)malloc(lEsub*sizeof(int));
    }
    while(indCom == NULL);

    for (k=0; k<lEsub; k++)
    {
        com[k] = ZHits[k] && ZArts[k];
        if (com[k]==1)
        {
            indCom[l]=k;
            l++;
        }
    }

    do
    {
        AA=(int*)malloc((int)(*lArts)*sizeof(int));
    }
    while(AA == NULL);
    do
    {
        BB=(int*)malloc((int)(*lArts)*sizeof(int));
    }
    while(BB == NULL);

    *NEmb=0;
    for (k=0; k<*lHits; k++)
    {
        l2=0;
        for (k2=0; k2<*lArts; k2++)
        {
            if (abs(PositionArts[k2]-PositionHits[k])<2*Tolerance)
            {
                AA[l2]=k2;
                l2++;
            }
        }
        l3=0;
        for (k2=0; k2<l; k2++)
        {
            if (abs(indCom[k2]-PositionHits[k])<2*Tolerance)
            {
                BB[l3]=k2;
                l3++;
            }
        }

        if ((l2==0) && (l3==0))
        {
            EmbPos[*NEmb]=PositionHits[k];
            EmbEner[*NEmb]=ESubNorm[PositionHits[k]]/ESubNormF[PositionHits[k]];
            (*NEmb)++;
        }
    }

    FREE(AA);
    FREE(BB);
    FREE(com);

    l2=0;
    a=0;
    for (k=0; k<(*NEmb); k++)
    {
        l2=0;
        for (k2=0; k2<*lHits; k2++)
        {
            if (PositionHits[k2]-EmbPos[k]==0)
            {
                a=k2;
            }
        }
        EmbDur[k]=DurationHits[a];
    }
}

// ---------------------------------------------------------------------------------------------------------------------------
// ---------------------------------------------------------------------------------------------------------------------------
// ---------------------------------------------------------------------------------------------------------------------------
// Detection complete : fonction a appeler depuis le main
// ---------------------------------------------------------------------------------------------------------------------------
// ---------------------------------------------------------------------------------------------------------------------------
// ---------------------------------------------------------------------------------------------------------------------------
void detectionEmbolus(const double* dataWave, const int Nmax, const int fs, const double Tobs, const double fc, const int tolerance, int *NEmbFinal, double *EmbPosFinal, double *EmbEnergFinal, int *EmbDuration, int *NArts)
{
    double *mpvApriori=NULL;
    do
    {
        mpvApriori = (double*) malloc(2*sizeof(double));
    }
    while(mpvApriori == NULL);

    pass1Detection(dataWave, Nmax, fs, fc, mpvApriori);

    pass2Detection(dataWave, Nmax, fs, Tobs, fc, tolerance, mpvApriori, NEmbFinal, EmbPosFinal, EmbEnergFinal, EmbDuration, NArts);
    FREE(mpvApriori);
}
