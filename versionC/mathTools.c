#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <string.h>
#include "mathTools.h"

#define FREE(_p) do{ \
	 free( (_p) ); \
     (_p) = NULL; \
}while( (_p) != NULL)

// fenetre de Hamming
double hamming(int n, int N)
{
    double temp;
    temp = 0.54-0.46*cos(2*M_PI*n/N);

    return temp;
}

// Histogramme
void hist(const int* data, int lData, int nClasse, int* histo, double* xo)
{
    int k1,k2;
    int mini=0, maxi=0;
    double *edges=NULL;

    minimum(data, lData, &mini, &k1);
    maximum(data, lData, &maxi, &k2);

    //printf("%f-%f\n",maxi,mini);

    do
    {
        edges=(double*)malloc((nClasse+1)*sizeof(double));
    }
    while(edges == NULL);
    edges[0]=mini;
    for (k1=0; k1<nClasse; k1++)
    {
        edges[k1+1] = edges[k1] + (maxi-mini)/(double)nClasse;
        if (k1 == nClasse-1)
        {
            edges[k1+1]=maxi;
        }
        xo[k1]=edges[k1] + (maxi-mini)/(double)nClasse/2.;
        histo[k1]=0;
        for(k2=0; k2<lData; k2++)
        {
            if ((k1>0) && (data[k2]>edges[k1]) && (data[k2]<=edges[k1+1]))
            {
                histo[k1]++;
            }
            if ((k1==0) && (data[k2]>=edges[k1]) && (data[k2]<=edges[k1+1]))
            {
                histo[k1]++;
            }
        }
    }
    FREE(edges);
}
void fhist(const double* data, int lData, int nClasse, int* histo, double* xo)
{
    int k1,k2;
    double mini=0, maxi=0;
    double *edges=NULL;

    fminimum(data, lData, &mini, &k1);
    fmaximum(data, lData, &maxi, &k2);

    do
    {
        edges=(double*)malloc((nClasse+1)*sizeof(double));
    }
    while(edges == NULL);

    edges[0]=mini;
    for (k1=0; k1<nClasse; k1++)
    {
        edges[k1+1] = edges[k1] + (maxi-mini)/(double)nClasse;
        if (k1 == nClasse-1)
        {
            edges[k1+1]=maxi;
        }
        xo[k1]=edges[k1] + (maxi-mini)/(double)nClasse/2.;
        histo[k1]=0;
        for(k2=0; k2<lData; k2++)
        {
            if ((k1>0) && (data[k2]>edges[k1]) && (data[k2]<=edges[k1+1]))
            {
                histo[k1]++;
            }
            if ((k1==0) && (data[k2]>=edges[k1]) && (data[k2]<=edges[k1+1]))
            {
                histo[k1]++;
            }
        }
//        printf("[%f\t%f]:\t%d\n", edges[k1], edges[k1+1], histo[k1]);
    }
    FREE(edges);
}

// minimum
void minimum(const int* data, int lData, int* mini, int* i_mini)
{
    int k;
    *mini = data[0];
    *i_mini = 0;
    for (k=0; k<lData; k++)
    {
        if (data[k]<*mini)
        {
            *mini = data[k];
            *i_mini = k;
        }
    }
}
void fminimum(const double* data, int lData, double* mini, int* i_mini)
{
    int k;
    *mini = data[0];
    *i_mini = 0;
    for (k=0; k<lData; k++)
    {
        if (data[k]<*mini)
        {
            *mini = data[k];
            *i_mini = k;
        }
    }
}

// maximum
void maximum(const int* data, int lData, int* maxi, int* i_maxi)
{
    int k;
    *maxi = data[0];
    *i_maxi = 0;
    for (k=0; k<lData; k++)
    {
        if (data[k]>*maxi)
        {
            *maxi = data[k];
            *i_maxi = k;
        }
    }
}
void fmaximum(const double* data, int lData, double* maxi, int* i_maxi)
{
    int k;
    *maxi = data[0];
    *i_maxi = 0;
    for (k=0; k<lData; k++)
    {
        if (data[k]>*maxi)
        {
            *maxi = data[k];
            *i_maxi = k;
        }
    }
}

// Tri rapide
void echanger(double* tableau, int a, int b)
{
    double temp = tableau[a];
    tableau[a] = tableau[b];
    tableau[b] = temp;
}

void quickSort(double* tableau, int debut, int fin)
{
    int gauche = debut-1;
    int droite = fin+1;
    const double pivot = tableau[debut];
    /* Si le tableau est de longueur nulle, il n'y a rien à faire. */
    if(debut >= fin)
        return;
    /* Sinon, on parcourt le tableau, une fois de droite à gauche, et une
       autre de gauche à droite, à la recherche d'éléments mal placés,
       que l'on permute. Si les deux parcours se croisent, on arrête. */
    while(1)
    {
        do droite--;
        while(tableau[droite] > pivot);
        do gauche++;
        while(tableau[gauche] < pivot);
        if(gauche < droite)
            echanger(tableau, gauche, droite);
        else break;
    }
    /* Maintenant, tous les éléments inférieurs au pivot sont avant ceux
       supérieurs au pivot. On a donc deux groupes de cases à trier. On utilise
       pour cela... la méthode quickSort elle-même ! */
    quickSort(tableau, debut, droite);
    quickSort(tableau, droite+1, fin);
}

// Mediane
void mediane(const double *data, int lData, double *med)
{
    double *dataTemp;
    do
    {
        dataTemp = (double*) malloc(lData*sizeof(double));
    }
    while(dataTemp == NULL);
    memcpy(dataTemp, data, lData*sizeof(double));

    quickSort(dataTemp, 0, lData);

    *med = dataTemp[(int)round((double)lData/2.)];

    FREE(dataTemp);
}

// Moyenne
double mean(const int *data, int lData)
{
    int k;
    double sum=0;
    for (k=0; k<lData; k++)
    {
        sum+=data[k];
    }

    return sum/lData;
}
double fmean(const double *data, int lData)
{
    int k;
    double sum=0;
    for (k=0; k<lData; k++)
    {
        sum+=data[k];
    }

    return sum/lData;
}

// Moyenne Complexe
void meanComplex(const double *data, int lData, double* moy)
{
    int k;
    moy[0]=0;
    moy[1]=0;
    for (k=0; k<lData; k++)
    {
        moy[0]+=data[k];
        moy[1]+=data[k+lData];
    }
    moy[0]/=lData;
    moy[1]/=lData;
}

// Ecart-Type
double stdDev(const double *data, int lData)
{
    int k;
    double sumCarre=0, moy=0;
    moy=fmean(data, lData);
    for (k=0; k<lData; k++)
    {
        sumCarre+=pow(data[k]-moy,2);
    }
    sumCarre/=(lData-1);

    return sqrt(sumCarre);
}

// Ecart-Type Complex
void stdDevComplex(const double *data, int lData, double *result)
{
    int k;
    double *moy=NULL;
    double sum=0;

    //moyenne complexe
    do
    {
        moy=(double*)malloc(2*sizeof(double));
    }
    while(moy == NULL);
    meanComplex(data, lData, moy);

    for (k=0; k<lData; k++)
    {
        sum+=pow(data[k]-moy[0],2)+pow(data[k+lData]-moy[1],2);
    }
    sum/=(lData-1);

    FREE(moy);

    *result=sqrt(sum);
}

// Smooth function
void smoothMean(const double* data, double *dataSmooth, int nMax, int lWind)
{
    int k1, k2;
    double *dataTemp=NULL;
    do
    {
        dataTemp=(double*)malloc(lWind*sizeof(double));
    }
    while(dataTemp == NULL);

    for(k1=0; k1<nMax; k1++)
    {
        for(k2=0; k2<lWind; k2++)
        {
            if ((k1-(int)round((double)lWind/2.)+k2+1 < 0) || (k1-(int)round((double)lWind/2.-1)+k2+1>nMax))
                dataTemp[k2]=0;
            else
                dataTemp[k2]=data[k1-(int)round((double)lWind/2.-1)+k2];

//            printf("%d\t%d\t%d\t%f\n",k2,k1-(int)round((double)lWind/2.)+k2+1, k1-(int)round((double)lWind/2.-1)+k2,dataTemp[k2]);
        }
//        printf("\n");
        dataSmooth[k1]=fmean(dataTemp, lWind);
    }
    FREE(dataTemp);
}

// extr
void extr(const double *data, int lData, int *indmin, int *Lmin, int *indmax, int *Lmax)
{
    int k;

    *Lmin=0;
    *Lmax=0;
    for (k=1; k<lData-1; k++)
    {
        if ((data[k]<data[k-1]) && (data[k]<data[k+1]))
        {
            indmin[*Lmin]=k;
            (*Lmin)++;
        }
        if ((data[k]>data[k-1]) && (data[k]>data[k+1]))
        {
            indmax[*Lmax]=k;
            (*Lmax)++;
        }
    }
}

