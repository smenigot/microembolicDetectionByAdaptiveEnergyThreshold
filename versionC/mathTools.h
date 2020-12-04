#ifndef MATHTOOLS_H_INCLUDED
#define MATHTOOLS_H_INCLUDED

double hamming(int n, int N);
void fhist(const double* data, int lData, int nClasse, int* histo, double* edges);
void hist(const int* data, int lData, int nClasse, int* histo, double* edges);
void minimum(const int* data, int lData, int* mini, int* i_mini);
void maximum(const int* data, int lData, int* maxi, int* i_maxi);
void fminimum(const double* data, int lData, double* mini, int* i_mini);
void fmaximum(const double* data, int lData, double* maxi, int* i_maxi);
void echanger(double* tableau, int a, int b);
void quickSort(double* tableau, int debut, int fin);
void mediane(const double *data, int lData, double *med);
double mean(const int *data, int lData);
double fmean(const double *data, int lData);
void meanComplex(const double *data, int lData, double* moy);
double stdDev(const double *data, int lData);
void stdDevComplex(const double *data, int lData, double *result);
void smoothMean(const double* data, double *dataSmooth, int nMax, int lWind);
void extr(const double *data, int lData, int *indmin, int *Lmin, int *indmax, int *Lmax);
#endif // MATHTOOLS_H_INCLUDED
