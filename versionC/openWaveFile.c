#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include "openWaveFile.h"

// lecture de l'en-tete wav
int readHeaderWave(FILE* fichier, Wavefile* header)
{
    if ( fread(header,sizeof(Wavefile),1,fichier) < 1 )
    {
        fprintf(stderr,"Can't read file header\n");
        return 0;
    }
    if ( (header->id[0] != 'R') || (header->id[1] != 'I') || (header->id[2] != 'F') || (header->id[3] != 'F') )
    {
        fprintf(stderr,"ERROR: Not wav format\n");
        return 0;
    }
    return 1;
}

// lecture des donnees
void readDataWave(FILE* fichier, Wavefile* header, double *data)
{
    short value=0;
    int Nmax=header->bytes_in_data/sizeof(short)/header->channels;
    int k=0,c=0;
    header->bytes_in_data=Nmax; //nombre de donnees
    k=0;
    c=0;
    while( fread(&value,sizeof(value),1,fichier) )
    {
        if (k%2 == 0)
        {
            data[c]=(double)value/(double)(pow(2,header->bits_per_sample-1)-1);//c + bytes_in_data*r
        }
        else
        {
            if (header->channels == 2)
                data[c+Nmax]=-(double)value/(double)(pow(2,header->bits_per_sample-1)-1);//c + bytes_in_data*r
            else
                data[c+Nmax]=0;
            c++;
        }
        k++;
    }
}

// ---------------------------------------------------------------------------------------------------------------------------
// ---------------------------------------------------------------------------------------------------------------------------
// ---------------------------------------------------------------------------------------------------------------------------
// Ecriture dans un fichier de sauvegarde
// ---------------------------------------------------------------------------------------------------------------------------
// ---------------------------------------------------------------------------------------------------------------------------
// ---------------------------------------------------------------------------------------------------------------------------
void saveEmbole(FILE* fichier, const Wavefile* header, const int *NEmbFinal, const double *EmbPosFinal, const double *EmbEnergFinal, const int *EmbDuration, const int *NArts)
{
    int k;

//    fprintf(fichier, "Wave File Information\n");
//    fprintf(fichier, "\tChannel Number: %d\n", header->channels);
//    fprintf(fichier, "\tSampling frequency: %2.2f kHz\n", header->frequency/1000.);
//    fprintf(fichier, "\tDuration: %d'%d''\n", (int)((double)header->bytes_in_data/(double)header->frequency)/60, (int)((double)header->bytes_in_data/(double)header->frequency)%60);
//    fprintf(fichier, "\tBits per sample: %d\n", header->bits_per_sample);
//    fprintf(fichier, "\tSample number: %d\n", header->bytes_in_data);
//    fprintf(fichier, "\n");
//
//    fprintf(fichier, "Emboli detection\n");
//    fprintf(fichier, "\tEmboli number: %d\n", *NEmbFinal);
//    fprintf(fichier, "\tArtefact number: %d\n", *NArts);
//    fprintf(fichier, "\n");
//
//    fprintf(fichier, "Details\n");
    fprintf(fichier, "Number\tPosition (s)\tEnergy\tDuration (ms)\n");
    for (k=0; k<*NEmbFinal; k++)
    {
        fprintf(fichier, "%d\t%f\t%f\t%f\n",k, EmbPosFinal[k], EmbEnergFinal[k], (double)EmbDuration[k]/(double)header->frequency*1000);
    }
}

