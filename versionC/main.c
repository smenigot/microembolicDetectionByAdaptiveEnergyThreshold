#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <string.h>

#include "openWaveFile.h"
#include "emboleDetection.h"
#include "mathTools.h"

#define FREE(_p) do{ \
	 free( (_p) ); \
     (_p) = NULL; \
}while( (_p) != NULL)

int main(int argc, char *argv[])
{
    // declarations de variables
    char nameWave[100];
    char nameFile[104];
    char nameSave[104];
    FILE* fichier = NULL;
    FILE* saveFileResults = NULL;
    Wavefile header;
    double *dataWave = NULL;
    int k=0, FIN=0;

    // reglage de la detection
    double Tobs=10;
    double fc=150; // [Hz] filtrage passe-haut pour eliminer les artefacts
    int tolerance=5;

    int nbEmboli=0, nbArtefacts=0;
    int *emboliDuration=NULL;
    double *emboliPosition=NULL, *emboliEnergy=NULL;


    do
    {
        // Ouverture fichier
        if (argc <= 1)
        {
            printf("Enter name file (without extension) : ");
            scanf("%s",nameWave);
            FIN=1;
        }
        else
        {
            k++;
            //strcpy(nameWave,"T17");
            strcpy(nameWave,argv[k]);

            if (k==argc-1)
                FIN=1;
        }


        strcpy(nameFile,nameWave);
        strcat(nameFile,".wav");
        fichier = fopen(nameFile, "rb");
        if ( fichier == NULL )
        {
            fprintf(stderr,"Can't open input file %s", nameFile);
            return EXIT_FAILURE;
        }
        else
        {
            // lecture du fichier
            printf("Open the file \"%s\"\n",nameFile);
            readHeaderWave(fichier, &header);
            if (header.channels != 2)
            {
                fprintf(stderr,"Wave file with %d channel", header.channels);
                return EXIT_FAILURE;
            }
            do
            {
                dataWave = (double*)malloc(header.channels*header.bytes_in_data*sizeof(double));
            }
            while(dataWave == NULL);

            readDataWave(fichier, &header, dataWave);

            // detection
            do
            {
                emboliPosition = (double*)malloc(sizeof(double)*header.bytes_in_data);
            }
            while(emboliPosition == NULL);
            do
            {
                emboliEnergy = (double*)malloc(sizeof(double)*header.bytes_in_data);
            }
            while(emboliEnergy == NULL);
            do
            {
                emboliDuration = (int*)malloc(sizeof(int)*header.bytes_in_data);
            }
            while(emboliDuration == NULL);

            detectionEmbolus(dataWave, header.bytes_in_data, header.frequency, Tobs, fc, tolerance, &nbEmboli, emboliPosition, emboliEnergy, emboliDuration, &nbArtefacts);

            //sauvegarde fichier texte
            strcpy(nameSave,nameWave);
            strcat(nameSave,".csv");
            saveFileResults = fopen(nameSave, "w+t");
            if ( fichier == NULL )
            {
                fprintf(stderr,"Can't save results %s", nameSave);
                return EXIT_FAILURE;
            }
            else
            {
                printf("Save results in the file \"%s\"\n",nameSave);
                saveEmbole(saveFileResults, &header, &nbEmboli, emboliPosition, emboliEnergy, emboliDuration, &nbArtefacts);
                fclose(saveFileResults);
            }

            FREE(dataWave);
            FREE(emboliPosition);
            FREE(emboliEnergy);
            FREE(emboliDuration);
            fclose(fichier);
        }
    }while (FIN!=1);


    return EXIT_SUCCESS;
}
