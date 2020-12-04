#ifndef OPENWAVEFILE_H_INCLUDED
#define OPENWAVEFILE_H_INCLUDED

typedef struct Wavefile Wavefile;
struct Wavefile
{
    char    id[4];          // should always contain "RIFF"
    int     totallength;    // total file length minus 8
    char    wavefmt[8];     // should be "WAVEfmt "
    int     format;         // 16 for PCM format
    short   pcm;            // 1 for PCM format
    short   channels;       // channels
    int     frequency;      // sampling frequency
    int     bytes_per_second;
    short   bytes_by_capture;
    short   bits_per_sample;
    char    data[4];        // should always contain "data"
    int     bytes_in_data;
};


int readHeaderWave(FILE* fichier, Wavefile* inHeader);
void readDataWave(FILE* fichier, Wavefile* inHeader, double *data);
void saveEmbole(FILE* fichier, const Wavefile* header, const int *NEmbFinal, const double *EmbPosFinal, const double *EmbEnergFinal, const int *EmbDuration, const int *NArts);



#endif // OPENWAVEFILE_H_INCLUDED
