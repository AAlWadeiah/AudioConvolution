#include <iostream>
#include <stdio.h>
#include <stdlib.h>
#include <stdint.h>
#include <string.h>
#include <math.h>
#include <fstream>
#include <cstdlib>
#include <vector>

// Sound file macros
#define DRUMS "DrumsDry.wav"
#define PIANO "Synthesized_Piano.wav"
#define HARPSICHORD "Synthesized_Harpsichord.wav"
#define CHANT "TibetanChant.wav"
#define GUITAR "GuitarDry.wav"
#define BIG_HALL_IR "BigHall_IR.wav"
#define TAJ_MAHAL_IR "TajMahal_IR.wav"
#define AUTO_PARK_IR "AutoPark_IR.wav"
#define OUTPUT_FILE "output.wav"

// Macros for calculations and writing headers
#define FREQUENCY         440.0
#define SAMPLE_RATE       44100.0
#define BITS_PER_SAMPLE   16
#define BYTES_PER_SAMPLE  (BITS_PER_SAMPLE/8)
#define MONOPHONIC        1
#define STEREOPHONIC      2

using namespace std;

struct wavHeader {
    char chunkID[4];
    unsigned int chunkSize;
    char format[4];
    char subChunk1ID[4];
    unsigned int subChunk1Size;
    unsigned short audioFormat;
    unsigned short numChannels;
    unsigned int sampleRate;			// sampleRate denotes the sampling rate.
    unsigned int byteRate;
    unsigned short blockAlign;
    unsigned short bitsPerSample;
} wavHeader;

struct dataHeader {
    char subChunk2ID[4];
    unsigned int subChunk2Size;			// subChunk2Size = NumSamples * NumChannels * BytesPerSample denotes the number of samples.
} dataHeader;

// Function prototypes
void readWav(const char* fileName, struct wavHeader* wHeader, struct dataHeader* dHeader, vector <double> &samples);
void writeWav(struct dataHeader* dHeader, vector <double> &samples);
void writeWavHeader(FILE* outfile, int numOfChannels, int soundSampleSize, double sampleRate);
size_t fwriteIntLSB(int data, FILE* fileStream);
size_t fwriteShortLSB(short int data, FILE* fileStream);
void displayWaveHeader(struct wavHeader* header);
void displayDataChunkHeader(struct dataHeader* header);
void timeConvolve(double x[], int N, double h[], int M, double y[], int P);

int main (int argc, char** argv) {
    // Sound file
    struct wavHeader soundFileHeader;
    struct dataHeader soundFileDataHeader;

    // Impulse response file
    struct wavHeader impulseFileHeader;
    struct dataHeader impulseFileDataHeader;

    vector<double> soundSamples;
    vector<double> impulseSamples;
    vector<double> outputSamples;

    cout << "Reading files\n";
    readWav(PIANO, &soundFileHeader, &soundFileDataHeader, soundSamples);
    readWav(AUTO_PARK_IR, &impulseFileHeader, &impulseFileDataHeader, impulseSamples);

    cout << "Convolving...\n";
    outputSamples.resize(impulseSamples.size() + soundSamples.size() - 1);
    timeConvolve(&soundSamples[0], soundSamples.size(), &impulseSamples[0], impulseSamples.size(), &outputSamples[0], outputSamples.size());
    
    cout << "Convolved. Writing to \"output.wav\"\n";
    writeWav(&soundFileDataHeader, outputSamples);

    return 0;
}

void timeConvolve(double x[], int N, double h[], int M, double y[], int P){
    int n, m;

    // Make sure the output buffer is the right size: P = N + M - 1
    if (P != (N + M - 1)) {
        printf("Output signal vector is the wrong size\n");
        printf("It is %-d, but should be %-d\n", P, (N + M - 1));
        printf("Aborting convolution\n");
        return;
    }

    // Clear the output buffer y[] to all zero values
    for (n = 0; n < P; n++) y[n] = 0.0;

    double max = 0;
    // Outer loop:  process each input value x[n] in turn
    for (n = 0; n < N; n++) {
        // Inner loop:  process x[n] with each sample of h[]
        for (m = 0; m < M; m++){
            y[n+m] += x[n] * h[m];
            if(abs(y[n+m]) > max) max = y[n+m];
        }
    }
    // Divide by max element in y to normalize y
    for(int i = 0; i < P; i++) y[i] /= max;
}

void writeWavHeader(FILE* outfile, int numOfChannels, int soundSampleSize, double sampleRate) {
    int numOfSamples = soundSampleSize / (numOfChannels * BYTES_PER_SAMPLE);
    int dataChunkSize = numOfChannels * numOfSamples * BYTES_PER_SAMPLE;
    int formSize = 36 + dataChunkSize;
    short int frameSize = numOfChannels * BYTES_PER_SAMPLE;
    int bytesPerSecond = (int)ceil(sampleRate * frameSize);
    
    // Writing values
    fputs("RIFF", outfile);
    fwriteIntLSB(formSize, outfile);
    fputs("WAVE", outfile);
    fputs("fmt ", outfile);
    fwriteIntLSB(16, outfile);
    fwriteShortLSB(1, outfile);
    fwriteShortLSB((short)numOfChannels, outfile);
    fwriteIntLSB((int)sampleRate, outfile);
    fwriteIntLSB(bytesPerSecond, outfile);
    fwriteShortLSB(frameSize, outfile);
    fwriteShortLSB(BITS_PER_SAMPLE, outfile);
    fputs("data", outfile);
    fwriteIntLSB(dataChunkSize, outfile);
}

size_t fwriteIntLSB(int data, FILE *stream) {
    unsigned char array[4];
    array[3] = (unsigned char)((data >> 24) & 0xFF);
    array[2] = (unsigned char)((data >> 16) & 0xFF);
    array[1] = (unsigned char)((data >> 8) & 0xFF);
    array[0] = (unsigned char)(data & 0xFF);
    return fwrite(array, sizeof(unsigned char), 4, stream);
}

size_t fwriteShortLSB(short int data, FILE *stream) {
    unsigned char array[2];
    array[1] = (unsigned char)((data >> 8) & 0xFF);
    array[0] = (unsigned char)(data & 0xFF);
    return fwrite(array, sizeof(unsigned char), 2, stream);
}

void writeWav(struct dataHeader* dHeader, vector <double> &samples) {
    FILE * outfile = fopen(OUTPUT_FILE, "wb+");
    int soundSampleSize = dHeader -> subChunk2Size;
    writeWavHeader(outfile, MONOPHONIC, soundSampleSize, SAMPLE_RATE);
    short sample_data[samples.size()];

    // write data to file
    for(int i = 0; i < samples.size(); i++) sample_data[i] = (short) (samples[i] * (double) (INT16_MAX - 1));

    fwrite(&sample_data, 1, samples.size(), outfile);
}

void readWav(const char* fileName, struct wavHeader * wHeader, struct dataHeader * dHeader, vector <double> &samples) {
    FILE * infile = fopen(fileName, "rb");

    // read header, skip 2 bytes if format chunk is 18 bytes
    fread(wHeader, sizeof(wavHeader), 1, infile);
    if(wHeader->subChunk1Size == 18) fseek(infile, 2, SEEK_CUR);
    fread(dHeader, sizeof(dataHeader), 1, infile);

    cout << "Making the array\n";
    int numOfSamples = dHeader -> subChunk2Size / (wHeader -> numChannels * BYTES_PER_SAMPLE);
    cout << "Number of samples " << numOfSamples << endl;
    short singleSample;
    
    for(int i = 0; i < dHeader -> subChunk2Size; i++) {
        fread(&singleSample, sizeof(short), 1, infile);
        samples.push_back((double) singleSample / (double) (INT16_MAX - 1));
    }
    displayWaveHeader(wHeader);
    displayDataChunkHeader(dHeader);
}

void displayWaveHeader(struct wavHeader * header) {
    printf("Size of header is %d bytes\n", sizeof(*header));
    printf("----------------------------------------------\n");
    printf("Chunk ID: %.4s\n", header -> chunkID);
    printf("Chunk Size: %d\n", header -> chunkSize);
    printf("Format: %.4s\n", header -> format);
    printf("Sub-Chunk 1 ID: %.4s\n", header -> subChunk1ID);
    printf("Sub-Chunk 1 Size: %d\n", header -> subChunk1Size);
    printf("Audio Format: %d\n", header -> audioFormat);
    printf("Number of Channels: %d\n", header -> numChannels);
    printf("Sample Rate: %d Hz\n", header -> sampleRate);
    printf("Byte Rate: %d\n", header -> byteRate);
    printf("Block Align: %d\n", header -> blockAlign);
    printf("Bits per Sample: %d\n", header -> bitsPerSample);
}

void displayDataChunkHeader(struct dataHeader * header) {
    printf("Size of data header is %d bytes\n", sizeof(*header));
    printf("----------------------------------------------\n");
    printf("Sub-Chunk 2 ID: %.4s\n", header -> subChunk2ID);
    printf("Sub-Chunk 2 Size: %d\n", header -> subChunk2Size);
}