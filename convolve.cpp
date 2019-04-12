#include <iostream>
#include <stdio.h>
#include <stdlib.h>
#include <stdint.h>
#include <string.h>
#include <math.h>
#include <fstream>
#include <cstdlib>
#include <vector>

#define FFT                 1
#define INVERSE_FFT         -1

// Macros for calculations and writing headers
#define FREQUENCY           440.0
#define SAMPLE_RATE         44100.0
#define BITS_PER_SAMPLE     16
#define BYTES_PER_SAMPLE    (BITS_PER_SAMPLE/8)
#define MONOPHONIC          1
#define STEREOPHONIC        2
#define PI                  3.141592653589793
#define TWO_PI              (2.0 * PI)
#define SWAP(a, b)          tempr = (a); (a) = (b); (b) = tempr

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
void writeWav(struct dataHeader* dHeader, vector <double> &samples, const char* outFileName);
void writeWavHeader(FILE* outfile, int numOfChannels, int soundSampleSize, double sampleRate);
size_t fwriteIntLSB(int data, FILE* fileStream);
size_t fwriteShortLSB(short int data, FILE* fileStream);
void displayWaveHeader(struct wavHeader* header);
void displayDataChunkHeader(struct dataHeader* header);
void timeConvolve(double x[], int N, double h[], int M, double y[], int P);
void fastFourierTransform(double data[], int nn, int isign);
void FFTConvolve(double x[], int N, double h[], int M, double y[], int P);
void timeToFrequency(double timeSignal[], unsigned int timeSize, double frequencySignal[], unsigned int frequencySize);
void multiplySignals(double drySignal[], double irSignal[], double outputSignal[], unsigned int complexSignalSize);

int main (int argc, char* argv[]) {
    if (argc != 4) {
        cout << "Usage: ./convolve [inputFile] [impulseResponseFile] [outputFile]" << endl;
        exit(0);
    }

    const char* inputFile = argv[1];
    const char* IRFile = argv[2];
    const char* outputFile = argv[3];

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
    readWav(inputFile, &soundFileHeader, &soundFileDataHeader, soundSamples);
    readWav(IRFile, &impulseFileHeader, &impulseFileDataHeader, impulseSamples);

    cout << "Convolving...\n";
    outputSamples.resize(impulseSamples.size() + soundSamples.size() - 1);
    FFTConvolve(&soundSamples[0], soundSamples.size(), &impulseSamples[0], impulseSamples.size(), &outputSamples[0], outputSamples.size());
    
    cout << "Convolved. Writing to \"" << outputFile << "\"\n";
    writeWav(&soundFileDataHeader, outputSamples, outputFile);
    cout << "written" << endl;
}

void timeToFrequency(double timeSignal[], unsigned int timeSize, double frequencySignal[], unsigned int frequencySize) {
    unsigned int i, j;
    for (i = 0, j = 0; i < timeSize; i+= 2, j += 4) {
        frequencySignal[j] = timeSignal[i];
        frequencySignal[j + 1] = 0.0;
        frequencySignal[j + 3] = timeSignal[i + 1];
        frequencySignal[j + 4] = 0.0;
    }
}

void multiplySignals(double drySignal[], double irSignal[], double outputSignal[], unsigned int complexSignalSize) {
    for (int i = 0; i < complexSignalSize; i += 2) { 
		outputSignal[i] = (drySignal[i] * irSignal[i]) - (drySignal[i + 1] * irSignal[i + 1]);
		outputSignal[i + 1] = (drySignal[i] * irSignal[i + 1]) + (drySignal[i + 1] * irSignal[i]);
	}
}

void FFTConvolve(double x[], int N, double h[], int M, double y[], int P) {
    unsigned int signalSize, complexSignalSize;

    signalSize = 1;
    while (signalSize < P) {
        signalSize <<= 1;
    }

    // unsigned int temp = static_cast<unsigned int>(log2(P));
    // signalSize = static_cast<unsigned int>(static_cast<int>(pow(2, temp)));

    // signalSize = 1;
    // for (int i=0; i<8*sizeof(unsigned int); i++) 
    // { 
    //     unsigned int curr = 1 << i; 
  
    //     // If current power is more than n, break 
    //     if (curr > P) 
    //        break; 
  
    //     signalSize = curr; 
    // }

    complexSignalSize = signalSize * 2;
    
    double* xComplexSignal = new double[complexSignalSize];
    double* hComplexSignal = new double[complexSignalSize];
    double* yComplexSignal = new double[complexSignalSize];

    timeToFrequency(x, N, xComplexSignal, signalSize);
    timeToFrequency(h, M, hComplexSignal, signalSize);

    fastFourierTransform((xComplexSignal - 1), signalSize, FFT);
    fastFourierTransform((hComplexSignal - 1), signalSize, FFT);

    multiplySignals(xComplexSignal, hComplexSignal, yComplexSignal, complexSignalSize);

    fastFourierTransform((yComplexSignal - 1), signalSize, INVERSE_FFT);

    for (int i = 0; i < complexSignalSize; i++)
        yComplexSignal[i] /= (double) complexSignalSize;
    
    double max = 0;
    for (int i = 0; i < P; i++) {
        y[i] = yComplexSignal[i*2];
        if(abs(y[i]) > max) max = y[i];
    }

    for(int i = 0; i < P; i++) y[i] /= max;
    
    delete[] xComplexSignal;
    delete[] hComplexSignal;
    delete[] yComplexSignal;
}

void fastFourierTransform(double data[], int nn, int isign) {
    unsigned long n, mmax, m, j, istep, i;
    double wtemp, wr, wpr, wpi, wi, theta, tempr, tempi;

    n = nn << 1;
    j = 1;

    for (i = 1; i < n; i += 2) {
        if (j > i) {
            SWAP(data[j], data[i]);
            SWAP(data[j+1], data[i+1]);
        }
        m = nn;
        while (m >= 2 && j > m) {
            j -= m;
            m >>= 1;
        }
        j += m;
    }

    mmax = 2;
    while (n > mmax) {
        istep = mmax << 1;
        theta = isign * (6.28318530717959 / mmax);
        wtemp = sin(0.5 * theta);
        wpr = -2.0 * wtemp * wtemp;
        wpi = sin(theta);
        wr = 1.0;
        wi = 0.0;
        for (m = 1; m < mmax; m += 2) {
            for (i = m; i <= n; i += istep) {
                j = i + mmax;
                tempr = wr * data[j] - wi * data[j+1];
                tempi = wr * data[j+1] + wi * data[j];
                data[j] = data[i] - tempr;
                data[j+1] = data[i+1] - tempi;
                data[i] += tempr;
                data[i+1] += tempi;
            }
            wr = (wtemp = wr) * wpr - wi * wpi + wr;
            wi = wi * wpr + wtemp * wpi + wi;
        }
        mmax = istep;
    }
}

void timeConvolve(double x[], int N, double h[], int M, double y[], int P) {
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
    // cout << "1" << endl;
    fwriteIntLSB(formSize, outfile);
    // cout << "2" << endl;
    fputs("WAVE", outfile);
    // cout << "3" << endl;
    fputs("fmt ", outfile);
    // cout << "4" << endl;
    fwriteIntLSB(16, outfile);
    // cout << "5" << endl;
    fwriteShortLSB(1, outfile);
    // cout << "6" << endl;
    fwriteShortLSB((short)numOfChannels, outfile);
    // cout << "7" << endl;
    fwriteIntLSB((int)sampleRate, outfile);
    // cout << "8" << endl;
    fwriteIntLSB(bytesPerSecond, outfile);
    // cout << "9" << endl;
    fwriteShortLSB(frameSize, outfile);
    // cout << "10" << endl;
    fwriteShortLSB(BITS_PER_SAMPLE, outfile);
    // cout << "11" << endl;
    fputs("data", outfile);
    // cout << "12" << endl;
    fwriteIntLSB(dataChunkSize, outfile);
    // cout << "13" << endl;
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

void writeWav(struct dataHeader* dHeader, vector <double> &samples, const char* outFileName) {
    FILE * outfile = fopen(outFileName, "wb+");
    int soundSampleSize = dHeader -> subChunk2Size;
    writeWavHeader(outfile, MONOPHONIC, soundSampleSize, SAMPLE_RATE);
    int sampleDataSize = samples.size();
    short sample_data[sampleDataSize];
    
    // write data to file
    for(int i = 0; i < sampleDataSize; i++) sample_data[i] = (short) (samples[i] * (double) (INT16_MAX - 1));

    fwrite(&sample_data, 1, sampleDataSize, outfile);
}

void readWav(const char* fileName, struct wavHeader * wHeader, struct dataHeader * dHeader, vector <double> &samples) {
    FILE * infile = fopen(fileName, "rb");

    // read header, skip 2 bytes if format chunk is 18 bytes
    fread(wHeader, sizeof(wavHeader), 1, infile);
    if(wHeader->subChunk1Size == 18) fseek(infile, 2, SEEK_CUR);
    fread(dHeader, sizeof(dataHeader), 1, infile);

    // int numOfSamples = dHeader -> subChunk2Size / (wHeader -> numChannels * BYTES_PER_SAMPLE);
    short singleSample;
    
    for(int i = 0; i < dHeader -> subChunk2Size; i++) {
        fread(&singleSample, sizeof(short), 1, infile);
        samples.push_back((double) singleSample / (double) (INT16_MAX - 1));
    }
    displayWaveHeader(wHeader);
    displayDataChunkHeader(dHeader);
}

void displayWaveHeader(struct wavHeader * header) {
    printf("^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^\n");
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
    printf("\nSize of data header is %d bytes\n", sizeof(*header));
    printf("----------------------------------------------\n");
    printf("Sub-Chunk 2 ID: %.4s\n", header -> subChunk2ID);
    printf("Sub-Chunk 2 Size: %d\n", header -> subChunk2Size);
}