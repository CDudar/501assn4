#include <iostream>
#include <stdio.h>
#include <string>
#include <fstream>
#include <unistd.h>
#include <stdint.h>
#include <math.h>


#define SWAP(a,b)  tempr=(a);(a)=(b);(b)=tempr

using namespace std;


class WAV_HEADER
{
    
    public:
    
        /* RIFF Chunk Descriptor */
        uint8_t         RIFF[4];        // RIFF Header Magic header
        uint32_t        ChunkSize;      // RIFF Chunk Size
        uint8_t         WAVE[4];        // WAVE Header
        /* "fmt" sub-chunk */
        uint8_t         fmt[4];         // FMT header
        uint32_t        Subchunk1Size;  // Size of the fmt chunk
        uint16_t        AudioFormat;    // Audio format 1=PCM
        uint16_t        NumOfChan;      // Number of channels 1=Mono 2=Sterio
        uint32_t        SamplesPerSec;  // Sampling Frequency in Hz
        uint32_t        bytesPerSec;    // bytes per second
        uint16_t        blockAlign;     // 2=16-bit mono, 4=16-bit stereo
        uint16_t        bitsPerSample;  // Number of bits per sample
        
        uint8_t*        extraData;
        /* "data" sub-chunk */
        uint8_t         Subchunk2ID[4]; // "data"  string
        uint32_t        Subchunk2Size;  // Sampled data length
        
        double*         dataArray;
        WAV_HEADER(string);
    
};


int resultLength;

//double* result;
void four1(double*, int, int);
void convolve(double*, int, double*, int, double*, int);
void writeWAVFile(string, double*, int, int);
void PrintWaveHeader(WAV_HEADER);


WAV_HEADER::WAV_HEADER(string fileName){
    ifstream input;
    
    input.open(fileName, ios::binary);
    
    input.read((char*)&RIFF, 4);
    
    input.read((char*)&ChunkSize, 4);
    
    input.read((char*)&WAVE, 4);
    
    input.read((char*)&fmt, 4);
    
    input.read((char*)&Subchunk1Size, 4);
    
    input.read((char*)&AudioFormat, 2);
    
    input.read((char*)&NumOfChan, 2);
    
    input.read((char*)&SamplesPerSec, 4);
    
    input.read((char*)&bytesPerSec, 4);
    
    input.read((char*)&blockAlign, 2);
    
    input.read((char*)&bitsPerSample, 2);
    
    extraData = new uint8_t[Subchunk1Size - 16];
    input.read((char*)&extraData, Subchunk1Size-16);
    
    input.read((char*)&Subchunk2ID, 4);
    
    input.read((char*)&Subchunk2Size, 4);
    
    
    int numSamples = Subchunk2Size / (bitsPerSample / 8);
    
    dataArray = new double[numSamples];
    
    int16_t sample;
    
    for(int i = 0; i < numSamples; i++){
        
        input.read((char*)&sample, 2);
        
        double converted = (double) sample / (double)INT16_MAX;
        
        if(converted < -1.0){
            converted = -1.0;
        }
        
        dataArray[i] = converted;
        
    }
    
    input.close();
}



/*Following algorithm taken from class handout */

//  The four1 FFT from Numerical Recipes in C,
//  p. 507 - 508.
//  Note:  changed float data types to double.
//  nn must be a power of 2, and use +1 for
//  isign for an FFT, and -1 for the Inverse FFT.
//  The data is complex, so the array size must be
//  nn*2. This code assumes the array starts
//  at index 1, not 0, so subtract 1 when
//  calling the routine (see main() below).

void four1(double data[], int nn, int isign)
{
    unsigned long n, mmax, m, j, istep, i;
    double wtemp, wr, wpr, wpi, wi, theta;
    double tempr, tempi;

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



void convolve(double x[], int N, double h[], int M, double y[], int P){
    

    
    cout << "entering convolve function" << endl;

    cout << "N: " << N << endl;
    cout << "M: " << M << endl;

    int longestInput = max(N, M);
    
    cout << "longest input: " << longestInput << endl;
    
    
    //powerOf2 AKA nn
    int powerOf2 = 1;
    while(longestInput > powerOf2){
        powerOf2 *= 2;
    }
    
    int actualSize = powerOf2 * 2;
    
    cout << "power of 2: " << powerOf2 << endl;
    
    cout << "actual size: " << actualSize << endl;
    
    double* inputDataPadded = new double[actualSize];
    
    double* impulseDataPadded = new double[actualSize];
    
    
    //Zero-padding audio and impulse files
    
    
    for (int i = 0; i < actualSize; i = i + 2){
        inputDataPadded[i] = 0;
        inputDataPadded[i + 1] = 0;
        impulseDataPadded[i] = 0;
        impulseDataPadded[i + 1] = 0;
    }
    cout << "zeroed input data" << endl;
    
    for(int i = 0; i < N; i++){
        inputDataPadded[i * 2] = x[i];
    }

    cout << "put input audio data into padded array" << endl;

    
    for(int i = 0; i < M; i++){
        impulseDataPadded[i * 2] = h[i];
    }
    

    
    cout << "entering fft algo" << endl;
    
    //running through fft
    four1(inputDataPadded - 1, powerOf2, 1);
    four1(impulseDataPadded - 1, powerOf2, 1);
    
    cout << "exiting fft algo" << endl;
    
    double* outputData = new double[actualSize];
    
    for(int i = 0; i < actualSize; i = i + 2){
        outputData[i] = (inputDataPadded[i] * impulseDataPadded[i]) - (inputDataPadded[i + 1] * impulseDataPadded[i + 1]);
        outputData[i + 1] = (inputDataPadded[i + 1] * impulseDataPadded[i]) + (inputDataPadded[i] * impulseDataPadded[i + 1]);
    }
    
    cout << "inversing fft" << endl;
    four1(outputData-1, powerOf2, -1);
    
    cout << "done inversing fft" << endl;


    double max = 0.0;
    for(int i = 0 ; i < resultLength; i++){
        y[i] = outputData[i * 2];
    
        if((abs(outputData[i * 2])) >  max)
            max = abs(outputData[i * 2]);
    }
    
    for(int i = 0; i < resultLength; i++){
        y[i] = y[i] / max;
    }
        
}




void writeWAVFile(string outputFile, double outputData[], int outputLength, int sampleRate){
    ofstream outputStream;
    outputStream.open(outputFile, ios::binary | ios::out);
    
    
    char* ChunkID = new char[4]{'R', 'I', 'F', 'F'};
    outputStream.write(ChunkID, 4);
    
    int dataBytesLength = outputLength * 2;
    
    int chunkSize = dataBytesLength + 36;
    
    outputStream.write((char*)&chunkSize, 4);
    
    char* format = new char[4]{'W', 'A', 'V', 'E'};
    outputStream.write(format, 4);
    
    char* subChunk1ID = new char[4]{'f', 'm', 't', ' '};
    outputStream.write(subChunk1ID, 4);
    
    int subChunk1Size = 16;
    outputStream.write((char*)&subChunk1Size, 4);
    
    short audioFmt = 1;
    outputStream.write((char*)&audioFmt, 2);

    short numChannels = 1;
    outputStream.write((char*)&numChannels, 2);
    
    outputStream.write((char*)&sampleRate, 4);
    
    
    short bitsPerSample = 16;
    
    int byteRate = sampleRate * numChannels * bitsPerSample/8;
    
    outputStream.write((char*)&byteRate, 4);
    
    short blockAlign = numChannels * (bitsPerSample / 8);
    
    outputStream.write((char*)&blockAlign, 2);

    outputStream.write((char*)&bitsPerSample, 2);
    
    char* subChunk2ID = new char[4]{'d', 'a', 't', 'a'};
    outputStream.write(subChunk2ID, 4);
    
    outputStream.write((char*)&dataBytesLength, 4);
    
    int16_t outData;
    for(int i = 0; i < outputLength; i++){
        
        outData = (int16_t) (outputData[i] * INT16_MAX);
                
        outputStream.write((char*)&outData, 2);
        
    }
    
    outputStream.close();
    
}


void PrintWaveHeader(WAV_HEADER wav_hdr){
    
        
    cout << "--- RIFF Chunk Descriptor ---" << endl;
    cout << "RIFF header                :" << wav_hdr.RIFF[0] << wav_hdr.RIFF[1] << wav_hdr.RIFF[2] << wav_hdr.RIFF[3] << endl;
    cout << "Chunk size                 :" << wav_hdr.ChunkSize << endl;
    cout << "Format                     :" << wav_hdr.WAVE[0] << wav_hdr.WAVE[1] << wav_hdr.WAVE[2] << wav_hdr.WAVE[3] << endl;
        
    cout << "--- fmt sub-chunk -----------" << endl;
    cout << "Subchunk1ID                :" << "'" << wav_hdr.fmt[0] << wav_hdr.fmt[1] << wav_hdr.fmt[2] << wav_hdr.fmt[3] << "'" << endl;
    cout << "Subchunk1-Size             :" << wav_hdr.Subchunk1Size << endl;
    cout << "Audio Format               :" << wav_hdr.AudioFormat << endl;
    cout << "Number of channels         :" << wav_hdr.NumOfChan << endl;
    cout << "Sampling Rate              :" << wav_hdr.SamplesPerSec << endl;    
    cout << "Byte Rate                  :" << wav_hdr.bytesPerSec << endl;
    cout << "Block Align                :" << wav_hdr.blockAlign << endl;
    cout << "Bits per sample            :" << wav_hdr.bitsPerSample << endl;
    
    cout << "-- data sub-chunk -----------" << endl;
    cout << "Subchunk2ID                :" << wav_hdr.Subchunk2ID[0] << wav_hdr.Subchunk2ID[1] << wav_hdr.Subchunk2ID[2] << wav_hdr.Subchunk2ID[3] << endl;
    cout << "Subchunk2Size (bytes)      :" << wav_hdr.Subchunk2Size << endl;
    
    
    
}

int main(int argc, char** argv) 
{
    cout << argc << " arguments entered" << endl;
    for(int i = 0; i < argc; i++){
        cout << "arg " << i << ": " << argv[i] << endl;
    }
    cout << endl;
            
    cout << "Reading in file: " << argv[1] << endl;
    WAV_HEADER inputAudio(argv[1]);
    
    cout << "Reading in file: " << argv[2] << endl;
    WAV_HEADER IRAudio(argv[2]);
    
    
    PrintWaveHeader(inputAudio);
    cout << "\n" << endl;
    
    PrintWaveHeader(IRAudio);
    cout << "\n" << endl;
    
    cout << "Done reading in waveFiles" << endl;    
    
    
    int inputDataLength = (inputAudio.Subchunk2Size / (inputAudio.bitsPerSample/8));
    
    int IRDataLength = (IRAudio.Subchunk2Size / (IRAudio.bitsPerSample/8));
    
    resultLength = inputDataLength + IRDataLength - 1;
    
    double* result = new double[resultLength];
    
    cout << "Convolving" << endl;
    convolve(inputAudio.dataArray, inputDataLength, IRAudio.dataArray, IRDataLength, result, resultLength);
    cout << "Finished convolving" << endl;
        
    
    cout << "writing wav file" << endl;
    writeWAVFile(argv[3], result, resultLength, inputAudio.SamplesPerSec);
    cout << "finished writing wav" << endl;
    
    cout << "Finished everything" << endl;
    
    return 0;
}

