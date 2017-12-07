#include <unistd.h>
#include <stdint.h>
#include <iostream>
#include <stdio.h>
#include <string>
#include <fstream>

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

void PrintWaveHeader(WAV_HEADER);
bool equalDataChunks(double*, double*); 

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
    
    dataArray = new double[Subchunk2Size / (bitsPerSample / 8)];
    
    int16_t sample;
    
    for(int i = 0; i < Subchunk2Size / (bitsPerSample / 8); i++){
        
        input.read((char*)&sample, 2);
        
        double converted = (double) sample / (double)INT16_MAX;
        
        if(converted < -1.0){
            converted = -1.0;
        }
        
        dataArray[i] = converted;
        
    }
    
    input.close();
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

bool equalDataChunks(double tdData[], double fftData[], int length){
    
    cout << "comparing..." << endl;
    
    bool equalDataValues = true;
    
    for(int i = 0; i < length; i++){
        if(tdData[i] != fftData[i]){
            equalDataValues = false;
        }
    }
    
    return equalDataValues;
}



int main(int argc, char** argv) 
{
    cout << argc << " arguments entered" << endl;
    for(int i = 0; i < argc; i++){
        cout << "arg " << i << ": " << argv[i] << endl;
    }
    cout << endl;
            
    cout << "Reading in file: " << argv[1] << endl;
    WAV_HEADER tdConvolve(argv[1]);
    
    cout << "Reading in file: " << argv[2] << endl;
    WAV_HEADER fftConvolve(argv[2]);
    
    
    PrintWaveHeader(tdConvolve);
    cout << "\n" << endl;
    
    PrintWaveHeader(fftConvolve);
    cout << "\n" << endl;
    
    cout << "Done reading in waveFiles" << endl;    
    
    
    cout << "Comparing file contents for equality" << endl;
    
    bool equal = true;
    
    cout << ".\n.\n." << endl;
    
    if(tdConvolve.Subchunk2Size == fftConvolve.Subchunk2Size){
        cout << "Data chunks lengths are equivalent" << endl;
    
        int dataLength = (tdConvolve.Subchunk2Size / (tdConvolve.bitsPerSample/8));
    
        if(equalDataChunks(tdConvolve.dataArray, fftConvolve.dataArray, dataLength)){
        cout << "All data chunk values are equal" << endl;
        }
        else{
        equal = false;
        }
    
    
    }
    else{
        equal = false;
        
    }
    
    
    if(equal){
        cout << "WAV files are identical" << endl;
    }
    else{
        cout << "WAV files are NOT identical" << endl;
    }
    
    return 0;
}

