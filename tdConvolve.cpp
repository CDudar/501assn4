#include <iostream>
#include <stdio.h>
using namespace std;


typedef struct  WAV_HEADER
{
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
    /* "data" sub-chunk */
    uint8_t         Subchunk2ID[4]; // "data"  string
    uint32_t        Subchunk2Size;  // Sampled data length
} wav_hdr;

// Function prototypes 
int getFileSize(FILE *inFile); 


void convolve(float x[], int N, float h[], int M, float y[], int P){
    
    int n, m;
    
    
    /*Clear output buffer y[]*/
    for(n = 0; n < P; n++){
        y[n] = 0.0;
    }
    
    /*Outer loop, process each input value x[n] in turn */
    
    for(n = 0; n < N; n++){
        /*Inner loop process x[n] with each sample of h[n]*/
        for(m = 0; m < M; m++){
            y[n + m] += x[n] * h[m];
        }
    }
        
}
    



int main(int argc, char** argv) 
{
    cout << argc << " arguments entered" << endl;
    for(int i = 0; i < argc; i++){
        cout << "arg " << i << ": " << argv[i] << endl;
    }
    cout << endl;
    
    /*Declare header struct and get header size */
    wav_hdr wavHeader;
    int headerSize = sizeof(wav_hdr);
    
    /* Open input wave file in read mode */
    cout << "Opening input wave file : " << argv[1] << endl;
    FILE * infile = fopen(argv[1], "rb");
    
    size_t bytesRead = fread(&wavHeader, 1, headerSize, infile);
    
    
    //Read the data
    uint16_t bytesPerSample = wavHeader.bitsPerSample / 8;      //Number of bytes per sample
    uint64_t numSamples = wavHeader.ChunkSize / bytesPerSample; //How many samples are in the wav file?
    static const uint16_t BUFFER_SIZE = 4096;
    int8_t* buffer = new int8_t[BUFFER_SIZE];
    while ((bytesRead = fread(buffer, sizeof buffer[0], BUFFER_SIZE / (sizeof buffer[0]), infile)) > 0)
    {
        cout << "Read " << bytesRead << " bytes." << endl;
    }
    delete [] buffer;
    buffer = nullptr;
    int filelength = getFileSize(infile);

    
    cout << "File is                    :" << filelength << " bytes." << endl;
    
    cout << "--- RIFF Chunk Descriptor ---" << endl;
    cout << "RIFF header                :" << wavHeader.RIFF[0] << wavHeader.RIFF[1] << wavHeader.RIFF[2] << wavHeader.RIFF[3] << endl;
    cout << "Chunk size                 :" << wavHeader.ChunkSize << endl;
    cout << "Format                     :" << wavHeader.WAVE[0] << wavHeader.WAVE[1] << wavHeader.WAVE[2] << wavHeader.WAVE[3] << endl;
    
    cout << "--- fmt sub-chunk -----------" << endl;
    cout << "Subchunk1ID                :" << "'" << wavHeader.fmt[0] << wavHeader.fmt[1] << wavHeader.fmt[2] << wavHeader.fmt[3] << "'" << endl;
    cout << "Subchunk1-Size             :" << wavHeader.Subchunk1Size << endl;
    cout << "Audio Format               :" << wavHeader.AudioFormat << endl;
    cout << "Number of channels         :" << wavHeader.NumOfChan << endl;
    cout << "Sampling Rate              :" << wavHeader.SamplesPerSec << endl;    
    cout << "Byte Rate                  :" << wavHeader.bytesPerSec << endl;
    cout << "Block Align                :" << wavHeader.blockAlign << endl;
    cout << "Bits per sample            :" << wavHeader.bitsPerSample << endl;
    
    cout << "-- data sub-chunk -----------" << endl;
    cout << "Subchunk2ID                :" << wavHeader.Subchunk2ID[0] << wavHeader.Subchunk2ID[1] << wavHeader.Subchunk2ID[2] << wavHeader.Subchunk2ID[3] << endl;
    cout << "Subchunk2Size (bytes)      :" << wavHeader.Subchunk2Size << endl;

    fclose(infile);

    return 0;
}

// find the file size
int getFileSize(FILE* inFile)
{
    int fileSize = 0;
    fseek(inFile, 0, SEEK_END);

    fileSize = ftell(inFile);

    fseek(inFile, 0, SEEK_SET);
    return fileSize;
}