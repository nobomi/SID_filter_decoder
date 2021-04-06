#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <string.h>
#include <stdint.h>
#include <complex.h>
#include "pocketfft/pocketfft.h"

typedef double complex cplx;

#define BLOCK   16384
#define WINDOW  16384
#define NWIN    (BLOCK/WINDOW)

#define NOISE_LEVEL 10.0    //-10dB

//////////////////////// WAV /////////////////////////

typedef struct __attribute__((packed)) {
    int8_t riff_header[4];
    uint32_t wav_size;
    int8_t wave_header[4];
    int8_t fmt_header[4];
    uint32_t fmt_chunk_size;
    short audio_format;
    short num_channels;
    uint32_t sample_rate;
    uint32_t byte_rate;
    short sample_alignment;
    short bit_depth;
    int8_t data_header[4];
    uint32_t data_bytes;
} t_wav_header;

#define IS_HEADER_ERR_RIFF  -1
#define IS_HEADER_ERR_MONO  -2
#define IS_HEADER_ERR_SPEED -3
#define IS_HEADER_ERR_BITS  -4

#define IS_HEADER_48MONO    1

static int leftright=-1;

static int is_header(char *name, t_wav_header *h) {
    if (memcmp(h->riff_header,"RIFF",4)) return IS_HEADER_ERR_RIFF;
    if (memcmp(h->wave_header,"WAVE",4)) return IS_HEADER_ERR_RIFF;
    if (memcmp(h->fmt_header, "fmt ",4)) return IS_HEADER_ERR_RIFF;
    if (memcmp(h->data_header,"data",4)) return IS_HEADER_ERR_RIFF;
    fprintf(stderr,"#WAV file %s sample rate %d/s depth %d bits channels %d\r\n",name,h->sample_rate,h->bit_depth,h->num_channels);
    if (h->num_channels==1) {
        leftright=-1;
    } else if ((h->num_channels!=2) || (leftright<0)) {
        if (h->num_channels==2) fprintf(stderr,"#WAV for a stereo file please specify left or right channel by -L or -R switch\r\n");
        return IS_HEADER_ERR_MONO;
    }
    if (h->sample_rate!=48000)  return IS_HEADER_ERR_SPEED;
    if (h->bit_depth!=16) return IS_HEADER_ERR_BITS;
    return IS_HEADER_48MONO;
}

size_t a_fread(void *buf, size_t vzorek, size_t len, FILE *f) {
    size_t res=0;
    size_t i,r;
    uint32_t u;
    if (leftright<0) {
        return fread(buf,vzorek,len,f);
    } else {
        if (vzorek!=2) {
            fprintf(stderr,"error: fread bad size\n");
            exit(-1);
        }
        for(i=0;i<len;i++) {
            r=fread(&u,4,1,f);
            if (r<0) return r;
            res+=r;
            memcpy(&((unsigned char *)buf)[i*2],&((unsigned char *)&u)[2*leftright],2);
        }
    }
    return res;
}

size_t addr_wav_header=sizeof(t_wav_header);
size_t a_fseek(FILE *f, long addr, int sset) {
    if (leftright<0) {
        return fseek(f, addr, sset);
    } else {
        if (sset!=SEEK_SET) {
            fprintf(stderr,"error: fseek not set\n");
            exit(-1);
        }
        if (addr>=addr_wav_header/2) return fseek(f, addr*2-addr_wav_header, sset)/2+addr_wav_header;
        else return fseek(f, addr*2, sset)/2;
    }
}

//////////////////////////////////////////////////////

static int last_n=0;
static rfft_plan plan;

static void FFT1(const double inputReal[], double outputreal[], const int n)
{
    memcpy (outputreal,inputReal,n*sizeof(double));
    if (last_n!=n) {
        if (last_n!=0) destroy_rfft_plan (plan);
        plan = make_rfft_plan(n);
    }
    rfft_forward (plan, outputreal, 1./n);
    for(int i=0;i<n/2;i++) outputreal[i]=sqrt(outputreal[i*2]*outputreal[i*2]+outputreal[i*2+1]*outputreal[i*2+1]);
}

void get_ref_kvadrat(const double inreal1[], const double inreal2[], double outreal[], int n)
{
    int i;
    for(i=0;i<n/2;i++) {
        if (inreal1[i]>inreal2[i]) outreal[i]=inreal1[i]*inreal1[i];
        else outreal[i]=inreal2[i]*inreal2[i];
    }
}

void Denoise2(double inreal1[], double inreal2[], const double refreal[], int n)
{
    int i,x=0;
    for(i=0;i<n/2;i++) {
        if (inreal2[i]*inreal2[i] < NOISE_LEVEL*refreal[i]) { // min. odstup
            inreal1[i]=0.0;x++;
        }
    }
}

void Filtr(double inreal[], double outreal[], int n)
{
    float K;
	int k;
	n/=2;
	int kleft,kright;
    double sum=0;
    double sum_w=0;
    K=1.414;
	k=1;kleft=1;kright=2;sum_w=1.0;sum=inreal[1]*inreal[1];outreal[1]=inreal[1];
	for (k = 2; k < n; k++) {  // For each output element
	    int i;
        if (n<=(float)k*K) K=(float)n/(float)k;
        if (256<=(float)k*(K-1.0)) K=(float)256/(float)k+1.0; // omezit sirku pasma nad cca 600Hz
	    for (i=kleft;i<k/K;i++) {
            sum_w-=1.0/(double)i;
            sum-=inreal[i]*inreal[i]/(double)i;
	    }
	    kleft=i;
	    for (i=kright;((i<=k*K) && (i<n));i++) {
            sum_w+=1.0/(double)i;
            sum+=inreal[i]*inreal[i]/(double)i;
	    }
	    kright=i;
        if (sum_w>0) outreal[k]=sqrt(sum/sum_w); else outreal[k]=0;
	}
}

float imin_mem=WINDOW/2048;

static float NajdiBod(const double x1[], const double x2[], const int n) {
    float imin=n/2048;
    int start1,start2,left1,left2,right1,right2,i;
    float max1,max2;
    // dej tam kvadrat
    static double k1[16384];
    static double k2[16384];
    for(i=0;i<n/2;i++) {
        k1[i]=x1[i];//*x1[i];
        k2[i]=x2[i];//*x2[i];
    }
    // najdi levy max
    start1=start2=n;
    max1=max2=0.0;
    left1=left2=0;
    for(i=(int)imin;i<n*15/32;i++) {
        if (k1[i]>k2[i]) {
            if (start1>i) start1=i;
            else if ((float)i/(float)start1 > max1) {
                max1=(float)i/(float)start1;
                left1=i;
            }
        } else start1=n;
        if (k1[i]>2.0*k2[i]) {
            if (start2>i) start2=i;
            else if ((float)i/(float)start2 > max2) {
                max2=(float)i/(float)start2;
                left2=i;
            }
        } else start2=n;
    }
    // najdi pravy min
    start1=start2=0;
    max1=max2=0.0;
    right1=right2=0;
    for(i=n*15/32-1;i>=(int)imin;i--) {
        if (k2[i]>k1[i]) {
            if (start1<i) start1=i;
            else if ((float)start1/(float)i > max1) {
                max1=(float)start1/(float)i;
                right1=i;
            }
        } else start1=0;
        if (k2[i]>2.0*k1[i]) {
            if (start2<i) start2=i;
            else if ((float)start2/(float)i > max2) {
                max2=(float)start2/(float)i;
                right2=i;
            }
        } else start2=0;
    }
    // prusecik nekde mezi ...
    if ((left1!=0)&&(right1!=0)&&(left1<right1)) { left2=left1; right2=right1; }
    if (right2==0) right2=n*15/32;
    if (left2==0) left2=imin;
    {
        double weight=0;
        double cnt_w=0;
        double cnt_sum=0;
        for(i=left2;i<=right2;i++) {
            double a=x1[i];
            double b=x2[i];
            if ((a==0.0) || (b==0.0)) continue;
            if (a>b) a=b/a;
            else a=a/b;
            if (a>0.5) {
                weight=(a-0.5)/0.5;
                weight=weight*weight;
                weight/=(float)i;
                if (x1[i]>x2[i]) cnt_sum+=(i/a)*weight;
                else cnt_sum+=(i*a)*weight;
                cnt_w+=weight;
            }
        }
        if (cnt_w>0) imin=cnt_sum/cnt_w;
        else imin=imin_mem;
    }
    imin_mem=imin;
    return imin*48000/n;
}

static void FFT(const double inputReal[], double outputreal[], const int n, const int N) {
    int i;
    for(i=0;i<n*N;i+=n) {
        FFT1(&inputReal[i], &outputreal[i], n);
        if (i>0) {
            int j;
            for(j=0;j<n;j++) outputreal[j]+=outputreal[i+j];
        }
    }
}

static void Normalize(double x[], double y[],double ref[],int n)
{
    int i=0;
    double r;
    for(i=0;i<n/2;i++) {
        r=sqrt(x[i]*x[i]+y[i]*y[i]);
        if (r>0.0) { x[i]*=ref[i]/r; y[i]*=ref[i]/r; }
    }
}

static short shortbuf[WINDOW];
static double inreal[WINDOW];
static double outreal1[WINDOW];
static double outreal2[WINDOW];
static double x1[WINDOW];
static double x2[WINDOW];
static double xref[WINDOW];

static double outreal[2][256][WINDOW/2];
static double outref[WINDOW/2];
static double x0[WINDOW/2];

int istone_and_not(double fr, double fr2) {
    int i,i2;
    int j;
    double sum=0;
    double sum1=0;
    double sum2=0;
    FFT(inreal,x1,1024,1);
    for (j=10;j<1024/2;j++) {
        sum+=x1[j]*x1[j];
    }
    sum=sqrt(sum);
    i=1024*fr/48000;
    i2=1024*fr2/48000;
    for (j=i-2;j<=i+2;j++) sum1+=x1[j]*x1[j];
    sum1=sqrt(sum1);
    for (j=i2-2;j<=i2+2;j++) sum2+=x1[j]*x1[j];
    sum2=sqrt(sum2);
    if ((sum1>500.0)&&(sum2<0.2*sum1)&&(sum1>0.8*sum)) return 1;
    return 0;
}

#define T0  ((int)((1.100)*48000))
#define T1  (((int)((0.340)*48000-8192)))
#define T4  ((int)((0.500)*48000))

int timeout(float t)
{
    return (int)(t*48000/1024);
}

int main(int argc, char *argv[])
{
    int seeky[256+1];
    int pokusy;
    FILE *f;
    FILE *fw;
    char *exe=argv[0];
    while (argc>3) {
        argv++;
        argc--;
        if (strcmp("-L",argv[0])==0) {
            leftright=0;
        } else if (strcmp("-R",argv[0])==0) {
            leftright=1;
        } else {
            argc=0;
        }
    }
    if (argc<3) {
        printf("usage: %s [-L] [-R] input_wav_audio_file output_text_file\n",exe);
        printf("  input_wav_audio_file = name of the audio sample file (48kHz/16bit wav)\n");
        printf("  output_text_file = name of the output text file\n");
        printf("  -L = for a stereo file, data in the left channel\n");
        printf("  -R = for a stereo file, data in the right channel\n");
        printf("  no extra parameter is needed for mono sampled files\n");
        return 0;
    }
    f=fopen(argv[1],"rb");
    if (!f) {
        fprintf(stderr,"error: input file can't open\n");
        return -1;
    } else {
        int i,j,rd,ad=0;
        t_wav_header header;

        rd=fread(&header,sizeof(header)-8,1,f);
        if (rd<1) {
            fclose(f);
            fprintf(stderr,"\r\nerror: file %s too short\n", argv[1]);
            return -1;
        }
        while (header.fmt_chunk_size>16) {
            int ch;
            addr_wav_header+=2;
            header.fmt_chunk_size-=2;
            ch=fgetc(f); if (ch!=EOF) ch=fgetc(f);
            ad++;
            if (ch==EOF) {
                fclose(f);
                fprintf(stderr,"\r\nerror: file %s format error, incomplete chunk\n", argv[1]);
                return -1;
            }
        }
        rd=fread(&header.data_header[0],8,1,f);
        if (rd<1) {
            fclose(f);
            fprintf(stderr,"\r\nerror: file %s too short\n", argv[1]);
            return -1;
        }
        ad+=sizeof(header)/2;
        switch (is_header(argv[1],&header)) {
            case IS_HEADER_ERR_RIFF:
                fprintf(stderr,"error: not a WAV file file %s\n", argv[1]);
                return -1;
            case IS_HEADER_ERR_SPEED:
                fprintf(stderr,"error: unsupported record speed of the audio file %s\n", argv[1]);
                return -1;
            case IS_HEADER_ERR_BITS:
                fprintf(stderr,"error: unsupported number of bits of the audio file %s\n", argv[1]);
                return -1;
            case IS_HEADER_48MONO:
                break;
            case IS_HEADER_ERR_MONO:
            default:
                fprintf(stderr,"error: unsupported audio file %s\n", argv[1]);
                return -1;
        }

        if (strcmp("-",argv[2])==0) {
            fw=stdout;
        } else {
            fw=fopen(argv[2],"r");
            if (fw) {
                char s[100];
                printf("File %s already exists\n",argv[2]);
                printf("Do you wish to overwrite? (yes or no): " );
                scanf("%s",s);
                if(strcmp(s, "yes")) {
                    printf("Canceled!\n");
                    return -1;
                }
                fclose(fw);
            }
            fw=fopen(argv[2],"w");
        }
        if (!fw) {
            fprintf(stderr,"error: output file %s can't open\n", argv[2]);
            return -1;
        }

        fprintf(stderr,"#looking for a header ... ");
znova:
        do {
            int rd=a_fread(shortbuf,2,1024,f);
            if (rd<1024) {
                fclose(f);
                fprintf(stderr,"\r\nerror: no header high tone\n");
                return -1;
            }
            ad+=rd;
            for (j=0;j<1024;j++) inreal[j]=(double)shortbuf[j];
        } while (!istone_and_not(2000,1000));
        pokusy=timeout(0.3);
        do {
            int rd=a_fread(shortbuf,2,1024,f);
            if (rd<1024) {
                fclose(f);
                fprintf(stderr,"\r\nerror: no header low tone\n");
                return -1;
            }
            ad+=rd;
            if (pokusy--<0) goto znova;
            for (j=0;j<1024;j++) inreal[j]=(double)shortbuf[j];
        } while (!istone_and_not(1000,2000));
        seeky[256]=ad+T1;
        ad+=T0;//((int)(1.80*48000));
        a_fseek(f,ad*2,SEEK_SET);

        fprintf(stderr,"done.\n");fflush(stderr);

        fprintf(stderr,"#searching patterns ... ");
        for (i=0;i<256;i++) {
//            fprintf(stderr,"searching pattern %d\n",i);
            pokusy=timeout(0.3);
            do {
                int rd=a_fread(shortbuf,2,1024,f);
                if ((rd<1024)||((pokusy--<0))) {
                    fclose(f);
                    fprintf(stderr,"\r\nerror: no data %d/256\n",i);
                    return 0;
                }
                ad+=rd;
                for (j=0;j<1024;j++) inreal[j]=(double)shortbuf[j];
            } while (!istone_and_not(1000,2000));
            pokusy=timeout(0.3);
            do {
                int rd=a_fread(shortbuf,2,1024,f);
                if ((rd<1024)||((pokusy--<0))) {
                    fclose(f);
                    fprintf(stderr,"\r\nerror: no data %d/256\n",i);
                    return 0;
                }
                ad+=rd;
                for (j=0;j<1024;j++) inreal[j]=(double)shortbuf[j];
            } while (!istone_and_not(2000,1000));
            seeky[i]=ad+T1;
            ad+=T0;
            a_fseek(f,ad*2,SEEK_SET);
        }
        {
            pokusy=timeout(0.5);
            do {
                int rd=a_fread(shortbuf,2,1024,f);
                if ((rd<1024)||((pokusy--<0))) {
                    fclose(f);
                    fprintf(stderr,"\r\nerror: no end tone\n");
                    return 0;
                }
                ad+=rd;
                for (j=0;j<1024;j++) inreal[j]=(double)shortbuf[j];
            } while (!istone_and_not(3000,1000));
        }
        fprintf(stderr,"done.\n");fflush(stderr);

        // ted nactem referenci

        ad=seeky[256];
        a_fseek(f,ad*2,SEEK_SET);
        ad+=a_fread(shortbuf,2,BLOCK,f);
        for (j=0;j<BLOCK;j++) inreal[j]=(double)shortbuf[j];
        FFT(inreal,outreal1,WINDOW,NWIN);
        ad=seeky[256]+T4;
        a_fseek(f,ad*2,SEEK_SET);
        ad+=a_fread(shortbuf,2,WINDOW,f);
        for (j=0;j<WINDOW;j++) inreal[j]=(double)shortbuf[j];
        FFT(inreal,outreal2,WINDOW,NWIN);
        get_ref_kvadrat(outreal1,outreal2,xref,WINDOW);

        // pokracujem 256ti vzorkama

        fprintf(stderr,"#reading patterns ... ");
        memset(outref,0,sizeof(outref));
        for (i=0;i<256;i++) {
            ad=seeky[i];
            a_fseek(f,ad*2,SEEK_SET);
            ad+=a_fread(shortbuf,2,BLOCK,f);
            for (j=0;j<BLOCK;j++) inreal[j]=(double)shortbuf[j];
            FFT(inreal,outreal1,WINDOW,NWIN);
            memcpy(&outreal[0][i][0],&outreal1[0],sizeof(double)*WINDOW/2);
            ad=seeky[i]+T4;
            a_fseek(f,ad*2,SEEK_SET);
            ad+=a_fread(shortbuf,2,BLOCK,f);
            for (j=0;j<BLOCK;j++) inreal[j]=(double)shortbuf[j];
            FFT(inreal,outreal2,WINDOW,NWIN);
            memcpy(&outreal[1][i][0],&outreal2[0],sizeof(double)*WINDOW/2);
            for (j=0;j<WINDOW/2;j++) outref[j]+=outreal1[j]*outreal1[j]+outreal2[j]*outreal2[j];
        }
        fclose(f);

        for (j=0;j<WINDOW/2;j++) x0[j]=outref[j]=sqrt(outref[j]/256.0);
        Denoise2(x0,outref,xref,WINDOW);
        Filtr(x0,outreal1,WINDOW);
        for (j=WINDOW/2048;j<WINDOW*15/32;j++) if (outreal1[j]<=0.0) {
            fprintf(stderr,"\r\nerror: too noisy !!! (frequency %dHz)\r\n",j*48000/WINDOW);
            return -1;
        }
        fprintf(stderr,"done.\n");fflush(stderr);

        fprintf(stderr,"#computing patterns ...\n");fflush(stderr);
        for (i=0;i<256;i++) {
            Denoise2(outreal1,outref,xref,WINDOW);
            Denoise2(outreal2,outref,xref,WINDOW);
            Normalize(outreal[0][i],outreal[1][i],x0,WINDOW);
            Filtr(outreal[0][i],x1,WINDOW);
            Filtr(outreal[1][i],x2,WINDOW);
            j=NajdiBod(x1,x2,WINDOW);
            fprintf(fw,"%d %d\n",i*8,j);fflush(fw);
        }
        fprintf(stderr,"#done.\n");fflush(stderr);
        fclose(fw);
    }
    return 0;
}

