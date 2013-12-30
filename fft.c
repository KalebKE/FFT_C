/* 
 * File:   FFT.c
 * Author: Kaleb
 *
 * Created on December 28, 2013, 5:41 PM
 */

#include <stdio.h>
#include <stdlib.h>
#include <math.h>

#define M_PI 3.14159265358979323846264338327

void calcImag(double radians_signal[], double signal[], int k);
void calcReal(double radians_signal[], double signal[], int k);
double sum(double c[]);


double signal[] = {2.28025, 1.32888, 0.39326, -0.49619, -1.31121, -2.02672, -2.62174,
    -3.08015, -3.39124, -3.55077, -3.55763, -3.42069, -3.15151,
    -2.76733, -2.28963, -1.74326, -1.15541, -0.55456, 0.03068,
    0.57271, 1.04606, 1.42835, 1.7122, 1.85105, 1.86948, 1.75376,
    1.50688, 1.13742, 0.65924, 0.09094, -0.54489, -1.22254,
    -1.91419, -2.59102, -3.22433, -3.78672, -4.25312, -4.60188,
    -4.81556, -4.8817, -4.79336, -4.54939, -4.15456, -3.61943,
    -2.95996, -2.19702, -1.35554, -0.46368, 0.44824, 1.34888,
    2.20699, 2.99264, 3.6783, 4.23987, 4.65761, 4.91685, 5.00855,
    4.92967, 4.68323, 4.27822, 3.72929, 3.05615, 2.28286, 1.43692,
    0.54824, -0.35197, -1.23239, -2.06272, -2.81485, -3.46385,
    -3.9889, -4.37404, -4.60875, -4.68828, -4.61378, -4.39223,
    -4.03611, -3.56286, -2.99417, -2.3551, -1.67308, -0.97682,
    -0.29515, 0.3441, 0.91522, 1.3956, 1.76664, 2.01448, 2.13052,
    2.1118, 1.96104, 1.68661, 1.30214, 0.826, 0.28055, -0.3087,
    -0.91416, -1.50721, -2.05935, -2.54337, -2.93442, -3.21099,
    -3.35583, -3.35666, -3.20667, -2.90487, -2.45619, -1.87128,
    -1.16626, -0.36208, 0.51622, 1.44039, 2.38004, 3.30375,
    4.18023, 4.97949, 5.67398, 6.23957, 6.65652, 6.91015, 6.99145,
    6.89738, 6.631, 6.20138, 5.62319, 4.91623, 4.10463, 3.216};

// size of table
int iN = 128;
// harmonics
int iK = 4;

double imc[4][128];
double rec[4][128];

/*
 * 
 */
int main()
{
    double amplitude[4];
    double phase[4];
    double degreesSignal[128];
    double imcSum[4];
    double radiansSignal[128];
    double recSum[4];

    int i = 0;

    for (i = 0; i < iN; i++)
    {
        degreesSignal[i] = (360.0 / (double) iN) * (double) i;

        //printf("Degree Signal: %f \n", degreesSignal[i]);

        radiansSignal[i] = degreesSignal[i] * (M_PI / 180.0);

        //printf("Radians Signal: %f \n", radiansSignal[i]);
    }

    // get the real coefficients
    for (i = 0; i < iK; i++)
    {
        calcReal(radiansSignal, signal, i + 1);
    }

    // sum the real coefficients
    for (i = 0; i < iK; i++)
    {
        recSum[i] = sum(rec[i+1]);

        //printf("Real Sum: %f \n", recSum[i]);
    }

    // get the imag coefficients
    for (i = 0; i < iK; i++)
    {
        calcImag(radiansSignal, signal, i + 1);
    }

    // sum the imag coefficients
    for (i = 0; i < iK; i++)
    {
        imcSum[i] = sum(imc[i + 1]);

        //printf("Imaginary Sum: %f \n", recSum[i]);
    }

    // calculate amplitude (Fourier coefficients)
    // a[i] = 2(abs(RECSum[i]*abs(IMCSum[i])
    for (i = 0; i < iK; i++)
    {
        amplitude[i] = 2.0 * (sqrt(pow(recSum[i], 2.0) + pow(imcSum[i], 2.0)));

        //printf("Amplitude: %f \n", amplitude[i]);
    }

    // calculate phase
    // atan(abs(b/a))
    for (i = 0; i < iK; i++)
    {
        double temp = (sqrt(pow(imcSum[i], 2)
                / pow(recSum[i], 2)));
        phase[i] = atan(temp);
        
        //printf("Phase: %f \n", phase[i]);
    }

    for (i = 0; i < iK; i++)
    {
        printf("REC %d: %f\n", (i + 1), recSum[i]);
        printf("IMC  %d: %f\n", (i + 1), imcSum[i]);
        printf("Amplitude(K=%d): %f\n", (i + 1), amplitude[i]);
        printf("Phase(K=%d): %f\n", (i + 1), phase[i]);
        printf("Degrees(K=%d): %f\n", (i + 1), phase[i]*(180.0/M_PI));
        printf("\n");
    }

    return (EXIT_SUCCESS);
}

void calcImag(double radians_signal[], double signal[], int k)
{
    imc[k][0] = signal[0] * sin(-radians_signal[0]);

    for (int i = 1; i < iN; i++)
    {
        // for when the modular value is 0
        // prevents array out of bounds
        // IMC[array.length - 0]...
        // use the first sin value instead
        if ((i * k) % iN == 0)
        {
            imc[k][i] = signal[i] * sin(radians_signal[0]);
        }
        else
        {
            imc[k][i] = (signal[i] * sin(radians_signal[iN - (i * k) % iN]));
        }

        //printf("Imaginary Signal: %f \n", imc[k][i]);
    }
}

void calcReal(double radians_signal[], double signal[], int k)
{
    rec[k][0] = signal[0] * cos(radians_signal[0]);

    for (int i = 1; i < iN; i++)
    {
        // for when the modular value is 0
        // prevents array out of bounds
        // REC[array.length - 0]
        // use the first sin value instead
        if (((i * k) % iN) == 0)
        {
            rec[k][i] = signal[i] * cos(radians_signal[0]);
        }
        else
        {
            rec[k][i] = +(signal[i] * cos(radians_signal[iN - (i * k) % iN]));
        }

        //printf("REC: %d (n=%d) %f \n", k, i, rec[k][i]);
    }
}

double sum(double c[])
{
    double dCSum = 0;

    for (int i = 0; i < iN; i++)
    {
        dCSum += c[i];
    }

    dCSum = dCSum / iN;

    return dCSum;
}

