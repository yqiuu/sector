#include<stdio.h>
#include<stdlib.h>
#include<string.h>
#include<math.h>


/* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * *
 *                                                                             *
 * Basic functions                                                             *
 *                                                                             *
 * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * */
FILE *open_file(char *fName, char *mode) {
/* Open the file with specific mode */
    FILE *fp;
    if ((fp = fopen(fName, mode)) == NULL) {
        printf("File open error: \"%s\"!\n", fName);
        exit(0);
    }
    printf("# File opened: \"%s\"!\n", fName);
    return fp;
}


double **malloc_2d_double(int nRow, int nCol) {
    int i;
    double **target = malloc(nRow*sizeof(double*));
    for(i = 0; i < nRow; ++i)
        target[i] = (double*)malloc(nCol*sizeof(double));
    return target;
}


double **memcpy_2d_double(double **source, int nRow, int nCol) {
    int i;
    double **target = malloc_2d_double(nRow, nCol);
    for(i = 0; i < nRow; ++i)
        memcpy(target[i], source[i], nCol*sizeof(double));
    return target;
}


void free_2d_double(double **p, int nRow) {
    int i;
    for(i = 0; i < nRow; ++i)
        free(p[i]);
    free(p);
}


int bisection_search(double a, double *x, int nX) {
    /* return idx such x[idx] <= a < x[idx + 1]
     * a must be x[0] <= a < x[nX - 1]
     */
    unsigned int idx0 = 0;
    unsigned int idx1 = nX - 1;
    unsigned int idxMid;
    while(idx1 - idx0 > 1) {
        idxMid = (idx0 + idx1)/2;
        if(a >= x[idxMid])
            idx0 = idxMid;
        else if(a < x[idxMid])
            idx1 = idxMid;
    }
    return idx0;
}


double interp(double xp, double *x, double *y, int nPts) {
    /* Interpolate a given points */
    int idx0, idx1;
    if((xp < x[0]) || (xp > x[nPts - 1])) {
        printf("Error: Point %10.5e is beyond the interpolation region\n", xp);
        exit(0);
    }
    if (xp == x[nPts - 1])
        return y[nPts - 1];
    else {
        idx0 = bisection_search(xp, x, nPts);
        if (x[idx0] == xp)
            return y[idx0];
        idx1 = idx0 + 1;
        return y[idx0] + (y[idx1] - y[idx0])*(xp - x[idx0])/(x[idx1] - x[idx0]);
    }
}


double trapz_table(double *y, double *x, int nPts, double a, double b) {
    /* Integrate tabular data from a to b */
    int i;
    int idx0, idx1;
    double ya, yb;
    double I;
    if (x[0] > a) {
        printf("Error: Integration range %10.5e is beyond the tabular data\n", a);
        exit(0);
    }
    if (x[nPts - 1] < b) {
        printf("Error: Integration range %10.5e is beyond the tabular data\n", b);
        exit(0);
    }
    if (a > b) {
        printf("Error: a must be smaller than b\n");
        exit(0);
    }
    idx0 = bisection_search(a, x, nPts);
    idx1 = idx0 + 1;

    ya = y[idx0] + (y[idx1] - y[idx0])*(a - x[idx0])/(x[idx1] - x[idx0]);
    if(b <= x[idx1]) {
        yb = y[idx0] + (y[idx1] - y[idx0])*(b - x[idx0])/(x[idx1] - x[idx0]);
        return (b - a)*(yb + ya)/2.;
    }
    else
        I = (x[idx1] - a)*(y[idx1] + ya)/2.;

    for(i = idx1; i < nPts - 1; ++i) {
        if (x[i + 1] < b)
            I += (x[i + 1] - x[i])*(y[i + 1] + y[i])/2.;
        else if (x[i] < b) {
            yb = y[i] + (y[i + 1] - y[i])*(b - x[i])/(x[i + 1] - x[i]);
            I += (b - x[i])*(yb + y[i])/2.;
        }
        else
            break;
    }
    return I;
}


double trapz_filter(double *filter, double *filterWaves, int nFilterWaves, 
                    double *flux, double *waves, int nWaves) {
    /* integrate the flux in a filter */
    int iFW;
    double f1 = interp(filterWaves[0], waves, flux, nWaves)*filter[0];
    double f2;
    double I = 0.;
    for(iFW = 1; iFW < nFilterWaves; ++iFW) {
        f2 = interp(filterWaves[iFW], waves, flux, nWaves)*filter[iFW];
        I += (filterWaves[iFW] - filterWaves[iFW - 1])*(f1 + f2);
        f1 = f2;
    }
    return I/2.;
}


struct linResult {
    double slope;
    double intercept;
    double R;
};


struct linResult linregress(double *x, double *y, int nPts) {
    int i;

    double xSum = 0.;
    for(i = 0; i < nPts; ++i)
        xSum += x[i];

    double ySum = 0.;
    for(i = 0; i < nPts; ++i)
        ySum += y[i];

    double xxSum = 0.;
    for(i = 0; i < nPts; ++i)
        xxSum += x[i]*x[i];

    double xySum = 0.;
    for(i = 0; i < nPts; ++i)
        xySum += x[i]*y[i];

    double denominator = nPts*xxSum - xSum*xSum;

    double slope = (nPts*xySum - xSum*ySum)/denominator;
    double intercept = (xxSum*ySum - xSum*xySum)/denominator;

    double yReg;
    double delta;
    double ssRes = 0.;
    for(i = 0; i < nPts; ++i) {
        yReg = slope*x[i] + intercept;
        delta = yReg - y[i];
        ssRes += delta*delta;
    }

    double yMean = ySum/nPts;
    double ssTot = 0.;
    for(i = 0; i < nPts; ++i) {
        delta = yMean - y[i];
        ssTot += delta*delta;
    }

    double R = sqrt(1. - ssRes/ssTot);

    struct linResult result;
    result.slope = slope;
    result.intercept = intercept;
    result.R = R;
    return result;
}
