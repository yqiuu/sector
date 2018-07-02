FILE *open_file(char *fName, char *mode);

double **malloc_2d_double(int nRow, int nCol);
double **memcpy_2d_double(double **source, int nRow, int nCol);
void free_2d_double(double **p, int nRow);

int bisection_search(double a, double *x, int nX);
double interp(double xp, double *x, double *y, int nPts);
double trapz_table(double *y, double *x, int nPts, double a, double b);
double trapz_filter(double *filter, double *flux, double *waves, int nWaves);

struct linResult {
    double slope;
    double intercept;
    double R;
};
struct linResult linregress(double *x, double *y, int nPts);
