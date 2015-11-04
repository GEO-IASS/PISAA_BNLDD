

int *ivector(long nl, long nh);
double *dvector(long nl, long nh);
double **dmatrix(long nrl, long nrh, long ncl, long nch);
int **imatrix(long nrl, long nrh, long ncl, long nch);
int ***i3tensor(long nrl, long nrh, long ncl, long nch, long ndl, long ndh);

void free_ivector(int *v, long nl, long nh);
void free_dmatrix(double **m, long nrl, long nrh, long ncl, long nch);
void free_dvector(double *v, long nl, long nh);
void free_imatrix(int **m, long nrl, long nrh, long ncl, long nch);
void free_i3tensor(int ***t, long nrl, long nrh, long ncl, long nch,
	long ndl, long ndh);

