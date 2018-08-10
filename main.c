#include <HsFFI.h>
#ifdef __GLASGOW_HASKELL__
#include "CMatrix_stub.h"
extern void __stginit_CMatrix(void);
#endif
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
enum
{
	BUF_LEN = 30000,
	BASE = 100,
	EXP = 2,
	MODE_SGAUSS   = 0,
	MODE_SGAUSSLE = 1,
	MODE_DGAUSS   = 2,
	MODE_SSOR     = 3,
	MODE_INV      = 4,
	MODE_CNNUM    = 5,
};

typedef struct {
	double **mtx;
	int n; 
} SqMtxWL;


SqMtxWL rdMtx(void)
{
	char    buf[BUF_LEN];
	int     len  = BASE * sizeof(double);
	double *frow = malloc(len);
	int     i    = 0;
	int     res  = 0;
	int     cnt  = 0;
	fgets(buf, sizeof(buf), stdin);
	while(sscanf(buf + cnt, "%lf%n", frow + i, &res) > 0) {
		cnt += res;
		i++;
		if(i == len) {
			len *= EXP;
			frow = realloc(frow, len);
		}
	}
	double **mtx = malloc(sizeof(double *) * i);
	mtx[0] = frow;
	for(int j = 1; j < i; j++) {
		mtx[j] = malloc(i * sizeof(double));
		for(int k = 0; k < i; k++) {
			scanf("%lf", &mtx[j][k]);
		}
	}
	SqMtxWL out = {mtx, i};
	return out;
}

double *rdVcr(int n)
{
	double *out = malloc(n * sizeof(double));
	for(int i = 0; i < n; i++) {
		scanf("%lf", out + i);
	}
	return out;
}

void printVector(double *v, int n) 
{
	for(int i = 0; i < n; i++) {
		printf("%.6e ", v[i]);
	}
	printf("\n");
}

void printMatrix(SqMtxWL m) 
{
	for(int i = 0; i < m.n; i++) {
		printVector(m.mtx[i], m.n);
	}
}

int main(int argc, char *argv[])
{
	hs_init(&argc, &argv);
#ifdef __GLASGOW_HASKELL__
	hs_add_root(__stginit_CMatrix);
#endif
	const char usage[] = "USAGE: mtxprog <MODE>\nPossible MODE values:\n  0 - Solve by Gauss method\n  1 - Solve by Gauss with leading element method\n  2 - Compute determinant with Gauss method\n  3 - Solve by Successive Over-Relaxation\n  4 - Compute invert matrix\n  5 - Compute condition number\n";
	if(argc == 1) {
		printf(usage);
		return 0;
	}
	SqMtxWL ipt;
	double *rhs;
	double *res;
	switch(atoi(argv[1])) {
		case MODE_SGAUSS:
			ipt = rdMtx();
			rhs = rdVcr(ipt.n);
			res = solveGauss_hs(ipt.mtx, rhs, ipt.n);
			printf("Result: ");
			printVector(res, ipt.n);
			break;
		case MODE_SGAUSSLE:
			ipt = rdMtx();
			rhs = rdVcr(ipt.n);
			res = solveGaussLE_hs(ipt.mtx, rhs, ipt.n);
			printf("Result: ");
			printVector(res, ipt.n);
			break;
		case MODE_DGAUSS:
			ipt = rdMtx();
			printf("Result: %.10e\n", detGauss_hs(ipt.mtx, ipt.n));
			break;
		case MODE_SSOR:
			ipt = rdMtx();
			rhs = rdVcr(ipt.n);
			double omega, eps;
			scanf("%lf%lf", &omega, &eps);
			res = sccOvRl_hs(ipt.mtx, rhs, ipt.n, omega, eps);
			printf("Result: ");
			printVector(res, ipt.n);
			break;
		case MODE_INV:
			ipt = rdMtx();
			SqMtxWL temp = {inv_hs(ipt.mtx, ipt.n), ipt.n};
			printf("Result: \n");
			printMatrix(temp);
			break;
		case MODE_CNNUM:
			ipt = rdMtx();
			printf("Result: %.10e\n", condNumber_hs(ipt.mtx, ipt.n));
			break;
		default:
			printf(usage);
			return 0;
	}
	hs_exit();
	return 0;
}
