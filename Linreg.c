#include<stdio.h>
#include<stdlib.h>
double ** multiply (double **,double**,int,int,int,int);
double ** transpose (double**,int,int);
double ** inverse (double **,int,int);
double * rowop (double *, double*, double, double,int);
void swapmax(double***,int,int);

// swap the maximum element in a column during gaussian elimination
void swapmax(double*** k, int b,int arow) {
	if (b == arow) return;
	int i,index,max;
	double ** a = *k;
	max = a[b][b];
	index = b;
	for (i = b; i < arow; i++) {
		if (a[i][b] > max) {
			max = a[i][b];
			index = i;
		} 
	}
	if (index == b) return;
	double * temp = rowop(a[index],NULL,1,0,2*arow);
	double * temp1 = rowop(a[b],NULL,1,0,2*arow);
	a[index] = temp1;
	a[b] = temp;
	return;
}

// perform a row operation of the form ax + by
double * rowop (double * a, double * b, double coefa, double coefb, int len) {
	if (a == NULL && b == NULL) {return NULL;}
	
	double * result = (double*) malloc(len*sizeof(double));
	int i;
	if (b == NULL || coefb == 0) {
		for (i = 0; i < len; i++) {
		result[i] = coefa*a[i];
		}
		return result;
	}
	if (coefa == 0) {
		for (i = 0; i < len; i++) {
		result[i] = coefb*b[i];
		}	
	return result;
	}
	else {
		for(i = 0; i < len; i++) {
		result [i] = coefa*a[i] + coefb*b[i];
		}		
		return result;
	}
	return NULL;

}
double ** inverse (double** a, int arow, int acol) {
	if (arow != acol) {return NULL;}
	double ** inv = (double**) malloc(arow*sizeof(double*));
	int i,j;
	for (i = 0; i < arow; i++) {
		inv[i] = (double*) malloc(2*acol*sizeof(double));
	}
	for (i = 0; i < arow; i++) {
		for (j = 0; j < 2*acol; j++) {
			if (j < acol) {
			inv[i][j] = a[i][j];
			}
			else if (j - acol == i) {
			inv[i][j] = 1;
			} 
		}
	}
	double * temp;
	for (i = 0; i < arow; i++) {
		swapmax(&inv,i,arow);
		temp = rowop(inv[i],NULL,1/inv[i][i],0,2*acol);
		free(inv[i]);
		inv[i] = temp;
		for (j = i+1; j < acol; j++) {
		temp = rowop(inv[i],inv[j],-1*inv[j][i],1,2*acol);
		free(inv[j]);
		inv[j] = temp;
		}
		
	}
	for (i = arow-1; i >= 0; i--) {
		for (j = i-1; j >= 0; j--) {
			temp = rowop(inv[i],inv[j],-1*inv[j][i],1,2*acol);
			free(inv[j]);
			inv[j] = temp;
			
		}
	}
	double ** result = (double**) malloc(acol*sizeof(double));
	for (i = 0; i < acol; i++) {
	result[i] = (double*) malloc(acol*sizeof(double));
	}
	for (i = 0; i < acol; i++) {
		for ( j = 0; j < acol; j++) {
			result[i][j] = inv[i][acol+j];
		}
	}
	return result;
}
double ** transpose(double** a, int arow, int acol) {
	if (a == NULL) {
		return NULL;
	}		
	double ** result = (double**) malloc(acol*sizeof(double*));
	int i,j;
	for (i = 0; i < acol; i++) {
		result[i] = (double*) malloc (arow*sizeof(double));
	}
	for (i = 0; i < arow; i++) {
		for (j = 0; j < acol; j++) {
			result[j][i] = a[i][j];
		}
	}
	return result;
}
double ** multiply(double** aa, double** bb, int arow, int acol, int brow, int bcol) {
	
	if(aa == NULL || bb == NULL || acol != brow) {
		return NULL;
	}
	int i,j,k;
	double temp = 0;
	double ** result = (double**) malloc(arow*sizeof(double*));
	for (i = 0; i < arow; i++) {
		result[i] = (double*) malloc(bcol*sizeof(double));
	}
	for (i = 0; i < arow; i++) {
		for (j = 0; j < bcol; j++) {
			for (k = 0; k < acol; k++) {
			temp = temp + aa[i][k] * bb[k][j];	
			}
			result[i][j] = temp;
			temp = 0;
		}
	}
	return result;
}
	
int main (int argc, char** argv) {
	if (argc != 3) {return 0;}
	FILE * fp = fopen(argv[1],"r");
	if (fp == NULL) {return 0;}
	int i,j,arow,acol;
	fscanf(fp,"%d",&acol);
	fscanf(fp,"%d",&arow);
	double ** matrix = (double**) malloc(arow*sizeof(double*));
	double ** yval = (double**) malloc(arow*sizeof(double*));
	for (i = 0; i < arow; i++) {
		matrix[i] = (double *) malloc((acol+1)*sizeof(double));
		yval[i] = (double*) malloc(sizeof(double));
	}
	char dummy;
	for (i = 0; i < arow; i++) {
		for (j = 0; j < acol+1; j++) {
		if (j != acol) {
		fscanf(fp,"%lf",&matrix[i][j+1]);
		fscanf(fp,"%c",&dummy);
		}
		else {
		fscanf(fp,"%lf",&yval[i][0]);
		fscanf(fp,"%c",&dummy);
		}
		}
	}
	for (j = 0; j < arow; j++) {
		matrix[j][0] = 1;
	}
	/*for (i = 0; i < arow; i++) {
		for (j = 0; j < acol+1; j++) {
		printf("%f ",matrix[i][j]);
		}
		printf("\n");
	}*/
	
	double ** p = inverse(multiply(transpose(matrix,arow,acol+1),matrix,acol+1,arow,arow,acol+1),acol+1,acol+1);  
	double ** alk = multiply(p,transpose(matrix,arow,acol+1),acol+1,acol+1,acol+1,arow);
	double ** res = multiply(alk,yval,acol+1,arow,arow,1);
	fclose(fp);
	fp = fopen(argv[2],"r");
	if (fp == NULL) return 0;
	
	int tol;
	fscanf(fp,"%d",&tol);
	double ** test = (double**) malloc(tol*sizeof(double*));
	for (i = 0; i < tol; i++)
	test[i] = (double*) malloc((acol+1)*sizeof(double));
	for (i = 0; i < tol; i++) {
		test[i][0] = 1;	
		for (j = 0; j < acol; j++) {
		fscanf(fp,"%lf",&test[i][j+1]);
		fscanf(fp,"%c",&dummy);
		}
	}
	double ** predict = multiply(test,res,tol,acol+1,acol+1,1);
	for (i = 0; i < tol; i++) {
		for (j = 0; j < 1; j++) {
		printf("%d\n",(int)predict[i][j]);
		}
	}

	return 0;
}
