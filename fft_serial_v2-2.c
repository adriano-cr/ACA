#include <stdio.h>
#include <complex.h>
#include <math.h>
#include <omp.h>
#include <stdlib.h>
#include <malloc.h>

#define VEC_SIZE 2
 
int check(int n) 
{
    return n > 0 && (n & (n - 1)) == 0;
}

void printvec(double complex* vec, int n)
{
	FILE *fptr;
	fptr = fopen("./results.txt","w");
	int i;
	for (i = 0; i < n; i++){
        if(creal(vec[i])>=0){
        	fprintf(fptr,"+%.2f", creal(vec[i]));
		}else{
			fprintf(fptr,"-%.2f", cabs(creal(vec[i])));
		}
		
		if(cimag(vec[i])>=0){
        	fprintf(fptr," +%.2fi\n", cimag(vec[i]));
		}else{
			fprintf(fptr," -%.2fi\n", cabs(cimag(vec[i])));
		}
	}
}

int readfile(double complex** vec, char* filepath)
{
	FILE* filePointer;
	char line[50] = {0};
	filePointer = fopen(filepath,"r");
	double complex* v;
	int n=0, i=0, vec_size=VEC_SIZE;						//number of elements in the vector
	
	v = malloc(VEC_SIZE*sizeof(complex));
	//read file and store elements in vec
    float tmp;
    
    for(i=0; i<vec_size; i++){
    	fgets(line, 50, filePointer);
    	if(n==(vec_size-1) && !feof(filePointer)){					//check if size is exceeded
			vec_size=2*vec_size;				
			v = realloc(v, vec_size*sizeof(complex));			//double initial allocated space
		}
		sscanf(line,"%f",&tmp);
		v[n] = tmp + 0*I;
		n++;
	}
    
    fclose(filePointer);
    *vec = v;
	return n;
}

int logint2(int N)    //compute log2 of an integer
{
    int k = N, i = 0;
    
	while (k) {
        k >>= 1;
        i++;
    }
    
	return i - 1;
}

int reverse(int N, int n)    //computes the inverse number of any integer n with respect to the max number N
{
    int j, p = 0;
    
	for (j = 1; j <= logint2(N); j++) {
        if (n & (1 << (logint2(N) - j)))
            p |= 1 << (j - 1);
    }
    
	return p;
}

void order(double complex* f1, int N)     //reverse order of vector's elements
{
    double complex f2[N];
    
    int i;
    int j;
	for (i = 0; i < N; i++)
        f2[i] = f1[reverse(N, i)];
    
	for (j = 0; j < N; j++)
        f1[j] = f2[j];
}

void transform(double complex* f, int N)     //performs actual FFT
{
    order(f, N);    						//vector is reverse-ordered
    double complex W[N/2]; 
    
    //prepare values for W
	float rew = 1*cos(-2.*M_PI/N);
	float imw = 1*sin(-2.*M_PI/N);
	double complex w = rew + imw * I;

	//initialize first values of W
	W[0] = 1;
	W[1] = w;
		
	int i,j;
	for (i = 2; i < N / 2; i++)
        W[i] = cpow(W[1], i);
    
	int n = 1;
    int a = N / 2;
    
	for (j = 0; j < logint2(N); j++) {
        for (i = 0; i < N; i++) {
            if (!(i & n)) {
                double complex temp = f[i];
                double complex Temp = W[(i * a) % (n * a)] * f[i + n];
                f[i] = temp + Temp;
                f[i + n] = temp - Temp;
            }
        }
        
		n *= 2;
        a = a / 2;
    }
}

int main(int argc, char **argv)
{
    double start=0;
    double end=0;
    start = omp_get_wtime();						//variables and functions to compute exec time
	
    double complex *vec;							//vector to be transformed
        
    if (argc < 2){
        printf("ERROR: no path");
		return 1;
	}
    
    int n = readfile(&vec, argv[1]);
        
    //check if n is a power of 2 before proceding
    if (!check(n)){
		printf("ERROR: power of 2");
		return 1;
	}
       
    //call function to perform the FFT
	transform(vec, n);
    
    //print result to file
	printvec(vec, n);
	
	end = omp_get_wtime();
	printf("\nComputed in %f s", end-start);
	
    return 0;
}