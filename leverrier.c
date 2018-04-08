#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include <time.h>

#define _MIN(a, b) (a<=b ? a : b)

// peut mettre des lignes et colonnes soit differentes soit meme
#define LIGNES					7
#define COLONNES				7

double matrix_values[] =
{
	8.0f, -1.0f, 3.0f, -1.0f, -1.0f, 6.0f, 2.0f, 0.0f, 3.0f, 2.0f, 9.0f, 1.0f, -1.0f, 0.0f, 1.0f, 7.0f
};

// TEST MATRICES
void matrice_bord(double **matrix, int lig);
void matrice_ding_dong(double **matrix, int lig);
void matrice_franc(double **matrix, int lig);
void matrice_hilbert(double **matrix, int lig);
void matrice_kms(double **matrix, int lig, double p);
void matrice_lehmer(double **matrix, int lig);
void matrice_lotkin(double **matrix, int lig);
void matrice_moler(double **matrix, int lig);
double calculate_time_in_seconds(clock_t start, clock_t end);
void print_time_elapsed(const char *matrix_name, double elapsed);

double *allocate_1d_array(int lignes)
{
	double *ret = (double*)calloc(sizeof(double), lignes);
	if( ret == NULL )
	{
		return NULL;
	}
	
	return ret;
}

double **allocate_2d_array(int lignes, int colonnes)
{
	// allocation de memoire pour les lignes
	double **a_retourner=(double**)malloc(sizeof(double*)*lignes);
	int i=0;

	for( i=0; i<lignes; ++i )
	{
		// allocation de memoire pour des colonnes
		a_retourner[i]=(double*)calloc(colonnes, sizeof(double));
	}

	return a_retourner;
}

void affiche_2d_array(double **matrix, int lig, int col)
{
	int i=0, j=0;

	for( i=0; i<lig; ++i )
	{
		printf("[");
		for( j=0; j<col; ++j )
		{
			printf(" %f ", matrix[i][j]);
		}
		printf("]\n");
	}
	printf("\n");

	return;
}

void fill_2d_array(double **matrix, double *src_array, int lig, int col)
{
	int i=0, j=0, k=0;

	for( i=0; i<lig; ++i )
	{
		for( j=0; j<col; ++j )
		{
			matrix[i][j] = src_array[k];
			k++;
		}
	}

	return;
}

void fill_1d_array(double *dest_array, double *src_array, int lignes)
{
	int i = 0;
	
	if( src_array == NULL )
	{
		for( i = 0; i < lignes; ++i )
		{
			dest_array[i] = 1.0f;
		}
	}
	else
	{
		for( i = 0; i < lignes; ++i )
		{
			dest_array[i] = src_array[i];
		}
	}
}

void free_2d_array(double ***matrix, int lig)
{
	int i=0;
	double **mat=*matrix;

	for( i=0; i<lig; ++i )
	{
		if( mat[i] )
		{
			free(mat[i]);
		}
	}

	if( mat )
	{
		free(mat);
		*matrix=NULL;
	}
}

void free_1d_array(double **matrix)
{
	double *mat=*matrix;

	if( mat )
	{
		free(mat);
		*matrix=NULL;
	}
}

void affiche_1d_array(double *matrix, int lig)
{
	int i=0;

	for( i=0; i<lig; ++i )
	{
		printf("[ p%d = %lf ]\n", i, matrix[i]);
	}
	printf("\n");
}

double **creer_matrice_identite(int lignes)
{
	int i = 0;
	double **mat = allocate_2d_array(lignes, lignes);
	for( i = 0; i < lignes; ++i )
	{
		mat[i][i] = 1.0f;
	}
	
	return mat;
}

// on supposera que les matrices peuvent etre multiplie sans
// verifier les lignes et colonnes regle
double **multiply_matrix(double **mat1, double **mat2, int lig, int col)
{
	int i = 0, j = 0, k = 0;
	double sum = 0.0f;
	double **vec_next = allocate_2d_array(lig, col);
	
	for( i = 0; i < lig; ++i )
	{
		for( j = 0; j < col; ++j )
		{
			sum = 0.0f;
			for( k = 0; k < col; ++k )
			{
				sum = sum + (mat1[i][k] * mat2[k][j]);
			}
			vec_next[i][j] = sum;
		}
	}
	
	return vec_next;
}

void multiply_matrix_by_scalar(double **matrix, double scalar, int lignes, int colonnes)
{
	int i = 0, j = 0;
	
	for( i = 0; i < lignes; ++i )
	{
		for( j = 0; j < colonnes; ++j )
		{
			matrix[i][j] = scalar * matrix[i][j];
		}
	}
}

// stocke result in mat2
double **subtract_matrix(double **mat1, double **mat2, int lignes, int colonnes)
{
	int i = 0, j = 0;
	double **ret_val = allocate_2d_array(lignes, colonnes);
	
	for( i = 0; i < lignes; ++i )
	{
		for( j = 0; j < colonnes; ++j )
		{
			ret_val[i][j] = mat1[i][j] - mat2[i][j];
		}
	}
	
	return ret_val;
}

double trace_matrix(double **mat, int lignes, int divisor)
{
	int i = 0;
	double res = 0.0f;
	
	for( i = 0; i < lignes; ++i )
	{
		res += mat[i][i];
	}
	
	return ( res / divisor );
}

void copy_2d_array(double **dest, double **src, int lignes, int colonnes)
{
	int i = 0, j = 0;
	
	for( i = 0; i < lignes; ++i )
	{
		for( j = 0; j < colonnes; ++j )
		{
			dest[i][j] = src[i][j];
		}
	}
	
	return;
}

void leverrier(double **A, int lignes, int colonnes)
{
	int i = 0, j = 0, index = 0;			// index for end of coefficients array
	double **B = allocate_2d_array(lignes, colonnes);
	double **idn = NULL;
	double *p_values = allocate_1d_array(lignes+1);
	
	// B1 = A
	copy_2d_array(B, A, lignes, colonnes);
	
	for( i = 0; i < lignes; ++i )
	{
		// creer matrice identite reinitialise
		idn = creer_matrice_identite(lignes);
		
		// fill coefficients beginning at the end of the array
		// p0 corresponds to p^0 where p is the coefficient
		index = lignes - i - 1;
		p_values[index] = trace_matrix(B, lignes, (i+1));
		
		// p(n-1) * idn
		// resultat stockee en idn donc il faut reinitialiser idn au debut du boucle
		multiply_matrix_by_scalar(idn, p_values[index], lignes, lignes);
		//printf("p * idn\n\n");
		//affiche_2d_array(idn, lignes, colonnes);
		
		// B(n-1) - p(n-1) * idn
		B = subtract_matrix(B, idn, lignes, colonnes);
		
		// A * B(n-1) - p(n-1) * idn
		B = multiply_matrix(A, B, lignes, colonnes);
		//printf("B %d\n\n", (i+1));
		//affiche_2d_array(B, lignes, colonnes);
		
		p_values[index] = p_values[index] * (-1.0f);
	}
	
	// last coefficient is always 1 according to the formula
	p_values[lignes] = 1.0f;
	
	affiche_1d_array(p_values, lignes+1);
	
	return;
}

int main()
{
	printf("#### les numero au devant de p correspond a la puissance dans le polynomial ####\np0 = p0 * x^0, p1 = p1 * x^1\n");
	double **matrix = NULL;
	
	double time_elapsed=0.0f, error_calculation=0.0f;
	clock_t start, end;
	
	/*
	// SI DECOMMENTE, REMPLACE LIGNES ET COLONNES PAR 4 POUR matrix_values
	matrix = allocate_2d_array(LIGNES, COLONNES);
	fill_2d_array(matrix, matrix_values, LIGNES, COLONNES);
	leverrier(matrix, LIGNES, COLONNES);
	free_2d_array(&matrix, LIGNES);
	*/
	
	start = clock();
	printf("MATRICE BORD\n\n");
	matrix = allocate_2d_array(LIGNES, COLONNES);
	matrice_bord(matrix, LIGNES);
	leverrier(matrix, LIGNES, COLONNES);
	free_2d_array(&matrix, LIGNES);
	end = clock();
	time_elapsed = calculate_time_in_seconds(start, end);
	print_time_elapsed("BORD time", time_elapsed);
	
	start = clock();
	printf("MATRICE DING DONG\n\n");
	matrix = allocate_2d_array(LIGNES, COLONNES);
	matrice_ding_dong(matrix, LIGNES);
	leverrier(matrix, LIGNES, COLONNES);
	free_2d_array(&matrix, LIGNES);
	end = clock();
	time_elapsed = calculate_time_in_seconds(start, end);
	print_time_elapsed("DING DONG time", time_elapsed);
	
	start = clock();
	printf("MATRICE FRANC\n\n");
	matrix = allocate_2d_array(LIGNES, COLONNES);
	matrice_franc(matrix, LIGNES);
	leverrier(matrix, LIGNES, COLONNES);
	free_2d_array(&matrix, LIGNES);
	end = clock();
	time_elapsed = calculate_time_in_seconds(start, end);
	print_time_elapsed("FRANC time", time_elapsed);
	
	start = clock();
	printf("MATRICE HILBERT\n\n");
	matrix = allocate_2d_array(LIGNES, COLONNES);
	matrice_hilbert(matrix, LIGNES);
	leverrier(matrix, LIGNES, COLONNES);
	free_2d_array(&matrix, LIGNES);
	end = clock();
	time_elapsed = calculate_time_in_seconds(start, end);
	print_time_elapsed("HILBERT time", time_elapsed);
	
	/*
	// PROBLEME ICI QLQPART
	start = clock();
	printf("MATRICE KMS\n\n");
	matrix = allocate_2d_array(LIGNES, COLONNES);
	matrice_kms(matrix, LIGNES, 0.005f);
	leverrier(matrix, LIGNES, COLONNES);
	free_2d_array(&matrix, LIGNES);
	end = clock();
	time_elapsed = calculate_time_in_seconds(start, end);
	print_time_elapsed("KMS time", time_elapsed);
	*/
	
	start = clock();
	printf("MATRICE LEHMER\n\n");
	matrix = allocate_2d_array(LIGNES, COLONNES);
	matrice_lehmer(matrix, LIGNES);
	leverrier(matrix, LIGNES, COLONNES);
	free_2d_array(&matrix, LIGNES);
	end = clock();
	time_elapsed = calculate_time_in_seconds(start, end);
	print_time_elapsed("LEHMER time", time_elapsed);
	
	start = clock();
	printf("MATRICE LOTKIN\n\n");
	matrix = allocate_2d_array(LIGNES, COLONNES);
	matrice_lotkin(matrix, LIGNES);
	leverrier(matrix, LIGNES, COLONNES);
	free_2d_array(&matrix, LIGNES);
	end = clock();
	time_elapsed = calculate_time_in_seconds(start, end);
	print_time_elapsed("LOTKIN time", time_elapsed);
	
	start = clock();
	printf("MATRICE MOLER\n\n");
	matrix = allocate_2d_array(LIGNES, COLONNES);
	matrice_moler(matrix, LIGNES);
	leverrier(matrix, LIGNES, COLONNES);
	free_2d_array(&matrix, LIGNES);
	end = clock();
	time_elapsed = calculate_time_in_seconds(start, end);
	print_time_elapsed("MOLER time", time_elapsed);
	
	return 0;
}

/*MATRICE TEST*/
void matrice_bord(double **matrix, int lig)
{
	int i=0;

	// je recopie ce code a partir de ton mail
	for( i=0; i<lig; ++i )
	{
		matrix[i][i]=1;
	}

	for( i=0; i<lig; ++i )
	{
		matrix[0][i]=pow(2.0f, (1-(i+1)));
		matrix[i][0]=matrix[0][i];
	}

	return;
}

void matrice_ding_dong(double **matrix, int lig)
{
	int i=0, j=0;

	for( i=0; i<lig; ++i )
	{
		for( j=0; j<lig; ++j )
		{
			matrix[i][j]=1/(2.0f*lig-(i+1)-(j+1)+1.5f);
		}
	}

	return;
}

void matrice_franc(double **matrix, int lig)
{
	int i=0, j=0;

	for( i=0; i<lig; ++i )
	{
		for( j=0; j<lig; ++j )
		{
			if( (i+1)>=((j+1)+2) )
			{
				matrix[i][j]=0;
			}
			else
			{
				matrix[i][j]=_MIN(i+1, j+1);
			}
		}
	}

	return;
}


void matrice_hilbert(double **matrix, int lig)
{
	int i=0, j=0;

	for( i=0; i<lig; ++i )
	{
		for( j=0; j<lig; ++j )
		{
			matrix[i][j]=(1.0f/((i+1)+(j+1)-1));
		}
	}

	return;
}

void matrice_kms(double **matrix, int lig, double p)
{
	if( p<= 0 || p>= 1 )
	{
		printf("error in matrice_kms\n");
	}
	else
	{
		int i=0, j=0;

		for( i=0; i<lig; ++i )
		{
			for( j=0; j<lig+1; ++j )
			{
				matrix[i][j]=pow(p, fabs((double)(i+1)-(j+1)));
			}
		}
	}

	return;
}

void matrice_lehmer(double **matrix, int lig)
{
	int i=0, j=0;

	for( i=0; i<lig; ++i)
	{
		for( j=0; j<lig; ++j)
		{
			if (i<=j)
			{
				matrix[i][j]=((double)(i+1)/(j+1));
			}
			else
			{
				matrix[i][j]=((double)(j+1)/(i+1));
			}
		}
	}

	return;
}

void matrice_lotkin(double **matrix, int lig)
{
	int i=0, j=0;

	for( i=0; i<lig; ++i)
	{
		for( j=0; j<lig; ++j)
		{
			if( i==1 )
			{
				matrix[i][j]=1;
			}
			else
			{
				matrix[i][j]=(1.0f/((i+1)+(j+1)-1));
			}
		}
	}

	return;
}

void matrice_moler(double **matrix, int lig)
{
	int i=0, j=0;

	for( i=0; i<lig; ++i)
	{
		for( j=0; j<lig; ++j )
		{
			if( i==j )
			{
				matrix[i][j]=((double)i+1);
			}
			else
			{
				matrix[i][j]=(_MIN(i+1 ,j+1)-2.0f);
			}
		}
	}

	return;
}

double calculate_time_in_seconds(clock_t start, clock_t end)
{
	return ((double)(end-start)/CLOCKS_PER_SEC);
}

void print_time_elapsed(const char *matrix_name, double elapsed)
{
	printf("matrix %s time elapsed : %f seconds\n\n", matrix_name, elapsed);

	return;
}

