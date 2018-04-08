#include <stdio.h>
#include <stdlib.h>
#include <string.h>

// peut mettre des lignes et colonnes soit differentes soit meme
#define LIGNES					4
#define COLONNES				4

double matrix_values[] =
{
	8.0f, -1.0f, 3.0f, -1.0f, -1.0f, 6.0f, 2.0f, 0.0f, 3.0f, 2.0f, 9.0f, 1.0f, -1.0f, 0.0f, 1.0f, 7.0f
};

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
	double **test_matrix = allocate_2d_array(LIGNES, COLONNES);
	
	fill_2d_array(test_matrix, matrix_values, LIGNES, COLONNES);
	leverrier(test_matrix, LIGNES, COLONNES);
	
	getchar();
	return 0;
}
