#include "matrix.h"
#include <stdio.h>
#include <math.h>


double _mahds_matrix_get(
	const mahds_matrix* matrix, 
	size_t             line, 
	size_t             column)
{
	return matrix->data[matrix->columns*line+column];
}

void _mahds_matrix_set(
	mahds_matrix* matrix, 
	size_t       line, 
	size_t       column, 
	double       value)
{
	matrix->data[matrix->columns*line + column] = value;
}

void mahds_matrix_set_all(mahds_matrix* matrix, double value) {
	for (size_t i = 0; i < (matrix->columns * matrix->lines); ++i) {
		matrix->data[i] = value;
	}
}

void mahds_matrix_free(mahds_matrix* matrix) {
	free(matrix->data);
	free(matrix);
}

void mahds_matrix_max(mahds_matrix* matrix, OUT double* max, OUT size_t* max_line, OUT size_t* max_column) {
	double lmax = mahds_matrix_get(matrix, 0, 0);
	size_t lmax_line = 0;
	size_t lmax_column = 0;
	for (size_t i = 0; i < matrix->lines; ++i) {
		for (size_t j = 0; j < matrix->columns; ++j) {
			double temp = mahds_matrix_get(matrix, i, j);
			if (temp > lmax) {
				lmax = temp;
				lmax_line = i;
				lmax_column = j;
			}
		}
	}
	if (max != NULL) {
		*max = lmax;
	}
	if (max_line != NULL) {
		*max_line = lmax_line;
	}
	if (max_column != NULL) {
		*max_column = lmax_column;
	}
}

mahds_matrix* mahds_matrix_init(size_t lines, size_t columns) {
	mahds_matrix* result = (mahds_matrix*)malloc(sizeof(mahds_matrix));
	result->lines = lines;
	result->columns = columns;
	result->data = (double*)calloc(sizeof(double*), lines*columns);
	return result;
}

void mahds_matrix_print(mahds_matrix* matrix) {
	for (size_t i = 0; i < matrix->lines; ++i) {
		for (size_t j = 0; j < matrix->columns; ++j) {
			printf("%f ", mahds_matrix_get(matrix, i, j));
		}
		printf("\n");
	}
}

void mahds_matrix_print_part(mahds_matrix* matrix, size_t line, size_t column) {
	for (size_t i = 0; i < line; ++i) {
		printf("[");
		for (size_t j = 0; j < column; ++j) {
			printf("%f, ", mahds_matrix_get(matrix, i, j));
		}
		printf("]\n");
	}
}

void mahds_matrix_distance(
	const mahds_matrix* m1,
	const mahds_matrix* m2,
	OUT double*        d) {
	*d = 0;
	for (size_t i = 0; i < m1->lines; ++i) {
		for (size_t j = 0; j < m1->columns; ++j) {
			double temp = mahds_matrix_get(m1, i, j) - mahds_matrix_get(m2, i, j);
			temp *= temp;
			*d += temp;
		}
	}
	*d = sqrt(*d);
}
