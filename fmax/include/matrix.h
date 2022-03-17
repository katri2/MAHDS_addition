#pragma once

#include <stdlib.h>

#define OUT

typedef struct {
	size_t  lines;
	size_t  columns;
	double* data;
} mahds_matrix;


mahds_matrix* mahds_matrix_init(
	size_t lines, 
	size_t columns);


#define mahds_matrix_get(matrix, line, column) ((matrix)->data[(matrix)->columns*(line)+(column)])

#define mahds_matrix_get_fast(matrix, position) ((matrix)->data[(position)])

#define mahds_matrix_set(matrix, line, column, value) ((matrix)->data[(matrix)->columns*(line)+(column)] = (value))

#define mahds_matrix_set_fast(matrix, position, value) ((matrix)->data[position] = (value))

void mahds_matrix_set_all(
	mahds_matrix* matrix, 
	double       value);

void mahds_matrix_free(
	mahds_matrix* matrix);

void mahds_matrix_max(
	mahds_matrix* matrix, 
	OUT double*  max, 
	OUT size_t*  max_line, 
	OUT size_t*  max_column);

void mahds_matrix_print(
	mahds_matrix* matrix);

void mahds_matrix_print_part(
	mahds_matrix* matrix, 
	size_t       line, 
	size_t       column);

void mahds_matrix_distance(
	const mahds_matrix* m1,
	const mahds_matrix* m2,
	OUT double*       d);
