#pragma once

#include "sequence.h"
#include "config.h"
#include "matrix.h"
#include <math.h>

#define mahds_max(a, b) (((a) <= (b)) ? (b) : (a))
#define mahds_min(a, b) (((a) <= (b)) ? (a) : (b))

double calculate_R_sqr (size_t const size_of_symbols_set, size_t const period_length, double const C);
void create_frequency_matrix(
	mahds_al_set* const aligned_orig_set,
	mahds_al_set* const aligned_artificial_set,
	OUT mahds_matrix*    V
);
void create_weight_matrix(
	mahds_matrix* const v, 
	OUT mahds_matrix*   m
);
void normalize_weight_matrix(
	mahds_matrix* const            m, 
	mahds_seq_set* const           seq_set, 
	double const                  R, 
	double const                  K_d, 
	size_t const                  period_length, 
	OUT mahds_matrix*              m_prime
);
