#include "utility.h"
#include <stdio.h>

double calculate_R_sqr (size_t const size_of_symbols_set, size_t const period_length, double const C){
    return C * period_length * size_of_symbols_set;
}

void create_frequency_matrix(
	mahds_al_set* const aligned_orig_set,
	mahds_al_set* const aligned_artificial_set,
	OUT mahds_matrix*    V
) {
	mahds_matrix_set_all(V, 0);

	for (size_t i = 0; i < aligned_artificial_set->seq_number; ++i) {
		for (size_t j = 0; j < aligned_artificial_set->sequences[i]->sequence.length; ++j) {
			mahds_symbol line = aligned_orig_set->sequences[i]->sequence.symbols[j];
			mahds_symbol column = aligned_artificial_set->sequences[i]->sequence.symbols[j];
			if (line >= V->lines || column >= V->columns) {
				continue;
			}
			mahds_matrix_set(V, line, column, 1 + mahds_matrix_get(V, line, column));
		}
	}
}

void create_weight_matrix(
	mahds_matrix* const v, 
	OUT mahds_matrix*   m
) {
	mahds_matrix_set_all(m, 0);
	double N = 0;
	double* x = calloc(sizeof(double), v->lines);
	double* y = calloc(sizeof(double), v->columns);
	for (size_t i = 0; i < v->lines; ++i) {
		for (size_t j = 0; j < v->columns; ++j) {
			double item = mahds_matrix_get(v, i, j);
			N += item;
			x[i] += item;
		}
	}
	for (size_t j = 0; j < v->columns; ++j) {
		for (size_t i = 0; i < v->lines; ++i) {
			double item = mahds_matrix_get(v, i, j);
			y[j] += item;
		}
	}
	for (size_t i = 0; i < v->lines; ++i) {
		for (size_t j = 0; j < v->columns; ++j) {
			double v_ij = mahds_matrix_get(v, i, j);
			double p_ij = x[i] * y[j] / (N*N) + 0.00001;
			mahds_matrix_set(m, i, j,
				(v_ij - N * p_ij)
				/
				sqrt(N*p_ij*(1 - p_ij))
			);
		}
	}
	free(x);
	free(y);
}

void normalize_weight_matrix(
	mahds_matrix* const            m, 
	mahds_seq_set* const           seq_set, 
	double const                  R, 
	double const                  K_d, 
	size_t const                  period_length, 
	OUT mahds_matrix*              m_prime
) {
	size_t* b_arr = calloc(sizeof(size_t), seq_set->sequences[0]->count_symbols);  // seq_set->sequences[k]->count_symbols should be equal for any k
	for (size_t i = 0; i < seq_set->seq_number; ++i) {
		for (size_t j = 0; j < seq_set->sequences[i]->length; ++j) {
			mahds_symbol sym = seq_set->sequences[i]->symbols[j];
			if (sym < seq_set->sequences[0]->count_symbols) {
				++b_arr[sym];
			}
		}
	}
	
	size_t set_length = 0;
	for (size_t i = 0; i < seq_set->seq_number; ++i) {
		set_length += seq_set->sequences[i]->length;
	}

	mahds_matrix* p = mahds_matrix_init(seq_set->sequences[0]->count_symbols, period_length);
	for (size_t i = 0; i < p->lines; ++i) {
		for (size_t j = 0; j < p->columns; ++j) {
			mahds_matrix_set(p, i, j,
				(1.0 / (double)period_length) * ((double)b_arr[i] / (double)set_length)
			);
		}
	}
	free(b_arr);

	double sum_p_2 = 0;
	double sum_p_m = 0;
	for (size_t i = 0; i < p->lines; ++i) {
		for (size_t j = 0; j < p->columns; ++j) {
			double p_ij = mahds_matrix_get(p, i, j);
			double m_ij = mahds_matrix_get(m, i, j);
			sum_p_m += p_ij * m_ij;
			sum_p_2 += p_ij * p_ij;
		}
	}
	
	double a = (K_d - sum_p_m) / sum_p_2;
	double b = K_d / sum_p_2;
	
	double sum_m_p_a_b = 0;
	for (size_t i = 0; i < p->lines; ++i) {
		for (size_t j = 0; j < p->columns; ++j) {
			double p_ij = mahds_matrix_get(p, i, j);
			double m_ij = mahds_matrix_get(m, i, j);
			double temp = m_ij + p_ij * (a - b);
			temp *= temp;
			sum_m_p_a_b += temp;
		}
	}
	
	double B = (2 * b*(sum_p_m + (a - b)*sum_p_2)) / sum_m_p_a_b;
	double C = (b*b*sum_p_2 - R * R) / sum_m_p_a_b;
	double D = B * B - 4 * C;
	double t_1 = 0.5*(-B - sqrt(D));
	double t_2 = 0.5*(-B + sqrt(D));
	double t = mahds_max(t_1, t_2);
	
	mahds_matrix_set_all(m_prime, 0);
	for (size_t i = 0; i < p->lines; ++i) {
		for (size_t j = 0; j < p->columns; ++j) {
			double p_ij = mahds_matrix_get(p, i, j);
			double m_ij = mahds_matrix_get(m, i, j);
			mahds_matrix_set(m_prime, i, j,
				p_ij*b + t * (m_ij + p_ij * (a - b))
			);
		}
	}
	mahds_matrix_free(p);
}
