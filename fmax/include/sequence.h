#pragma once

#include <stdlib.h>
#include <limits.h>
#include <string.h>

#define SYMBOL_MAX USHRT_MAX
typedef unsigned short mahds_symbol;

typedef struct {
	size_t count_symbols;
	size_t length;
	mahds_symbol* symbols;  //=chars 0-x.     (Input symbols are mappend on 0-x)
} mahds_symbolic_sequence;

typedef struct {
	mahds_symbolic_sequence sequence;
	size_t orig_begin; // 0th element
	size_t orig_end; // just like length so idx < orig_end
} mahds_alignment_sequence;

typedef struct {
	size_t                   seq_number;
	mahds_symbolic_sequence** sequences;
} mahds_seq_set;

typedef struct {
	size_t                    seq_number;
	mahds_alignment_sequence** sequences;
} mahds_al_set;


void ins_symbols(mahds_symbolic_sequence* seq, size_t ind, size_t num, mahds_symbol *sym);
void ins_symbol(mahds_symbolic_sequence* seq, size_t ind, mahds_symbol sym);
void del_symbol(mahds_symbolic_sequence* seq, size_t ind);
void del_symbols(mahds_symbolic_sequence* seq, size_t from, size_t to); // [from, to)
double calc_mean_length(size_t num, mahds_symbolic_sequence** sequences);
mahds_symbolic_sequence* copy_mahds_seq (mahds_symbolic_sequence* seq);
