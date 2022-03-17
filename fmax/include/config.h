#pragma once

#include "sequence.h"
#include <stdio.h>
#include <string.h>
#include <stdlib.h>

#define CONSTANT
#define OUT
#define _CRT_SECURE_NO_WARNINGS

#define PROT_SYM_COUNT         20
#define DNA_SYM_COUNT          4
#define K_D_CONST              -0.1
#define R_MULT_CONST           5.
#define D_PRICE                40.
#define E_PRICE                D_PRICE/10.
#define STABLE_DIST_C          1.12515
#define BORDERS                60
#define BORDERS_C              2
#define CHECKING_SET           100
#define SELECTION_SET_SIZE     1000000
#define PROT_SYMBOLS           20
#define PERIODS_SCATTER        5
#define SYMBOLS_IN_CFG_STR     256
#define THREADS_NUM            7
#define DEFAULT_MATR_PATH      "./init_matrix_sets/"
#define DEFAULT_SYS_CFG_PATH   "./config/system.cfg"
#define DEFAULT_AL_CFG_PATH    "./config/alignment.cfg"

#define MAHDS_ALIGNMENT_TYPE_LOCAL  ((char)0)
#define MAHDS_ALIGNMENT_TYPE_GLOBAL ((char)1)
#define MAHDS_ALIGNMENT_TYPE_MIXED  ((char)2)

typedef struct {
	size_t  CONSTANT size_of_symbols_set;
	double  CONSTANT d_price;
	double  CONSTANT e_price;
	double  CONSTANT K_d;
	double  CONSTANT R_sqr_mult_const;
	size_t  CONSTANT alignment_matrix_borders;
	size_t  CONSTANT size_of_checking_matrix_set;
	char    CONSTANT alignment_type;
	size_t           spectre_min;
	size_t           spectre_max;
	size_t  CONSTANT threads_num;
	char*   CONSTANT matrix_path;
} mahds_config;


void load_config(
	char*           system_settings_path,
	char*           alignment_settings_path,
	size_t          symbol_types,
	size_t          spectre_middle,
	OUT mahds_config* config
); // already got seq
