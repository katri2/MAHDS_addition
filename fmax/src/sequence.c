#include "sequence.h"


void ins_symbols(mahds_symbolic_sequence* seq, size_t ind, size_t num, mahds_symbol *sym){
    seq->length+=num;
    seq->symbols = (mahds_symbol*)realloc((void*)seq->symbols, seq->length*sizeof(mahds_symbol));
    for (long i = seq->length-1; i > ind + num - 1; --i){
        seq->symbols[i] = seq->symbols[i-num];
    }
    for (size_t i = ind; i<ind+num; ++i){
        seq->symbols[i] = sym[i-ind];
    }
}

void ins_symbol(mahds_symbolic_sequence* seq, size_t ind, mahds_symbol sym){
    seq->length++;
    seq->symbols = (mahds_symbol*)realloc((void*)seq->symbols, seq->length*sizeof(mahds_symbol));
    for (long i = seq->length-1; i > ind; --i){
        seq->symbols[i] = seq->symbols[i-1];
    }
    seq->symbols[ind] = sym;
}

void del_symbol(mahds_symbolic_sequence* seq, size_t ind){
    for (size_t i = ind; i<seq->length-1; i++){
        seq->symbols[i] = seq->symbols[i+1];
    }
    seq->length--;
    seq->symbols = (mahds_symbol*)realloc((void*)seq->symbols, seq->length*sizeof(mahds_symbol));
}

void del_symbols(mahds_symbolic_sequence* seq, size_t from, size_t to){ // [from, to)
    size_t delta = to - from;
    for (size_t i = from; i < seq->length - delta; ++i) {
        seq->symbols[i] = seq->symbols[i+delta];
    }
    seq->length-=delta;
    seq->symbols = (mahds_symbol*)realloc((void*)seq->symbols, seq->length*sizeof(mahds_symbol));
}

double calc_mean_length(size_t num, mahds_symbolic_sequence** sequences) {
    double res = 0;
    for (size_t i = 0; i < num; ++i) {
        res += sequences[i]->length;
    }
    return res/(double)num;
}

mahds_symbolic_sequence* copy_mahds_seq (mahds_symbolic_sequence* seq) {
    mahds_symbolic_sequence* res = (mahds_symbolic_sequence*)malloc(sizeof(mahds_symbolic_sequence));
    res->length = seq->length;
    res->count_symbols = seq->count_symbols;
    res->symbols = (mahds_symbol*)malloc(res->length * sizeof(mahds_symbol));
    memcpy(res->symbols, seq->symbols, res->length * sizeof(mahds_symbol));
    return res;
}
