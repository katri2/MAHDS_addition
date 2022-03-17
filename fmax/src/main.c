#include "fmax.h"


int main(int argc, char **argv) {
    double fmax;
    if (argc == 4) {
        fmax = get_fmax(argv[1], MAPPING_ARTIFICIAL, argv[2], argv[3]);
    } else {
        fmax = get_fmax(argv[1], MAPPING_ARTIFICIAL, "./cfgdir/system.cfg", "./cfgdir/alignment.cfg");
    }
    printf("%lf\n", fmax);
    return 0;
}
