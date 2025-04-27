/* xmlstub.c — stubs for XML functions (GCG‐only build) */
#include <stdio.h>
#include <stdlib.h>
#include "score.h"   /* for ALN */

/* Always returns zero sequences */
int count_xml_seqs(FILE *fin) {
    (void)fin;
    return 0;
}

/* Never actually used under -DGCG; just abort if it somehow gets called */
ALN read_xml(FILE *fin, int first_seq) {
    fprintf(stderr, "Error: XML reader not supported in GCG‐only build.\n");
    exit(1);
}
