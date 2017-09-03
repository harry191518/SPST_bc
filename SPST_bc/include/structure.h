struct ENTRY3 {
    unsigned short int Port[4], type, n;
    unsigned char proto;
};

struct level1 { 
    unsigned int *endpoint;
    struct level2 *lv2;
    int n, *n1; // n = total endpoint of level 1, n1 = number of each interval
};

struct level2 {
    unsigned int *endpoint;
    struct bucket *b;
    int n, *n2, *rule, r; // n2 = number of each interval
};

struct bucket {
    int *rule, r, set, *rule2, r2;
    unsigned int BV;
};

extern struct level1 gp[4][32];
extern int thres[4], thres2[4], thres3[4], count2[4], numcombine;
extern int group[4][32], pnp[4], rnr[4], groupp[4], num_bucket[4], uni_bucket[4];
extern struct ENTRY3 *table3;

void groupping();
void first_level();
void second_level();
void convert();
void software_compress();
