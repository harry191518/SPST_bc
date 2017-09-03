#include <stdlib.h>
#include <stdio.h>
#include <string.h>
#include "header.h"
#include "structure.h"
#include "function.h"
#include "result.h"

double onm1[4], onm2[4], offm[4];

void size() {
    int i, j, combine[4] = {0};
    double on1, on2, off, a[2] = {0}, b[2] = {0}, c[2] = {0}, soft, hard;
    int numofmb = 0, countcombine[4][2500] = {0};

    for(i = 0; i < num_entry; i++) {
        countcombine[table[i].group][table[i].rule]++;
    }
    for(i = 0; i < 4; i++) {
        for(j = 0; j < 2500; j++) {
            if(countcombine[i][j]) combine[i]++;
        }
    }

    for(i = 0; i < 4; i++) {
        numofmb += rnr[i];
        on1 = uni_bucket[i] * (count_bit(rnr[i]) + thres3[i]);
        on2 = rnr[i] * thres3[i] * 73;
        off = thres3[i] * uni_bucket[i] * count_bit(groupp[i]);

        onm1[i] = on1;
        onm2[i] = on2;
        offm[i] = off;
    }
}

static double level_node(int segtable[][2048], int level, int next0) {
    int i, j, node[10] = {0}, t = 0;

    for(i = 0; i < 32; i++) {
        for(j = 0; j < 2048; j++) {
            int c = 1, l = 0;
            if(segtable[i][j]) {
                while(segtable[i][j] > c) {
                    node[l++] += c;
                    segtable[i][j] -= c;
                    c *= 2;
                } 
                node[l] += segtable[i][j];
            }
        }
    }

    double m1 = 0;

    for(i = 0; i < level - 1; i++) {
        int b = (node[i + 1] > next0) ? node[i + 1] : next0;
        m1 += (1 + count_bit(b) * 2 + 16) * node[i];
        t += node[i];
    }
    m1 += (1 + count_bit(next0) * 2 + 16) * node[level - 1];
    t += node[level - 1];

    printf("    # of field 1 nodes : %d\n", t);
    printf("    memory(field 1 node information): %f KB\n", m1 / 8096);

    return m1;
}

void level() {
    int i, j, m, max, l1, l2, t1, t2, t1null, na, max1[4] = {0}, max2, tol[8] = {0}, numtol[8] = {0};
    int segtable[4][32][2048] = {0}, segtable2[4][32][2048] = {0}, numroot[4] = {0}, root5[4][32] = {0};

    for(m = 0; m < 4; m++) {
        for(na = 0; na < 32; na++) {
            for(i = 1; i < gp[m][na].n; i++) {
                unsigned int p = (gp[m][na].endpoint[i] >> 16) & 2047;
                segtable2[m][na][p]++;
            }
        }

        for(na = 0; na < 32; na++) {
            for(i = 0; i < 2048; i++) {
                segtable[m][na][i] = count_bit(segtable2[m][na][i]);
                if(segtable[m][na][i] > max1[m]) {
                    max1[m] = segtable[m][na][i];
                }

                if(segtable[m][na][i]) {
                    numroot[m]++;
                    root5[m][na] = 1;
                }
            }

            for(i = 0; i < 2048; i++) {
                if(segtable[m][na][i]) {
                    numtol[m]++;
                    tol[m] += segtable[m][na][i];
                }
            }
        }
    }

    for(m = 0; m < 4; m++) {
        double tolmem = 0;
        int SA5 = 0;
        for(i = 0; i < 32; i++) {
            if(root5[m][i]) {
                SA5++;
            }

            int SA11 = 0;
            for(j = 0; j < 2048; j++) {
                if(segtable[m][i][j]) {
                    SA11++;
                }
            }
        }

        tolmem += 32 * count_bit(SA5);
        tolmem += SA5 * 2048 * count_bit(numroot[m]);

        int node[10] = {0};
        t1     = 0;
        t1null = 0;
        t2     = 0;
        max    = 0;
        max2   = 0;
        for(na = 0; na < 32; na++) {
            t1 += gp[m][na].n - 1;

            for(i = 1; i < gp[m][na].n; i++) {
                if(gp[m][na].n1[i] == 0) {
                    t1null++;
                    continue;
                }

                int tmp = gp[m][na].lv2[i].n - 1;
                t2 += gp[m][na].lv2[i].n - 1;
                l2 = count_bit(gp[m][na].lv2[i].n - 1);
                unsigned int p = (gp[m][na].endpoint[i] >> 16) & 2047;

                if(max < segtable[m][na][p] + l2) {
                    max = segtable[m][na][p] + l2;
                }
                if(max2 < l2) max2 = l2;

                numtol[m + 4]++;
                tol[m + 4] += l2;

                int c = 1, l = 0;
                while(tmp > c) {
                    node[l++] += c;
                    tmp -= c;
                    c *= 2;
                } 

                node[l] += tmp;
            }
        }

        double lev1 = tol[m], lev2 = tol[m + 4];
        printf("group %c\n", m + 65);
        printf("    # of rules : %d\n", groupp[m]);
        printf("    # of 5-bit / 11-bit segmentation table root: %d / %d\n", SA5, numroot[m]);
        printf("    memory(segmentation table pointer): %f KB\n", tolmem / 8192);

        int f2n = 0;
        double m1, m2 = 0;

        m1 = level_node(segtable2[m], max1[m], node[0]);

        for(i = 0; i < max2 - 1; i++) {
            int b = (node[i + 1] > uni_bucket[m]) ? node[i + 1] : uni_bucket[m];
            m2  += (1 + count_bit(b) * 2 + 32) * node[i];
            f2n += node[i];
        }
        m2 += (1 + count_bit(uni_bucket[m]) * 2 + 32) * node[max2 - 1];
        f2n += node[max2 - 1];

        printf("    # of field 2 nodes : %d\n", f2n);
        printf("    memory(field 2 node information): %f KB\n", m2 / 8192);
        printf("    # of buckets / unique buckets: %d / %d\n", num_bucket[m], uni_bucket[m]);
        printf("    # of buckets after merged / bucket size : %d / %d\n", rnr[m], thres3[m]);
        printf("    memory(bucket information): %f KB\n", onm1[m] / 8192);
        printf("    memory(merged bucket set): %f KB\n", onm2[m] / 8192);
        printf("    total memory(on-chip): %f KB\n", (tolmem + m1 + m2 + onm1[m] + onm2[m]) / 8192);
        printf("    total memory(off-chip bitmap): %f KB\n", offm[m] / 8192);
        printf("\n");
    }
}

void f1_duplicate_count() {
    int i, j, k, m, T, N, na; // nor = num of rule
    int *R[250000];
    int rr[250000];
    
    printf("duplication of field 1 bucket\n");
    for(m = 0; m < 4; m++) {
        int uni = 0, total = 0;
        T = thres[m];

        for(na = 0; na < 32; na++) {
            N = gp[m][na].n;

            for(i = 1; i < N; i++) {
                if(gp[m][na].lv2[i].r == 0)
                    continue;

                total++;

                for(k = 0; k < uni; k++) {
                    if(gp[m][na].lv2[i].r != rr[k])
                        continue;

                    if(rule_check(R[k], gp[m][na].lv2[i].rule, T, T, 2) != 2)
                        break;
                }

                if(k == uni) {
                    rr[uni]  = gp[m][na].lv2[i].r;
                    R[uni++] = gp[m][na].lv2[i].rule;
                }
            }
        }

        float per = uni;
        per /= total;
        printf("    group %c uni/ total/ percent : %5d/ %5d/ %f\n", 'A' + m, uni, total, per);
    }
}

void f2_tree_duplicate() {
    int i, j, k, m, T, N, na; // nor = num of rule
    unsigned int *R[450000];
    int rr[450000];
    
    printf("\nduplication of field 2 tree\n");
    for(m = 0; m < 4; m++) {
        int max = 0;
        for(na = 0; na < 32; na++) {
            N = gp[m][na].n;

            for(i = 1; i < N; i++) {
                if(gp[m][na].lv2[i].n > max)
                    max = gp[m][na].lv2[i].n;
            }
        }

        int num = 0;
        for(na = 0; na < 32; na++) {
            N = gp[m][na].n;

            for(i = 1; i < N; i++) {
                if(gp[m][na].lv2[i].n == 0)
                    continue;

                rr[num]  = gp[m][na].lv2[i].n;
                R[num++] = (int *) malloc (max * sizeof(int));
                for(k = 0; k < max; k++)
                    R[num - 1][k] = gp[m][na].lv2[i].endpoint[k];
            }
        }

        int uni = 0;
        for(i = 0; i < num; i++) {
            for(j = i + 1; j < num; j++) {
                if(rr[i] > rr[j])
                    continue;

                if(rule_check(R[j], R[i], max, max, 2) != 2)
                    break;
            }

            if(j == num) {
                uni++;
            }
        }

        float per = uni;
        per /= num;
        printf("    group %c uni/ total/ percent : %5d/ %5d/ %f\n", 'A' + m, uni, num, per);
    }
}
