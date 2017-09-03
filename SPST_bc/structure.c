#include <stdlib.h>
#include <stdio.h>
#include <string.h>
#include "header.h"
#include "structure.h"
#include "function.h"

struct level1 gp[4][32];
int thres[4], thres2[4], thres3[4], count2[4], numcombine;
int group[4][32], pnp[4], rnr[4], groupp[4], num_bucket[4], uni_bucket[4];
struct ENTRY3 *table3;

void groupping() {
    int i, na, max[4] = {0};

    char s[] = "start groupping ...";
    printf("%-40s", s);

    for(i = 0; i < num_entry; i++) {
        group[table[i].group][table[i].srcIP >> 27]++;
        groupp[table[i].group]++;
    }

    for(i = 0; i < 4; i++){
        for(na = 0; na < 32; na++) {
            if(max[i] < group[i][na]) {
                max[i] = group[i][na];
            }
        }
        //printf("group %c number of rules(max / total): %d / %d\n", i + 65, max[i], groupp[i]);
    }

    printf("finish\n");
}

void first_level() { 
    int i, j, g, N, na;

    char s[] = "start computing first level ...";
    printf("%-40s", s);

    for(i = 0; i < 4; i++) {
        for(na = 0; na < 32; na++) {
            gp[i][na].endpoint = (unsigned int *) malloc ((group[i][na] * 2 + 1) * sizeof(unsigned int));
            gp[i][na].lv2      = (struct level2 *) malloc ((group[i][na] * 2 + 1) * sizeof(struct level2));
            gp[i][na].n1       = (int *) malloc ((group[i][na] * 2 + 1) * sizeof(int));
            gp[i][na].n        = 1;
        }
    }

    unsigned int l, r, ip;
    int len;

    for(i = 0; i < num_entry; i++) {
        g = table[i].group;

        ip = table[i].srcIP;
        len = table[i].srclen;

        l = (ip == 0)  ?  0 : ip - 1;
        r = (len == 0) ? -1 : (((ip >> (32 - len)) + 1) << (32 - len)) - 1;

        unsigned int addr1 = table[i].srcIP >> 27;
        add_endpoint(gp[g][addr1].endpoint, &gp[g][addr1].n, l, r);
    }

    for(i = 0; i < num_entry; i++) {
        g = table[i].group;

        unsigned int addr1 = table[i].srcIP >> 27;
        interval_operation(1, i, table[i].srcIP, table[i].srclen, 0, 0, gp[g][addr1].n, gp[g][addr1].endpoint, gp[g][addr1].n1, NULL, 0);
    }

    for(i = 0; i < 4; i++) {
        for(na = 0; na < 32; na++) {
            for(j = 1; j < gp[i][na].n; j++) {
                if(gp[i][na].n1[j] > thres[i]) thres[i] = gp[i][na].n1[j];
            }
        }

        //printf("group %c bucket size(lv1): %d\n", i + 65, thres[i]);
    }

    for(i = 0; i < 4; i++) {
        for(na = 0; na < 32; na++) {
            N = gp[i][na].n;
            
            for(j = 0; j < N; j++) {
                if(gp[i][na].n1[j] != 0) {
                    gp[i][na].lv2[j].endpoint = (unsigned int *) malloc ((gp[i][na].n1[j] * 2 + 1) * sizeof(unsigned int));
                    gp[i][na].lv2[j].n2       = (int *) malloc ((gp[i][na].n1[j] * 2 + 1) * sizeof(int));
                    gp[i][na].lv2[j].b        = (struct bucket *) malloc ((gp[i][na].n1[j] * 2 + 1) * sizeof(struct bucket));
                    gp[i][na].lv2[j].n        = 1;
                    gp[i][na].lv2[j].rule     = (int *) malloc ((gp[i][na].n1[j] * 2 + 1) * sizeof(int));
                    gp[i][na].lv2[j].r        = 0;
                }
            }
        }
    }

    printf("finish\n");
}

void second_level() {
    int i, j, k, count, g, N, T, na;

    char s[] = "start computing second level ...";
    printf("%-40s", s);

    unsigned int l, r, ip, l2, r2, ip2;
    int len, len2;
    for(i = 0; i < num_entry; i++) {
        g = table[i].group;

        unsigned int addr1 = table[i].srcIP >> 27;
        interval_operation(2, i, table[i].srcIP, table[i].srclen, table[i].dstIP, table[i].dstlen, gp[g][addr1].n, gp[g][addr1].endpoint, NULL, gp[g][addr1].lv2, 0);
    }

    int rr, m;

    for(i = 0; i < 4; i++) {
        for(na = 0; na < 32; na++) {
            N = gp[i][na].n;

            for(j = 1; j < N; j++) {
                for(k = 0; k < gp[i][na].lv2[j].r; k++) {
                    rr = gp[i][na].lv2[j].rule[k];

                    unsigned int addr1 = table[rr].srcIP >> 27;
                    interval_operation(1, 0, table[rr].dstIP, table[rr].dstlen, 0, 0, gp[i][addr1].lv2[j].n, gp[i][addr1].lv2[j].endpoint, gp[i][addr1].lv2[j].n2, NULL, 0);
                }
            }

            for(j = 1; j < N; j++) {
                for(k = 1; k < gp[i][na].lv2[j].n; k++) {
                    if(gp[i][na].lv2[j].n2[k] > thres2[i]) thres2[i] = gp[i][na].lv2[j].n2[k];
                }
            }

            for(j = 1; j < N; j++) {
                for(k = 1; k < gp[i][na].lv2[j].n; k++){
                    gp[i][na].lv2[j].b[k].rule  = (int *) malloc (thres[i] * sizeof(int));
                    gp[i][na].lv2[j].b[k].rule2 = (int *) malloc (thres[i] * sizeof(int));
                    gp[i][na].lv2[j].b[k].set   = 0;
                    gp[i][na].lv2[j].b[k].r     = 0;
                    gp[i][na].lv2[j].b[k].r2    = 0;
                    gp[i][na].lv2[j].b[k].BV    = 0;
                }
            }

            for(j = 1; j < N; j++) {
                for(k = 0; k < gp[i][na].lv2[j].r; k++) {
                    rr = gp[i][na].lv2[j].rule[k];

                    unsigned int addr1 = table[rr].srcIP >> 27;
                    interval_operation(3, rr, table[rr].dstIP, table[rr].dstlen, 0, 0, gp[i][addr1].lv2[j].n, gp[i][addr1].lv2[j].endpoint, NULL, gp[i][addr1].lv2, j);
                }
            }
        }

        //printf("group %c bucket size(lv2): %d\n", i + 65, thres2[i]);
    }

    printf("finish\n");
}

void convert() {
    int i, j, k, l, m, n = 0, N, nn;

    char s[] = "start converting to new rule ID ...";
    printf("%-40s", s);

    table3 = (struct ENTRY3 *) malloc (2500 * sizeof(struct ENTRY3));
    
    for(i = 0; i < num_entry; i++) {
        if(table[i].srcPort[0] == table[i].srcPort[1])
            table[i].type += 4;
        else if(table[i].srcPort[0] != 0 || table[i].srcPort[1] != 65535)
            table[i].type += 8;

        if(table[i].dstPort[0] == table[i].dstPort[1])
            table[i].type += 1;
        else if(table[i].dstPort[0] != 0 || table[i].dstPort[1] != 65535)
            table[i].type += 2;

        for(j = 0; j < n; j++) {
            if(table[i].srcPort[0] == table3[j].Port[0] && table[i].srcPort[1] == table3[j].Port[1] && table[i].dstPort[0] == table3[j].Port[2] && table[i].dstPort[1] == table3[j].Port[3] && table[i].proto == table3[j].proto) {
                table[i].rule = j + 1;
                break;
            }
        }

        if(j == n) {
            table[i].rule = n + 1;
            table3[n].n = 1;
            table3[n].type = table[i].type;
            table3[n].Port[0] = table[i].srcPort[0];
            table3[n].Port[1] = table[i].srcPort[1];
            table3[n].Port[2] = table[i].dstPort[0];
            table3[n].Port[3] = table[i].dstPort[1];
            table3[n++].proto = table[i].proto;
        }
    }

    for(i = 0; i < n; i++) {
        switch(table3[i].type) {
            case  1:
                table3[i].Port[0] = table3[i].Port[2];  
                break;
            case  2:
                table3[i].Port[0] = table3[i].Port[2];
                table3[i].Port[1] = table3[i].Port[3];
                break;
            case  5:
                table3[i].Port[1] = table3[i].Port[2];
                break;
            case  6:
                table3[i].Port[0] = table3[i].Port[2];
                table3[i].Port[1] = table3[i].Port[3];
                table3[i].Port[2] = table3[i].Port[0];
                break;
            default:
                break;
        }
    }

    numcombine = n;
    //printf("number of combination of port & protocol: %d\n\n", n);
    nn = n;

    int *tmp, na;
    for(m = 0; m < 4; m++){
        for(na = 0; na < 32; na++) {
            N = gp[m][na].n;

            for(i = 1; i < N; i++) {
                for(j = 1; j < gp[m][na].lv2[i].n; j++) {
                    n = 0;

                    for(k = 0; k < gp[m][na].lv2[i].b[j].r; k++) {
                        for(l = 0; l < n; l++) {
                            if(gp[m][na].lv2[i].b[j].rule2[l] == table[gp[m][na].lv2[i].b[j].rule[k]].rule) break;
                        }

                        if(l == n) gp[m][na].lv2[i].b[j].rule2[n++] = table[gp[m][na].lv2[i].b[j].rule[k]].rule;
                    }
                    
                    gp[m][na].lv2[i].b[j].r2 = gp[m][na].lv2[i].b[j].r;
                    gp[m][na].lv2[i].b[j].r = n;
                    tmp = gp[m][na].lv2[i].b[j].rule;
                    gp[m][na].lv2[i].b[j].rule = gp[m][na].lv2[i].b[j].rule2;
                    gp[m][na].lv2[i].b[j].rule2 = tmp;

                    if(n > thres3[m]) thres3[m] = n;
                }
            }
        }
        //printf("group %c bucket size: %-4d\n", m + 65, thres3[m]);
    }

    printf("finish\n");
}

void software_compress() {
    int i, j, k, l, m, T, *pn, *rn, nor1, nor2, N; // nor = num of rule
    int aa, bb, zz, xx;
    unsigned int b;
    struct bucket *p[500000];
    int *R[250000];
    int na;

    char s[] = "start compressing bucket ...";
    printf("%-40s", s);
    
    for(m = 0; m < 4; m++) {
        pn = &pnp[m]; // wait to combine num
        rn = &rnr[m]; // bucket num
        T = thres3[m];
        //printf("compress group %c ...\n", m + 65);

        int same = 0;

        for(na = 0; na < 32; na++) {

            N = gp[m][na].n;

            for(i = 1; i < N; i++) {
                for(j = 1; j < gp[m][na].lv2[i].n; j++) {
                    nor1 = gp[m][na].lv2[i].b[j].r;
                    if(nor1 == 0) continue;

                    num_bucket[m]++;

                    if(nor1 == T) {
                        gp[m][na].lv2[i].b[j].set = 1;
                        for(k = 0; k < *rn; k++)
                            if(rule_check(R[k], gp[m][na].lv2[i].b[j].rule, T, T, 2) != 2) break;

                        if(k == *rn) {
                            R[(*rn)++] = gp[m][na].lv2[i].b[j].rule;
                        }
                    }

                    for(k = 0; k < *pn; k++) {
                        if(gp[m][na].lv2[i].b[j].r != (*p[k]).r) {
                            continue;
                        }

                        if(rule_check((*p[k]).rule, gp[m][na].lv2[i].b[j].rule, (*p[k]).r, (*p[k]).r, 2) != 2) {
                            gp[m][na].lv2[i].b[j].set = 1;
                            same++;
                            break;
                        }
                    }

                    if(k == *pn) p[(*pn)++] = &(gp[m][na].lv2[i].b[j]);
                }
            }
        }

        uni_bucket[m] = num_bucket[m] - same;
        //printf("same buckets: %d\n", same);
        //printf("remain %-6d buckets to compress ...\n", *pn);

        for(i = 0; i < *pn; i++) {
            if((*p[i]).set != 0)
                continue;
        
            for(j = 0; j < *pn; j++) {
                if((*p[j]).set == 2 || i == j)
                    continue;

                if(rule_check((*p[j]).rule, (*p[i]).rule, (*p[j]).r, (*p[i]).r, 2) != 2) {
                    (*p[i]).set = 2;
                    break;
                }
            }
        }

        int count = 0;
        for(i = 0; i < *pn; i++) {
            if((*p[i]).set == 0) count++;
        }
        //printf("remain %-6d buckets to compress ...\n", count);

        int max, max2, maxn, sum;
        for(i = 0; i < *pn; i++) {
            if((*p[i]).set != 0) continue;
        
            maxn = -1;
            max = 0;
            max2 = 0;
            for(j = 0; j < *pn; j++) {
                if((*p[j]).set != 0 || i == j) continue;

                sum = rule_check((*p[j]).rule, (*p[i]).rule, (*p[j]).r, (*p[i]).r, 3); 
                if(sum <= T) {
                    if((*p[i]).r + (*p[j]).r - sum > max2) {
                        max2 = (*p[i]).r + (*p[j]).r - sum;
                        max = sum;
                        maxn = j;
                    }
                    if((*p[i]).r + (*p[j]).r - sum == max2) {
                        if(sum > max) {
                            max = sum;
                            maxn = j;
                        }
                    }
                }
            }

            if(maxn == -1) {
                (*p[i]).set = 3;
                continue;
            }
                
            (*p[i]).set = 2;    
            for(j = 0; j < (*p[i]).r; j++) {
                l = 0;

                for(k = 0; k < (*p[maxn]).r; k++) {
                    if((*p[i]).rule[j] == (*p[maxn]).rule[k]) {
                        l = 1;
                        break;
                    }
                }

                if(l == 0) (*p[maxn]).rule[(*p[maxn]).r++] = (*p[i]).rule[j];
            }

            if((*p[maxn]).r == T) {
                (*p[maxn]).set = 1;
                R[(*rn)++] = (*p[maxn]).rule;
            }
        }

        for(i = 0; i < *pn; i++) {
            if((*p[i]).set == 0 || (*p[i]).set == 3) R[(*rn)++] = (*p[i]).rule;
        }

        count = 0;

        for(i = 0; i < *pn; i++) {
            if((*p[i]).set == 0 || (*p[i]).set == 3) {
                count++;
                count2[m] += T - (*p[i]).r;
            }
        }

        //printf("total buckets after compress: %d\n", *rn);
        //printf("bucket do not well compress: %d %d\n", count, count2[m]);
    }

    printf("finish\n\n");
}
