#include <stdlib.h>
#include <stdio.h>
#include <string.h>
#include "header.h"

int num_entry;
struct ENTRY *table;

static void read_table(char *str, int n) {
    int i = 2;
	char tok[] = "./@ :\t\n";
	char buf[100], *str1;
	unsigned int ip = 0, len;

    while(i--) {
        if(i == 1)
	        sprintf(buf, "%s\0", strtok(str, tok));
        if(i == 0)
	        sprintf(buf, "%s\0", strtok(NULL, tok));
	    ip = atoi(buf);
        ip <<= 8;
	    sprintf(buf, "%s\0", strtok(NULL, tok));
	    ip += atoi(buf);
        ip <<= 8;
	    sprintf(buf, "%s\0", strtok(NULL, tok));
	    ip += atoi(buf);
        ip <<= 8;
    	sprintf(buf, "%s\0", strtok(NULL, tok));
	    ip += atoi(buf);
	    str1 = (char *) strtok(NULL, tok);

	    if(str1 != NULL) {
		    sprintf(buf, "%s\0", str1);
	    	len = atoi(buf);
    	}

        if(i == 1) {
            table[n].srcIP  = ip;
            table[n].srclen = len;
        }
        else {
            table[n].dstIP  = ip;
            table[n].dstlen = len;
        }
    }

    if(table[n].srclen <= 2 && table[n].dstlen <= 2)
        table[n].group = 0;
    else if(table[n].srclen <= 2 && table[n].dstlen > 2)
        table[n].group = 1;
    else if(table[n].srclen > 2 && table[n].dstlen <= 2)
        table[n].group = 2;
    else if(table[n].srclen > 2 && table[n].dstlen > 2)
        table[n].group = 3;

	sprintf(buf, "%s\0", strtok(NULL, tok));
	table[n].srcPort[0] = atoi(buf);
	sprintf(buf, "%s\0", strtok(NULL, tok));
	table[n].srcPort[1] = atoi(buf);
	sprintf(buf, "%s\0", strtok(NULL, tok));
	table[n].dstPort[0] = atoi(buf);
	sprintf(buf, "%s\0", strtok(NULL, tok));
	table[n].dstPort[1] = atoi(buf);

	str1 = (char *) strtok(NULL, tok);
    if((str1[2] > 57 && str1[2] < 97) || (str1[3] > 57 && str1[3] < 97))
        printf("reading table error!\n");
    if(str1[3] > 57)
        table[n].proto = (str1[2] - 48) * 16 + (str1[3] - 87);
    else
        table[n].proto = (str1[2] - 48) * 16 + (str1[3] - 48);

    if(table[n].group == 1) {
        ip = table[n].srcIP;
        table[n].srcIP  = table[n].dstIP;
        table[n].dstIP  = ip;
        len             = table[n].srclen;
        table[n].srclen = table[n].dstlen;
        table[n].dstlen = len;
    }
}

void set_table(char *file_name) {
	FILE *fp;
	char string[100];
	fp = fopen(file_name, "r");

	while(fgets(string, 100, fp) != NULL) {
		num_entry++;
	}

	rewind(fp);
	table = (struct ENTRY *) malloc (num_entry * sizeof(struct ENTRY));
	num_entry = 0;

	while(fgets(string, 100, fp) != NULL) {
        table[num_entry].rule = num_entry + 1;
        read_table(string, num_entry++);
	}
}
