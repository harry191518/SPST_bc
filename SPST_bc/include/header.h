struct ENTRY {
    unsigned int srcIP, dstIP;
	unsigned char srclen, dstlen, proto;
	unsigned short int srcPort[2], dstPort[2], type;
	int rule, group;
};

extern int num_entry;
extern struct ENTRY *table;

void set_table(char *);
