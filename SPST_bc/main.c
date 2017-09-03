#include <stdio.h>
#include "header.h"
#include "structure.h"
#include "function.h"
#include "result.h"

int main(int argc, char *argv[]) {
    set_table(argv[1]);

    groupping();
    first_level();
    second_level();
    convert();
    software_compress();

    size();
    level();

    //f1_duplicate_count();
    //f2_tree_duplicate();
    
    return 0;
}
