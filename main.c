#define _GNU_SOURCE
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include "find_prime.h"

int main(int argc, char **argv)
{
    if (argc < 2) {
        printf("Usage: %s decimal_digits\n", argv[0]);
        return 0;
    }

    char *s = find_prime(atoi(argv[1]));
    printf("%s\n",s);

    char *path;
    asprintf(&path,"results/%lu.txt", strlen(s));
    FILE *fout = fopen(path, "w");
    if (!fout) {
        fprintf(stderr, "Failed to open output file\n");
    } else {
        fputs(s, fout);
        fclose(fout);
    }

    free(s);
    free(path);
    return 0;
}
