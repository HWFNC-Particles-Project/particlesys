#define _GNU_SOURCE

#include <stdlib.h>
#include <stdio.h>

#include <sys/mman.h>

#include <codegen.h>

void free_code(void *fun) {
    if(fun != NULL) {
        binary_header *chunk = ((binary_header*)fun)-1;
        munmap(chunk, chunk->chunksize);
    }
}

void dump_code(void *fun) {
    binary_header *chunk = ((binary_header*)fun)-1;
    uint8_t *buffer = fun;

    FILE *fp = fopen("out.bin", "wb");
    fwrite(buffer, 1, chunk->size, fp);
    fclose(fp);
    if(system("ndisasm -b 64 out.bin")) {
        printf("ndisasm not available?\n");
    }
}
