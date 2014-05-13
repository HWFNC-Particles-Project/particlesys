#ifndef MALLOC_ALIGN_H_INCLUDED
#define MALLOC_ALIGN_H_INCLUDED

#include <stdint.h>
#include <stddef.h>

void *malloc_align(size_t size, size_t align, void **mem);

#endif // MALLOC_ALIGN_H_INCLUDED
