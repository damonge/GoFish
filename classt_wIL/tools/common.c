#include "common.h"

void class_protect_sprintf(char* dest, char* tpl,...) {
  va_list args;
  va_start(args,tpl);
  vsnprintf(dest, 2048,tpl,args);
  va_end(args);
}

void class_protect_fprintf(FILE* stream, char* tpl,...) {
  va_list args;
  char dest[6000];
  va_start(args,tpl);
  vsnprintf(dest, 2048,tpl,args);
  va_end(args);
  fprintf(stream,dest);
}

void* class_protect_memcpy(void* dest, void* from, size_t sz) {
  return memcpy(dest, from,sz);
}

void my_abort(const int errtyp,const char *fmt, ...)
{
    va_list argp;
    
    va_start(argp,fmt);
    vfprintf(stderr,fmt,argp);
    va_end(argp);

    exit(errtyp);
}
