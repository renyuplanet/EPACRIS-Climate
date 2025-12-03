#include <stdio.h>
#include <time.h>

void datetime() 
{
        time_t mytime;
        time(&mytime);
        printf("%s", ctime(&mytime));
}

void pexit(void) {printf("\n%s\n","=========================\n==== EXITING PROGRAM ====\n=========================\n");} 