#include "GNFS.h"

bool NFS(const fmpz_t n);

void init()
{
	srand((unsigned)time(0));
}

int main(int argc, char *argv[])
{
	init();
	printf("linear\n");
	fmpz_t n; fmpz_init(n);
	char *ptr = NULL;
	char num[] = "268435457";
	if(argc > 1) ptr = argv[1];
	else ptr = num;
	fmpz_set_str(n,ptr,10);
	printf("n = "); fmpz_print(n); printf("\n");
	int Tries = 50;
	for(int i = 0; i < Tries; i++)
	{
		printf("--------------------------------------------------------------------------------\n");
		printf("Trying for the %d time...\n",i+1);
		if(NFS(n)) break;
		printf("--------------------------------------------------------------------------------\n");
	}
	fmpz_clear(n);
	return 0;
}
