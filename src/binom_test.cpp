#include "binom.h"

#include <stdio.h>
#include <stdlib.h>

using namespace binom;
using namespace binom::detail;

static bool test(int maxn, int maxk)
{
	for (int n = 0; n <= maxn; n++)
	{
		if (fastbinomial(n, 1) != n)
		{
			printf("problem with k = 1\n");
			return false;
		}
	}

	for (int k = 2; k <= maxk; k++)
	{
		for (int n = k; n <= maxn; n++)
		{
			if (fastbinomial(n, k) != fastbinomial(n - 1, k) + fastbinomial(n - 1, k - 1))
			{
				printf("problem with Pascal's rule at k = %d n = %d\n", k, n);
				return false;
			}
		}
	}

	return true;
}

int main(int argc, char* argv[])
{
	for (int k = 2; k <= 64; k++)
	{
		if (!test(safen[k], k))
		{
			printf("bug\n");
			abort();
		}
	}

	if (!test(100, 10))
	{
		printf("bug\n");
		abort();
	}

	printf("OK\n");

	return EXIT_SUCCESS;
}

