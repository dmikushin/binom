#include "binom.h"

#include <stdio.h>

#include "gtest/gtest.h"

using namespace binom;
using namespace binom::detail;

static bool test(int maxn, int maxk)
{
	for (int n = 0; n <= maxn; n++)
	{
		uint64_t result;
		bool overflow = fastbinomial(n, 1, result);
		if ((result != n) && !overflow)
		{
			printf("problem with k = 1\n");
			return false;
		}
	}

	for (int k = 2; k <= maxk; k++)
	{
		for (int n = k; n <= maxn; n++)
		{
			uint64_t result1, result2, result3;
			bool overflow =
				fastbinomial(n, k, result1) |
				fastbinomial(n - 1, k, result2) |
				fastbinomial(n - 1, k - 1, result3);
			if ((result1 != result2 + result3) && !overflow)
			{
				printf("problem with Pascal's rule at k = %d n = %d\n", k, n);
				return false;
			}
		}
	}

	return true;
}

TEST(binom, test)
{
	for (int k = 2; k <= 64; k++)
		ASSERT_TRUE(test(safen[k], k));

	ASSERT_TRUE(test(100, 10));
}

int main(int argc, char ** argv)
{
	::testing::InitGoogleTest( & argc, argv);
	return RUN_ALL_TESTS();
}

