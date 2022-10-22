#ifndef BINOM_H
#define BINOM_H

#include <limits>
#include <stdint.h>

namespace binom {

namespace detail {

struct fastdiv_t
{
	int shift;
	uint64_t inverse;
};

constexpr const fastdiv_t precomputed[] =
{
	{ 0, 0 },
	{ 0, 0x1 },
	{ 1, 0x1 },
	{ 0, 0xaaaaaaaaaaaaaaab },
	{ 2, 0x1 },
	{ 0, 0xcccccccccccccccd },
	{ 1, 0xaaaaaaaaaaaaaaab },
	{ 0, 0x6db6db6db6db6db7 },
	{ 3, 0x1 },
	{ 0, 0x8e38e38e38e38e39 },
	{ 1, 0xcccccccccccccccd },
	{ 0, 0x2e8ba2e8ba2e8ba3 },
	{ 2, 0xaaaaaaaaaaaaaaab },
	{ 0, 0x4ec4ec4ec4ec4ec5 },
	{ 1, 0x6db6db6db6db6db7 },
	{ 0, 0xeeeeeeeeeeeeeeef },
	{ 4, 0x1 },
	{ 0, 0xf0f0f0f0f0f0f0f1 },
	{ 1, 0x8e38e38e38e38e39 },
	{ 0, 0x86bca1af286bca1b },
	{ 2, 0xcccccccccccccccd },
	{ 0, 0xcf3cf3cf3cf3cf3d },
	{ 1, 0x2e8ba2e8ba2e8ba3 },
	{ 0, 0xd37a6f4de9bd37a7 },
	{ 3, 0xaaaaaaaaaaaaaaab },
	{ 0, 0x8f5c28f5c28f5c29 },
	{ 1, 0x4ec4ec4ec4ec4ec5 },
	{ 0, 0x84bda12f684bda13 },
	{ 2, 0x6db6db6db6db6db7 },
	{ 0, 0x34f72c234f72c235 },
	{ 1, 0xeeeeeeeeeeeeeeef },
	{ 0, 0xef7bdef7bdef7bdf },
	{ 5, 0x1 },
	{ 0, 0xf83e0f83e0f83e1 },
	{ 1, 0xf0f0f0f0f0f0f0f1 },
	{ 0, 0xaf8af8af8af8af8b },
	{ 2, 0x8e38e38e38e38e39 },
	{ 0, 0x14c1bacf914c1bad },
	{ 1, 0x86bca1af286bca1b },
	{ 0, 0x6f96f96f96f96f97 },
	{ 3, 0xcccccccccccccccd },
	{ 0, 0x8f9c18f9c18f9c19 },
	{ 1, 0xcf3cf3cf3cf3cf3d },
	{ 0, 0x82fa0be82fa0be83 },
	{ 2, 0x2e8ba2e8ba2e8ba3 },
	{ 0, 0x4fa4fa4fa4fa4fa5 },
	{ 1, 0xd37a6f4de9bd37a7 },
	{ 0, 0x51b3bea3677d46cf },
	{ 4, 0xaaaaaaaaaaaaaaab },
	{ 0, 0x7d6343eb1a1f58d1 },
	{ 1, 0x8f5c28f5c28f5c29 },
	{ 0, 0xfafafafafafafafb },
	{ 2, 0x4ec4ec4ec4ec4ec5 },
	{ 0, 0x21cfb2b78c13521d },
	{ 1, 0x84bda12f684bda13 },
	{ 0, 0x6fb586fb586fb587 },
	{ 3, 0x6db6db6db6db6db7 },
	{ 0, 0x823ee08fb823ee09 },
	{ 1, 0x34f72c234f72c235 },
	{ 0, 0xcbeea4e1a08ad8f3 },
	{ 2, 0xeeeeeeeeeeeeeeef },
	{ 0, 0x4fbcda3ac10c9715 },
	{ 1, 0xef7bdef7bdef7bdf },
	{ 0, 0xefbefbefbefbefbf },
};

constexpr const int safen[] =
{
	0,   std::numeric_limits<int>::max(), 2642246, 77936, 10206, 2762, 1122, 585, 359, 247,
	184, 146,                             121,     104,   92,    83,   77,   72,  68,  65,
	63,  61,                              59,      58,    58,    57,   57,   56,  56,  56,
	56,  56,                              57,      57,    57,    58,   58,   58,  59,  59,
	60,  61,                              61,      62,    62,    62,   62,   62,  62,  62,
	62,  62,                              62,      62,    62,    62,   62,   62,  62,  62,
	62,  62,                              62,      62,    62
};

} // namespace detail

// correct for n <= 100, k <= 10
template<
	int n,
	int k
>
inline bool fastbinomial(uint64_t& result)
{
	if (0 == k || n == k)
	{
		result = 1;
		return false;
	}

	if (k > n)
	{
		result = 0;
		return false;
	}

	if (1 == k)
	{
		result = n;
		return false;
	}

	result = std::numeric_limits<uint64_t>::max();

	if ((k <= 0) || (n <= 0))
		return true;

	if (k > 64)
		return true;

	if (n > detail::safen[k])
		return true;

	// Make k as small as possible.
	if (k > n / 2) k = n - k;

	uint64_t np = n - k;
	result = np + 1;
	for (uint64_t z = 2; z <= (uint64_t)k; z++)
	{
		result =
			result *
			(np + z); // this could overflow! but it won't for our
				  // range of values for n and k: require that n <= safen[k]

		auto f = detail::precomputed[z];
		result >>= f.shift;
		result *= f.inverse;
	}

	return false;
}

// correct for n <= 100, k <= 10
inline bool fastbinomial(int n, int k, uint64_t& result)
{
	if (0 == k || n == k)
	{
		result = 1;
		return false;
	}

	if (k > n)
	{
		result = 0;
		return false;
	}

	if (1 == k)
	{
		result = n;
		return false;
	}

	result = std::numeric_limits<uint64_t>::max();

	if ((k <= 0) || (n <= 0))
		return true;

	if (k > 64)
		return true;

	if (n > detail::safen[k])
		return true;

	// Make k as small as possible.
	if (k > n / 2) k = n - k;

	uint64_t np = n - k;
	result = np + 1;
	for (uint64_t z = 2; z <= (uint64_t)k; z++)
	{
		result =
			result *
			(np + z); // this could overflow! but it won't for our
				  // range of values for n and k: require that n <= safen[k]

		auto f = detail::precomputed[z];
		result >>= f.shift;
		result *= f.inverse;
	}

	return false;
}

} // namespace binom

#endif // BINOM_H

