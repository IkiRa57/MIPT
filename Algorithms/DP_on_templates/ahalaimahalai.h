#include <bits/stdc++.h>

typedef long long ll;
const ll MOD = MODULUS;

template<size_t M, size_t N, size_t mask = 0, size_t pos = 0, bool correct = 1>
class AhalaiMahalai {
public:
    static constexpr ll value = (AhalaiMahalai<M, N, mask ^ (1 << pos), pos + 1, (mask >> pos) & 1>::value +
                                 AhalaiMahalai<M, N, mask ^ (1 << pos), pos + 1, !((mask >> pos) & 1)>::value +
                                 AhalaiMahalai<M, N,
                                         mask ^ (1 << (pos + 1)),
                                         pos + 1, !((mask >> pos) & 1) && !((mask >> (pos + 1)) & 1)>::value) % MOD;
};

template<size_t M>
class AhalaiMahalai<M, 0, 0, 0, 1> {
public:
    static constexpr ll value = 1;
};

template<size_t M, size_t mask, size_t pos>
class AhalaiMahalai<M, 0, mask, pos, 1> {
public:
    static constexpr ll value = 0;
};

template<size_t M, size_t N, size_t mask>
class AhalaiMahalai<M, N, mask, M, 1> {
public:
    static constexpr ll value = AhalaiMahalai<M, N - 1, mask, 0, 1>::value;
};

template<size_t M, size_t N, size_t mask, size_t pos>
class AhalaiMahalai<M, N, mask, pos, 0> {
public:
    static constexpr ll value = 0;
};

