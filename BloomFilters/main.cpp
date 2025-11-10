#include <bits/stdc++.h>
using namespace std;

struct BloomFilter {
    // Bitset stored in a byte vector for compactness and speed.
    vector<uint8_t> bits;
    uint64_t m; // number of bits
    uint32_t k; // number of hash functions

    explicit BloomFilter(uint64_t m_bits, uint32_t k_hashes)
        : bits((m_bits + 7) / 8, 0), m(m_bits), k(k_hashes) {}

    // Construct from expected items n and target false positive rate eps in (0,1).
    static BloomFilter from_n_eps(uint64_t n, double eps) {
        if (n == 0) n = 1;
        if (!(eps > 0.0 && eps < 1.0)) eps = 0.01;
        // m ≈ -n ln ε / (ln 2)^2, k ≈ (m/n) ln 2
        long double ln2 = 0.693147180559945309417232121458176568L;
        long double m_ld = (- (long double)n * log(eps)) / (ln2 * ln2);
        uint64_t m_bits = (uint64_t)ceil(m_ld);
        uint32_t k_hashes = (uint32_t)ceil((m_ld / (long double)n) * ln2);
        if (k_hashes == 0) k_hashes = 1;
        return BloomFilter(m_bits, k_hashes);
    }

    // Simple 64-bit splitmix64 hash; good distribution and speed.
    static inline uint64_t splitmix64(uint64_t x) {
        x += 0x9E3779B97F4A7C15ull;
        x = (x ^ (x >> 30)) * 0xBF58476D1CE4E5B9ull;
        x = (x ^ (x >> 27)) * 0x94D049BB133111EBull;
        return x ^ (x >> 31);
    }

    // Mix a byte span into a 64-bit seed using a simple streaming mixer.
    static uint64_t mix_bytes(const uint8_t* data, size_t len, uint64_t seed) {
        // Feed 8 bytes at a time; tail handled naturally.
        const uint64_t mul = 0x9ddfea08eb382d69ull; // from wyhash finalizer inspiration
        uint64_t h = seed ^ (len * mul);
        const uint8_t* p = data;
        while (len >= 8) {
            uint64_t v;
            memcpy(&v, p, 8);
            p += 8; len -= 8;
            h ^= splitmix64(v);
            h *= mul;
            h ^= (h >> 47);
        }
        // Tail
        uint64_t t = 0;
        switch (len) {
            case 7: t ^= (uint64_t)p[6] << 48;
            case 6: t ^= (uint64_t)p[5] << 40;
            case 5: t ^= (uint64_t)p[4] << 32;
            case 4: t ^= (uint64_t)p[3] << 24;
            case 3: t ^= (uint64_t)p[2] << 16;
            case 2: t ^= (uint64_t)p[1] << 8;
            case 1: t ^= (uint64_t)p[0];
            default: break;
        }
        if (len) {
            h ^= splitmix64(t);
            h *= mul;
            h ^= (h >> 47);
        }
        // Final avalanche
        h = splitmix64(h);
        return h;
    }

    // Double hashing: generate k indices using h1 + i*h2 mod m (Kirsch–Mitzenmacher).
    // See: “Less Hashing, Same Performance: Building a Better Bloom Filter”.
    inline pair<uint64_t,uint64_t> base_hashes(const uint8_t* data, size_t len) const {
        // Two independent seeds
        uint64_t h1 = mix_bytes(data, len, 0x1234567890abcdefull);
        uint64_t h2 = mix_bytes(data, len, 0xfedcba0987654321ull);
        // Avoid h2 being multiple of m causing cycles; force odd.
        if ((h2 & 1ull) == 0) h2 ^= 1ull;
        return {h1, h2};
    }

    inline void set_bit(uint64_t idx) {
        uint64_t byte = idx >> 3;
        uint8_t mask = 1u << (idx & 7u);
        bits[byte] |= mask;
    }

    inline bool get_bit(uint64_t idx) const {
        uint64_t byte = idx >> 3;
        uint8_t mask = 1u << (idx & 7u);
        return (bits[byte] & mask) != 0;
    }

    void add_bytes(const uint8_t* data, size_t len) {
        auto [h1, h2] = base_hashes(data, len);
        // Compute k indices: (h1 + i*h2) % m
        // Use 128-bit for product to reduce bias; then modulo m.
        for (uint32_t i = 0; i < k; ++i) {
            __uint128_t idx = ( (__uint128_t)h1 + (__uint128_t)i * (__uint128_t)h2 );
            uint64_t bit = (uint64_t)(idx % m);
            set_bit(bit);
        }
    }

    bool possibly_contains_bytes(const uint8_t* data, size_t len) const {
        auto [h1, h2] = base_hashes(data, len);
        for (uint32_t i = 0; i < k; ++i) {
            __uint128_t idx = ( (__uint128_t)h1 + (__uint128_t)i * (__uint128_t)h2 );
            uint64_t bit = (uint64_t)(idx % m);
            if (!get_bit(bit)) return false;
        }
        return true;
    }

    // String helpers
    inline void add(const string& s) { add_bytes(reinterpret_cast<const uint8_t*>(s.data()), s.size()); }
    inline bool possibly_contains(const string& s) const { return possibly_contains_bytes(reinterpret_cast<const uint8_t*>(s.data()), s.size()); }

    // Int helpers
    inline void add(uint64_t x) { add_bytes(reinterpret_cast<const uint8_t*>(&x), sizeof(x)); }
    inline bool possibly_contains(uint64_t x) const { return possibly_contains_bytes(reinterpret_cast<const uint8_t*>(&x), sizeof(x)); }

    // Estimate current false-positive probability after inserting n elements:
    // p ≈ (1 - e^{-k n / m})^k
    double fp_rate_estimate(uint64_t n_inserted) const {
        long double kd_m = ((long double)k * (long double)n_inserted) / (long double)m;
        long double p = pow((long double)1.0 - expl(-kd_m), (long double)k);
        return (double)p;
    }

    // Convenience: return bit density (fraction of 1s)
    double fill_ratio() const {
        uint64_t ones = 0;
        for (uint8_t b : bits) ones += __builtin_popcount(b);
        return (double)ones / (double)m;
    }
};

int main(int argc, char** argv) {
    ios::sync_with_stdio(false);
    cin.tie(nullptr);
    // Example usage:
    // Choose by target items and false positive rate
    uint64_t expected_n = 1'000;      // expected insertions
    double target_fp = 0.01;          // ~1% false positive
    BloomFilter bf = BloomFilter::from_n_eps(expected_n, target_fp);

    // Insert some keys
    vector<string> keys = {"alice", "bob", "carol", "dave", "erin", "mallory"};
    for (const auto& s : keys) bf.add(s);

    // Queries
    vector<string> tests = {"alice", "bob", "trent", "mallory", "oscar", "carol", "peggy"};
    cout << "m(bits)=" << bf.m << " k=" << bf.k << "\n";
    cout << "fill_ratio=" << fixed << setprecision(4) << bf.fill_ratio() << "\n";
    cout << "fp_rate_estimate(n=" << keys.size() << ")=" << bf.fp_rate_estimate((uint64_t)keys.size()) << "\n";
    for (const auto& q : tests) {
        bool maybe = bf.possibly_contains(q);
        cout << q << ": " << (maybe ? "possibly present" : "definitely absent") << "\n";
    }

    // Optional: command-line mode
    // Usage: ./a.out n eps; then read words from stdin to insert and query.
    if (argc >= 3) {
        uint64_t n = strtoull(argv[1], nullptr, 10);
        double eps = strtod(argv[2], nullptr);
        BloomFilter cli_bf = BloomFilter::from_n_eps(n, eps);
        vector<string> inserted;
        string line;
        cerr << "Insert words, empty line to end:\n";
        while (true) {
            if (!getline(cin, line)) break;
            if (line.empty()) break;
            cli_bf.add(line);
            inserted.push_back(line);
        }
        cerr << "Query words, Ctrl-D to end:\n";
        while (getline(cin, line)) {
            bool maybe = cli_bf.possibly_contains(line);
            cout << (maybe ? "maybe " : "no ") << line << "\n";
        }
        cerr << "m=" << cli_bf.m << " k=" << cli_bf.k
            << " fill=" << cli_bf.fill_ratio()
            << " fp_est=" << cli_bf.fp_rate_estimate((uint64_t)inserted.size()) << "\n";
    }

    // Print separator
    cout << "--------------------------\n";

    expected_n = 5000;
    target_fp = 0.01;
    BloomFilter bf2 = BloomFilter::from_n_eps(expected_n, target_fp);

    // Create 1000 random integers and insert half of them
    random_device rd;
    mt19937_64 rng(rd());
    vector<uint64_t> random_ints(20000);
    for (size_t i = 0; i < random_ints.size(); ++i) {
        random_ints[i] = rng();
    }
    for (size_t i = 0; i < random_ints.size(); i += 2) {
        bf2.add(random_ints[i]);
    }

    // Query all integers and report results
    int false_positives = 0;
    for (size_t i = 0; i < random_ints.size(); ++i) {
        bool maybe = bf2.possibly_contains(random_ints[i]);
        bool actually_present = (i % 2 == 0);
        if (maybe && !actually_present) {
            false_positives++;
        }
        if(!maybe && actually_present) {
            cout << "Error: false negative for " << random_ints[i] << "\n";
        }
    }

    cout << "Integer test:\n";
    cout << "m(bits)=" << bf2.m << " k=" << bf2.k << "\n";
    cout << "fill_ratio=" << fixed << setprecision(4) << bf2.fill_ratio() << "\n";
    cout << "fp_rate_estimate(n=" << random_ints.size() / 2 << ")=" << bf2.fp_rate_estimate(random_ints.size() / 2) << "\n";
    cout << "False positives: " << false_positives << " out of " << (random_ints.size() / 2) << "\n";
    return 0;
}