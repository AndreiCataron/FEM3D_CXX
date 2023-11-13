#ifndef FEM_UTILS_HPP
#define FEM_UTILS_HPP

#include <vector>
#include <unordered_map>
#include <set>
#include "params.hpp"

// hash used for hashing unordered_set<size_t>
struct MyHash
{
    std::size_t operator()(const std::set<std::size_t>& s) const noexcept
    {
        std::size_t h1 = std::hash<int>{}(int(*s.begin()));
        std::size_t h2 = std::hash<int>{}(int(*next(s.begin(), 1)));
        std::size_t h3 = std::hash<int>{}(int(*next(s.begin(), 2)));
        return h1 ^ (h2 << 1) ^ (h3 << 2);
    }
};

namespace utils {
    int binomialCoefficient(int, int);

    void checkParamsLE(ParamsLE&);

    void generateCombinations(std::vector<std::vector<int> >&, int, int);

    std::vector<double> midpoint(std::vector<double>&, std::vector<double>&);

    template <typename T> void deleteDuplicatesFromVector(std::vector<T> &vec) {
        std::sort(vec.begin(), vec.end());
        vec.erase(std::unique(vec.begin(), vec.end()), vec.end());
    }

    template<typename K, typename V>
    std::unordered_map<V, K> inverseMap(std::unordered_map<K, V> &map)
    {
        std::unordered_map<V, K> inv;
        std::for_each(map.begin(), map.end(),
                      [&inv] (const std::pair<K, V> &p) {
                          inv.insert(std::make_pair(p.second, p.first));
                      });
        return inv;
    }
}

#endif