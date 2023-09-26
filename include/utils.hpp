#ifndef FEM_UTILS_HPP
#define FEM_UTILS_HPP

#include <vector>
#include <unordered_map>
#include "LinearElasticity3D.hpp"

namespace utils {
    int binom(int, int);

    void checkParamsLE(LinearElasticity3D::ParamsLE&);

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