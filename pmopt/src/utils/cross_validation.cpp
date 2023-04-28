#include "pmopt/utils/cross_validation.hpp"

#include <unordered_map>
#include <numeric>
#include <random>

using namespace PMOpt;


unsigned int Utils::make_leave_group_out(
    std::vector<unsigned int> & indices,
    const std::vector<unsigned int> & group
)
noexcept
{
    std::unordered_map<unsigned int, unsigned int> group2fold;

    indices.clear();
    indices.reserve(group.size());
    for (auto g : group)
    {
        if (group2fold.find(g) == group2fold.end())
        {
            // fold_id must be >= 1
            // fold_id = 0 for entire dataset fitting
            group2fold[g] = group2fold.size() + 1; 
        }
        indices.push_back(group2fold[g]);
    }
    return group2fold.size() + 1;
}


unsigned int Utils::make_leave_one_out(
    std::vector<unsigned int> & indices,
    unsigned int size,
    unsigned int seed,
    bool shuffle
)
noexcept
{
    indices.resize(size);
    std::iota(indices.begin(), indices.end(), 1);
    
    if (shuffle)
    {
        std::mt19937 random_generator(seed);
        std::shuffle(indices.begin(), indices.end(), random_generator);
    }
    
    return size + 1;
}


unsigned int Utils::make_k_fold(
    std::vector<unsigned int> & indices,
    unsigned int size,
    unsigned int n_folds,
    unsigned int seed
)
noexcept
{
    (void) size;
    (void) seed;
    (void) indices;
    return n_folds + 1;
}


void Utils::make_test_mask(
    std::vector<unsigned int> & mask,
    unsigned int & n_training,
    const std::vector<unsigned int> & fold_ids,
    unsigned int fold_id
)
noexcept
{
    // mask[i] is zero if fold id of the i-th instance equals to `fold_id`.
    // This means the i-th is a test instance.
    auto n_instances = fold_ids.size();
    mask.resize(n_instances);
    n_training = 0;
    for (unsigned int i_instance = 0; i_instance < n_instances; ++i_instance)
    {
        mask[i_instance] = fold_ids[i_instance] != fold_id;
        n_training += mask[i_instance];
    }
}