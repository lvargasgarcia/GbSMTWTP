#ifndef NEIGHBORHOOD_HPP
#define NEIGHBORHOOD_HPP

#include <memory>
#include <vector>
#include <algorithm> // Required for std::next_permutation
#include <string>
#include <random>
#include <iostream>

using permutation_t = std::vector<int>;
using perms_set_t = std::vector<permutation_t>;

long factorial (int n) {
    if (n <= 1) return 1;
    return n * factorial(n - 1);
}

permutation_t compose(const permutation_t& a, const permutation_t& b, int n) {
    permutation_t result(n);
    for (int i = 0; i < n; ++i) {
        result[i] = a[b[i]];
    }
    return result;
}

permutation_t inverse(const permutation_t& a, int n) {
    permutation_t result(n);
    for (int i = 0; i < n; ++i) {
        result[a[i]] = i;
    }
    return result;
}

std::vector<std::vector<int>> permutations(int k) {
    
    std::vector<std::vector<int>> result;
    std::vector<int> perm(k);
    
    for (int i = 0; i < k; ++i) {
        perm[i] = i;
    }

    std::next_permutation(perm.begin(), perm.end());

    do {
        std::vector<int> new_perm(k);
        std::copy(perm.begin(), perm.end(), new_perm.begin());
        result.push_back(new_perm);
    } while(std::next_permutation(perm.begin(), perm.end()));


    return result;
}

// void generate_permutations_set(int n, int position, int k, perms_set_t& perms_set, long perms_count) {
        
//     auto perms = permutations(k);
    
//     for(size_t j = 0; j < perms_count; j++){
        
//         perms_set[position * perms_count + j] = permutation_t(n);
//         std::copy(perms[j].begin(), perms[j].end(), perms_set[position * perms_count + j].begin());
//     }
    

// }

class Neighborhood {

    public:
        
        perms_set_t perms_set;
        int n;
        int k;
        long size;
        int categories;
        long perms_count;

        std::vector<std::pair<int,int>> improving_moves_set;
        std::vector<int> improving_moves_indexes;
        std::vector<long> scores;
        std::vector<long> common_tardiness_set;
        // std::vector<std::vector<int>> independent_moves_matrix;

        Neighborhood() = default;
        Neighborhood(const Neighborhood&) = default;

        Neighborhood(int n, int k) {
            
            this->n = n;
            this->k = k;
            this->categories = n - k + 1;
            this->perms_count = (factorial(k) - 1);
            this->size = this->perms_count * this->categories;
            this->perms_set = permutations(k);
            this->common_tardiness_set = std::vector<long>(this->size);
            this->improving_moves_indexes = std::vector<int>(this->size);
            this->scores = std::vector<long>(this->size);


            for (size_t i = 0; i < this->size; i++) {
                this->improving_moves_indexes[i] = -1;
                this->scores[i] = 0;
            }

            // generate_permutations_set(n,i,k,this->perms_set, this->perms_count);

            // this->independent_moves_matrix = std::vector<std::vector<int>>(this->size);

            // for(size_t i = 0; i < this->size; i++){
                
            //     int lower_bound = std::max(0, static_cast<int>((i/this->perms_count) - this->k + 1));
            //     int upper_bound = std::min(this->categories - 1, static_cast<int>((i/this->perms_count) + this->k - 1));
                
            //     for(size_t j = lower_bound*this->perms_count; j < (upper_bound + 1)*this->perms_count; j++){
            //         if(!this->independent(i,j)){
            //             this->independent_moves_matrix[i].push_back(j);
            //         }
            //     }
            
            // }

        
        }

        // bool independent(int i, int j) {
            
        //     int category_i = i / this->perms_count;
        //     int category_j = j / this->perms_count;

        //     int position_i = i % this->perms_count;
        //     int position_j = j % this->perms_count;

        //     auto perm_i = this->perms_set[position_i];
        //     auto perm_j = this->perms_set[position_j];

        //     auto interval_i = std::make_pair(-1,-1);
        //     auto intrerval_j = std::make_pair(-1,-1);

        //     int k = 0;

        //     while(k < perm_i.size() && interval_i.first == -1) {
        //         if(perm_i[k] != k) {
        //             interval_i.first = k + category_i;
        //         }
        //         k++;
        //     }

        //     k = perm_i.size() - 1;

        //     while(k >= 0 && interval_i.second == -1) {
        //         if(perm_i[k] != k) {
        //             interval_i.second = k + category_i;
        //         }
        //         k--;
        //     }

        //     k = 0;

        //     while(k < perm_j.size() && intrerval_j.first == -1) {
        //         if(perm_j[k] != k) {
        //             intrerval_j.first = k + category_j;
        //         }
        //         k++;
        //     }

        //     k = perm_j.size() - 1;

        //     while(k >= 0 && intrerval_j.second == -1) {
        //         if(perm_j[k] != k) {
        //             intrerval_j.second = k + category_j;
        //         }
        //         k--;
        //     }

        //     return (interval_i.first > intrerval_j.second || interval_i.second < intrerval_j.first);

        // }

        void insert_move(const int position, const int index) {
            
            if(this->improving_moves_indexes[index] == -1) {
                
                this->improving_moves_set.push_back(std::make_pair(position, index));
                this->improving_moves_indexes[index] = this->improving_moves_set.size() - 1;
            
            }

        }


        void remove_move(const int index) {

            int idx = this->improving_moves_indexes[index];

            if(idx != -1) {

                auto moved_h = this->improving_moves_set.back();
                this->improving_moves_indexes[moved_h.second] = idx;

                std::swap(this->improving_moves_set[idx], this->improving_moves_set.back());
                this->improving_moves_set.pop_back();
                this->improving_moves_indexes[index] = -1;
                                
            
            }

        }

        std::pair<int,int> select_improving_move(std::mt19937& gen) {
            
            std::uniform_int_distribution<> dis(0, this->improving_moves_set.size() - 1);
            int random_index = dis(gen);
            auto selected_move = this->improving_moves_set[random_index];
            
            return selected_move;

        }

};

#endif
