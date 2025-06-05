#ifndef LOCAL_SEARCH_HPP
#define LOCAL_SEARCH_HPP

#include "Neighborhood.hpp"
#include "smwtp.hpp"
#include <memory>
#include <algorithm>
#include <random>
#include <ctime>

permutation_t comp (permutation_t p, permutation_t move, int position) {
    auto result = permutation_t(p.size());
    for (size_t i = 0; i < position; i++) {
        result[i] = p[i];
    }
    for (size_t i = position; i < position + move.size(); i++) {
        result[i] = p[move[i - position] + position];
    }
    for (size_t i = position + move.size(); i < p.size(); i++) {
        result[i] = p[i];
    }
    return result;
}

std::vector<int> get_filtered_list(const int x, const int k, const int n) {

    int start = x - k + 1;
    int end = x + k - 1;

    int first = std::max(0, start);
    int last = std::min(n - k, end);

    int result_size = std::max(0, last - first + 1);

    auto result = std::vector<int>(result_size);

    for (int i = 0, val = first; val <= last; ++i, ++val) {
        result[i] = val;
    }

    return result;
}

std::string to_string(const permutation_t& perm, int n) {
    std::string result = "[";
    for (int i = 0; i < n; ++i) {
        result += std::to_string(perm[i]);
        if (i < n - 1) {
            result += ", ";
        }
    }
    result += "]";
    return result;
}

void compute_scores(const permutation_t& pi, Neighborhood& neighborhood, SMWTP& instance) {
    
    auto identity = permutation_t(pi.size());
    for (size_t i = 0; i < pi.size(); i++) {
        identity[i] = i;
    }

    int k = neighborhood.k;
    
    //#pragma omp parallel for num_threads(6)
    for (size_t i = 0; i < neighborhood.categories; i++) {
        for(size_t j = 0; j < neighborhood.perms_count; j++){
            
            auto current = i*neighborhood.perms_count + j;
            
            //Compute initial common tardiness
            
            long common_tardiness = 0;
            auto move = neighborhood.perms_set[j];
            int l = 0;

            int count = 0;

            while(count < k && move[count] == count){
                count++;
            }
            
            for(int l = 0; l < i + count; l++){
                common_tardiness += instance.times[pi[l]];
            }

            neighborhood.common_tardiness_set[current] = common_tardiness;

            auto change = instance.initial_delta(pi, move, i, count, k, common_tardiness);

            neighborhood.scores[current] = change;

            // auto change2 = -instance.evaluate(pi) + instance.evaluate(comp(pi, move, i));

            // if(change != change2){
            //     std::cout << "initial Error: " << change << " != " << change2 << "\n";
            //     change2 = -instance.evaluate(pi) + instance.evaluate(comp(pi, move, i));
            //     change = instance.initial_delta(pi, move, i, count, k, common_tardiness);
            // }

            if (change < 0) {
                //#pragma omp critical
                neighborhood.insert_move(i, current);
                
            }

        }
    }

}

void update_scores(permutation_t& pi, Neighborhood& neighborhood, SMWTP& instance, int selected_move_position,
                    int selected_move_index) {                              
    
        auto selected_move = neighborhood.perms_set[selected_move_index % neighborhood.perms_count];
        int k = neighborhood.k;

        int lower_bound = std::max(0, selected_move_position - neighborhood.k + 1)*neighborhood.perms_count;
        int upper_bound = (std::min(neighborhood.categories - 1, selected_move_position + neighborhood.k - 1) + 1)*neighborhood.perms_count;

        //#pragma omp parallel for num_threads(6)
        for (int idx = lower_bound; idx < upper_bound; ++idx) {
            
            int i = idx / neighborhood.perms_count;
            int j = idx % neighborhood.perms_count;
            auto move = neighborhood.perms_set[j];

            auto common_tardiness = neighborhood.common_tardiness_set[idx];

            int l = selected_move_position;

            int count = 0;
            while(count < k && move[count] == count){
                count++;
            }

            for (int l = selected_move_position; l < i + count && (l  - selected_move_position < k); ++l) {
                common_tardiness = common_tardiness - instance.times[pi[l]] + instance.times[pi[selected_move[l - selected_move_position] + selected_move_position]];
            }

            // for(int l = i; l < selected_move_position; ++l){
            //     common_tardiness -= instance.times[pi[l]];
            // }

            neighborhood.common_tardiness_set[idx] = common_tardiness;

            auto change = instance.delta(pi, selected_move, move, selected_move_position, i, count, k, common_tardiness);

            neighborhood.scores[idx] = change;

            // auto change2 = instance.evaluate(comp(comp(pi,selected_move,selected_move_position), move, i)) - instance.evaluate(comp(pi,selected_move,selected_move_position));

            // if(change != change2){
            //     std::cout << "Error: " << change << " != " << change2 << "\n";
                
            //     auto p_1 = comp(comp(pi,selected_move,selected_move_position), move, i);
            //     auto p_2 = comp(pi,selected_move,selected_move_position);
                
            //     change2 = instance.evaluate(p_1) - instance.evaluate(p_2);

            //     change = instance.delta(pi, selected_move, move, selected_move_position, i, std::max(0,l-i), k, common_tardiness);
            // }

            if(change < 0){
                //#pragma omp critical
                neighborhood.insert_move(i, idx);
            } else {
                //#pragma omp critical
                neighborhood.remove_move(idx);
            }

        }

        auto aux = std::make_unique<int[]>(k);

        for(int i = selected_move_position; i < selected_move_position + k; i++){
            aux[i - selected_move_position] = pi[selected_move[i - selected_move_position] + selected_move_position];
        }

        for(int i = selected_move_position; i < selected_move_position + k; i++){
            pi[i] = aux[i - selected_move_position];
        }

}

std::pair<permutation_t,long> local_search(Neighborhood neighborhood, SMWTP& instance, const permutation_t& pi, std::mt19937& device){

    int n = instance.getN();

    permutation_t p = permutation_t(n);
    std::copy(pi.begin(), pi.end(), p.begin());

    long current_fit = instance.evaluate(p);

    // std::cout << "Initial permutation: " << to_string(p, n) << "\n";

    // std::cout << "Initial fitness: " << current_fit << "\n";

    // std::cout << neighborhood.to_string() << "\n";

    compute_scores(p, neighborhood, instance);

    while (neighborhood.improving_moves_set.size() > 0) {

        //std::cout << neighborhood.improving_moves_set.size() << "\n";
        
        auto selected_move = neighborhood.select_improving_move(device);
        current_fit += neighborhood.scores[selected_move.second];

        //std::cout << "New fitness: " << current_fit << "\n";

        update_scores(p, neighborhood, instance, selected_move.first, selected_move.second);

        // std::cout << "New permutation: " << to_string(p, n) << "\n";   

    }

    return std::make_pair(p, current_fit);

}

#endif 

