#ifndef DRILS_HPP
#define DRILS_HPP

#include "Neighborhood.hpp"
#include "smwtp.hpp"
#include <memory>
#include <algorithm>
#include <random>
#include <ctime>
#include "local_search.hpp"
#include "Neighborhood.hpp"
#include <chrono>


permutation_t partition_crossover(const permutation_t& sigma_1, const permutation_t& sigma_2, SMWTP& instance) {

    int n = sigma_1.size();
    std::vector<int> child = sigma_1;
    std::vector<int> pi = compose(inverse(sigma_1, n), sigma_2, n);

    int i = 0;
    long common_tardiness = 0;

    permutation_t identity = permutation_t(n);

    for (int j = 0; j < n; ++j) {
        identity[j] = j;
    }

    while (i < n) {

        permutation_t h(n, 0);

        for (int k = 0; k < n; ++k) {
            h[k] = k;
        }

        int l = i;
        int initial_pos = l;

        long aux_common_tardiness = instance.times[sigma_1[l]];

        h[l] = pi[l];
        int j = h[l];

        while (l < j) {
            l += 1;
            aux_common_tardiness += instance.times[sigma_1[l]];
            h[l] = pi[l];
            j = std::max(j, h[l]);
        }

        int final_pos = l;
        int incr = final_pos - initial_pos;

        double delta_1 = instance.delta_px(identity, sigma_1, h, initial_pos, incr + 1, common_tardiness);
        common_tardiness += aux_common_tardiness;

        // double delta_2 = instance.evaluate(compose(sigma_1, h, n)) - instance.evaluate(sigma_1);

        // if(delta_1 != delta_2){
        //     std::cout << "Error: " << delta_1 << " != " << delta_2 << "\n";
        //     delta_2 = instance.evaluate(compose(sigma_1, h, n)) - instance.evaluate(sigma_1);
        //     delta_1 = instance.delta_px(identity, sigma_1, h, initial_pos, incr, common_tardiness);
        // }

        if (delta_1 < 0) {

            permutation_t aux(incr + 1);

            for(int k = initial_pos; k <= final_pos; k++){
                aux[k - initial_pos] = sigma_1[h[k]];
            }

            for(int k = initial_pos; k <= final_pos; k++){
                child[k] = aux[k - initial_pos];
            }

        }

        i = l + 1;
    }

    return child;
}

permutation_t random_permutation(const permutation_t& elems, std::mt19937& gen) {
    permutation_t vec(elems.begin(), elems.end());
    std::shuffle(vec.begin(), vec.end(), gen);
    return vec;
}

permutation_t perturbation_function(const permutation_t& p, std::mt19937& gen, int initial_position, int final_position) {

    int n = p.size();
    std::uniform_int_distribution<> dist;

    int r = std::uniform_int_distribution<>(1, final_position - initial_position)(gen);
    int total_elems = r + initial_position;

    permutation_t current_permutation(n);
    current_permutation.reserve(n);

    for(int i = 0; i < initial_position; ++i){
        current_permutation[i] = p[i];
    }
    
    permutation_t first_chunk(p.begin() + initial_position, p.begin() + total_elems);
    permutation_t new_perm = random_permutation(first_chunk, gen);

    for(int i = initial_position; i < total_elems; ++i){
        current_permutation[i] = new_perm[i - initial_position];
    }

    while (total_elems < final_position) {
        
        int remaining = final_position - total_elems;
        r = std::uniform_int_distribution<>(1, remaining)(gen);

        permutation_t next_chunk(p.begin() + total_elems, p.begin() + total_elems + r);
        permutation_t next_permutation = random_permutation(next_chunk, gen);

        for(int i = total_elems; i < total_elems + r; ++i){
            current_permutation[i] = next_permutation[i - total_elems];
        }

        total_elems += r;
    }

    for(int i = total_elems; i < n; ++i){
        current_permutation[i] = p[i];
    }

    // for (int i = 0; i < n; i++){
    //     std::cout << current_permutation[i] << " ";
    // }

    // std::cout << "\n";

    // auto dist = std::uniform_int_distribution<>(0, p.size() - 1);

    // int position = dist(gen);

    // permutation_t perturb(p.begin() + position, position + 150 >= p.size() ? p.end() : p.begin() + position + 150);

    // permutation_t random_perm = random_permutation(perturb, gen);

    // permutation_t current_permutation(p.size());

    // for (int i = 0; i < position; ++i) {
    //     current_permutation[i] = p[i];
    // }
    // for (int i = position; i < position + random_perm.size(); ++i) {
    //     current_permutation[i] = random_perm[i - position];
    // }
    // for (int i = position + random_perm.size(); i < p.size(); ++i) {
    //     current_permutation[i] = p[i];
    // }

    // int n = p.size();
    // std::uniform_int_distribution<> dist;
    // permutation_t current_permutation(n);
    // current_permutation.reserve(n);

    // for(int i = 0; i < initial_position; ++i){
    //     current_permutation[i] = p[i];
    // }

    // std::shuffle(current_permutation.begin() + initial_position, current_permutation.end(), gen);

    return current_permutation;
}

permutation_t perturb(const permutation_t& p, std::mt19937& gen, int initial_position, int final_position, double perturbation_factor) {

    int n = p.size();

    int elements_to_perturb = (int)(n*perturbation_factor);

    // std::cout << elements_to_perturb << "\n";

    std::uniform_int_distribution<> dist(0, final_position-2);

    permutation_t resp(n);

    for(size_t i = 0; i < n; i++){
        resp[i] = p[i];
    }

    
    for(size_t i = 0; i < elements_to_perturb; i++){
        int k = dist(gen);
        int l = dist(gen);

        // if(k % 2 == 0){
        //     k += (int)(final_position/2);
        //     l += (int)(final_position/2);
        // }

        std::swap(resp[k], resp[l]);

    }

    return resp;

}


// std::vector<std::pair<time_t,long>> DRILS(SMWTP& instance, int neighborhood_size, std::mt19937& device, long time_to_run) {

//     std::vector<std::pair<time_t, long>> resp;
    
//     long t_0 = std::time(nullptr);
//     //auto pi = partition_crossover(partition_crossover(instance.greedy_solution_wspt(), instance.greedy_solution_due_date(), instance), instance.greedy_solution_weights(), instance);
//     auto n = instance.getN();

//     auto identity = permutation_t(n);

//     for(int i = 0; i < n; i++){
//         identity[i] = i;
//     }

//     //seeding

//     auto pi = random_permutation(identity, device);

//     // int counter = 400;
//     auto neighborhood = Neighborhood(n, neighborhood_size);
//     //auto bigNeighborhood = Neighborhood(n, neighborhood_size + 1);
//     // int old_quotient = 4;

//     std::vector<permutation_t> population = {
//         local_search(neighborhood, instance, pi, device).first,
//         local_search(neighborhood, instance, instance.greedy_solution_atc(), device).first,
//         local_search(neighborhood, instance, instance.greedy_solution_wspt(), device).first,
//         local_search(neighborhood, instance, instance.greedy_solution_due_date(), device).first
//     };

//     auto curr = partition_crossover(population[0], population[1], instance);

//     for(int i = 2; i < population.size(); i++){
//         curr = partition_crossover(curr, population[i], instance);
//     }

//     auto current = local_search(neighborhood, instance, curr, device);

//     auto best = current.first;
//     long best_value = current.second;

//     resp.push_back(std::make_pair(std::time(nullptr) - t_0, best_value));

//     std::cerr << "Initial value: " << best_value << "\n";
//     int k = 0;

//     while(std::time(nullptr) - t_0 < time_to_run) {
        
//         auto next = local_search(neighborhood, instance, perturbation_function(current.first, device, ((k % 4) + 4)*(n/8), n), device);
//         k++;

//         auto child = partition_crossover(current.first, next.first, instance);


//         if(child != current.first && child != next.first){
//             std::cerr << "perturbation pxSuccess" << "\n";
//             current = local_search(neighborhood, instance, child, device);
//         }else if(child == next.first){
//             // if(current.second < instance.evaluate(child)){
//             //     std::cout << "Error" << "\n";
//             //     partition_crossover(current.first, next.first, instance);
//             //     return;
//             // }
//             current = next;
//         }

//         // auto next2 = local_search(neighborhood, instance, perturbation_function(current.first, device, 0, ((k%4) + 5)*(n/8)), device);

//         // child = partition_crossover(current.first, next2.first, instance);

//         // if(child != current.first && child != next2.first){
//         //     std::cerr << "perturbation pxSuccess" << "\n";
//         //     current = local_search(neighborhood, instance, child, device);
//         // }else if(child == next2.first){
//         //     // if(current.second < instance.evaluate(child)){
//         //     //     std::cout << "Error" << "\n";
//         //     //     partition_crossover(current.first, next2.first, instance);
//         //     //     return;
//         //     // }
//         //     current = next2;
//         // }

//         // for(int i = 0; i < 125; i++){
            
//         //     int pos =  (i*(n/125));
//         //     std::cerr << "Pos: " << pos << "\n";
            
//         //     child = partition_crossover(current.first, perturbation_function(current.first, device, pos, n), instance);
//         //     if(child != current.first){
//         //         // if(current.second < instance.evaluate(child)){
//         //         //     std::cout << "Error" << "\n";
//         //         //     partition_crossover(current.first, child, instance);
//         //         //     return;
//         //         // }
//         //         current = local_search(neighborhood, instance, child, device);
//         //         std::cerr << "pxSuccess" << "\n";
//         //     }

//         // }

//         // if(std::time(nullptr) - t_0 > 120){
//         auto next3 = local_search(neighborhood, instance, perturb(current.first, device, 0, n), device);

//         child = partition_crossover(current.first, next3.first, instance);

//         if(child != current.first && child != next3.first){
//             std::cerr << "perturbation pxSuccess" << "\n";
//             current = local_search(neighborhood, instance, child, device);
//         }else if(child == next3.first){
//             // if(current.second < instance.evaluate(child)){
//             //     std::cout << "Error" << "\n";
//             //     partition_crossover(current.first, next3.first, instance);
//             //     return;
//             // }
//             current = next3;
//         }
//         // }


//         if(current.second < best_value){
//             best = current.first;
//             best_value = current.second;
//             resp.push_back(std::make_pair(std::time(nullptr) - t_0, best_value));
//         }

//         std::cerr << "Current value: " << current.second << "\n";
//         std::cerr << "Best value: " << best_value << "\n";

//         // std::cerr << "Counter: " << counter << "\n";

//         // counter++;
//         // if(counter / 100 != old_quotient){
//         //     neighborhood = Neighborhood(n, counter / 100);
//         // }
//         // old_quotient = counter / 100;

//     }

//     return resp;

// }

std::vector<std::tuple<long,long, double>> basic_DRILS(SMWTP& instance, int neighborhood_size, std::mt19937& device, long time_to_run, double perturbation_factor, int tolerance) {

    auto t_0 = std::chrono::steady_clock::now();
    
    std::vector<std::tuple<long, long, double>> resp; // Use long for time in seconds
    
    int iters = 0;
    double px_counter = 0;

    auto n = instance.getN();
    auto identity = permutation_t(n);

    for(int i = 0; i < n; i++){
        identity[i] = i;
    }

    auto pi = random_permutation(identity, device);
    auto neighborhood = Neighborhood(n, neighborhood_size);

    auto current = local_search(neighborhood, instance, pi, device);

    auto best = current.first;
    long best_value = current.second;

    resp.push_back(std::make_tuple(
        std::chrono::duration_cast<std::chrono::milliseconds>(std::chrono::steady_clock::now() - t_0).count(),
        best_value,
        0
    ));

    int chunk_size = n/10;
    int chunk = 0;

    //std::cerr << "Initial value: " << best_value << "\n";
    int k = 0;
    // bool changed = false;

    while(std::chrono::duration_cast<std::chrono::seconds>(std::chrono::steady_clock::now() - t_0).count() < time_to_run) {
        
        std::pair<permutation_t, long> next;
        iters++;

        next = local_search(neighborhood, instance, perturb(current.first, device, 0, n, perturbation_factor), device);

        auto child = partition_crossover(current.first, next.first, instance);

        if(child == current.first){
            k++;
            if(k > tolerance){
                current = next;
                k = 0;
            }
        }else if(child == next.first){
            current = next;
            k = 0;
        }else{
            current = local_search(neighborhood, instance, child, device);
            px_counter++;
            k = 0;
        }

        if(current.second < best_value){
            best = current.first;
            best_value = current.second;
            long time = std::chrono::duration_cast<std::chrono::milliseconds>(std::chrono::steady_clock::now() - t_0).count();
            resp.push_back(std::make_tuple(
                time,
                best_value,
                (px_counter/iters)*100
            ));

            // if(!changed && time > 300000){
            //     neighborhood = Neighborhood(n, 4);
            //     changed = true;
            //     tolerance = 10000;
            // }

            // std::cerr << "Current value: " << current.second << "\n";
            // std::cerr << "Best value: " << best_value << "\n";
        }

        // std::cerr << "Current value: " << current.second << "\n";
        // std::cerr << "Best value: " << best_value << "\n";
    }
    // std:cerr << iters << "\n";

    return resp;
}

#endif
