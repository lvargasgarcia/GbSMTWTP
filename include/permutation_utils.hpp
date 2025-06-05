#ifndef PERMUTATION_UTILS_HPP
#define PERMUTATION_UTILS_HPP

#include <vector>
#include <iostream>
#include <unordered_set>
#include <unordered_map>
#include <algorithm>
#include <string>
#include <cstdlib>
#include <ctime>
#include <algorithm>
#include <random>
using namespace std;


vector<int> compose(vector<int> a, vector<int> b) {
    vector<int> result(a.size());
    for (int i = 0; i < a.size(); i++) {
        result[i] = a[b[i]-1];
    }
    return result;
}

vector<int> inverse(vector<int> a) {
    vector<int> result(a.size());
    for (int i = 0; i < a.size(); i++) {
        result[a[i]-1] = i+1;
    }
    return result;
}

std::string to_string(const std::vector<int>& vec) {
    std::ostringstream oss;
    for (size_t i = 0; i < vec.size(); ++i) {
        oss << vec[i];
        if (i != vec.size() - 1) {
            oss << ",";
        }
    }
    return oss.str();
}

vector<int> to_vec(string str) {
    vector<int> result;
    std::istringstream iss(str);
    std::string token;
    while (std::getline(iss, token, ',')) {
        result.push_back(std::stoi(token));
    }
    return result;
}

string print_vec(vector<int> a) {
    string result = "[";
    for (int i = 0; i < a.size(); i++) {
        result += to_string(a[i]);
        if (i < a.size()-1) {
            result += ", ";
        }
    }
    result += "]";
    return result;
}

bool williamsCondition(vector<int> p, int n){
    int index = -1;
    int j = 0;
    while(index==-1){
        if(p[j] == n){
            index = j;
        }
        j++;
    }
    int r = p[(index % (n-1)) + 1];
    return r % (n-1) + 1 == p[0];
}

string williams_sequence(int n){
    
    string sequence = "";

    vector<int> tau;
    tau.push_back(2);
    tau.push_back(1);
    for(size_t i = 3; i <= n; i++){
        tau.push_back(i);
    }

    vector<int> sigma;
    for(size_t i = 1; i < n; i++){
        sigma.push_back(i+1);
    }
    sigma.push_back(1);

    vector<int> q;
    for(size_t i = 1; i <= n; i++){
        q.push_back(n-i+1);
    }

    auto qtau = compose(q, tau);
    auto qsigma = compose(q, sigma);
    auto qsigmatau = compose(qsigma, tau);
    auto invSigma = inverse(sigma);

    auto p = compose(qsigma, tau);

    while(p != qtau){
        if(p != qsigmatau){
            if(williamsCondition(p,n) && p != qsigma){
                p = compose(p, tau);
                sequence += "t";
            }else{
                p = compose(p, invSigma);
                sequence += "i";
            }
        }else{
            p = compose(p, invSigma);
            sequence += "i";
        }
    }
    return sequence;
}

// vector<vector<int>> permutations(int n){
    
//     vector<vector<int>> perms;

//     vector<int> tau;
//     tau.push_back(2);
//     tau.push_back(1);
//     for(size_t i = 3; i <= n; i++){
//         tau.push_back(i);
//     }

//     vector<int> sigma;
//     for(size_t i = 1; i < n; i++){
//         sigma.push_back(i+1);
//     }
//     sigma.push_back(1);

//     vector<int> q;
//     for(size_t i = 1; i <= n; i++){
//         q.push_back(n-i+1);
//     }

//     auto qtau = compose(q, tau);
//     auto qsigma = compose(q, sigma);
//     auto qsigmatau = compose(qsigma, tau);
//     auto invSigma = inverse(sigma);

//     auto p = compose(qsigma, tau);
//     perms.push_back(p);

//     while(p != qtau){
//         if(p != qsigmatau){
//             if(williamsCondition(p,n) && p != qsigma){
//                 p = compose(p, tau);
//             }else{
//                 p = compose(p, invSigma);
//             }
//         }else{
//             p = compose(p, invSigma);
//         }
//         perms.push_back(p);
//     }

//     std::sort(perms.begin(), perms.end(), [](const vector<int>& a, const vector<int>& b) {
//         return to_string(a) < to_string(b);
//     });

//     return perms;
// }

// std::vector<std::vector<int>> permutations(int n) {
//     std::vector<std::vector<int>> result;
//     std::vector<int> perm(n);
//     for (int i = 0; i < n; ++i) {
//         perm[i] = i; // Initialize with 1, 2, ..., n
//     }
    
//     do {
//         result.push_back(perm);
//     } while (std::next_permutation(perm.begin(), perm.end()));
    
//     return result;
// }

int to_int(vector<int> vec) {
    int number = 0;
    for (int i = 0; i < vec.size(); ++i) {
        number = number * 10 + vec[i];  // Shift current number left by one digit and add the next digit
    }
    return number;
}

vector<vector<int>> integer_partitions(int n) { 
    int p[n]; // An array to store a partition
    for(size_t i = 0; i < n; i++){
        p[i] = 0;
    }
    int k = 0; // Index of last element in a partition 
    p[k] = n; // Initialize first partition as number itself 
    vector<vector<int>> result;
    // This loop first prints current partition then generates next 
    // partition. The loop stops when the current partition has all 1s 
    while (true) 
    { 
        // print current partition 
        vector<int> partition;
        for (int i = 0; i <= k; i++){
            if(p[i] != 0){
                partition.push_back(p[i]);
            }
        }
        result.push_back(partition);
        // Generate next partition 
 
        // Find the rightmost non-one value in p[]. Also, update the 
        // rem_val so that we know how much value can be accommodated 
        int rem_val = 0; 
        while (k >= 0 && p[k] == 1) 
        { 
            rem_val += p[k]; 
            k--; 
        } 
 
        // if k < 0, all the values are 1 so there are no more partitions 
        if (k < 0) break; 
 
        // Decrease the p[k] found above and adjust the rem_val 
        p[k]--; 
        rem_val++; 
 
 
        // If rem_val is more, then the sorted order is violated. Divide 
        // rem_val in different values of size p[k] and copy these values at 
        // different positions after p[k] 
        while (rem_val > p[k]) 
        { 
            p[k+1] = p[k]; 
            rem_val = rem_val - p[k]; 
            k++; 
        } 
 
        // Copy rem_val to next position and increment position 
        p[k+1] = rem_val; 
        k++; 
    }
    return result;
}

std::vector<int> random_permutation(std::unordered_set<int> set) {
    // Convert the unordered_set to a vector
    std::vector<int> resp(set.begin(), set.end());

    // Seed the random number generator (only once per program execution)
    static bool seeded = false;
    if (!seeded) {
        std::srand(std::time(nullptr));
        seeded = true;
    }

    // Shuffle the vector randomly
    std::shuffle(resp.begin(), resp.end(), std::mt19937(std::random_device{}()));

    return resp;
}

std::unordered_map<std::string, int> permutations_set(int n, int pos, int k) {
    
    std::unordered_map<std::string, int> perms;

    // Generate all permutations of size k
    auto perms_k = permutations(k);

    perms_k.erase(perms_k.begin()); // Remove the identity permutation

    // Create move_set_beginning
    std::vector<int> move_set_beginning;
    for (int i = pos; i < pos + k; ++i) {
        move_set_beginning.push_back(i);
    }
    for (int i = 1; i < pos; ++i) {
        move_set_beginning.push_back(i);
    }
    for (int i = pos + k; i <= n; ++i) {
        move_set_beginning.push_back(i);
    }

    // Create move_set_end (inverse of move_set_beginning)
    std::vector<int> move_set_end = inverse(move_set_beginning);

    // Generate permutations
    for (const auto& tau : perms_k) {
        std::vector<int> tau_extended = tau;
        for (int i = k + 1; i <= n; ++i) {
            tau_extended.push_back(i);
        }

        std::vector<int> composed = compose(move_set_beginning, compose(tau_extended, move_set_end));
        perms[to_string(composed)] = 0;
    }

    return perms;
}

// Function to generate movements for every set of k consecutive numbers
std::unordered_map<int, std::unordered_map<std::string, int>> generate_movements(int n, int k) {
    std::unordered_map<int, std::unordered_map<std::string, int>> resp;

    for (int i = 1; i <= n - k + 1; ++i) {
        resp[i] = permutations_set(n, i, k);
    }

    return resp;
}

#endif