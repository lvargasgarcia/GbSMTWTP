#ifndef SMWTP_HPP
#define SMWTP_HPP

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <vector>
#include <fstream>
#include <sstream>
#include <iostream>
#include <stdexcept>
#include <map>
#include <bits/stl_numeric.h>
#include <bits/algorithmfwd.h>
#include "permutation_utils.hpp"
#include "Neighborhood.hpp"

class SMWTP {
    
    public:
        
        std::vector<int> times;
        std::vector<int> weights;
        std::vector<int> due;
        permutation_t optima;
        double globalMin;
        double globalMax;

        SMWTP(const std::string& file) {
            
            globalMin = -1;
            globalMax = -1;
    
            std::ifstream infile(file);
            if (!infile.is_open()) {
                throw std::runtime_error("Error: Unable to open file " + file);
            }
    
            std::string line;
            int n = 0;
    
            // Read the first line for the number of elements
            if (std::getline(infile, line)) {
                std::istringstream iss(line);
                iss >> n;
            } else {
                throw std::runtime_error("Error: File is empty or invalid format");
            }
    
            // Read the remaining lines for times, weights, and due dates
            while (std::getline(infile, line)) {
                std::istringstream iss(line);
                double time, weight, dueDate;
                if (iss >> time >> weight >> dueDate) {
                    times.push_back(time);
                    weights.push_back(weight);
                    due.push_back(dueDate);
                } else {
                    throw std::runtime_error("Error: Invalid line format in file");
                }
            }
    
            if (n != static_cast<int>(times.size())) {
                std::cerr << "Warning: Elements count does not match" << std::endl;
            }
        }

        double evaluate(const permutation_t& permutation) {
            
            int t = 0;
            double fit = 0.0;
    
            for (size_t i = 0; i < this->times.size(); ++i) {
                int job = permutation[i];
                t += times[job];
                if (t > due[job]) {
                    fit += weights[job] * (t - due[job]);
                }
            }
    
            return fit;
        }

        double initial_delta(const permutation_t& pi, const permutation_t& h, int position, int offset, int k, long common_tardiness) {
                
            double result = 0.0;
    
            int t_pi = common_tardiness;
            int t_hpi = common_tardiness;
    
            for (int i = position + offset; i < position + k; i++) {
                
                int job_pi = pi[i];
                int job_hpi = pi[h[i - position] + position];
    
                t_pi += times[job_pi];
                t_hpi += times[job_hpi];
    
                if (t_pi > due[job_pi]) {
                    result -= weights[job_pi] * (t_pi - due[job_pi]);
                }
    
                if (t_hpi > due[job_hpi]) {
                    result += weights[job_hpi] * (t_hpi - due[job_hpi]);
                }
            }
    
            return result;
        }

       // Delta function
       double delta_px(const permutation_t& pi, const permutation_t& selected_move, const permutation_t& h, int position, int k, long common_tardiness) {
                
        double result = 0.0;

        int t_pi = common_tardiness;
        int t_hpi = common_tardiness;

        for (int i = position; i < position + k; i++) {
            
            int job_pi = pi[selected_move[i]];
            int job_hpi = pi[selected_move[h[i]]];

            t_pi += times[job_pi];
            t_hpi += times[job_hpi];

            if (t_pi > due[job_pi]) {
                result -= weights[job_pi] * (t_pi - due[job_pi]);
            }

            if (t_hpi > due[job_hpi]) {
                result += weights[job_hpi] * (t_hpi - due[job_hpi]);
            }
        }

        return result;
    }
    
        // Delta function
        double delta(const permutation_t& pi, const permutation_t& selected_move, const permutation_t& h, int position_move, int position_h, int offset, int k, long common_tardiness) {
                
            double result = 0.0;
    
            int t_pi = common_tardiness;
            int t_hpi = common_tardiness;

            // int ct = 0;

            // for(int i = 0; i < position_move; i++){
            //     ct += times[pi[i]];
            // }

            // for(int i = position_h; i < position_move; i++){
            //     ct -= times[pi[i]];
            // }

            // for(int i = position_move; i < position_h; i++){
            //     ct += times[pi[selected_move[i - position_move] + position_move]];
            // }

            // t_pi = ct;
            // t_hpi = ct;
    
            for (int i = position_h + offset; i < position_h + k; i++) {
                
                int job_pi, job_hpi;
                
                if(i >= position_move && i < (position_move + k)){
                    job_pi =  pi[selected_move[i - position_move] + position_move];
                }else{
                    job_pi = pi[i];
                }

                if(h[i - position_h] + position_h >= position_move && h[i - position_h] + position_h < (position_move + k)){
                    job_hpi = pi[selected_move[h[i - position_h] + position_h - position_move] + position_move];
                }else{
                    job_hpi = pi[h[i - position_h] + position_h];
                }
    
                t_pi += times[job_pi];
                t_hpi += times[job_hpi];
    
                if (t_pi > due[job_pi]) {
                    result -= weights[job_pi] * (t_pi - due[job_pi]);
                }
    
                if (t_hpi > due[job_hpi]) {
                    result += weights[job_hpi] * (t_hpi - due[job_hpi]);
                }
            }
    
            return result;
        }
    
        // Get the number of jobs
        int getN() const {
            return times.size();
        }
    
        // Get the function values for all permutations
        std::map<std::string, double> getFunction() {
            
            int n = times.size();
            std::map<std::string, double> dict;

            vector<vector<int>> perms = permutations(n);
    
            for (const auto& perm : perms) {
                std::string key = to_string(perm);
                double value = evaluate(perm);
                dict[key] = value;

                if (globalMin == -1 || value < globalMin) {
                    globalMin = value;
                    optima = perm;
                }

                if (globalMax == -1 || value > globalMax) {
                    globalMax = value;
                }
            }

            return dict;
        }
   
        // 1. Heaviest Jobs First (Weight-Based Priority)
        permutation_t greedy_solution_weights() {
            permutation_t solution(times.size());
            for (int i = 0; i < times.size(); ++i)
                solution[i] = i;
    
            std::sort(solution.begin(), solution.end(), [this](int a, int b) {
                return weights[a] > weights[b]; // higher weight first
            });
    
            return solution;
        }
    
        // // 2. Shortest Processing Time First
        // permutation_t greedy_solution_processing_time() {
        //     permutation_t solution(times.size());
        //     for (int i = 0; i < times.size(); ++i)
        //         solution[i] = i;
    
        //     std::sort(solution.begin(), solution.end(), [this](int a, int b) {
        //         return times[a] < times[b]; // shorter processing time first
        //     });
    
        //     return solution;
        // }
    
        // 3. Weighted Earliest Due Date First (WEDD)
        permutation_t greedy_solution_due_date() {
            permutation_t solution(times.size());
            for (int i = 0; i < times.size(); ++i)
                solution[i] = i;

            std::sort(solution.begin(), solution.end(), [this](int a, int b) {
                return static_cast<double>(due[a]) / weights[a] < static_cast<double>(due[b]) / weights[b];
            });

            return solution;
        }
    
        // 4. WSPT: Weighted Shortest Processing Time First
        permutation_t greedy_solution_wspt() {
            permutation_t solution(times.size());
            for (int i = 0; i < times.size(); ++i)
                solution[i] = i;
    
            std::sort(solution.begin(), solution.end(), [this](int a, int b) {
                double ratio_a = static_cast<double>(weights[a]) / times[a];
                double ratio_b = static_cast<double>(weights[b]) / times[b];
                return ratio_a > ratio_b; // prioritize higher weight per time
            });
    
            return solution;
        }

        permutation_t greedy_solution_atc(double k = 4.0) {
            permutation_t solution;
            std::vector<bool> scheduled(times.size(), false);
            double avg_p = std::accumulate(times.begin(), times.end(), 0.0) / times.size();
            double current_time = 0.0;
        
            for (size_t step = 0; step < times.size(); ++step) {
                int best_idx = -1;
                double best_priority = -1.0;
        
                for (size_t i = 0; i < times.size(); ++i) {
                    if (scheduled[i]) continue;
        
                    double slack = std::max(0.0, static_cast<double>(due[i]) - current_time - times[i]);
                    double priority = (static_cast<double>(weights[i]) / times[i]) *
                                      std::exp(-slack / (k * avg_p));
        
                    if (priority > best_priority) {
                        best_priority = priority;
                        best_idx = i;
                    }
                }
        
                scheduled[best_idx] = true;
                solution.push_back(best_idx);
                current_time += times[best_idx];
            }
        
            return solution;
        }
            
};


#endif