#include "solve_gridgraph.h"
#include <omp.h>
#include <filesystem>
#include <iostream>
#include <stdexcept>
#include <string>

void run_single_instance(int N, int run_id, const std::string& folder, double mean, double std) {
    try {
        // Construct output filename: ../folder/network_N_run_k.txt
        //std::string filename = folder + "/network_N" + std::to_string(N) + "_run" + std::to_string(run_id) + "_" + std::to_string(std)+ ".txt";

        std::ostringstream oss;
        oss << folder << "/network_N" << N << "_run" << run_id << "_" << std::fixed << std::setprecision(1) << std << ".txt";
        std::string filename = oss.str();
        
        // Build graph and solve max-flow
        double max_flow = buildAndSolveMaxFlow(N, run_id, filename, mean, std);
        
        std::cout << "Run " << run_id << ": Completed, Max Flow = " << max_flow << ", Output written to " << filename << std::endl;
    } catch (const std::exception& e) {
        std::cerr << "Run " << run_id << ": Error: " << e.what() << std::endl;
    }
}
int main(int argc, char* argv[]) {
if (argc != 6) {
        std::cerr << "Usage: " << argv[0] << " <N> <num_tries> <folder> <mean> <std>" << std::endl;
        return 1;
    }

    int N, num_tries;
    double mean, std;
    try {
        N = std::stoi(argv[1]);
        num_tries = std::stoi(argv[2]);
        mean = std::stod(argv[4]);
        std = std::stod(argv[5]);
    } catch (const std::exception& e) {
        std::cerr << "Error: Invalid N, num_tries, mean, or std. Must be valid numbers." << std::endl;
        return 1;
    }

    if (N <= 0 || num_tries <= 0) {
        std::cerr << "Error: N and num_tries must be positive." << std::endl;
        return 1;
    }
    if (std <= 0) {
        std::cerr << "Error: Standard deviation must be positive." << std::endl;
        return 1;
    }

    std::string folder = argv[3];

    // Create output directory if it doesn't exist
    try {
        std::filesystem::create_directories("../" + folder);
    } catch (const std::exception& e) {
        std::cerr << "Error: Unable to create directory ../" + folder << ": " << e.what() << std::endl;
        return 1;
    }

    // Parallel loop over num_tries
    #pragma omp parallel for schedule(dynamic)
    for (int run_id = 0; run_id < num_tries; ++run_id) {
        run_single_instance(N, run_id, folder, mean, std);
    }

    std::cout << "All runs completed." << std::endl;
    return 0;
}
