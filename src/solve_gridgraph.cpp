// src/solve_gridgraph.cpp
#include "solve_gridgraph.h"
#include <lemon/list_graph.h>
#include <lemon/preflow.h>
#include <lemon/edmonds_karp.h>
#include <random>
#include <fstream>
#include <iostream>
#include <stdexcept>
#include <boost/math/distributions/normal.hpp>

using namespace lemon;
using boost::math::normal;

// Helper function to generate truncated normal random variable in [0, 1]
double truncated_normal(double mean, double std, std::mt19937& gen) {
    if (std <= 0) {
        throw std::invalid_argument("Standard deviation must be positive");
    }

    normal dist(0.0, 1.0); 
    std::uniform_real_distribution<> uniform(0.0, 1.0);

    double z0 = (0.0 - mean) / std; 
    double z1 = (1.0 - mean) / std; 
    double phi0 = cdf(dist, z0);    
    double phi1 = cdf(dist, z1);    
    if (phi1 <= phi0) {
        throw std::runtime_error("Invalid CDF bounds for truncated normal");
    }

    double u = uniform(gen);
    double p = phi0 + u * (phi1 - phi0); 
    double z = quantile(dist, p);        
    return mean + std * z;               
}

double buildAndSolveMaxFlow(int N, int run, const std::string& output_filename, double mean, double std) {

    bool periodicna = false;
    ListDigraph graph; // ustvari usmerjen grag iz knjiznice <lemon>

    ListDigraph::ArcMap<double> capacity(graph);
    float bound = 0.;

    std::random_device rd;
    std::mt19937 gen(rd() + run);  //razlicni seedi v primeru parallelnega pogona
    //std::mt19937 gen(42 + run);  //razlicni seedi v primeru parallelnega pogona
    //std::uniform_real_distribution<> dis(0.0, 1.0);

    int num_nodes = N * N;
    std::vector<ListDigraph::Node> nodes(num_nodes);
    for (int i = 0; i < num_nodes; ++i) {
        nodes[i] = graph.addNode();
    }

    auto idx = [=](int v, int s) {
        return v * N + s; // primer kako pretvorimo v(rstice) in s(tolpce) v enolicno vozlisce
    };

    for (int i = 0; i < N; ++i) {
        for (int j = 0; j < N; ++j) {
            int u_idx = idx(i, j);
            ListDigraph::Node u = nodes[u_idx];

            // povezave gor-dol
            if (i > 0) {
                if (i == 1) {
                    int v_idx = idx(i - 1, j);
                    //double cap = truncated_normal(mean, std, gen);
                    ListDigraph::Node v = nodes[v_idx];
                    //double cap = dis(gen);
                    ListDigraph::Arc a1 = graph.addArc(u, v);
                    //ListDigraph::Arc a2 = graph.addArc(v, u);
                    capacity[a1] = 1;
                    //capacity[a1] = capacity[a2] = 1.;
                }
                else if (i == N-1) {
                    int v_idx = idx(i - 1, j);
                    ListDigraph::Node v = nodes[v_idx];
                    //double cap = dis(gen);
                    ListDigraph::Arc a1 = graph.addArc(u, v);
                    //ListDigraph::Arc a2 = graph.addArc(v, u);
                    //double cap = truncated_normal(mean, std, gen);
                    capacity[a1] = 1;
                    //capacity [a1] = capacity[a2] = 1.;
                }
                else {
                    int v_idx = idx (i - 1, j);
                    ListDigraph::Node v = nodes[v_idx];
                    //double cap = dis(gen);
                    double cap = truncated_normal(mean, 1000, gen);
                    if (cap > bound) {
                        ListDigraph::Arc a1 = graph.addArc(u, v);
                        ListDigraph::Arc a2 = graph.addArc(v, u);
                        capacity[a1] = capacity[a2] = cap;
                    }
                }
            }
            
            // povezve levo-desno brez s_i in o_i
            if (j > 0 ) {
                if (i > 0 && i < N-1) {
                    int v_idx = idx(i, j - 1);
                    ListDigraph::Node v = nodes[v_idx];
                    //double cap = dis(gen);
                    double cap = truncated_normal(mean, 1000, gen);
                    if (cap > bound) {
                        ListDigraph::Arc a1 = graph.addArc(u, v);
                        ListDigraph::Arc a2 = graph.addArc(v, u);
                        capacity[a1] = capacity[a2] = cap;
                    }
                }

             }
            // periodicna povezav
            if (periodicna == true) {
                if (j == 0) {
                    if (i > 0 && i < N-1) {
                        int v_idx = idx(i, N-1);
                        ListDigraph::Node v = nodes[v_idx];
                        //double cap = dis(gen);
                        double cap = truncated_normal(mean, std, gen);
                        if (cap > bound) {
                            ListDigraph::Arc a1 = graph.addArc(u, v);
                            ListDigraph::Arc a2 = graph.addArc(v, u);
                            capacity[a1] = capacity[a2] = 1;
                        }
                    }
                }
           }
         }
    }

    double cap_inf = 1;

    ListDigraph::Node supersource = graph.addNode(); // S_idx = N^2  
    ListDigraph::Node supersink = graph.addNode(); // O_idx = N^2 + 1

    for (int j = 0; j < N; ++j) {
        ListDigraph::Node s_j = nodes[idx(N-1, j)];
        ListDigraph::Arc a = graph.addArc(supersource, s_j);

        capacity[a] =  cap_inf;
    }

    for (int j = 0; j < N; ++j) {
        ListDigraph::Node o_j = nodes[idx(0, j)];
        ListDigraph::Arc a = graph.addArc(o_j, supersink);
        capacity[a] =  cap_inf;
    }

    ListDigraph::ArcMap<double> capacity_copy(graph);
    for (ListDigraph::ArcIt a(graph); a != INVALID; ++a) {
        capacity_copy[a] = capacity[a];
    }

    //Run the Preflow algorithm on the temporary copy
    Preflow<ListDigraph, ListDigraph::ArcMap<double>> preflow_solver(graph, capacity_copy, supersource, supersink);
    preflow_solver.run();

    //Get the max-flow value
    double max_flow = preflow_solver.flowValue();
    ListDigraph::NodeMap<bool> min_cut_map(graph);
    preflow_solver.minCutMap(min_cut_map);
    //
    //EdmondsKarp<ListDigraph, ListDigraph::ArcMap<double>> ek_solver(graph, capacity, supersource, supersink);
    //ek_solver.run();

    //ListDigraph::NodeMap<bool> min_cut_map(graph);
    //ek_solver.minCutMap(min_cut_map);
    //double max_flow = ek_solver.flowValue();

    int min_cut_edges_count = 0;
    for (ListDigraph::ArcIt a(graph); a != INVALID; ++a) {
         ListDigraph::Node u = graph.source(a);
        ListDigraph::Node v = graph.target(a);

        if (min_cut_map[u] && !min_cut_map[v]) {
            min_cut_edges_count++;
          }
    }

    // Write results to file in CSV format
    std::ofstream out_file(output_filename);
    if (!out_file) {
        throw std::runtime_error("Unable to open output file: " + output_filename);
    }
    out_file << "Total Max Flow: " << max_flow << "\n\n";
    out_file << "Total Min Cut Edges: " << min_cut_edges_count << "\n\n";
    out_file << "SourceID,TargetID,FlowValue,Capacity\n";
    for (ListDigraph::ArcIt a(graph); a != INVALID; ++a) {
        int u_id = graph.id(graph.source(a));
        int v_id = graph.id(graph.target(a));
        double flow_val = preflow_solver.flow(a);
        //double flow_val = ek_solver.flow(a);
        double cap = capacity[a];
        out_file << u_id << "," << v_id << "," << flow_val << "," << cap << "\n";
    }
    out_file.close();

    return max_flow;
}
    
