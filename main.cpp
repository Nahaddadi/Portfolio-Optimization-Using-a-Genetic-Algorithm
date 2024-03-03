#include <iostream>
#include <vector>
#include <random>
#include <algorithm>
#include <numeric>
#include <cmath>

// Include the contents of Statistics.cpp directly
#include "Statistics.cpp"

// Include the contents of Stock_simulator.cpp directly
#include "Stock_simulator.cpp"

// Include the contents of Genetic_algorithm.cpp directly
#include "Genetic_algorithm.cpp"

int main() {
    // User input for parameters
    int num_assets, population_size, generations;
    double risk_free_rate, mutation_rate, expected_return, volatility;

    std::cout << "Enter the number of assets in the portfolio: ";
    std::cin >> num_assets;

    std::cout << "Enter the size of the population: ";
    std::cin >> population_size;

    std::cout << "Enter the number of generations for the genetic algorithm: ";
    std::cin >> generations;

    std::cout << "Enter the annual risk-free rate for Sharpe Ratio calculation (as a decimal): ";
    std::cin >> risk_free_rate;

    std::cout << "Enter the mutation rate for the genetic algorithm (as a decimal): ";
    std::cin >> mutation_rate;

    std::cout << "Enter the expected return for the simulated stock data (as a decimal): ";
    std::cin >> expected_return;

    std::cout << "Enter the volatility for the simulated stock data (as a decimal): ";
    std::cin >> volatility;

    // Simulate some historical stock price data
    std::vector<std::vector<double>> stock_data = simulateStockData(num_assets, 100, expected_return, volatility);

    // Print the simulated stock data
    std::cout << "Simulated Stock Data:" << std::endl;
    printStockData(stock_data);

    // Initialize a population of chromosomes
    Population population = initializePopulation(num_assets, population_size);

    // Run the genetic algorithm for the specified number of generations
    for (int i = 0; i < generations; ++i) {
        // Evaluate fitness of the population
        evaluateFitness(population, stock_data, risk_free_rate);

        // Select the best chromosomes
        selection(population);

        // Generate a new population through crossover and mutation
        Population new_population;
        while (new_population.size() < population_size) {
            for (int j = 0; j < population_size / 2; ++j) {
                Chromosome offspring1(num_assets), offspring2(num_assets);
                int parent1_index = j;
                int parent2_index = (j + 1) % (population_size / 2);

                crossover(population[parent1_index], population[parent2_index], offspring1, offspring2);
                mutation(offspring1, mutation_rate);
                mutation(offspring2, mutation_rate);

                evaluateFitness(new_population, stock_data, risk_free_rate);

                new_population.push_back(offspring1);
                new_population.push_back(offspring2);
            }
        }

        population = new_population; // Replace the old population with the new one
    }

    // After the final generation, select the best chromosome as the optimized portfolio
    selection(population);
    Chromosome best_portfolio = population.front();

    // Display the weights of the optimized portfolio
    std::cout << "Optimized portfolio weights: ";
    for (double weight : best_portfolio.weights) {
        std::cout << weight << " ";
    }
    std::cout << std::endl;
    
    return 0;
}