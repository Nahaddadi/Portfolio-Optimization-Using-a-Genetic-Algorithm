#include <vector>
#include <algorithm>
#include <numeric>
#include <random>
#include <cmath>

// Define a chromosome as a vector of doubles representing portfolio weights
struct Chromosome {
    std::vector<double> weights; // Portfolio weights
    double fitness; // Fitness score (Sharpe Ratio)

    Chromosome(int num_assets) : weights(num_assets), fitness(0.0) {}
};

// Define a population as a vector of chromosomes
using Population = std::vector<Chromosome>;

// Function to initialize a population with randomly generated chromosomes
Population initializePopulation(int num_assets, int population_size) {
    Population population;

    // Random number generation setup
    std::random_device rd;
    std::mt19937 gen(rd());
    std::uniform_real_distribution<> dis(0.0, 1.0);

    for (int i = 0; i < population_size; ++i) {
        Chromosome chromosome(num_assets);
        double sum = 0.0;

        // Generate random weights and normalize them
        for (double &weight : chromosome.weights) {
            weight = dis(gen);
            sum += weight;
        }
        for (double &weight : chromosome.weights) {
            weight /= sum;
        }

        population.push_back(chromosome);
    }

    return population;
}

// Function to calculate returns from portfolio weights and historical returns
double calculatePortfolioReturn(const Chromosome& chromosome, const std::vector<double>& asset_returns) {
    // Calculate portfolio return based on the weights and asset returns
    return std::inner_product(chromosome.weights.begin(), chromosome.weights.end(), asset_returns.begin(), 0.0);
}

// Function to calculate the actual variance of portfolio returns
double calculatePortfolioVariance(const Chromosome& chromosome, const std::vector<std::vector<double>>& stock_data) {
    size_t num_days = stock_data[0].size();
    std::vector<double> portfolio_returns(num_days);

    // Calculate daily return for the portfolio
    for (size_t day = 1; day < num_days; ++day) {
        double total_return = 0.0;
        for (size_t asset = 0; asset < stock_data.size(); ++asset) {
            total_return += chromosome.weights[asset] * (stock_data[asset][day] - stock_data[asset][day - 1]) /
                            stock_data[asset][day - 1];
        }
        portfolio_returns[day] = total_return;
    }

    // Return the standard deviation of the portfolio returns
    return stdDev(portfolio_returns);
}

// Fitness function to calculate the Sharpe Ratio for a chromosome
double calculateSharpeRatio(const Chromosome& chromosome, const std::vector<std::vector<double>>& stock_data, double risk_free_rate) {
    // Calculate portfolio return and variance
    std::vector<double> asset_returns = calculateReturns(stock_data);
    double portfolio_return = calculatePortfolioReturn(chromosome, asset_returns);
    double portfolio_variance = calculatePortfolioVariance(chromosome, stock_data);
    double portfolio_std_dev = std::sqrt(portfolio_variance);

    // Calculate and return the Sharpe Ratio
    return (portfolio_return - risk_free_rate) / portfolio_std_dev;
}

// Function to evaluate the fitness of each chromosome in the population
void evaluateFitness(Population& population, const std::vector<std::vector<double>>& stock_data, double risk_free_rate) {
    // Calculate fitness (Sharpe Ratio) for each chromosome
    for (Chromosome &chromosome : population) {
        chromosome.fitness = calculateSharpeRatio(chromosome, stock_data, risk_free_rate);
    }
}

// Selection: Sort the population based on fitness and select the best
void selection(Population& population) {
    // Sorting the population in descending order of fitness
    std::sort(population.begin(), population.end(), [](const Chromosome& a, Chromosome& b) {
        return a.fitness > b.fitness; // Descending order
    });
}

// Crossover: Combine two chromosomes to produce offspring
void crossover(const Chromosome& parent1,  Chromosome& parent2, Chromosome& offspring1, Chromosome& offspring2) {
    // Implementing a simple one-point crossover
    int size = parent1.weights.size();
    int crossover_point = std::rand() % size;

    for (int i = 0; i < size; ++i) {
        if (i < crossover_point) {
            offspring1.weights[i] = parent1.weights[i];
            offspring2.weights[i] = parent2.weights[i];
        } else {
            offspring1.weights[i] = parent2.weights[i];
            offspring2.weights[i] = parent1.weights[i];
        }
    }
}

// Mutation: Randomly alter the chromosome's weights
void mutation(Chromosome& chromosome, double mutation_rate) {
    // Random number generation setup
    std::random_device rd;
    std::mt19937 gen(rd());
    std::uniform_real_distribution<> dis(0.0, 1.0);

    // Mutation implementation
    for (double &weight : chromosome.weights) {
        if (dis(gen) < mutation_rate) {
            weight = dis(gen); // Assign a new weight randomly
        }
    }

    // Normalize the weights to ensure they sum up to 1
    double sum = std::accumulate(chromosome.weights.begin(), chromosome.weights.end(), 0.0);
        if (sum != 0.0) {
        for (double &weight : chromosome.weights) {
            weight /= sum;
        }
    } else {
        // Handle the case where the sum is zero after mutation
        for (double &weight : chromosome.weights) {
            weight = dis(gen); // Assign random weights
            sum += weight;
        }
        for (double &weight : chromosome.weights) {
            weight /= sum; // Normalize the new weights
        }
    }
}
