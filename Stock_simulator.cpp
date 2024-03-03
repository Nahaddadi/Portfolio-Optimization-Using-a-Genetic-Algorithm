#include <iostream>
#include <vector>
#include <random>
#include <cmath>

// Function to simulate stock price data using Geometric Brownian Motion
std::vector<std::vector<double>> simulateStockData(int num_assets, int num_days, double mu, double sigma) {
    // Creating a 2D vector to store stock data for each asset over the days
    std::vector<std::vector<double>> stock_data(num_assets, std::vector<double>(num_days));

    // Setting up random number generation
    std::random_device rd;
    std::mt19937 gen(rd());
    std::normal_distribution<> d(0, 1);  // Standard normal distribution

    double dt = 1.0;  // Assuming each step is one time unit

    // Simulating the stock price data for each asset
    for (auto &asset_data : stock_data) {
        asset_data[0] = 40;  // Starting price of the stock
        for (size_t day = 1; day < num_days; ++day) {
            // Generate daily returns using GBM formula and update the price
            double daily_return = d(gen);
            asset_data[day] = asset_data[day - 1] * exp((mu - 0.5 * sigma * sigma) * dt + sigma * sqrt(dt) * daily_return);
        }
    }

    return stock_data;
}

// Function to print stock data
void printStockData(const std::vector<std::vector<double>>& stock_data) {
    // Iterating through each asset and day to print stock prices
    for (size_t asset = 0; asset < stock_data.size(); ++asset) {
        std::cout << "Asset " << asset + 1 << ": ";
        for (size_t day = 0; day < stock_data[asset].size(); ++day) {
            std::cout << stock_data[asset][day] << " ";
        }
        std::cout << std::endl;
    }
}