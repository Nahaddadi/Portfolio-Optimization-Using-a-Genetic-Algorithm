#include <vector>
#include <numeric>
#include <cmath>

// Function to calculate mean of a dataset
double mean(const std::vector<double>& vec) {
    // Calculating the mean of the given vector
    return std::accumulate(vec.begin(), vec.end(), 0.0) / vec.size();
}

// Function to calculate standard deviation of a dataset
double stdDev(const std::vector<double>& vec) {
    // Calculating the standard deviation
    double mu = mean(vec);
    double sq_sum = std::inner_product(vec.begin(), vec.end(), vec.begin(), 0.0,
                                       [](double a, double b) { return a + b; },
                                       [mu](double a, double b) { return (a - mu) * (b - mu); });
    return std::sqrt(sq_sum / vec.size());
}

// Function to calculate covariance between two datasets
double covariance(const std::vector<double>& x, const std::vector<double>& y) {
    // Calculating covariance between two datasets
    double mean_x = mean(x);
    double mean_y = mean(y);
    double sum = 0.0;
    for (size_t i = 0; i < x.size(); ++i) {
        sum += (x[i] - mean_x) * (y[i] - mean_y);
    }
    return sum / (x.size() - 1);
}

// Function to calculate returns from historical stock price data
std::vector<double> calculateReturns(const std::vector<std::vector<double>>& stock_data) {
    // Vector to store average daily returns
    std::vector<double> returns(stock_data[0].size(), 0.0);

    // Calculating average returns for each day
    for (size_t day = 1; day < stock_data[0].size(); ++day) {
        double total_return = 0.0;
        for (const auto &asset_data : stock_data) {
            total_return += (asset_data[day] - asset_data[day - 1]) / asset_data[day - 1];
        }
        returns[day] = total_return / stock_data.size();
    }

    return returns;
}