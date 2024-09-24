#include <iostream>
#include <vector>
#include <cmath>
#include <algorithm>
#include <string>
#include <cstdint>
#include <fstream>
#include <sstream>
#include <limits>
#include <unordered_set>
#include <chrono>
#include <iomanip>
#include "MurmurHash3.h"

// Función hash utilizando MurmurHash3
uint32_t hash_function(const std::string& str) {
    uint32_t hash;
    MurmurHash3_x86_32(str.c_str(), str.length(), 0, &hash);
    return hash;
}

class HyperLogLog {
private:
    std::vector<uint8_t> M;
    int m;  // Número de registros
    double alpha;

public:
    HyperLogLog(int b) : m(1 << b) {
        M.resize(m, 0);
        switch (m) {
            case 16: alpha = 0.673; break;
            case 32: alpha = 0.697; break;
            case 64: alpha = 0.709; break;
            default: alpha = 0.7213 / (1 + 1.079 / m); break;
        }
    }

    void add(const std::string& element) {
        uint32_t hash = hash_function(element);
        int j = hash >> (32 - log2(m));
        uint8_t w = hash & ((1 << (32 - log2(m))) - 1);
        M[j] = std::max(M[j], (uint8_t)(__builtin_clz(w) + 1));
    }

    double estimate() const {
        double sum = 0;
        for (int i = 0; i < m; i++) {
            sum += 1.0 / (1 << M[i]);
        }
        double E = alpha * m * m / sum;
        
        // Aplicar correcciones para valores pequeños y grandes
        if (E <= 2.5 * m) {
            int V = std::count(M.begin(), M.end(), 0);
            if (V != 0) {
                E = m * log(static_cast<double>(m) / V);
            }
        } else if (E > (1.0 / 30.0) * pow(2, 32)) {
            E = -pow(2, 32) * log(1 - E / pow(2, 32));
        }
        
        return E;
    }

    HyperLogLog merge(const HyperLogLog& other) const {
        if (m != other.m) {
            throw std::runtime_error("Cannot merge HyperLogLog with different number of registers");
        }
        HyperLogLog result(log2(m));
        for (int i = 0; i < m; i++) {
            result.M[i] = std::max(M[i], other.M[i]);
        }
        return result;
    }
};

double jaccard_similarity(const HyperLogLog& hll1, const HyperLogLog& hll2) {
    HyperLogLog union_hll = hll1.merge(hll2);
    double union_estimate = union_hll.estimate();
    double intersection_estimate = hll1.estimate() + hll2.estimate() - union_estimate;
    return intersection_estimate / union_estimate;
}

std::vector<std::string> generate_kmers(const std::string& sequence, int k) {
    std::vector<std::string> kmers;
    for (size_t i = 0; i <= sequence.length() - k; ++i) {
        kmers.push_back(sequence.substr(i, k));
    }
    return kmers;
}

std::vector<std::string> generate_minimizers(const std::string& sequence, int k, int w) {
    std::vector<std::string> minimizers;
    for (size_t i = 0; i <= sequence.length() - w; ++i) {
        std::string min_kmer = sequence.substr(i, k);
        for (size_t j = i + 1; j < i + w - k + 1; ++j) {
            std::string kmer = sequence.substr(j, k);
            if (kmer < min_kmer) {
                min_kmer = kmer;
            }
        }
        minimizers.push_back(min_kmer);
    }
    return minimizers;
}

std::string read_genome(const std::string& filename) {
    std::ifstream file(filename);
    if (!file.is_open()) {
        throw std::runtime_error("Unable to open file: " + filename);
    }
    std::stringstream buffer;
    buffer << file.rdbuf();
    return buffer.str();
}

double calculate_true_jaccard(const std::vector<std::string>& set1, const std::vector<std::string>& set2) {
    std::unordered_set<std::string> set1_unique(set1.begin(), set1.end());
    std::unordered_set<std::string> set2_unique(set2.begin(), set2.end());

    std::unordered_set<std::string> intersection;
    for (const auto& item : set1_unique) {
        if (set2_unique.count(item) > 0) {
            intersection.insert(item);
        }
    }

    std::unordered_set<std::string> union_set = set1_unique;
    union_set.insert(set2_unique.begin(), set2_unique.end());

    return static_cast<double>(intersection.size()) / union_set.size();
}

void evaluate_performance(const std::vector<std::string>& genomes, int k, int w, int b) {
    std::vector<double> true_jaccards, estimated_jaccards;
    std::vector<double> kmer_times, minimizer_times;

    for (size_t i = 0; i < genomes.size(); ++i) {
        for (size_t j = i + 1; j < genomes.size(); ++j) {
            std::cout << "Processing genomes " << i+1 << " and " << j+1 << std::endl;

            // K-mers approach
            auto start = std::chrono::high_resolution_clock::now();
            
            auto kmers1 = generate_kmers(genomes[i], k);
            auto kmers2 = generate_kmers(genomes[j], k);
            double true_jaccard = calculate_true_jaccard(kmers1, kmers2);

            HyperLogLog hll1(b), hll2(b);
            for (const auto& kmer : kmers1) hll1.add(kmer);
            for (const auto& kmer : kmers2) hll2.add(kmer);
            double estimated_jaccard = jaccard_similarity(hll1, hll2);

            auto end = std::chrono::high_resolution_clock::now();
            std::chrono::duration<double> diff = end - start;
            kmer_times.push_back(diff.count());

            true_jaccards.push_back(true_jaccard);
            estimated_jaccards.push_back(estimated_jaccard);

            std::cout << "  K-mers - True Jaccard: " << true_jaccard 
                      << ", Estimated: " << estimated_jaccard 
                      << ", Time: " << diff.count() << " seconds" << std::endl;

            // Minimizers approach
            start = std::chrono::high_resolution_clock::now();

            auto minimizers1 = generate_minimizers(genomes[i], k, w);
            auto minimizers2 = generate_minimizers(genomes[j], k, w);
            true_jaccard = calculate_true_jaccard(minimizers1, minimizers2);

            hll1 = HyperLogLog(b);
            hll2 = HyperLogLog(b);
            for (const auto& minimizer : minimizers1) hll1.add(minimizer);
            for (const auto& minimizer : minimizers2) hll2.add(minimizer);
            estimated_jaccard = jaccard_similarity(hll1, hll2);

            end = std::chrono::high_resolution_clock::now();
            diff = end - start;
            minimizer_times.push_back(diff.count());

            true_jaccards.push_back(true_jaccard);
            estimated_jaccards.push_back(estimated_jaccard);

            std::cout << "  Minimizers - True Jaccard: " << true_jaccard 
                      << ", Estimated: " << estimated_jaccard 
                      << ", Time: " << diff.count() << " seconds" << std::endl;
        }
    }

    // Calculate ERM and EAM
    double erm = 0, eam = 0;
    for (size_t i = 0; i < true_jaccards.size(); ++i) {
        erm += std::abs(estimated_jaccards[i] - true_jaccards[i]) / true_jaccards[i];
        eam += std::abs(estimated_jaccards[i] - true_jaccards[i]);
    }
    erm /= true_jaccards.size();
    eam /= true_jaccards.size();

    double avg_kmer_time = std::accumulate(kmer_times.begin(), kmer_times.end(), 0.0) / kmer_times.size();
    double avg_minimizer_time = std::accumulate(minimizer_times.begin(), minimizer_times.end(), 0.0) / minimizer_times.size();

    std::cout << std::fixed << std::setprecision(6);
    std::cout << "Error Relativo Medio (ERM): " << erm << std::endl;
    std::cout << "Error Absoluto Medio (EAM): " << eam << std::endl;
    std::cout << "Average K-mer processing time: " << avg_kmer_time << " seconds" << std::endl;
    std::cout << "Average Minimizer processing time: " << avg_minimizer_time << " seconds" << std::endl;
}

int main() {
    std::vector<std::string> genome_files = {"genome1.txt", "genome2.txt", "genome3.txt", "genome4.txt", "genome5.txt"};
    std::vector<std::string> genomes;

    std::cout << "Reading genome files..." << std::endl;
    for (const auto& file : genome_files) {
        try {
            genomes.push_back(read_genome(file));
            std::cout << "Successfully read: " << file << std::endl;
        } catch (const std::exception& e) {
            std::cerr << "Error reading file " << file << ": " << e.what() << std::endl;
        }
    }

    if (genomes.empty()) {
        std::cerr << "No genomes were successfully read. Exiting." << std::endl;
        return 1;
    }

    std::vector<int> k_values = {20, 25};
    std::vector<int> w_values = {50, 100};
    std::vector<int> b_values = {10, 14};  // 2^10 = 1024, 2^14 = 16384 registros

    for (int k : k_values) {
        for (int w : w_values) {
            for (int b : b_values) {
                std::cout << "\nEvaluating with k=" << k << ", w=" << w << ", b=" << b << std::endl;
                std::cout << "----------------------------------------" << std::endl;
                evaluate_performance(genomes, k, w, b);
            }
        }
    }

    return 0;
}