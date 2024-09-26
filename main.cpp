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
#include <numeric>
#include "MurmurHash3.h"

class HyperLogLog {
private:
    std::vector<uint8_t> M;
    int m;  // Número de registros
    double alpha;
    int b;  // Número de bits para indexación

public:
    HyperLogLog(int b) : m(1 << b), b(b) {
        M.resize(m, 0);
        switch (m) {
            case 16: alpha = 0.673; break;
            case 32: alpha = 0.697; break;
            case 64: alpha = 0.709; break;
            default: alpha = 0.7213 / (1 + 1.079 / m); break;
        }
    }

    void add(const std::string& element) {
        uint32_t hash[4];
        MurmurHash3_x86_128(element.c_str(), element.length(), 0, hash);
        uint32_t index = hash[0] >> (32 - b);
        uint8_t w = __builtin_clz((hash[0] << b) | (1 << (b - 1))) + 1;
        M[index] = std::max(M[index], w);
    }

    double estimate() const {
        double sum = 0;
        for (int i = 0; i < m; i++) {
            sum += std::pow(2.0, -M[i]);
        }
        double E = alpha * m * m / sum;
        
        if (E <= 2.5 * m) {
            int V = std::count(M.begin(), M.end(), 0);
            if (V != 0) {
                E = m * std::log(static_cast<double>(m) / V);
            }
        } else if (E > (1.0 / 30.0) * std::pow(2, 32)) {
            E = -std::pow(2, 32) * std::log(1 - E / std::pow(2, 32));
        }
        
        return E;
    }

    HyperLogLog merge(const HyperLogLog& other) const {
        if (m != other.m) {
            throw std::runtime_error("Cannot merge HyperLogLog with different number of registers");
        }
        HyperLogLog result(b);
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
    
    intersection_estimate = std::max(0.0, intersection_estimate);
    
    if (union_estimate == 0) return 0;
    
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
        for (size_t j = i + 1; j < i + w - k + 1 && j + k <= sequence.length(); ++j) {
            std::string kmer = sequence.substr(j, k);
            if (kmer < min_kmer) {
                min_kmer = kmer;
            }
        }
        minimizers.push_back(min_kmer);
    }
    return minimizers;
}

double jaccard_similarity_kmers(const std::string& seq1, const std::string& seq2, int k, int b) {
    HyperLogLog hll1(b), hll2(b);
    auto kmers1 = generate_kmers(seq1, k);
    auto kmers2 = generate_kmers(seq2, k);
    
    for (const auto& kmer : kmers1) hll1.add(kmer);
    for (const auto& kmer : kmers2) hll2.add(kmer);
    
    return jaccard_similarity(hll1, hll2);
}

double jaccard_similarity_minimizers(const std::string& seq1, const std::string& seq2, int k, int w, int b) {
    HyperLogLog hll1(b), hll2(b);
    auto minimizers1 = generate_minimizers(seq1, k, w);
    auto minimizers2 = generate_minimizers(seq2, k, w);
    
    for (const auto& minimizer : minimizers1) hll1.add(minimizer);
    for (const auto& minimizer : minimizers2) hll2.add(minimizer);
    
    return jaccard_similarity(hll1, hll2);
}

double calculate_true_jaccard(const std::vector<std::string>& set1, const std::vector<std::string>& set2) {
    std::unordered_set<std::string> set1_unique(set1.begin(), set1.end());
    std::unordered_set<std::string> set2_unique(set2.begin(), set2.end());

    size_t intersection_size = 0;
    for (const auto& item : set1_unique) {
        if (set2_unique.count(item) > 0) {
            intersection_size++;
        }
    }

    size_t union_size = set1_unique.size() + set2_unique.size() - intersection_size;

    if (union_size == 0) return 0;

    return static_cast<double>(intersection_size) / union_size;
}

struct EvaluationResults {
    double ERM;
    double EAM;
};

EvaluationResults calculate_errors(const std::vector<double>& true_values, const std::vector<double>& estimated_values) {
    double sum_relative_error = 0.0;
    double sum_absolute_error = 0.0;
    int valid_comparisons = 0;

    for (size_t i = 0; i < true_values.size(); ++i) {
        if (true_values[i] > 0) {
            sum_relative_error += std::abs(estimated_values[i] - true_values[i]) / true_values[i];
            valid_comparisons++;
        }
        sum_absolute_error += std::abs(estimated_values[i] - true_values[i]);
    }

    EvaluationResults results;
    results.ERM = valid_comparisons > 0 ? sum_relative_error / valid_comparisons : 0;
    results.EAM = sum_absolute_error / true_values.size();

    return results;
}

std::string read_first_genome(const std::string& filename) {
    std::ifstream file(filename);
    if (!file.is_open()) {
        throw std::runtime_error("Unable to open file: " + filename);
    }

    std::string line, genome;
    bool first_genome_found = false;

    while (std::getline(file, line)) {
        if (line[0] == '>') {
            if (first_genome_found) {
                // Si ya encontramos el primer genoma y vemos otro encabezado, terminamos
                break;
            }
            first_genome_found = true;
        } else if (first_genome_found) {
            genome += line;
        }
    }

    if (genome.empty()) {
        throw std::runtime_error("No genome found in file: " + filename);
    }

    return genome;
}

int main() {
    std::vector<std::string> genome_files = {
        "instances/GCF_000331305.1_ASM33130v1_genomic.fna",
        "instances/GCF_000373685.1_ASM37368v1_genomic.fna",
        "instances/GCF_000583735.1_ASM58373v1_genomic.fna",
        "instances/GCF_000716715.1_ASM71671v1_genomic.fna",
        "instances/GCF_000959725.1_ASM95972v1_genomic.fna"
    };
    std::vector<std::string> genomes;

    std::cout << "Reading genome files..." << std::endl;
    for (const auto& file : genome_files) {
        try {
            std::string genome = read_first_genome(file);
            genomes.push_back(genome);
            std::cout << "Read first genome from " << file << " (length: " << genome.length() << ")" << std::endl;
        } catch (const std::exception& e) {
            std::cerr << "Error reading file " << file << ": " << e.what() << std::endl;
        }
    }

    std::vector<int> k_values = {20, 25};
    std::vector<int> w_values = {50, 100};
    std::vector<int> b_values = {10, 14};  // 2^10 = 1024, 2^14 = 16384 registros

    for (int k : k_values) {
        for (int w : w_values) {
            for (int b : b_values) {
                std::cout << "\nEvaluating with k=" << k << ", w=" << w << ", b=" << b << std::endl;
                std::cout << "----------------------------------------" << std::endl;

                std::vector<double> true_jaccards, estimated_jaccards_kmers, estimated_jaccards_minimizers;
                std::vector<double> kmer_times, minimizer_times;

                for (size_t i = 0; i < genomes.size(); ++i) {
                    for (size_t j = i + 1; j < genomes.size(); ++j) {
                        std::cout << "Processing genomes " << i+1 << " and " << j+1 << std::endl;

                        // K-mers approach
                        auto start = std::chrono::high_resolution_clock::now();
                        double jaccard_kmers = jaccard_similarity_kmers(genomes[i], genomes[j], k, b);
                        auto end = std::chrono::high_resolution_clock::now();
                        std::chrono::duration<double> diff = end - start;
                        kmer_times.push_back(diff.count());

                        // Minimizers approach
                        start = std::chrono::high_resolution_clock::now();
                        double jaccard_minimizers = jaccard_similarity_minimizers(genomes[i], genomes[j], k, w, b);
                        end = std::chrono::high_resolution_clock::now();
                        diff = end - start;
                        minimizer_times.push_back(diff.count());

                        // Calculate true Jaccard (using k-mers for simplicity)
                        auto kmers1 = generate_kmers(genomes[i], k);
                        auto kmers2 = generate_kmers(genomes[j], k);
                        double true_jaccard = calculate_true_jaccard(kmers1, kmers2);

                        true_jaccards.push_back(true_jaccard);
                        estimated_jaccards_kmers.push_back(jaccard_kmers);
                        estimated_jaccards_minimizers.push_back(jaccard_minimizers);

                        std::cout << "  True Jaccard: " << true_jaccard << std::endl;
                        std::cout << "  K-mers Estimated: " << jaccard_kmers << ", Time: " << kmer_times.back() << " seconds" << std::endl;
                        std::cout << "  Minimizers Estimated: " << jaccard_minimizers << ", Time: " << minimizer_times.back() << " seconds" << std::endl;
                    }
                }

                // Calculate errors
                auto errors_kmers = calculate_errors(true_jaccards, estimated_jaccards_kmers);
                auto errors_minimizers = calculate_errors(true_jaccards, estimated_jaccards_minimizers);

                double avg_kmer_time = std::accumulate(kmer_times.begin(), kmer_times.end(), 0.0) / kmer_times.size();
                double avg_minimizer_time = std::accumulate(minimizer_times.begin(), minimizer_times.end(), 0.0) / minimizer_times.size();

                std::cout << std::fixed << std::setprecision(6);
                std::cout << "K-mers results:" << std::endl;
                std::cout << "  Error Relativo Medio (ERM): " << errors_kmers.ERM << std::endl;
                std::cout << "  Error Absoluto Medio (EAM): " << errors_kmers.EAM << std::endl;
                std::cout << "  Average processing time: " << avg_kmer_time << " seconds" << std::endl;

                std::cout << "Minimizers results:" << std::endl;
                std::cout << "  Error Relativo Medio (ERM): " << errors_minimizers.ERM << std::endl;
                std::cout << "  Error Absoluto Medio (EAM): " << errors_minimizers.EAM << std::endl;
                std::cout << "  Average processing time: " << avg_minimizer_time << " seconds" << std::endl;
            }
        }
    }

    return 0;
}