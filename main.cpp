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
#include <cstring>
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

struct ExperimentResult {
    int k, w, b;
    double true_jaccard;
    double estimated_jaccard_kmers;
    double estimated_jaccard_minimizers;
    double time_kmers;
    double time_minimizers;
    double erm_kmers;
    double eam_kmers;
    double erm_minimizers;
    double eam_minimizers;
};

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

void write_csv_header(std::ofstream& file) {
    file << "k,w,b,true_jaccard,estimated_jaccard_kmers,estimated_jaccard_minimizers,"
         << "time_kmers,time_minimizers,erm_kmers,eam_kmers,erm_minimizers,eam_minimizers\n";
}

void write_csv_row(std::ofstream& file, const ExperimentResult& result) {
    file << result.k << "," << result.w << "," << result.b << ","
         << result.true_jaccard << "," << result.estimated_jaccard_kmers << ","
         << result.estimated_jaccard_minimizers << "," << result.time_kmers << ","
         << result.time_minimizers << "," << result.erm_kmers << ","
         << result.eam_kmers << "," << result.erm_minimizers << ","
         << result.eam_minimizers << "\n";
}

int main(int argc, char* argv[]) {
    bool generate_log = false;
    for (int i = 1; i < argc; ++i) {
        if (strcmp(argv[i], "--log") == 0) {
            generate_log = true;
            break;
        }
    }

    std::ofstream log_file;
    if (generate_log) {
        log_file.open("experiment_results.csv");
        if (!log_file.is_open()) {
            std::cerr << "Error: Unable to create log file." << std::endl;
            return 1;
        }
        write_csv_header(log_file);
    }

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

    std::vector<ExperimentResult> all_results;

    for (int k : k_values) {
        for (int w : w_values) {
            for (int b : b_values) {
                std::cout << "\nEvaluating with k=" << k << ", w=" << w << ", b=" << b << std::endl;
                std::cout << "----------------------------------------" << std::endl;

                for (size_t i = 0; i < genomes.size(); ++i) {
                    for (size_t j = i + 1; j < genomes.size(); ++j) {
                        std::cout << "Processing genomes " << i+1 << " and " << j+1 << std::endl;

                        ExperimentResult result;
                        result.k = k;
                        result.w = w;
                        result.b = b;

                        // K-mers approach
                        auto start = std::chrono::high_resolution_clock::now();
                        result.estimated_jaccard_kmers = jaccard_similarity_kmers(genomes[i], genomes[j], k, b);
                        auto end = std::chrono::high_resolution_clock::now();
                        result.time_kmers = std::chrono::duration<double>(end - start).count();

                        // Minimizers approach
                        start = std::chrono::high_resolution_clock::now();
                        result.estimated_jaccard_minimizers = jaccard_similarity_minimizers(genomes[i], genomes[j], k, w, b);
                        end = std::chrono::high_resolution_clock::now();
                        result.time_minimizers = std::chrono::duration<double>(end - start).count();

                        // Calculate true Jaccard
                        auto kmers1 = generate_kmers(genomes[i], k);
                        auto kmers2 = generate_kmers(genomes[j], k);
                        result.true_jaccard = calculate_true_jaccard(kmers1, kmers2);

                        // Calculate ERM and EAM
                        if (result.true_jaccard == 0) {
                            // Si true_jaccard es 0, usamos el valor estimado como error
                            result.erm_kmers = result.estimated_jaccard_kmers;
                            result.erm_minimizers = result.estimated_jaccard_minimizers;
                        } else {
                            result.erm_kmers = std::abs(result.estimated_jaccard_kmers - result.true_jaccard) / result.true_jaccard;
                            result.erm_minimizers = std::abs(result.estimated_jaccard_minimizers - result.true_jaccard) / result.true_jaccard;
                        }
                        result.eam_kmers = std::abs(result.estimated_jaccard_kmers - result.true_jaccard);
                        result.eam_minimizers = std::abs(result.estimated_jaccard_minimizers - result.true_jaccard);

                        all_results.push_back(result);

                        if (generate_log) {
                            write_csv_row(log_file, result);
                        }

                        std::cout << "  True Jaccard: " << result.true_jaccard << std::endl;
                        std::cout << "  K-mers Estimated: " << result.estimated_jaccard_kmers 
                                  << ", Time: " << result.time_kmers << " seconds" << std::endl;
                        std::cout << "  Minimizers Estimated: " << result.estimated_jaccard_minimizers 
                                  << ", Time: " << result.time_minimizers << " seconds" << std::endl;
                    }
                }

                // Calculate and display ERM and EAM for this configuration
                double erm_kmers = 0, eam_kmers = 0, erm_minimizers = 0, eam_minimizers = 0;
                int valid_comparisons = 0;
                for (const auto& result : all_results) {
                    if (result.k == k && result.w == w && result.b == b) {
                        erm_kmers += result.erm_kmers;
                        eam_kmers += result.eam_kmers;
                        erm_minimizers += result.erm_minimizers;
                        eam_minimizers += result.eam_minimizers;
                        valid_comparisons++;
                    }
                }

                if (valid_comparisons > 0) {
                    erm_kmers /= valid_comparisons;
                    eam_kmers /= valid_comparisons;
                    erm_minimizers /= valid_comparisons;
                    eam_minimizers /= valid_comparisons;
                }

                std::cout << std::fixed << std::setprecision(6);
                std::cout << "K-mers results:" << std::endl;
                std::cout << "  Error Relativo Medio (ERM): " << erm_kmers << std::endl;
                std::cout << "  Error Absoluto Medio (EAM): " << eam_kmers << std::endl;

                std::cout << "Minimizers results:" << std::endl;
                std::cout << "  Error Relativo Medio (ERM): " << erm_minimizers << std::endl;
                std::cout << "  Error Absoluto Medio (EAM): " << eam_minimizers << std::endl;
            }
        }
    }

    if (generate_log) {
        // Calcular y escribir estadísticas globales al final del archivo CSV
        double avg_true_jaccard = 0, min_true_jaccard = std::numeric_limits<double>::max(), max_true_jaccard = 0;
        double avg_erm_kmers = 0, avg_eam_kmers = 0, avg_erm_minimizers = 0, avg_eam_minimizers = 0;
        double avg_time_kmers = 0, avg_time_minimizers = 0;
        
        for (const auto& result : all_results) {
            avg_true_jaccard += result.true_jaccard;
            min_true_jaccard = std::min(min_true_jaccard, result.true_jaccard);
            max_true_jaccard = std::max(max_true_jaccard, result.true_jaccard);
            avg_erm_kmers += result.erm_kmers;
            avg_eam_kmers += result.eam_kmers;
            avg_erm_minimizers += result.erm_minimizers;
            avg_eam_minimizers += result.eam_minimizers;
            avg_time_kmers += result.time_kmers;
            avg_time_minimizers += result.time_minimizers;
        }
        
        int n = all_results.size();
        avg_true_jaccard /= n;
        avg_erm_kmers /= n;
        avg_eam_kmers /= n;
        avg_erm_minimizers /= n;
        avg_eam_minimizers /= n;
        avg_time_kmers /= n;
        avg_time_minimizers /= n;
        
        log_file << "\nGlobal Statistics\n";
        log_file << "Average True Jaccard," << avg_true_jaccard << "\n";
        log_file << "Min True Jaccard," << min_true_jaccard << "\n";
        log_file << "Max True Jaccard," << max_true_jaccard << "\n";
        log_file << "Average ERM K-mers," << avg_erm_kmers << "\n";
        log_file << "Average EAM K-mers," << avg_eam_kmers << "\n";
        log_file << "Average ERM Minimizers," << avg_erm_minimizers << "\n";
        log_file << "Average EAM Minimizers," << avg_eam_minimizers << "\n";
        log_file << "Average Time K-mers," << avg_time_kmers << "\n";
        log_file << "Average Time Minimizers," << avg_time_minimizers << "\n";

        log_file.close();
        std::cout << "\nExperiment results have been written to experiment_results.csv" << std::endl;
    }

    return 0;
}