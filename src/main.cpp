/**
 * File: main.cpp
 * Project: foldcomp
 * Created: 2021-12-23 17:44:53
 * Author: Hyunbin Kim (khb7840@gmail.com)
 * Contributor: Milot Mirdita (milot@mirdita.de), Martin Steinegger (themartinsteinegger@gmail.com)
 * Description:
 *     This code contains main function for "foldcomp".
 *     Foldcomp is a fast lossy compression algorithm for protein structure.
 *     It encodes torsion angles with optimal number of bits and reconstruct
 *     3D coordinates from the encoded angles.
 * Usage:
 *    foldcomp compress input.pdb output.fcz
 *    foldcomp decompress input.fcz output.pdb
 * ---
 * Last Modified: 2026-02-10 20:09:37
 * Modified By: Hyunbin Kim (khb7840@gmail.com)
 * ---
 * Copyright © 2021 Hyunbin Kim, All rights reserved
 */
// Headers in the project
#include "atom_coordinate.h"
#include "foldcomp.h"
#include "structure_reader.h"
#include "utility.h"
#include "database_writer.h"
#include "tcbspan.h"
#include "execution_timer.h"
#include "input_processor.h"

// Standard libraries
#include <algorithm>
#include <cstdint>
#include <cstring>
#include <fstream> // IWYU pragma: keep
#include <iomanip>
#ifdef _WIN32
#include <direct.h>
#include "windows/getopt.h"
#include "windows/dirent.h"
#else
#include <getopt.h>
#endif
#include <iostream>
#include <sstream> // IWYU pragma: keep
#include <set>
#include <string>
#include <vector>

#include <sys/stat.h>

// OpenMP for parallelization
#ifdef OPENMP
#include <omp.h>
#endif

static int use_alt_order = 0;
static int anchor_residue_threshold = DEFAULT_ANCHOR_THRESHOLD;
static int save_as_tar = 0;
static int ext_mode = 0;
static int ext_plddt_digits = 1;
static int ext_merge = 1;
static int ext_use_title = 0;
static int overwrite = 0;

// version
#define FOLDCOMP_VERSION "1.0.0"
static constexpr uint8_t FOLDCOMP_FORMAT_VERSION_CONTAINER = 1;
static constexpr uint8_t FOLDCOMP_FORMAT_FLAG_CONTAINER = 1;

struct ContainerFragment {
    int16_t model = 1;
    std::string chain;
    std::string payload;
};

struct DecompressedSegment {
    int model = 1;
    std::vector<AtomCoordinate> atoms;
};

template <typename T>
static bool readBinaryValue(const char* data, size_t size, size_t& offset, T& out) {
    if (offset + sizeof(T) > size) {
        return false;
    }
    memcpy(&out, data + offset, sizeof(T));
    offset += sizeof(T);
    return true;
}

static bool hasFoldcompMagic(const char* data, size_t size) {
    if (size < MAGICNUMBER_LENGTH) {
        return false;
    }
    return memcmp(data, MAGICNUMBER, MAGICNUMBER_LENGTH) == 0;
}

static bool isContainerHeader(const CompressedFileHeader& header) {
    return header.version == FOLDCOMP_FORMAT_VERSION_CONTAINER &&
           (header.flags & FOLDCOMP_FORMAT_FLAG_CONTAINER) != 0 &&
           header.nResidue == 0 &&
           header.nAtom == 0;
}

static bool writeContainer(
    std::ostream& os, const std::string& title,
    const std::vector<ContainerFragment>& fragments
) {
    CompressedFileHeader header = {};
    header.version = FOLDCOMP_FORMAT_VERSION_CONTAINER;
    header.flags = FOLDCOMP_FORMAT_FLAG_CONTAINER;
    header.lenTitle = static_cast<uint32_t>(title.size());

    os.write(MAGICNUMBER, MAGICNUMBER_LENGTH);
    os.write(reinterpret_cast<const char*>(&header), sizeof(header));
    uint32_t nFragments = static_cast<uint32_t>(fragments.size());
    os.write(reinterpret_cast<const char*>(&nFragments), sizeof(nFragments));
    if (!title.empty()) {
        os.write(title.data(), static_cast<std::streamsize>(title.size()));
    }

    for (const auto& fragment : fragments) {
        uint8_t chainLen = static_cast<uint8_t>(std::min<size_t>(255, fragment.chain.size()));
        uint32_t payloadSize = static_cast<uint32_t>(fragment.payload.size());
        os.write(reinterpret_cast<const char*>(&fragment.model), sizeof(fragment.model));
        os.write(reinterpret_cast<const char*>(&chainLen), sizeof(chainLen));
        if (chainLen > 0) {
            os.write(fragment.chain.data(), chainLen);
        }
        os.write(reinterpret_cast<const char*>(&payloadSize), sizeof(payloadSize));
        if (payloadSize > 0) {
            os.write(fragment.payload.data(), static_cast<std::streamsize>(payloadSize));
        }
    }

    return static_cast<bool>(os);
}

static bool readContainer(
    const char* data, size_t size, std::string& title,
    std::vector<ContainerFragment>& fragments
) {
    title.clear();
    fragments.clear();
    if (!hasFoldcompMagic(data, size)) {
        return false;
    }
    if (size < MAGICNUMBER_LENGTH + sizeof(CompressedFileHeader) + sizeof(uint32_t)) {
        return false;
    }

    size_t offset = MAGICNUMBER_LENGTH;
    CompressedFileHeader header = {};
    memcpy(&header, data + offset, sizeof(header));
    offset += sizeof(header);
    if (!isContainerHeader(header)) {
        return false;
    }

    uint32_t nFragments = 0;
    if (!readBinaryValue(data, size, offset, nFragments)) {
        return false;
    }

    if (offset + header.lenTitle > size) {
        return false;
    }
    if (header.lenTitle > 0) {
        title.assign(data + offset, data + offset + header.lenTitle);
        offset += header.lenTitle;
    }

    fragments.reserve(nFragments);
    for (uint32_t i = 0; i < nFragments; i++) {
        ContainerFragment fragment;
        uint8_t chainLen = 0;
        uint32_t payloadSize = 0;
        if (!readBinaryValue(data, size, offset, fragment.model) ||
            !readBinaryValue(data, size, offset, chainLen)) {
            return false;
        }
        if (offset + chainLen > size) {
            return false;
        }
        if (chainLen > 0) {
            fragment.chain.assign(data + offset, data + offset + chainLen);
            offset += chainLen;
        }
        if (!readBinaryValue(data, size, offset, payloadSize)) {
            return false;
        }
        if (offset + payloadSize > size) {
            return false;
        }
        fragment.payload.assign(data + offset, data + offset + payloadSize);
        offset += payloadSize;
        fragments.push_back(std::move(fragment));
    }
    return true;
}

static void writeTitleLines(const std::string& title, std::ostream& pdb_stream) {
    if (title.empty()) {
        return;
    }
    const char* headerData = title.c_str();
    size_t headerLen = title.length();
    int remainingHeader = static_cast<int>(headerLen);
    char buffer[128];
    int written = snprintf(
        buffer, sizeof(buffer), "TITLE     %.*s\n",
        std::min(70, remainingHeader), headerData
    );
    if (written >= 0 && written < static_cast<int>(sizeof(buffer))) {
        pdb_stream << buffer;
    }
    remainingHeader -= 70;
    int continuation = 2;
    while (remainingHeader > 0) {
        written = snprintf(
            buffer, sizeof(buffer), "TITLE  % 3d%.*s\n",
            continuation, std::min(70, remainingHeader),
            headerData + (headerLen - remainingHeader)
        );
        if (written >= 0 && written < static_cast<int>(sizeof(buffer))) {
            pdb_stream << buffer;
        }
        remainingHeader -= 70;
        continuation++;
    }
}

static void writeSegmentsToPDB(
    const std::vector<DecompressedSegment>& segments,
    const std::string& title,
    std::ostream& pdb_stream
) {
    if (segments.empty()) {
        return;
    }
    writeTitleLines(title, pdb_stream);
    std::set<int> modelSet;
    for (const auto& segment : segments) {
        modelSet.insert(segment.model);
    }
    bool writeModels = modelSet.size() > 1;
    int currModel = -1;
    for (const auto& segment : segments) {
        if (writeModels && segment.model != currModel) {
            if (currModel != -1) {
                pdb_stream << "ENDMDL\n";
            }
            pdb_stream << "MODEL     " << std::setw(4) << segment.model << '\n';
            currModel = segment.model;
        }
        std::vector<AtomCoordinate> atoms = segment.atoms;
        writeAtomCoordinatesToPDB(atoms, "", pdb_stream);
    }
    if (writeModels) {
        pdb_stream << "ENDMDL\n";
    }
}

static std::string makeFragmentSuffix(const ContainerFragment& fragment, size_t index) {
    std::ostringstream oss;
    oss << "|m" << fragment.model;
    if (!fragment.chain.empty()) {
        oss << "|c" << fragment.chain;
    }
    oss << "|f" << index;
    return oss.str();
}

int print_usage(void) {
    std::cout << "Usage: foldcomp compress <pdb|cif> [<fcz>]" << std::endl;
    std::cout << "       foldcomp compress [-t number] <dir|tar(.gz)> [<dir|tar|db>]" << std::endl;
    std::cout << "       foldcomp decompress <fcz|tar> [<pdb>]" << std::endl;
    std::cout << "       foldcomp decompress [-t number] <dir|tar(.gz)|db> [<dir|tar>]" << std::endl;
    std::cout << "       foldcomp extract [--plddt|--amino-acid] <fcz> [<fasta>]" << std::endl;
    std::cout << "       foldcomp extract [--plddt|--amino-acid] [-t number] <dir|tar(.gz)|db> [<fasta_out>]" << std::endl;
    std::cout << "       foldcomp check <fcz>" << std::endl;
    std::cout << "       foldcomp check [-t number] <dir|tar(.gz)|db>" << std::endl;
    std::cout << "       foldcomp rmsd <pdb|cif> <pdb|cif>" << std::endl;
    std::cout << " -h, --help               print this help message" << std::endl;
    std::cout << " -v, --version            print version" << std::endl;
    std::cout << " -t, --threads            threads for (de)compression of folders/tar files [default=1]" << std::endl;
    std::cout << " -r, --recursive          recursively look for files in directory [default=0]" << std::endl;
    std::cout << " -f, --file               input is a list of files [default=0]" << std::endl;
    std::cout << " -a, --alt                use alternative atom order [default=false]" << std::endl;
    std::cout << " -b, --break              interval size to save absolute atom coordinates [default=" << anchor_residue_threshold << "]" << std::endl;
    std::cout << " -z, --tar                save as tar file [default=false]" << std::endl;
    std::cout << " -d, --db                 save as database [default=false]" << std::endl;
    std::cout << " -y, --overwrite          overwrite existing files [default=false]" << std::endl;
    std::cout << " -l, --id-list            a file of id list to be processed (only for database input)" << std::endl;
    std::cout << " -m, --id-mode            id mode for database input. 0: database keys, 1: names (.lookup) [default=1]" << std::endl;
    std::cout << " --skip-discontinuous     skip PDB with with discontinuous residues (only batch compression)" << std::endl;
    std::cout << " --check                  check FCZ before and skip entries with error (only for batch decompression)" << std::endl;
    std::cout << " --plddt                  extract pLDDT score (only for extraction mode)" << std::endl;
    std::cout << " -p, --plddt-digits       extract pLDDT score with specified number of digits (only for extraction mode)" << std::endl;
    std::cout << "                          - 1: single digit (fasta-like format), 2: 2-digit(00-99; tsv), 3: 3-digit, 4: 4-digit (max)" << std::endl;
    std::cout << " --fasta, --amino-acid    extract amino acid sequence (only for extraction mode)" << std::endl;
    std::cout << " --no-merge               do not merge output files (only for extraction mode)" << std::endl;
    std::cout << " --use-title              use TITLE as the output file name (only for extraction mode)" << std::endl;
    std::cout << " --time                   measure time for compression/decompression" << std::endl;
    std::cout << " --use-cache              use cached index for database input [default=false]" << std::endl;
    return 0;
}

inline int print_version(void) {
    std::cout << "foldcomp " << FOLDCOMP_VERSION << std::endl;
    return 0;
}

int rmsd(const std::string& pdb1, const std::string& pdb2) {
    // RMSD calculation between two PDB/mmCIF files
    StructureReader reader;
    reader.load(pdb1);
    std::vector<AtomCoordinate> atomCoordinates1;
    reader.readAllAtoms(atomCoordinates1);
    reader.load(pdb2);
    std::vector<AtomCoordinate> atomCoordinates2;
    reader.readAllAtoms(atomCoordinates2);
    // Check
    if (atomCoordinates1.size() == 0) {
        std::cerr << "[Error] No atoms found in the input file: " << pdb1 << std::endl;
        return 1;
    }
    if (atomCoordinates2.size() == 0) {
        std::cerr << "[Error] No atoms found in the input file: " << pdb2 << std::endl;
        return 1;
    }
    if (atomCoordinates1.size() != atomCoordinates2.size()) {
        std::cerr << "[Error] The number of atoms in the two files are different." << std::endl;
        return 1;
    }
    std::vector<AtomCoordinate> backbone1 = filterBackbone(atomCoordinates1);
    std::vector<AtomCoordinate> backbone2 = filterBackbone(atomCoordinates2);
    // Print
    std::cout << pdb1 << '\t' << pdb2 << '\t';
    std::cout << backbone1.size() / 3 << '\t' << atomCoordinates1.size() << '\t';
    std::cout << RMSD(backbone1, backbone2) << '\t';
    std::cout << RMSD(atomCoordinates1, atomCoordinates2) << std::endl;
    return 0;
}

int main(int argc, char* const *argv) {
    if (argc < 3) {
        // Check if version is requested
        if (argc == 2 && (strcmp(argv[1], "--version") == 0 || strcmp(argv[1], "-v") == 0)) {
            return print_version();
        }
        return print_usage();
    }

    int option_index = 0;
    int num_threads = 1;
    int has_output = 0;
    int recursive = 0;
    int file_input = 0;
    int db_output = 0;
    int measure_time = 0;
    int skip_discontinuous = 0;
    int check_before_decompression = 0;
    int id_mode = 1;
    int use_cache = 0;
    std::string user_id_file = "";
    std::vector<std::string> user_names;
    std::vector<uint32_t> user_ids;

    // Mode - non-optional argument
    enum {
        COMPRESS,
        DECOMPRESS,
        EXTRACT,
        CHECK,
        RMSD
    } mode = COMPRESS;

    // Define command line options
    static struct option long_options[] = {
            {"help",               no_argument,                           0, 'h'},
            {"alt",                no_argument,                           0, 'a'},
            {"tar",                no_argument,                           0, 'z'},
            {"recursive",          no_argument,                           0, 'r'},
            {"file",               no_argument,                           0, 'f'},
            {"plddt",              no_argument,                   &ext_mode,  0 },
            {"fasta",              no_argument,                   &ext_mode,  1 },
            {"amino-acid",         no_argument,                   &ext_mode,  1 },
            {"no-merge",           no_argument,                  &ext_merge,  0 },
            {"overwrite",          no_argument,                  &overwrite,  1 },
            {"time",               no_argument,               &measure_time,  1 },
            {"skip-discontinuous", no_argument,         &skip_discontinuous,  1 },
            {"check",              no_argument, &check_before_decompression,  1 },
            {"use-title",          no_argument,              &ext_use_title,  1 },
            {"db",                 no_argument,                           0, 'd'},
            {"version",            no_argument,                           0, 'v'},
            {"threads",      required_argument,                           0, 't'},
            {"break",        required_argument,                           0, 'b'},
            {"id-list",      required_argument,                           0, 'l'},
            {"id-mode",      required_argument,                           0, 'm'},
            {"plddt-digits", required_argument,                           0, 'p'},
            {"use-cache",      no_argument,                          &use_cache,  1 },
            {0,                              0,                           0,  0 }
    };

    // Parse command line options with getopt_long
    int flag = getopt_long(argc, argv, "hadzrfyvt:b:l:p:", long_options, &option_index);
    while (flag != -1) {
        switch (flag) {
            case 'h':
                return print_usage();
            case 't':
                num_threads = atoi(optarg);
                break;
            case 'a':
                use_alt_order = 1;
                break;
            case 'z':
                save_as_tar = 1;
                break;
            case 'r':
                recursive = 1;
                break;
            case 'f':
                file_input = 1;
                break;
            case 'y':
                overwrite = 1;
                break;
            case 'b':
                anchor_residue_threshold = atoi(optarg);
                break;
            case 'l':
                user_id_file = std::string(optarg);
                break;
            case 'm':
                id_mode = atoi(optarg);
                if (id_mode != 0 && id_mode != 1) {
                    std::cerr << "[Error] Invalid id mode. Please use 0 or 1." << std::endl;
                    return print_usage();
                }
                break;
            case 'd':
                db_output = 1;
                break;
            case 'p':
                ext_plddt_digits = atoi(optarg);
                break;
            case 'v':
                return print_version();
            case '?':
                return print_usage();
            default:
                break;
        }
        flag = getopt_long(argc, argv, "hadzrfyt:b:l:p:", long_options, &option_index);
    }

    // Parse non-option arguments
    // argv[optind]: MODE
    // argv[optind + 1]: INPUT
    // argv[optind + 2]: OUTPUT (optional)
    if ((optind + 1) >= argc) {
        std::cerr << "[Error] Not enough arguments." << std::endl;
        return print_usage();
    }

    struct stat inputStat;
    int inputExists = stat(argv[optind + 1], &inputStat);
    const char* outputSuffix = "";
    bool mayHaveOutput = false;
    // get mode from command line
    if (strcmp(argv[optind], "compress") == 0) {
        mode = COMPRESS;
#ifdef HAVE_GCS
        if ((optind + 1) < argc && stringStartsWith("gcs://", argv[optind + 1])) {
            fileExists = 1;
        }
#endif
        mayHaveOutput = true;
        outputSuffix = "fcz";
    } else if (strcmp(argv[optind], "decompress") == 0) {
        mode = DECOMPRESS;
        mayHaveOutput = true;
        outputSuffix = "pdb";
    } else if (strcmp(argv[optind], "extract") == 0) {
        mode = EXTRACT;
        mayHaveOutput = true;
        if (ext_mode == 0){
            outputSuffix = "plddt";
            if (ext_plddt_digits != 1) {
                outputSuffix = "plddt.tsv";
            }
        } else if (ext_mode == 1) {
            outputSuffix = "fasta";
        }
    } else if (strcmp(argv[optind], "check") == 0){
        mode = CHECK;
    } else if (strcmp(argv[optind], "rmsd") == 0) {
        mode = RMSD;
    } else {
        return print_usage();
    }

    // Error if no input file given
    if (inputExists == -1) {
        std::cerr << "[Error] " << argv[optind + 1] << " does not exist." << std::endl;
        return EXIT_FAILURE;
    }

    std::string input = argv[optind + 1];
    while (input.back() == '/') {
        input.pop_back();
    }
    std::vector<std::string> inputs;
    std::vector<std::string> single_file_inputs;
    if (file_input) {
        std::ifstream inputFile(input);
        if (!inputFile) {
            std::cerr << "[Error] Could not open file " << input << std::endl;
            return EXIT_FAILURE;
        }
        std::string line;
        while (std::getline(inputFile, line)) {
            // If the file ends with .pdb, .pdb.gz, .cif, .cif.gz, .fcz, assume it is a single file input
            if (stringEndsWith(".pdb", line) || stringEndsWith(".pdb.gz", line) ||
                stringEndsWith(".cif", line) || stringEndsWith(".cif.gz", line) ||
                stringEndsWith(".fcz", line)) {
                single_file_inputs.push_back(line);
            } else {
                inputs.push_back(line);
            }
        }
    } else {
        inputs.push_back(input);
    }

    std::string output;
    if (argc == optind + 3) {
        has_output = 1;
        output = argv[optind + 2];
        if (stringEndsWith(".tar", output)) {
            save_as_tar = 1;
        }
    }

    bool isSingleFileInput = !file_input && (S_ISREG(inputStat.st_mode) || S_ISLNK(inputStat.st_mode));
    if (!file_input) {
        struct stat st;
        if (stringEndsWith(".tar", input) || stringEndsWith(".tar.gz", input) || stringEndsWith(".tgz", input)) {
            isSingleFileInput = false;
        } else if (stat((input + ".dbtype").c_str(), &st) == 0) {
            isSingleFileInput = false;
        }
#ifdef HAVE_GCS
        else if (stringStartsWith("gcs://", input)) {
            isSingleFileInput = false;
        }
#endif
    }

    while (has_output && output.back() == '/') {
        output.pop_back();
    }
    if (mayHaveOutput && !has_output) {
        if (db_output) {
            output = input + "_db";
        } else if (save_as_tar) {
            output = input + "." +  outputSuffix + ".tar";
        } else {
            if (isSingleFileInput) {
                output = getFileWithoutExt(input) + "." + outputSuffix;
            } else {
                output = input + "_" + outputSuffix + "/";
            }
        }
    }

    if (mode == RMSD) {
        // Calculate RMSD between two PDB files
        rmsd(input, output);
    } else if (mode == COMPRESS) {
        // output variants
        void* handle;
        mtar_t tar_out;
        if (save_as_tar) {
            mtar_open(&tar_out, output.c_str(), "w");
        } else if (db_output) {
            handle = make_writer(output.c_str(), (output + ".index").c_str());
        } else if (!isSingleFileInput) {
            struct stat st;
            if (stat(output.c_str(), &st) == -1) {
#ifdef _WIN32
                _mkdir(output.c_str());
#else
                mkdir(output.c_str(), 0755);
#endif
            }
        }

        if (isSingleFileInput) {
            std::cout << "Compressing " << input << " to " << output << std::endl;
        } else {
            std::cout << "Compressing files in " << input;
            std::cout << " using " << num_threads << " threads" << std::endl;
            if (db_output) {
                std::cout << "Output database: " << output << std::endl;
            } else if (save_as_tar) {
                std::cout << "Output tar file: " << output << std::endl;
            } else {
                std::cout << "Output directory: " << output << std::endl;
            }
        }

        unsigned int key = 0;
        for (size_t i = 0; i < inputs.size() + 1; i++) {
            const std::string& input = (i == inputs.size()) ? "" : inputs[i];
            Processor* processor;
            struct stat st;
            if (i != inputs.size()) {
                if (stringEndsWith(".tar", input) || stringEndsWith(".tar.gz", input) || stringEndsWith(".tgz", input)) {
                    processor = new TarProcessor(input);
                }
                else if (stat((input + ".dbtype").c_str(), &st) == 0) {
                    if (user_id_file.size() > 0) {
                        processor = new DatabaseProcessor(input, user_id_file, id_mode, use_cache);
                    }
                    else {
                        processor = new DatabaseProcessor(input);
                    }
                }
#ifdef HAVE_GCS
                else if (stringStartsWith("gcs://", input)) {
                    processor = new GcsProcessor(input);
                }
#endif
                else {
                    processor = new DirectoryProcessor(input, recursive);
                }
            } else {
                if (single_file_inputs.size() > 0) {
                    processor = new DirectoryProcessor(single_file_inputs);
                } else {
                    continue;
                }
            }
            process_entry_func func = [&](const char* name, const char* dataBuffer, size_t size) -> bool {
                TimerGuard guard(name, measure_time);
                std::vector<AtomCoordinate> atomCoordinates;
                std::vector<BackboneChain> compData;
                StructureReader reader;

                std::string base = baseName(name);
                std::pair<std::string, std::string> outputParts = getFileParts(base);
                std::string outputFile;
                if (save_as_tar) {
                    outputFile = outputParts.first;
                } else if (db_output) {
                    outputFile = outputParts.first;
                } else if (isSingleFileInput) {
                    outputParts = getFileParts(output);
                    outputFile = outputParts.first;
                } else {
                    outputFile = output + "/" + outputParts.first;
                }
                reader.loadFromBuffer(dataBuffer, size, base);
                reader.readAllAtoms(atomCoordinates);
                if (atomCoordinates.size() == 0) {
                    std::cerr << "[Error] No atoms found in the input file: " << base << std::endl;
                    return false;
                }

                // replace the title with only the base name if it was the same as the file name
                std::string title = reader.title == base ? outputParts.first : reader.title;

                removeAlternativePosition(atomCoordinates);

                std::vector<ContainerFragment> encodedFragments;
                std::vector<std::pair<size_t, size_t>> chain_indices = identifyChains(atomCoordinates);
                // Check if there are multiple chains or regions with discontinous residue indices
                for (size_t i = 0; i < chain_indices.size(); i++) {
                    std::vector<std::pair<size_t, size_t>> frag_indices = identifyDiscontinousResInd(
                        atomCoordinates, chain_indices[i].first, chain_indices[i].second
                    );
                    if (skip_discontinuous && frag_indices.size() > 1) {
                        std::string message = "Skipping discontinuous chain: " + base + "\n";
                        std::cerr << message;
                        continue;
                    }
                    for (size_t j = 0; j < frag_indices.size(); j++) {
                        tcb::span<AtomCoordinate> frag_span = tcb::span<AtomCoordinate>(
                            &atomCoordinates[frag_indices[j].first],
                            atomCoordinates.data() + frag_indices[j].second
                        );
                        Foldcomp compRes;
                        compRes.strTitle = title;
                        compRes.anchorThreshold = anchor_residue_threshold;
                        compData = compRes.compress(frag_span);
                        if (compData.empty()) {
                            std::cerr << "[Warning] Skipping fragment with incomplete backbone: " << base << std::endl;
                            continue;
                        }
                        std::ostringstream oss;
                        compRes.writeStream(oss);
                        ContainerFragment containerFragment;
                        containerFragment.model = atomCoordinates[frag_indices[j].first].model;
                        containerFragment.chain = atomCoordinates[frag_indices[j].first].chain;
                        containerFragment.payload = oss.str();
                        encodedFragments.push_back(std::move(containerFragment));
                        compData.clear();
                    }
                }

                if (encodedFragments.empty()) {
                    std::cerr << "[Warning] No encodable fragments found in input file: " << base << std::endl;
                    return true;
                }

                bool useContainer = encodedFragments.size() > 1;
                std::string encoded;
                if (useContainer) {
                    std::ostringstream oss;
                    if (!writeContainer(oss, title, encodedFragments)) {
                        std::cerr << "[Error] Failed to write container payload: " << base << std::endl;
                        return false;
                    }
                    encoded = oss.str();
                } else {
                    encoded = encodedFragments[0].payload;
                }

                std::string filename = outputFile;
                if (!db_output) {
                    if (isCompressible(outputParts)) {
                        filename += ".fcz";
                    } else {
                        filename += "." + outputParts.second;
                    }
                }

                if (db_output) {
                    std::string dbKey = baseName(filename);
                    std::replace(dbKey.begin(), dbKey.end(), '.', '_');
#pragma omp critical
                    {
                        writer_append(handle, encoded.c_str(), encoded.size(), key, dbKey.c_str());
                        key++;
                    }
                } else if (save_as_tar) {
#pragma omp critical
                    {
                        mtar_write_file_header(&tar_out, baseName(filename).c_str(), encoded.size());
                        mtar_write_data(&tar_out, encoded.c_str(), encoded.size());
                    }
                } else {
                    if (stat(filename.c_str(), &st) == 0 && !overwrite) {
                        std::cerr << "[Error] Output file already exists: " << baseName(outputFile) << std::endl;
                        return false;
                    }
                    std::ofstream outFile(filename, std::ios::binary);
                    if (!outFile) {
                        std::cerr << "[Error] Failed to open output file: " << filename << std::endl;
                        return false;
                    }
                    outFile.write(encoded.data(), encoded.size());
                    outFile.close();
                }
                atomCoordinates.clear();
                return true;
            };
            processor->run(func, num_threads);
            delete processor;
        }
        if (db_output) {
            free_writer(handle);
        } else if (save_as_tar) {
            mtar_write_finalize(&tar_out);
            mtar_close(&tar_out);
        }
    } else if (mode == DECOMPRESS) {
        void* handle;
        mtar_t tar_out;
        if (save_as_tar) {
            mtar_open(&tar_out, output.c_str(), "w");
        } else if (db_output) {
            handle = make_writer(output.c_str(), (output + ".index").c_str());
        } else if (!isSingleFileInput) {
            struct stat st;
            if (stat(output.c_str(), &st) == -1) {
#ifdef _WIN32
                _mkdir(output.c_str());
#else
                mkdir(output.c_str(), 0755);
#endif
            }
        }

        if (isSingleFileInput) {
            std::cout << "Decompressing " << input << " to " << output << std::endl;
        } else {
            std::cout << "Decompressing files in " << input;
            std::cout << " using " << num_threads << " threads" << std::endl;
            if (db_output) {
                std::cout << "Output database: " << output << std::endl;
            } else if (save_as_tar) {
                std::cout << "Output tar file: " << output << std::endl;
            } else {
                std::cout << "Output directory: " << output << std::endl;
            }
        }

        unsigned int key = 0;
        for (size_t i = 0; i < inputs.size() + 1; i++) {
            const std::string& input = (i == inputs.size()) ? "" : inputs[i];
            Processor* processor;
            struct stat st;
            if (i != inputs.size()) {
                if (stringEndsWith(".tar", input) || stringEndsWith(".tar.gz", input) || stringEndsWith(".tgz", input)) {
                    processor = new TarProcessor(input);
                }
                else if (stat((input + ".dbtype").c_str(), &st) == 0) {
                    if (user_id_file.size() > 0) {
                        processor = new DatabaseProcessor(input, user_id_file, id_mode, use_cache);
                    }
                    else {
                        processor = new DatabaseProcessor(input);
                    }
                }
#ifdef HAVE_GCS
                else if (stringStartsWith("gcs://", input)) {
                    processor = new GcsProcessor(input);
                }
#endif
                else {
                    processor = new DirectoryProcessor(input, recursive);
                }
            } else {
                if (single_file_inputs.size() > 0) {
                    processor = new DirectoryProcessor(single_file_inputs);
                }
                else {
                    continue;
                }
            }

            process_entry_func func = [&](const char* name, const char* dataBuffer, size_t size) -> bool {
                TimerGuard guard(name, measure_time);
                std::vector<DecompressedSegment> segments;
                std::string structureTitle;

                std::string containerTitle;
                std::vector<ContainerFragment> containerFragments;
                bool isContainer = readContainer(dataBuffer, size, containerTitle, containerFragments);
                if (isContainer) {
                    structureTitle = containerTitle;
                    for (const auto& fragment : containerFragments) {
                        Foldcomp compRes;
                        std::istringstream fragmentInput(fragment.payload);
                        int flag = compRes.read(fragmentInput);
                        if (flag != 0) {
                            std::cerr << "[Error] Failed to read a container fragment." << std::endl;
                            return false;
                        }
                        compRes.useAltAtomOrder = use_alt_order;
                        if (check_before_decompression) {
                            ValidityError err = compRes.checkValidity();
                            if (err != ValidityError::SUCCESS) {
                                printValidityError(err, compRes.strTitle);
                                return true;
                            }
                        }
                        std::vector<AtomCoordinate> atomCoordinates;
                        flag = compRes.decompress(atomCoordinates);
                        if (flag != 0) {
                            std::cerr << "[Error] decompressing container fragment." << std::endl;
                            return false;
                        }
                        for (auto& atom : atomCoordinates) {
                            atom.model = fragment.model;
                        }
                        segments.push_back({fragment.model, std::move(atomCoordinates)});
                        if (structureTitle.empty()) {
                            structureTitle = compRes.strTitle;
                        }
                    }
                } else {
                    Foldcomp compRes;
                    std::istringstream input(std::string(dataBuffer, size));
                    int flag = compRes.read(input);
                    if (flag != 0) {
                        if (flag == -1) {
                            std::cerr << "[Error] File is not a valid fcz file" << std::endl;
                        } else if (flag == -2) {
                            std::cerr << "[Error] Could not restore prevAtoms" << std::endl;
                        } else {
                            std::cerr << "[Error] Unknown read error" << std::endl;
                        }
                        return false;
                    }
                    compRes.useAltAtomOrder = use_alt_order;
                    if (check_before_decompression) {
                        ValidityError err = compRes.checkValidity();
                        if (err != ValidityError::SUCCESS) {
                            printValidityError(err, compRes.strTitle);
                            return true;
                        }
                    }
                    std::vector<AtomCoordinate> atomCoordinates;
                    flag = compRes.decompress(atomCoordinates);
                    if (flag != 0) {
                        std::cerr << "[Error] decompressing compressed data." << std::endl;
                        return false;
                    }
                    segments.push_back({1, std::move(atomCoordinates)});
                    structureTitle = compRes.strTitle;
                }

                if (segments.empty()) {
                    std::cerr << "[Error] No segments were decompressed." << std::endl;
                    return false;
                }

                std::string base = baseName(name);
                std::pair<std::string, std::string> outputParts = getFileParts(base);
                std::string outputFile;
                if (save_as_tar) {
                    outputFile = outputParts.first + "." + outputSuffix;
                } else if (db_output) {
                    outputFile = outputParts.first;
                } else if (isSingleFileInput) {
                    outputFile = output;
                } else {
                    outputFile = output + "/" + outputParts.first + "." + outputSuffix;
                }

                std::ostringstream pdbStream;
                writeSegmentsToPDB(segments, structureTitle, pdbStream);
                std::string pdbText = pdbStream.str();

                if (db_output) {
                    pdbText.push_back('\0');
#pragma omp critical
                    {
                        writer_append(handle, pdbText.c_str(), pdbText.size(), key, outputFile.c_str());
                        key++;
                    }
                } else if (save_as_tar) {
#pragma omp critical
                    {
                        mtar_write_file_header(&tar_out, outputFile.c_str(), pdbText.size());
                        mtar_write_data(&tar_out, pdbText.c_str(), pdbText.size());
                    }
                } else {
                    // Write decompressed data to file
                    // Check output file exists
                    if (stat(outputFile.c_str(), &st) == 0 && !overwrite) {
                        std::cerr << "[Error] Output file already exists: " << baseName(outputFile) << std::endl;
                        return false;
                    }
                    std::ofstream out(outputFile, std::ios::binary);
                    if (!out) {
                        std::cerr << "[Error] Writing decompressed data to file: " << output << std::endl;
                        return false;
                    }
                    out.write(pdbText.data(), pdbText.size());
                }
                return true;
            };
            processor->run(func, num_threads);
            delete processor;
        }
        if (db_output) {
            free_writer(handle);
        } else if (save_as_tar) {
            mtar_write_finalize(&tar_out);
            mtar_close(&tar_out);
        }
    } else if (mode == EXTRACT) {
        void* handle;
        mtar_t tar_out;
        if (save_as_tar) {
            mtar_open(&tar_out, output.c_str(), "w");
        } else if (db_output) {
            handle = make_writer(output.c_str(), (output + ".index").c_str());
        } else {
            struct stat st;
            if (stat(output.c_str(), &st) == -1 && !ext_merge) {
#ifdef _WIN32
                _mkdir(output.c_str());
#else
                mkdir(output.c_str(), 0755);
#endif
            }
        }

        if (isSingleFileInput) {
            std::cout << "Extracting " << input << " to " << output << std::endl;
        } else {
            std::cout << "Extracting files in " << input;
            std::cout << " using " << num_threads << " threads" << std::endl;
            if (db_output) {
                std::cout << "Output database: " << output << std::endl;
            } else if (save_as_tar) {
                std::cout << "Output tar file: " << output << std::endl;
            } else if (!ext_merge) {
                std::cout << "Output directory: " << output << std::endl;
            } else {
                // Single file output. Remove "/" from end of output
                if (output[output.size() - 1] == '/') {
                    output = output.substr(0, output.size() - 1);
                }
                std::cout << "Output: " << output << std::endl;
            }
        }

        bool isMergedOutput = false;
        std::ofstream default_out;
        if (!save_as_tar && !db_output && !isSingleFileInput && ext_merge) {
            default_out.open(output, std::ios::out);
            isMergedOutput = true;
        }

        unsigned int key = 0;
        for (size_t i = 0; i < inputs.size() + 1; i++) {
            const std::string& input = (i == inputs.size()) ? "" : inputs[i];
            Processor* processor;
            struct stat st;
            if (i != inputs.size()) {
                if (stringEndsWith(".tar", input) || stringEndsWith(".tar.gz", input) || stringEndsWith(".tgz", input)) {
                    processor = new TarProcessor(input);
                }
                else if (stat((input + ".dbtype").c_str(), &st) == 0) {
                    if (user_id_file.size() > 0) {
                        processor = new DatabaseProcessor(input, user_id_file, id_mode, use_cache);
                    }
                    else {
                        processor = new DatabaseProcessor(input);
                    }
                }
#ifdef HAVE_GCS
                else if (stringStartsWith("gcs://", input)) {
                    processor = new GcsProcessor(input);
                }
#endif
                else {
                    processor = new DirectoryProcessor(input, recursive);
                }
            } else {
                if (single_file_inputs.size() > 0) {
                    processor = new DirectoryProcessor(single_file_inputs);
                }
                else {
                    continue;
                }
            }

            process_entry_func func = [&](const char* name, const char* dataBuffer, size_t size) -> bool {
                std::string strName(name);
                std::vector<ContainerFragment> containerFragments;
                std::string containerTitle;
                bool isContainer = readContainer(dataBuffer, size, containerTitle, containerFragments);

                std::ostringstream extractedStream;
                if (!isContainer) {
                    std::string outputTitle = strName;
                    std::string extractedData;
                    int totalResidue = 0;
                    std::istringstream input(std::string(dataBuffer, size));
                    Foldcomp compRes;
                    int flag = compRes.read(input);
                    if (flag != 0) {
                        if (flag == -1) {
                            std::cerr << "[Error] File is not a valid fcz file" << std::endl;
                        } else if (flag == -2) {
                            std::cerr << "[Error] Could not restore prevAtoms" << std::endl;
                        } else {
                            std::cerr << "[Error] Unknown read error" << std::endl;
                        }
                        return false;
                    }
                    if (ext_use_title && !compRes.strTitle.empty()) {
                        outputTitle = compRes.strTitle;
                    }
                    compRes.extract(extractedData, ext_mode, ext_plddt_digits);
                    totalResidue = compRes.nResidue;

                    if (ext_mode == 0 && ext_plddt_digits > 1) {
                        extractedStream << outputTitle << "\t" << totalResidue << "\t" << extractedData << "\n";
                    } else {
                        extractedStream << ">" << outputTitle << "\n" << extractedData << "\n";
                    }
                } else {
                    struct ExtractGroup {
                        int16_t model = 1;
                        std::string chain;
                        int totalResidue = 0;
                        std::string data;
                    };
                    std::vector<ExtractGroup> groups;
                    std::string baseTitle = ext_use_title ? containerTitle : strName;

                    for (size_t fragmentIndex = 0; fragmentIndex < containerFragments.size(); fragmentIndex++) {
                        Foldcomp compRes;
                        std::istringstream fragmentInput(containerFragments[fragmentIndex].payload);
                        int flag = compRes.read(fragmentInput);
                        if (flag != 0) {
                            std::cerr << "[Error] Failed to read a container fragment during extraction." << std::endl;
                            return false;
                        }
                        if (baseTitle.empty() && ext_use_title && !compRes.strTitle.empty()) {
                            baseTitle = compRes.strTitle;
                        }

                        size_t groupIdx = groups.size();
                        for (size_t i = 0; i < groups.size(); i++) {
                            if (groups[i].model == containerFragments[fragmentIndex].model &&
                                groups[i].chain == containerFragments[fragmentIndex].chain) {
                                groupIdx = i;
                                break;
                            }
                        }
                        if (groupIdx == groups.size()) {
                            ExtractGroup g;
                            g.model = containerFragments[fragmentIndex].model;
                            g.chain = containerFragments[fragmentIndex].chain;
                            groups.push_back(std::move(g));
                        }

                        std::string fragmentData;
                        compRes.extract(fragmentData, ext_mode, ext_plddt_digits);
                        if (ext_mode == 0 && ext_plddt_digits > 1 && !groups[groupIdx].data.empty()) {
                            groups[groupIdx].data.push_back(',');
                        }
                        groups[groupIdx].data += fragmentData;
                        groups[groupIdx].totalResidue += compRes.nResidue;
                    }

                    if (baseTitle.empty()) {
                        baseTitle = strName;
                    }
                    bool annotateGroup = groups.size() > 1;
                    for (const auto& group : groups) {
                        std::string groupTitle = baseTitle;
                        if (annotateGroup) {
                            std::ostringstream suffix;
                            suffix << "|m" << group.model;
                            if (!group.chain.empty()) {
                                suffix << "|c" << group.chain;
                            }
                            groupTitle += suffix.str();
                        }
                        if (ext_mode == 0 && ext_plddt_digits > 1) {
                            extractedStream << groupTitle << "\t" << group.totalResidue << "\t" << group.data << "\n";
                        } else {
                            extractedStream << ">" << groupTitle << "\n" << group.data << "\n";
                        }
                    }
                }
                std::string extractedText = extractedStream.str();
                std::string base = baseName(name);
                std::pair<std::string, std::string> outputParts = getFileParts(base);
                std::string outputFile;
                if (save_as_tar) {
                    outputFile = outputParts.first + "." + outputSuffix;
                } else if (db_output) {
                    outputFile = outputParts.first;
                } else if (isSingleFileInput) {
                    outputFile = output;
                } else {
                    outputFile = output + "/" + outputParts.first + "." + outputSuffix;
                }

                if (isMergedOutput) {
#pragma omp critical
                    {
                        default_out << extractedText;
                    }
                } else if (db_output) {
                    extractedText.push_back('\0');
#pragma omp critical
                    {
                        writer_append(handle, extractedText.c_str(), extractedText.size(), key, outputFile.c_str());
                        key++;
                    }
                } else if (save_as_tar) {
#pragma omp critical
                    {
                        mtar_write_file_header(&tar_out, outputFile.c_str(), extractedText.size());
                        mtar_write_data(&tar_out, extractedText.c_str(), extractedText.size());
                    }
                } else {
                    std::ofstream output(outputFile);
                    if (!output) {
                        std::cerr << "[Error] Could not open file " << outputFile << std::endl;
                        return false;
                    }
                    output << extractedText;
                }

                return true;
            };
            processor->run(func, num_threads);
            delete processor;
        }
        if (db_output) {
            free_writer(handle);
        } else if (save_as_tar) {
            mtar_write_finalize(&tar_out);
            mtar_close(&tar_out);
        }
    } else if (mode == CHECK) {
        if (inputs.size() == 1) {
            std::cout << "Checking " << input << std::endl;
        } else {
            std::cout << "Checking files in " << input;
            std::cout << " using " << num_threads << " threads" << std::endl;
        }

        for (size_t i = 0; i < inputs.size() + 1; i++) {
            const std::string& input = (i == inputs.size()) ? "" : inputs[i];
            Processor* processor;
            struct stat st;
            if (i != inputs.size()) {
                if (stringEndsWith(".tar", input) || stringEndsWith(".tar.gz", input) || stringEndsWith(".tgz", input)) {
                    processor = new TarProcessor(input);
                }
                else if (stat((input + ".dbtype").c_str(), &st) == 0) {
                    if (user_id_file.size() > 0) {
                        processor = new DatabaseProcessor(input, user_id_file, id_mode, use_cache);
                    }
                    else {
                        processor = new DatabaseProcessor(input);
                    }
                }
#ifdef HAVE_GCS
                else if (stringStartsWith("gcs://", input)) {
                    processor = new GcsProcessor(input);
                }
#endif
                else {
                    processor = new DirectoryProcessor(input, recursive);
                }
            } else {
                if (single_file_inputs.size() > 0) {
                    processor = new DirectoryProcessor(single_file_inputs);
                }
                else {
                    continue;
                }
            }

            process_entry_func func = [&](const char* name, const char* dataBuffer, size_t size) -> bool {
                std::vector<ContainerFragment> containerFragments;
                std::string containerTitle;
                bool isContainer = readContainer(dataBuffer, size, containerTitle, containerFragments);
                if (isContainer) {
                    for (size_t fragmentIndex = 0; fragmentIndex < containerFragments.size(); fragmentIndex++) {
                        Foldcomp compRes;
                        std::istringstream fragmentInput(containerFragments[fragmentIndex].payload);
                        int flag = compRes.read(fragmentInput);
                        if (flag != 0) {
                            std::cerr << "[Error] Failed to read a container fragment during check." << std::endl;
                            return false;
                        }
                        ValidityError err = compRes.checkValidity();
                        std::string sname(name);
                        sname += makeFragmentSuffix(containerFragments[fragmentIndex], fragmentIndex);
                        printValidityError(err, sname);
                    }
                    return true;
                } else {
                    std::istringstream input(std::string(dataBuffer, size));
                    Foldcomp compRes;
                    int flag = compRes.read(input);
                    if (flag != 0) {
                        if (flag == -1) {
                            std::cerr << "[Error] File is not a valid fcz file" << std::endl;
                        } else if (flag == -2) {
                            std::cerr << "[Error] Could not restore prevAtoms" << std::endl;
                        } else {
                            std::cerr << "[Error] Unknown read error" << std::endl;
                        }
                        return false;
                    }
                    ValidityError err = compRes.checkValidity();
                    std::string sname(name);
                    printValidityError(err, sname);
                    return true;
                }
            };
            processor->run(func, num_threads);
            delete processor;
        }
    }
    return EXIT_SUCCESS;
}
