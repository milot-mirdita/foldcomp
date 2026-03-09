#include "amino_acid.h"
#include "atom_coordinate.h"
#include "database_reader.h"
#include "foldcomp.h"
#include "structure_reader.h"
#include "utility.h"

#include <algorithm>
#include <cstring>
#include <filesystem>
#include <iostream>
#include <map>
#include <set>
#include <sstream>
#include <string>
#include <tuple>
#include <vector>

namespace {

struct ResidueStats {
    std::set<std::string> atoms;
    std::string residue;
};

struct Summary {
    size_t atomCount = 0;
    size_t residueCount = 0;
    size_t completeBackboneResidueCount = 0;
    size_t caOnlyResidueCount = 0;
};

using ResidueKey = std::tuple<int, std::string, int, int>;
namespace fs = std::filesystem;

Summary summarize(const std::vector<AtomCoordinate>& atoms);
std::map<ResidueKey, ResidueStats> groupResidues(const std::vector<AtomCoordinate>& atoms);

static constexpr uint8_t FOLDCOMP_FORMAT_VERSION_CONTAINER = 2;
static constexpr uint8_t FOLDCOMP_FORMAT_FLAG_CONTAINER = 1;
static constexpr uint8_t CONTAINER_FRAGMENT_KIND_FCZ = 0;
static constexpr uint8_t CONTAINER_FRAGMENT_KIND_RAW_ATOMS = 2;

struct ContainerFragment {
    uint8_t kind = CONTAINER_FRAGMENT_KIND_FCZ;
    int16_t model = 1;
    std::string chain;
    std::string payload;
};

enum class RawReason {
    NonProtein = 0,
    CaOnlyResidue,
    IncompleteBackboneResidue,
    IsolatedCompleteResidue,
    UnknownResidue,
    PartialCanonicalAtomSet,
    ExtraNonCanonicalAtoms,
};

struct FallbackStats {
    size_t fileCount = 0;
    size_t proteinFileCount = 0;
    size_t rawFileCount = 0;
    size_t residueCount = 0;
    std::map<RawReason, size_t> rawResiduesByReason;
    std::map<RawReason, size_t> rawFilesByReason;
};

const char* rawReasonName(RawReason reason) {
    switch (reason) {
        case RawReason::NonProtein: return "non_protein";
        case RawReason::CaOnlyResidue: return "ca_only_residue";
        case RawReason::IncompleteBackboneResidue: return "incomplete_backbone_residue";
        case RawReason::IsolatedCompleteResidue: return "isolated_complete_residue";
        case RawReason::UnknownResidue: return "unknown_residue";
        case RawReason::PartialCanonicalAtomSet: return "partial_canonical_atom_set";
        case RawReason::ExtraNonCanonicalAtoms: return "extra_noncanonical_atoms";
    }
    return "unknown";
}

template <typename T>
bool readBinaryValue(const char* data, size_t size, size_t& offset, T& out) {
    if (offset + sizeof(T) > size) {
        return false;
    }
    memcpy(&out, data + offset, sizeof(T));
    offset += sizeof(T);
    return true;
}

bool hasFoldcompMagic(const char* data, size_t size) {
    if (size < MAGICNUMBER_LENGTH) {
        return false;
    }
    return memcmp(data, MAGICNUMBER, MAGICNUMBER_LENGTH) == 0;
}

bool isContainerHeader(const CompressedFileHeader& header) {
    bool supportedVersion = header.version == 1 || header.version == FOLDCOMP_FORMAT_VERSION_CONTAINER;
    return supportedVersion &&
           (header.flags & FOLDCOMP_FORMAT_FLAG_CONTAINER) != 0 &&
           header.nResidue == 0 &&
           header.nAtom == 0;
}

bool readContainer(
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
        if (header.version >= FOLDCOMP_FORMAT_VERSION_CONTAINER) {
            if (!readBinaryValue(data, size, offset, fragment.kind)) {
                return false;
            }
        }
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

bool decodeEntryToAtoms(const char* dataBuffer, size_t size, std::vector<AtomCoordinate>& atoms) {
    atoms.clear();
    std::string containerTitle;
    std::vector<ContainerFragment> containerFragments;
    bool isContainer = readContainer(dataBuffer, size, containerTitle, containerFragments);
    if (isContainer) {
        for (const auto& fragment : containerFragments) {
            std::vector<AtomCoordinate> fragmentAtoms;
            if (fragment.kind == CONTAINER_FRAGMENT_KIND_RAW_ATOMS) {
                if (!deserializeAtomCoordinates(
                        fragment.payload.data(), fragment.payload.size(), fragmentAtoms)) {
                    return false;
                }
            } else if (fragment.kind == CONTAINER_FRAGMENT_KIND_FCZ) {
                Foldcomp compRes;
                std::istringstream fragmentInput(fragment.payload);
                int flag = compRes.read(fragmentInput);
                if (flag != 0) {
                    return false;
                }
                flag = compRes.decompress(fragmentAtoms);
                if (flag != 0) {
                    return false;
                }
            } else {
                return false;
            }
            for (auto& atom : fragmentAtoms) {
                atom.model = fragment.model;
                atom.chain = fragment.chain;
                atoms.push_back(atom);
            }
        }
        return true;
    }

    Foldcomp compRes;
    std::istringstream input(std::string(dataBuffer, size));
    int flag = compRes.read(input);
    if (flag != 0) {
        return false;
    }
    return compRes.decompress(atoms) == 0;
}

std::string sourcePathToDbKey(const std::string& path) {
    std::string base = baseName(path);
    std::pair<std::string, std::string> outputParts = getFileParts(base);
    std::string dbKey = baseName(outputParts.first);
    std::replace(dbKey.begin(), dbKey.end(), '.', '_');
    return dbKey;
}

std::vector<RawReason> classifyResidue(
    const std::vector<AtomCoordinate>& residueAtoms,
    bool isolatedCompleteResidue
) {
    std::vector<RawReason> reasons;
    std::set<std::string> atomNames;
    for (const auto& atom : residueAtoms) {
        atomNames.insert(atom.atom);
    }

    bool hasN = atomNames.count("N") != 0;
    bool hasCA = atomNames.count("CA") != 0;
    bool hasC = atomNames.count("C") != 0;
    if (atomNames.size() == 1 && hasCA) {
        reasons.push_back(RawReason::CaOnlyResidue);
    } else if (!(hasN && hasCA && hasC)) {
        reasons.push_back(RawReason::IncompleteBackboneResidue);
    }
    if (isolatedCompleteResidue) {
        reasons.push_back(RawReason::IsolatedCompleteResidue);
    }

    auto aaIt = Foldcomp::AAS.find(residueAtoms[0].residue);
    if (aaIt == Foldcomp::AAS.end()) {
        reasons.push_back(RawReason::UnknownResidue);
        return reasons;
    }
    std::set<std::string> canonicalAtoms(aaIt->second.atoms.begin(), aaIt->second.atoms.end());
    std::set<std::string> canonicalAtomsWithOxt = canonicalAtoms;
    canonicalAtomsWithOxt.insert("OXT");
    bool hasExtraNonCanonicalAtoms = false;
    for (const auto& atom : atomNames) {
        if (canonicalAtomsWithOxt.count(atom) == 0) {
            hasExtraNonCanonicalAtoms = true;
        }
    }
    if (hasExtraNonCanonicalAtoms) {
        reasons.push_back(RawReason::ExtraNonCanonicalAtoms);
    } else if (atomNames != canonicalAtoms && atomNames != canonicalAtomsWithOxt) {
        reasons.push_back(RawReason::PartialCanonicalAtomSet);
    }
    return reasons;
}

std::map<ResidueKey, std::vector<RawReason>> classifyRawResidues(const std::vector<AtomCoordinate>& atoms) {
    std::map<ResidueKey, std::vector<RawReason>> output;
    if (atoms.empty()) {
        return output;
    }
    std::map<std::tuple<int, std::string, int>, int> occurrences;

    std::vector<std::pair<size_t, size_t>> chainRegions = identifyChains(atoms);
    for (const auto& chainRegion : chainRegions) {
        std::vector<std::pair<size_t, size_t>> discontinuous = identifyDiscontinousResInd(
            atoms, chainRegion.first, chainRegion.second
        );
        if (discontinuous.empty()) {
            discontinuous.push_back(chainRegion);
        }
        for (const auto& fragment : discontinuous) {
            tcb::span<AtomCoordinate> fragmentSpan(
                const_cast<AtomCoordinate*>(&atoms[fragment.first]),
                const_cast<AtomCoordinate*>(atoms.data()) + fragment.second
            );
            std::vector<BackboneRegion> regions = identifyBackboneRegions(fragmentSpan);
            for (const auto& region : regions) {
                if (region.encodable) {
                    bool needsRaw = false;
                    std::vector<std::vector<AtomCoordinate>> residues = splitAtomByResidue(
                        tcb::span<AtomCoordinate>(fragmentSpan.data() + region.start,
                                                  fragmentSpan.data() + region.end)
                    );
                    for (const auto& residue : residues) {
                        auto reasons = classifyResidue(residue, false);
                        if (!reasons.empty()) {
                            needsRaw = true;
                            break;
                        }
                    }
                    if (!needsRaw) {
                        continue;
                    }
                }

                std::vector<std::vector<AtomCoordinate>> residues = splitAtomByResidue(
                    tcb::span<AtomCoordinate>(fragmentSpan.data() + region.start,
                                              fragmentSpan.data() + region.end)
                );
                for (const auto& residue : residues) {
                    bool atomicallyComplete = false;
                    {
                        std::set<std::string> names;
                        for (const auto& atom : residue) {
                            names.insert(atom.atom);
                        }
                        atomicallyComplete = names.count("N") && names.count("CA") && names.count("C");
                    }
                    bool isolatedComplete = region.encodable == false && atomicallyComplete;
                    auto reasons = classifyResidue(residue, isolatedComplete);
                    if (reasons.empty()) {
                        reasons.push_back(RawReason::PartialCanonicalAtomSet);
                    }
                    std::tuple<int, std::string, int> baseKey = std::make_tuple(
                        residue[0].model, residue[0].chain, residue[0].residue_index
                    );
                    ResidueKey key = std::make_tuple(
                        residue[0].model, residue[0].chain, residue[0].residue_index,
                        occurrences[baseKey]++
                    );
                    output[key] = reasons;
                }
            }
        }
    }
    return output;
}

int commandFallbackReport(const std::string& input) {
    FallbackStats stats;
    std::vector<std::string> files;
    if (fs::is_directory(input)) {
        for (const auto& entry : fs::directory_iterator(input)) {
            if (entry.is_regular_file()) {
                files.push_back(entry.path().string());
            }
        }
        std::sort(files.begin(), files.end());
    } else {
        files.push_back(input);
    }

    for (const auto& file : files) {
        stats.fileCount++;
        StructureReader reader;
        if (!reader.load(file)) {
            continue;
        }
        std::vector<AtomCoordinate> atoms;
        reader.readAllAtoms(atoms);
        removeAlternativePosition(atoms);
        if (atoms.empty()) {
            stats.rawFilesByReason[RawReason::NonProtein]++;
            continue;
        }
        stats.proteinFileCount++;
        auto rawResidues = classifyRawResidues(atoms);
        auto groupedResidues = groupResidues(atoms);
        stats.residueCount += groupedResidues.size();
        if (!rawResidues.empty()) {
            stats.rawFileCount++;
        }
        std::set<RawReason> fileReasons;
        for (const auto& kv : rawResidues) {
            for (RawReason reason : kv.second) {
                stats.rawResiduesByReason[reason]++;
                fileReasons.insert(reason);
            }
        }
        for (RawReason reason : fileReasons) {
            stats.rawFilesByReason[reason]++;
        }
    }

    std::cout << "files_total=" << stats.fileCount
              << " protein_files=" << stats.proteinFileCount
              << " raw_files=" << stats.rawFileCount
              << " residues_total=" << stats.residueCount
              << "\n";
    for (RawReason reason : {
             RawReason::NonProtein,
             RawReason::CaOnlyResidue,
             RawReason::IncompleteBackboneResidue,
             RawReason::IsolatedCompleteResidue,
             RawReason::UnknownResidue,
             RawReason::PartialCanonicalAtomSet,
             RawReason::ExtraNonCanonicalAtoms
         }) {
        size_t rawResidues = stats.rawResiduesByReason[reason];
        size_t rawFiles = stats.rawFilesByReason[reason];
        double residueRatio = stats.residueCount == 0 ? 0.0 : static_cast<double>(rawResidues) / static_cast<double>(stats.residueCount);
        double fileRatio = stats.proteinFileCount == 0 ? 0.0 : static_cast<double>(rawFiles) / static_cast<double>(stats.proteinFileCount);
        std::cout << rawReasonName(reason)
                  << " raw_files=" << rawFiles
                  << " file_ratio=" << fileRatio
                  << " raw_residues=" << rawResidues
                  << " residue_ratio=" << residueRatio
                  << "\n";
    }
    return 0;
}

Summary summarize(const std::vector<AtomCoordinate>& atoms) {
    std::map<ResidueKey, ResidueStats> residues = groupResidues(atoms);

    Summary summary;
    summary.atomCount = atoms.size();
    summary.residueCount = residues.size();
    for (const auto& kv : residues) {
        const auto& residue = kv.second;
        bool hasN = residue.atoms.count("N") != 0;
        bool hasCA = residue.atoms.count("CA") != 0;
        bool hasC = residue.atoms.count("C") != 0;
        if (hasN && hasCA && hasC) {
            summary.completeBackboneResidueCount++;
        }
        if (residue.atoms.size() == 1 && hasCA) {
            summary.caOnlyResidueCount++;
        }
    }
    return summary;
}

void printSummary(const std::string& label, const Summary& summary) {
    std::cout << label
              << " atoms=" << summary.atomCount
              << " residues=" << summary.residueCount
              << " complete_backbone=" << summary.completeBackboneResidueCount
              << " ca_only=" << summary.caOnlyResidueCount
              << "\n";
}

std::map<ResidueKey, ResidueStats> groupResidues(const std::vector<AtomCoordinate>& atoms) {
    std::map<ResidueKey, ResidueStats> residues;
    if (atoms.empty()) {
        return residues;
    }

    std::map<std::tuple<int, std::string, int>, int> occurrences;
    size_t residueStart = 0;
    for (size_t i = 1; i <= atoms.size(); i++) {
        bool endOfResidue = (i == atoms.size()) ||
                            atoms[i].model != atoms[i - 1].model ||
                            atoms[i].chain != atoms[i - 1].chain ||
                            atoms[i].residue_index != atoms[i - 1].residue_index;
        if (!endOfResidue) {
            continue;
        }

        const AtomCoordinate& firstAtom = atoms[residueStart];
        std::tuple<int, std::string, int> baseKey = std::make_tuple(
            firstAtom.model, firstAtom.chain, firstAtom.residue_index
        );
        ResidueKey key = std::make_tuple(
            firstAtom.model, firstAtom.chain, firstAtom.residue_index,
            occurrences[baseKey]++
        );
        ResidueStats& stats = residues[key];
        stats.residue = firstAtom.residue;
        for (size_t j = residueStart; j < i; j++) {
            stats.atoms.insert(atoms[j].atom);
        }
        residueStart = i;
    }
    return residues;
}

std::string formatResidueKey(const ResidueKey& key, const ResidueStats& residue) {
    std::ostringstream oss;
    oss << "model=" << std::get<0>(key)
        << " chain=" << std::get<1>(key)
        << " resid=" << std::get<2>(key)
        << " occ=" << std::get<3>(key)
        << " resn=" << residue.residue;
    return oss.str();
}

std::string joinAtomNames(const std::set<std::string>& atoms) {
    std::ostringstream oss;
    bool first = true;
    for (const auto& atom : atoms) {
        if (!first) {
            oss << ",";
        }
        oss << atom;
        first = false;
    }
    return oss.str();
}

void printResidueStream(const std::vector<AtomCoordinate>& atoms) {
    if (atoms.empty()) {
        return;
    }
    std::map<std::tuple<int, std::string, int>, int> occurrences;
    size_t residueStart = 0;
    for (size_t i = 1; i <= atoms.size(); i++) {
        bool endOfResidue = (i == atoms.size()) ||
                            atoms[i].model != atoms[i - 1].model ||
                            atoms[i].chain != atoms[i - 1].chain ||
                            atoms[i].residue_index != atoms[i - 1].residue_index;
        if (!endOfResidue) {
            continue;
        }
        const AtomCoordinate& atom = atoms[residueStart];
        std::tuple<int, std::string, int> baseKey = std::make_tuple(atom.model, atom.chain, atom.residue_index);
        int occ = occurrences[baseKey]++;
        std::cout << "residue model=" << atom.model
                  << " chain=" << atom.chain
                  << " resid=" << atom.residue_index
                  << " occ=" << occ
                  << " resn=" << atom.residue
                  << "\n";
        residueStart = i;
    }
}

int commandSummarize(const std::string& input) {
    StructureReader reader;
    if (!reader.load(input)) {
        std::cerr << "[Error] Failed to read " << input << "\n";
        return 1;
    }
    std::vector<AtomCoordinate> atoms;
    reader.readAllAtoms(atoms);
    removeAlternativePosition(atoms);
    printSummary(input, summarize(atoms));
    return 0;
}

int compareAtomVectors(const std::vector<AtomCoordinate>& sourceAtoms,
                       const std::vector<AtomCoordinate>& roundtripAtoms) {
    Summary sourceSummary = summarize(sourceAtoms);
    Summary roundtripSummary = summarize(roundtripAtoms);
    printSummary("source", sourceSummary);
    printSummary("roundtrip", roundtripSummary);

    auto sourceResidues = groupResidues(sourceAtoms);
    auto roundtripResidues = groupResidues(roundtripAtoms);

    size_t missingResidues = 0;
    for (const auto& kv : sourceResidues) {
        if (roundtripResidues.find(kv.first) == roundtripResidues.end()) {
            std::cout << "missing_residue " << formatResidueKey(kv.first, kv.second) << "\n";
            missingResidues++;
        }
    }
    size_t extraResidues = 0;
    for (const auto& kv : roundtripResidues) {
        if (sourceResidues.find(kv.first) == sourceResidues.end()) {
            std::cout << "extra_residue " << formatResidueKey(kv.first, kv.second) << "\n";
            extraResidues++;
        }
    }

    size_t changedResidues = 0;
    for (const auto& kv : sourceResidues) {
        auto it = roundtripResidues.find(kv.first);
        if (it == roundtripResidues.end()) {
            continue;
        }
        if (kv.second.atoms != it->second.atoms) {
            std::set<std::string> missingAtoms;
            std::set<std::string> extraAtoms;
            std::set_difference(kv.second.atoms.begin(), kv.second.atoms.end(),
                                it->second.atoms.begin(), it->second.atoms.end(),
                                std::inserter(missingAtoms, missingAtoms.begin()));
            std::set_difference(it->second.atoms.begin(), it->second.atoms.end(),
                                kv.second.atoms.begin(), kv.second.atoms.end(),
                                std::inserter(extraAtoms, extraAtoms.begin()));
            std::cout << "changed_residue " << formatResidueKey(kv.first, kv.second)
                      << " source_atoms=" << kv.second.atoms.size()
                      << " roundtrip_atoms=" << it->second.atoms.size()
                      << " roundtrip_resn=" << it->second.residue
                      << " missing_atom_names=" << joinAtomNames(missingAtoms)
                      << " extra_atom_names=" << joinAtomNames(extraAtoms)
                      << "\n";
            changedResidues++;
        }
    }

    std::cout << "missing_residue_count=" << missingResidues
              << " extra_residue_count=" << extraResidues
              << " changed_residue_count=" << changedResidues
              << "\n";
    return 0;
}

int commandCompare(const std::string& source, const std::string& roundtrip) {
    StructureReader sourceReader;
    StructureReader roundtripReader;
    if (!sourceReader.load(source)) {
        std::cerr << "[Error] Failed to read source " << source << "\n";
        return 1;
    }
    if (!roundtripReader.load(roundtrip)) {
        std::cerr << "[Error] Failed to read roundtrip " << roundtrip << "\n";
        return 1;
    }

    std::vector<AtomCoordinate> sourceAtoms;
    std::vector<AtomCoordinate> roundtripAtoms;
    sourceReader.readAllAtoms(sourceAtoms);
    roundtripReader.readAllAtoms(roundtripAtoms);
    removeAlternativePosition(sourceAtoms);
    removeAlternativePosition(roundtripAtoms);
    return compareAtomVectors(sourceAtoms, roundtripAtoms);
}

int commandDbEntryCompare(const std::string& sourcePath, const std::string& dbPath) {
    StructureReader sourceReader;
    if (!sourceReader.load(sourcePath)) {
        std::cerr << "[Error] Failed to read source " << sourcePath << "\n";
        return 1;
    }

    std::vector<AtomCoordinate> sourceAtoms;
    sourceReader.readAllAtoms(sourceAtoms);
    removeAlternativePosition(sourceAtoms);
    if (sourceAtoms.empty()) {
        std::cerr << "[Error] No protein atoms found in source " << sourcePath << "\n";
        return 1;
    }

    void* reader = make_reader(dbPath.c_str(), (dbPath + ".index").c_str(),
                               DB_READER_USE_DATA | DB_READER_USE_LOOKUP);
    if (reader == NULL) {
        std::cerr << "[Error] Failed to open database " << dbPath << "\n";
        return 1;
    }

    std::string dbKeyStr = sourcePathToDbKey(sourcePath);
    uint32_t key = reader_lookup_entry(reader, dbKeyStr.c_str());
    if (key == UINT32_MAX) {
        std::cerr << "[Error] Missing db entry for key " << dbKeyStr << "\n";
        free_reader(reader);
        return 1;
    }
    int64_t id = reader_get_id(reader, key);
    if (id == -1) {
        std::cerr << "[Error] Missing db id for key " << dbKeyStr << "\n";
        free_reader(reader);
        return 1;
    }

    std::vector<AtomCoordinate> dbAtoms;
    const char* data = reader_get_data(reader, id);
    int64_t length = reader_get_length(reader, id);
    if (data == NULL || length <= 0 || !decodeEntryToAtoms(data, static_cast<size_t>(length), dbAtoms)) {
        std::cerr << "[Error] Failed to decode db entry for key " << dbKeyStr << "\n";
        free_reader(reader);
        return 1;
    }
    free_reader(reader);

    removeAlternativePosition(dbAtoms);
    return compareAtomVectors(sourceAtoms, dbAtoms);
}

int commandDumpResidues(const std::string& sourcePath) {
    StructureReader reader;
    if (!reader.load(sourcePath)) {
        std::cerr << "[Error] Failed to read source " << sourcePath << "\n";
        return 1;
    }
    std::vector<AtomCoordinate> atoms;
    reader.readAllAtoms(atoms);
    removeAlternativePosition(atoms);
    printResidueStream(atoms);
    return 0;
}

int commandDbEntryDumpResidues(const std::string& sourcePath, const std::string& dbPath) {
    void* reader = make_reader(dbPath.c_str(), (dbPath + ".index").c_str(),
                               DB_READER_USE_DATA | DB_READER_USE_LOOKUP);
    if (reader == NULL) {
        std::cerr << "[Error] Failed to open database " << dbPath << "\n";
        return 1;
    }
    std::string dbKeyStr = sourcePathToDbKey(sourcePath);
    uint32_t key = reader_lookup_entry(reader, dbKeyStr.c_str());
    int64_t id = (key == UINT32_MAX) ? -1 : reader_get_id(reader, key);
    std::vector<AtomCoordinate> atoms;
    const char* data = (id == -1) ? NULL : reader_get_data(reader, id);
    int64_t length = (id == -1) ? -1 : reader_get_length(reader, id);
    if (id == -1 || data == NULL || length <= 0 ||
        !decodeEntryToAtoms(data, static_cast<size_t>(length), atoms)) {
        std::cerr << "[Error] Failed to decode db entry for key " << dbKeyStr << "\n";
        free_reader(reader);
        return 1;
    }
    free_reader(reader);
    removeAlternativePosition(atoms);
    printResidueStream(atoms);
    return 0;
}

int commandFragmentReport(const std::string& sourcePath) {
    StructureReader reader;
    if (!reader.load(sourcePath)) {
        std::cerr << "[Error] Failed to read source " << sourcePath << "\n";
        return 1;
    }
    std::vector<AtomCoordinate> atoms;
    reader.readAllAtoms(atoms);
    removeAlternativePosition(atoms);
    if (atoms.empty()) {
        std::cerr << "[Error] No protein atoms found in " << sourcePath << "\n";
        return 1;
    }

    std::vector<std::pair<size_t, size_t>> chainRegions = identifyChains(atoms);
    for (size_t chainIdx = 0; chainIdx < chainRegions.size(); chainIdx++) {
        const auto& chainRegion = chainRegions[chainIdx];
        std::cout << "chain_region idx=" << chainIdx
                  << " start_atom=" << chainRegion.first
                  << " end_atom=" << chainRegion.second
                  << " start_resid=" << atoms[chainRegion.first].residue_index
                  << " end_resid=" << atoms[chainRegion.second - 1].residue_index
                  << " chain=" << atoms[chainRegion.first].chain
                  << "\n";
        std::vector<std::pair<size_t, size_t>> fragments = identifyDiscontinousResInd(
            atoms, chainRegion.first, chainRegion.second
        );
        if (fragments.empty()) {
            fragments.push_back(chainRegion);
        }
        for (size_t fragIdx = 0; fragIdx < fragments.size(); fragIdx++) {
            const auto& fragment = fragments[fragIdx];
            tcb::span<AtomCoordinate> chainSpan(
                &atoms[fragment.first], atoms.data() + fragment.second
            );
            std::vector<BackboneRegion> regions = identifyBackboneRegions(chainSpan);
            std::cout << "fragment idx=" << fragIdx
                      << " start_atom=" << fragment.first
                      << " end_atom=" << fragment.second
                      << " start_resid=" << atoms[fragment.first].residue_index
                      << " end_resid=" << atoms[fragment.second - 1].residue_index
                      << " regions=" << regions.size()
                      << "\n";
            for (size_t regionIdx = 0; regionIdx < regions.size(); regionIdx++) {
                const auto& region = regions[regionIdx];
                size_t startAtom = fragment.first + region.start;
                size_t endAtom = fragment.first + region.end - 1;
                std::cout << "region idx=" << regionIdx
                          << " encodable=" << (region.encodable ? 1 : 0)
                          << " start_atom=" << startAtom
                          << " end_atom=" << endAtom
                          << " start_resid=" << atoms[startAtom].residue_index
                          << " end_resid=" << atoms[endAtom].residue_index
                          << "\n";
            }
        }
    }
    return 0;
}

int commandContainerReport(const std::string& encodedPath) {
    std::ifstream input(encodedPath, std::ios::binary);
    if (!input) {
        std::cerr << "[Error] Failed to open " << encodedPath << "\n";
        return 1;
    }
    std::string data((std::istreambuf_iterator<char>(input)), std::istreambuf_iterator<char>());
    std::string title;
    std::vector<ContainerFragment> fragments;
    bool isContainer = readContainer(data.data(), data.size(), title, fragments);
    if (!isContainer) {
        Foldcomp compRes;
        std::istringstream singleInput(data);
        if (compRes.read(singleInput) != 0) {
            std::cerr << "[Error] Failed to read FCZ " << encodedPath << "\n";
            return 1;
        }
        std::cout << "single_fcz"
                  << " idxResidue=" << compRes.header.idxResidue
                  << " nResidue=" << compRes.header.nResidue
                  << " idxAtom=" << compRes.header.idxAtom
                  << " nAtom=" << compRes.header.nAtom
                  << " nAnchor=" << compRes.header.nAnchor
                  << " nSideChainTorsion=" << compRes.header.nSideChainTorsion
                  << " chain=" << getChainName(compRes.header)
                  << "\n";
        std::cout << "anchor_indices";
        for (int anchorIndex : compRes.anchorIndices) {
            std::cout << " " << anchorIndex;
        }
        std::cout << "\n";
        return 0;
    }

    std::cout << "container title=" << title
              << " fragments=" << fragments.size()
              << "\n";
    for (size_t i = 0; i < fragments.size(); i++) {
        const auto& fragment = fragments[i];
        std::cout << "fragment idx=" << i
                  << " kind=" << static_cast<int>(fragment.kind)
                  << " model=" << fragment.model
                  << " chain=" << fragment.chain
                  << " payload_size=" << fragment.payload.size();
        if (fragment.kind == CONTAINER_FRAGMENT_KIND_RAW_ATOMS) {
            std::vector<AtomCoordinate> atoms;
            if (!deserializeAtomCoordinates(fragment.payload.data(), fragment.payload.size(), atoms)) {
                std::cout << " raw_decode_error=1\n";
                continue;
            }
            std::cout << " raw_atoms=" << atoms.size();
            if (!atoms.empty()) {
                std::cout << " start_resid=" << atoms.front().residue_index
                          << " end_resid=" << atoms.back().residue_index;
            }
            std::cout << "\n";
            continue;
        }
        if (fragment.kind != CONTAINER_FRAGMENT_KIND_FCZ) {
            std::cout << " unsupported_kind=" << static_cast<int>(fragment.kind) << "\n";
            continue;
        }

        Foldcomp compRes;
        std::istringstream fragmentInput(fragment.payload);
        if (compRes.read(fragmentInput) != 0) {
            std::cout << " read_error=1\n";
            continue;
        }
        std::cout << " idxResidue=" << compRes.header.idxResidue
                  << " nResidue=" << compRes.header.nResidue
                  << " idxAtom=" << compRes.header.idxAtom
                  << " nAtom=" << compRes.header.nAtom
                  << " nAnchor=" << compRes.header.nAnchor
                  << " nSideChainTorsion=" << compRes.header.nSideChainTorsion
                  << " fcz_chain=" << getChainName(compRes.header)
                  << "\n";
        std::cout << "anchor_indices idx=" << i;
        for (int anchorIndex : compRes.anchorIndices) {
            std::cout << " " << anchorIndex;
        }
        std::cout << "\n";
    }
    return 0;
}

int commandDbCompare(const std::string& sourceDir, const std::string& dbPath) {
    std::vector<std::string> files = getFilesInDirectory(sourceDir, false);
    std::sort(files.begin(), files.end());

    void* reader = make_reader(dbPath.c_str(), (dbPath + ".index").c_str(),
                               DB_READER_USE_DATA | DB_READER_USE_LOOKUP);
    if (reader == NULL) {
        std::cerr << "[Error] Failed to open database " << dbPath << "\n";
        return 1;
    }

    size_t totalFiles = 0;
    size_t proteinFiles = 0;
    size_t matchedFiles = 0;
    size_t missingEntries = 0;
    size_t decodeFailures = 0;
    size_t filesWithDiffs = 0;
    size_t totalMissingResidues = 0;
    size_t totalExtraResidues = 0;
    size_t totalChangedResidues = 0;

    for (const auto& file : files) {
        totalFiles++;
        StructureReader sourceReader;
        if (!sourceReader.load(file)) {
            continue;
        }
        std::vector<AtomCoordinate> sourceAtoms;
        sourceReader.readAllAtoms(sourceAtoms);
        removeAlternativePosition(sourceAtoms);
        if (sourceAtoms.empty()) {
            continue;
        }
        proteinFiles++;

        std::string dbKeyStr = sourcePathToDbKey(file);
        uint32_t key = reader_lookup_entry(reader, dbKeyStr.c_str());
        if (key == UINT32_MAX) {
            std::cout << "missing_db_entry " << baseName(file) << " key=" << dbKeyStr << "\n";
            missingEntries++;
            continue;
        }
        int64_t id = reader_get_id(reader, key);
        if (id == -1) {
            std::cout << "missing_db_id " << baseName(file) << " key=" << dbKeyStr << "\n";
            missingEntries++;
            continue;
        }

        std::vector<AtomCoordinate> dbAtoms;
        const char* data = reader_get_data(reader, id);
        int64_t length = reader_get_length(reader, id);
        if (data == NULL || length <= 0 || !decodeEntryToAtoms(data, static_cast<size_t>(length), dbAtoms)) {
            std::cout << "decode_failure " << baseName(file) << " key=" << dbKeyStr << "\n";
            decodeFailures++;
            continue;
        }
        removeAlternativePosition(dbAtoms);
        matchedFiles++;

        auto sourceResidues = groupResidues(sourceAtoms);
        auto dbResidues = groupResidues(dbAtoms);
        size_t missingResidues = 0;
        size_t extraResidues = 0;
        size_t changedResidues = 0;
        for (const auto& kv : sourceResidues) {
            auto it = dbResidues.find(kv.first);
            if (it == dbResidues.end()) {
                missingResidues++;
            } else if (kv.second.atoms != it->second.atoms) {
                changedResidues++;
            }
        }
        for (const auto& kv : dbResidues) {
            if (sourceResidues.find(kv.first) == sourceResidues.end()) {
                extraResidues++;
            }
        }

        totalMissingResidues += missingResidues;
        totalExtraResidues += extraResidues;
        totalChangedResidues += changedResidues;
        if (missingResidues != 0 || extraResidues != 0 || changedResidues != 0) {
            filesWithDiffs++;
            std::cout << "diff_file " << baseName(file)
                      << " missing_residue_count=" << missingResidues
                      << " extra_residue_count=" << extraResidues
                      << " changed_residue_count=" << changedResidues
                      << "\n";
        }
    }

    free_reader(reader);
    std::cout << "files_total=" << totalFiles
              << " protein_files=" << proteinFiles
              << " matched_files=" << matchedFiles
              << " missing_entries=" << missingEntries
              << " decode_failures=" << decodeFailures
              << " files_with_diffs=" << filesWithDiffs
              << " total_missing_residues=" << totalMissingResidues
              << " total_extra_residues=" << totalExtraResidues
              << " total_changed_residues=" << totalChangedResidues
              << "\n";
    return 0;
}

}  // namespace

int main(int argc, char** argv) {
    if (argc < 3) {
        std::cerr << "Usage:\n"
                  << "  foldcomp_audit summarize <structure>\n"
                  << "  foldcomp_audit compare <source_structure> <roundtrip_structure>\n"
                  << "  foldcomp_audit fallback-report <structure_or_directory>\n"
                  << "  foldcomp_audit db-compare <source_directory> <db>\n"
                  << "  foldcomp_audit db-entry-compare <source_structure> <db>\n"
                  << "  foldcomp_audit dump-residues <structure>\n"
                  << "  foldcomp_audit db-entry-dump-residues <source_structure> <db>\n"
                  << "  foldcomp_audit fragment-report <structure>\n"
                  << "  foldcomp_audit container-report <encoded_structure>\n";
        return 1;
    }

    std::string command = argv[1];
    if (command == "summarize") {
        return commandSummarize(argv[2]);
    }
    if (command == "compare") {
        if (argc < 4) {
            std::cerr << "[Error] compare requires two input files\n";
            return 1;
        }
        return commandCompare(argv[2], argv[3]);
    }
    if (command == "fallback-report") {
        return commandFallbackReport(argv[2]);
    }
    if (command == "db-compare") {
        if (argc < 4) {
            std::cerr << "[Error] db-compare requires source_directory and db path\n";
            return 1;
        }
        return commandDbCompare(argv[2], argv[3]);
    }
    if (command == "db-entry-compare") {
        if (argc < 4) {
            std::cerr << "[Error] db-entry-compare requires source_structure and db path\n";
            return 1;
        }
        return commandDbEntryCompare(argv[2], argv[3]);
    }
    if (command == "dump-residues") {
        return commandDumpResidues(argv[2]);
    }
    if (command == "db-entry-dump-residues") {
        if (argc < 4) {
            std::cerr << "[Error] db-entry-dump-residues requires source_structure and db path\n";
            return 1;
        }
        return commandDbEntryDumpResidues(argv[2], argv[3]);
    }
    if (command == "fragment-report") {
        return commandFragmentReport(argv[2]);
    }
    if (command == "container-report") {
        return commandContainerReport(argv[2]);
    }

    std::cerr << "[Error] Unknown command: " << command << "\n";
    return 1;
}
