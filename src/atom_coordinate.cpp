/**
 * File: atom_coordinate.cpp
 * Project: foldcomp
 * Created: 2021-01-18 12:53:34
 * Author: Hyunbin Kim (khb7840@gmail.com)
 * Description:
 *     The data type to handle atom coordinate comes here.
 * ---
 * Last Modified: Fri Mar 03 2023
 * Modified By: Hyunbin Kim
 * ---
 * Copyright © 2021 Hyunbin Kim, All rights reserved
 */
#include "atom_coordinate.h"

#include <cstring>
#include <iomanip>
#include <iostream>
#include <sstream> // IWYU pragma: keep

namespace {

bool startsNewResidue(const AtomCoordinate& current, const AtomCoordinate& previous) {
    if (current.chain != previous.chain) {
        return true;
    }
    if (current.model != previous.model) {
        return true;
    }
    if (current.residue_index != previous.residue_index) {
        return true;
    }
    if (current.atom == "N" && previous.atom != "N") {
        // Reused residue indices can still denote a new residue in author numbering.
        return true;
    }
    return false;
}

}

/**
 * @brief Construct a new Atom Coordinate:: Atom Coordinate object
 *
 * @param a A string for atom name
 * @param r A string for residue name
 * @param ai An integer for atom index
 * @param ri An integer for residue index
 * @param x A float for x coordinate
 * @param y A float for y coordinate
 * @param z A float for z coordinate
 */
AtomCoordinate::AtomCoordinate(
    std::string a, std::string r, std::string c,
    int ai, int ri, float x, float y, float z,
    float occupancy, float tempFactor, int model
): atom(a), residue(r), chain(c), atom_index(ai), residue_index(ri), occupancy(occupancy), tempFactor(tempFactor), model(model) {
    this->coordinate = {x, y, z};
}

/**
 * @brief Construct a new Atom Coordinate:: Atom Coordinate object
 *
 * @param a A string for atom name
 * @param r A string for residue name
 * @param ai An integer for atom index
 * @param ri An integer for residue index
 * @param coord A float vector for x,y,z coordinates.
 */
AtomCoordinate::AtomCoordinate(
    std::string a, std::string r, std::string c,
    int ai, int ri, float3d coord,
    float occupancy, float tempFactor, int model
): atom(a), residue(r), chain(c), atom_index(ai), residue_index(ri), coordinate(coord), occupancy(occupancy), tempFactor(tempFactor), model(model) {
}

bool AtomCoordinate::operator==(const AtomCoordinate& other) const {
    return (
        (this->atom == other.atom) &&
        (this->atom_index == other.atom_index) &&
        (this->residue == other.residue) &&
        (this->residue_index == other.residue_index) &&
        (this->chain == other.chain) &&
        (this->model == other.model) &&
        (this->coordinate.x == other.coordinate.x) &&
        (this->coordinate.y == other.coordinate.y) &&
        (this->coordinate.z == other.coordinate.z)
    );
}
bool AtomCoordinate::operator!=(const AtomCoordinate& other) const {
    return !(*this == other);
}

bool AtomCoordinate::isBackbone() const {
    return ((this->atom == "N") ||(this->atom == "CA") ||(this->atom == "C"));
}

void AtomCoordinate::print(int option) const {
    std::cout << "Atom: " << this->atom << std::endl;
    if (option != 0) {
        std::cout << "Residue: " << this->residue << std::endl;
        std::cout << "Chain: " << this->chain << std::endl;
        std::cout << "Atom Index: " << this->atom_index << std::endl;
        std::cout << "Residue Index: " << this->residue_index << std::endl;
        if (option == 2) {
            std::cout << "Coordinate: ";
            std::cout << this->coordinate.x << " ";
            std::cout << this->coordinate.y << " ";
            std::cout << this->coordinate.z << " ";
            std::cout << std::endl;
        }
    }
}

/**
 * @brief Extracts coordinates from AtomCoordinate vector
 *
 * @param atoms A vector of AtomCoordinate
 * @return std::vector< std::vector<float> >
 */
std::vector<float3d> extractCoordinates(
    const std::vector<AtomCoordinate>& atoms
) {
    std::vector<float3d> output(atoms.size());
    for (size_t i = 0; i < atoms.size(); i++) {
        output[i] = atoms[i].coordinate;
    }
    return output;
}

std::vector<AtomCoordinate> extractChain(
    std::vector<AtomCoordinate>& atoms, std::string chain
) {
    std::vector<AtomCoordinate> output;
    int total = atoms.size();
    for (int i = 0; i < total; i++) {
        const AtomCoordinate& curr_atm = atoms[i];
        if (i < (total-1)) {
            const AtomCoordinate& next_atm = atoms[i + 1];
            if (next_atm.atom == curr_atm.atom) {
                continue;
            }
        }
        if (curr_atm.chain == chain) {
            output.push_back(curr_atm);
        }
    }
    return output;
}

void printAtomCoordinateVector(std::vector<AtomCoordinate>& atoms, int option) {
    for (const AtomCoordinate& curr_atm : atoms) {
        curr_atm.print(option);
    }
}

std::vector<AtomCoordinate> filterBackbone(const tcb::span<AtomCoordinate>& atoms) {
    std::vector<AtomCoordinate> output;
    std::vector<std::vector<AtomCoordinate>> atomByResidue = splitAtomByResidue(atoms);
    for (const auto& residueAtoms : atomByResidue) {
        const AtomCoordinate* n = nullptr;
        const AtomCoordinate* ca = nullptr;
        const AtomCoordinate* c = nullptr;
        for (const auto& atom : residueAtoms) {
            if (atom.atom == "N" && n == nullptr) {
                n = &atom;
            } else if (atom.atom == "CA" && ca == nullptr) {
                ca = &atom;
            } else if (atom.atom == "C" && c == nullptr) {
                c = &atom;
            }
        }
        if (n != nullptr && ca != nullptr && c != nullptr) {
            output.emplace_back(*n);
            output.emplace_back(*ca);
            output.emplace_back(*c);
        }
    }
    return output;
}

std::vector<AtomCoordinate> weightedAverage(
    const std::vector<AtomCoordinate>& origAtoms, const std::vector<AtomCoordinate>& revAtoms
) {
    std::vector<AtomCoordinate> output;
    output.reserve(origAtoms.size());
    int total = origAtoms.size();
    for (int i = 0; i < total; i++) {
        const AtomCoordinate& curr_atm = origAtoms[i];
        const AtomCoordinate&  rev_atm = revAtoms[i];
        output.emplace_back(
            curr_atm.atom, curr_atm.residue, curr_atm.chain,
            curr_atm.atom_index, curr_atm.residue_index,
            ((curr_atm.coordinate.x * (float)(total - i)) + (rev_atm.coordinate.x * (float)i)) / (float)total,
            ((curr_atm.coordinate.y * (float)(total - i)) + (rev_atm.coordinate.y * (float)i)) / (float)total,
            ((curr_atm.coordinate.z * (float)(total - i)) + (rev_atm.coordinate.z * (float)i)) / (float)total
        );
    }
    return output;
}

void reverse(char* s) {
    for (int i = 0, j = strlen(s)-1; i < j; i++, j--) {
        char c = s[i];
        s[i] = s[j];
        s[j] = c;
    }
}

void itoa_pos_only(int n, char* s) {
    int i = 0;
    do {
        // generate digits in reverse order
        // get next digit
        s[i++] = n % 10 + '0';
    // shift to next
    } while ((n /= 10) > 0);
    s[i] = '\0';
    reverse(s);
}

template <int32_t T, int32_t P>
void fast_ftoa(float n, char* s) {
    float rounded = n + ((n < 0) ? -(0.5f / T) : (0.5f / T));
    int32_t integer = (int32_t)rounded;
    int32_t decimal = (int32_t)((rounded - (float)integer) * (float)T);
    char* data = s;
    if (n < 0) {
        integer = std::abs(integer);
        decimal = std::abs(decimal);
        *data = '-';
        data++;
    }
    itoa_pos_only(integer, data);
    data += strlen(data);
    *data = '.';
    data++;
    char buffer[10];
    itoa_pos_only(decimal, buffer);
    int32_t len = strlen(buffer);
    for (int32_t i = 0; i < (P - len); i++) {
        *data = '0';
        data++;
    }
    memcpy(data, buffer, len);
    // add a null terminator
    data += len;
    *data = '\0';
    // std::string check(s);
    // std::ostringstream ss;
    // ss << std::fixed << std::setprecision(P) << n;
    // if (ss.str() != check) {
    //     std::cout << "ERROR: " << ss.str() << " != " << check << " ORIG: " << std::fixed << std::setprecision(10) << n << std::endl;
    // }
}

void writeAtomCoordinatesToPDB(
    std::vector<AtomCoordinate>& atoms, std::string title, std::ostream& pdb_stream
) {
    // Write title
    // Check if title is too long and if so, write the title in multiple lines
    if (title != "") {
        const char* headerData = title.c_str();
        size_t headerLen = title.length();
        int remainingHeader = headerLen;
        char buffer[128];
        int written = snprintf(buffer, sizeof(buffer), "TITLE     %.*s\n",  std::min(70, (int)remainingHeader), headerData);
        if (written >= 0 && written < (int)sizeof(buffer)) {
            pdb_stream << buffer;
        }
        remainingHeader -= 70;
        int continuation = 2;
        while (remainingHeader > 0) {
            written = snprintf(buffer, sizeof(buffer), "TITLE  % 3d%.*s\n", continuation, std::min(70, (int)remainingHeader), headerData + (headerLen - remainingHeader));
            if (written >= 0 && written < (int)sizeof(buffer)) {
                pdb_stream << buffer;
            }
            remainingHeader -= 70;
            continuation++;
        }
    }

    int total = atoms.size();
    std::string residue;
    for (int i = 0; i < total; i++) {
        pdb_stream << "ATOM  "; // 1-4 ATOM
        pdb_stream << std::setw(5) << atoms[i].atom_index; // 7-11
        pdb_stream << " "; // 12
        if (atoms[i].atom.size() == 4) {
            pdb_stream << std::setw(4) << std::left << atoms[i].atom; // 13-16
        } else {
            pdb_stream << " ";
            pdb_stream << std::setw(3) << std::left << atoms[i].atom; // 13-16
        }
        pdb_stream << " "; // 17
        pdb_stream << std::setw(3) << std::right << atoms[i].residue; // 18-20
        pdb_stream << " "; // 21
        char chainId = atoms[i].chain.empty() ? ' ' : atoms[i].chain[0];
        pdb_stream << chainId; // 22
        pdb_stream << std::setw(4) << atoms[i].residue_index; // 23-26
        pdb_stream << "    "; // 27-30
        char buffer[16];
        fast_ftoa<1000, 3>(atoms[i].coordinate.x, buffer);
        pdb_stream << std::setw(8) << buffer; // 31-38
        fast_ftoa<1000, 3>(atoms[i].coordinate.y, buffer);
        pdb_stream << std::setw(8) << buffer; // 39-46
        fast_ftoa<1000, 3>(atoms[i].coordinate.z, buffer);
        pdb_stream << std::setw(8) << buffer; // 47-54
        pdb_stream << "  1.00"; // 55-60
        fast_ftoa<100, 2>(atoms[i].tempFactor, buffer);
        pdb_stream << std::setw(6) << buffer; // 61-66
        pdb_stream << "          "; // 67-76
        // First one character from atom
        pdb_stream << std::setw(2) << atoms[i].atom[0]; // 77-78
        pdb_stream << "  \n"; // 79-80
        if (i == (total-1)) {
            // TER
            // 1-6 Record name "TER   "
            // 7-11 Atom serial number.
            // 18-20 Residue name.
            // 22 Chain identifier.
            // 23-26 Residue sequence number.
            pdb_stream << "TER   " << std::setw(5) << atoms[i].atom_index + 1 << "      ";
            pdb_stream << std::setw(3) << std::right << atoms[i].residue;
            pdb_stream << " " << chainId;
            pdb_stream << std::setw(4) << atoms[i].residue_index << std::endl;
        }
    }
}

int writeAtomCoordinatesToPDBFile(
    std::vector<AtomCoordinate>& atoms, std::string title, std::string pdb_path
) {
    std::ofstream pdb_file(pdb_path);
    if (!pdb_file) {
        return 1;
    }
    writeAtomCoordinatesToPDB(atoms, title, pdb_file);
    return 0;
}

std::vector< std::vector<AtomCoordinate> > splitAtomByResidue(
    const tcb::span<AtomCoordinate>& atomCoordinates
) {
    std::vector< std::vector<AtomCoordinate> > output;
    if (atomCoordinates.empty()) {
        return output;
    }
    std::vector<AtomCoordinate> currentResidue;

    for (size_t i = 0; i < atomCoordinates.size(); i++) {
        if (i == 0) {
            currentResidue.push_back(atomCoordinates[i]);
            continue;
        }

        if (startsNewResidue(atomCoordinates[i], atomCoordinates[i - 1])) {
            output.push_back(currentResidue);
            currentResidue.clear();
        }
        currentResidue.push_back(atomCoordinates[i]);
    }
    output.push_back(currentResidue);

    return output;
}

std::vector<std::string> getResidueNameVector(
    const tcb::span<AtomCoordinate>& atomCoordinates
) {
    std::vector<std::string> output;
    std::vector<std::vector<AtomCoordinate>> atomByResidue = splitAtomByResidue(atomCoordinates);
    output.reserve(atomByResidue.size());
    for (const auto& residueAtoms : atomByResidue) {
        if (!residueAtoms.empty()) {
            output.push_back(residueAtoms[0].residue);
        }
    }
    return output;
}

AtomCoordinate findFirstAtom(const std::vector<AtomCoordinate>& atoms, std::string atom_name) {
    for (const AtomCoordinate& curr_atm : atoms) {
        if (curr_atm.atom == atom_name) {
            return curr_atm;
        }
    }
    return AtomCoordinate();
}

void setAtomIndexSequentially(std::vector<AtomCoordinate>& atoms, int start) {
    for (size_t i = 0; i < atoms.size(); i++) {
        atoms[i].atom_index = start + i;
    }
}

void removeAlternativePosition(std::vector<AtomCoordinate>& atoms) {
    // If there is an alternative position, remove it
    for (size_t i = 1; i < atoms.size(); i++) {
        if (atoms[i].atom == atoms[i - 1].atom &&
            atoms[i].residue == atoms[i - 1].residue &&
            atoms[i].residue_index == atoms[i - 1].residue_index &&
            atoms[i].chain == atoms[i - 1].chain) {
            atoms.erase(atoms.begin() + i);
            i--;
        }
    }
}

std::vector<AtomCoordinate> getAtomsWithResidueIndex(
    std::vector<AtomCoordinate>& atoms, int residue_index
) {
    std::vector<AtomCoordinate> output;
    for (const AtomCoordinate& curr_atm : atoms) {
        if (curr_atm.residue_index == residue_index) {
            output.emplace_back(curr_atm);
        }
    }
    return output;
}

std::vector<AtomCoordinate> getAtomsWithResidueIndiceRange(
    std::vector<AtomCoordinate>& atoms, int start, int end
) {
    std::vector<AtomCoordinate> output;
    for (const AtomCoordinate& curr_atm : atoms) {
        if (curr_atm.residue_index >= start && curr_atm.residue_index < end) {
            output.emplace_back(curr_atm);
        }
    }
    return output;
}

std::vector<AtomCoordinate> getAtomsWithResidueIndex(
    const tcb::span<AtomCoordinate>& atoms, int residue_index,
    std::vector<std::string> atomNames
) {
    std::vector<AtomCoordinate> output;
    for (const AtomCoordinate& curr_atm : atoms) {
        if (curr_atm.residue_index == residue_index) {
            for (const std::string& atom_name : atomNames) {
                if (curr_atm.atom == atom_name) {
                    output.emplace_back(curr_atm);
                }
            }
        }
    }
    return output;
}

std::vector< std::vector<AtomCoordinate> > getAtomsWithResidueIndex(
    const tcb::span<AtomCoordinate>& atoms, std::vector<int> residue_index,
    std::vector<std::string> atomNames
) {
    std::vector<std::vector<AtomCoordinate>> output;
    for (int curr_index : residue_index) {
        output.emplace_back(getAtomsWithResidueIndex(atoms, curr_index, atomNames));
    }
    return output;
}

float RMSD(std::vector<AtomCoordinate>& atoms1, std::vector<AtomCoordinate>& atoms2) {
    // RMSD: Root Mean Square Deviation
    float sum = 0;
    // Sum of square of distance
    for (size_t i = 0; i < atoms1.size(); i++) {
        sum += pow(atoms1[i].coordinate.x - atoms2[i].coordinate.x, 2);
        sum += pow(atoms1[i].coordinate.y - atoms2[i].coordinate.y, 2);
        sum += pow(atoms1[i].coordinate.z - atoms2[i].coordinate.z, 2);
    }
    return sqrt(sum / atoms1.size());
}

std::vector<AtomCoordinate> _subsetAtomVectorWithIndices(
    std::vector<AtomCoordinate>& atoms,
    std::pair<size_t, size_t>& indices
) {
    std::vector<AtomCoordinate> output;
    for (size_t i = indices.first; i < indices.second; i++) {
        output.push_back(atoms[i]);
    }
    return output;
}

void _splitAtomVectorWithIndices(
    std::vector<AtomCoordinate>& atoms,
    std::vector< std::pair<size_t, size_t> >& indices,
    std::vector< std::vector<AtomCoordinate> >& output
) {
    if (indices.size() == 0) {
        output.push_back(atoms);
    } else {
        for (size_t i = 0; i < indices.size(); i++) {
            output.push_back(_subsetAtomVectorWithIndices(atoms, indices[i]));
        }
    }
}

/**
 * @brief Identify discontinuous regions in atom coordinate vector and return
 *        vector of indices of the start and end of each region.
 *        start: inclusive, end: exclusive [start, end)
 * @param atoms
 * @param mode
 * @return std::vector<std::pair<size_t, size_t>>
 */
std::vector< std::pair<size_t, size_t> > identifyChains(const std::vector<AtomCoordinate>& atoms) {
    std::vector< std::pair<size_t, size_t> > output;
    if (atoms.empty()) {
        return output;
    }
    size_t start = 0;
    for (size_t i = 1; i < atoms.size(); i++) {
        if (atoms[i].model != atoms[i - 1].model || atoms[i].chain != atoms[i - 1].chain) {
            output.emplace_back(start, i);
            start = i;
        }
    }
    if (start < atoms.size()) {
        output.emplace_back(start, atoms.size());
    }
    return output;
}

/**
 * @brief Identify discontinuous residue indices in atom coordinate vector and return
 *        coordinates that have the same chain
 * @param atoms
 * @return std::vector< std::pair<size_t, size_t> >
 */
std::vector<std::pair<size_t, size_t>> identifyDiscontinousResInd(
    const std::vector<AtomCoordinate>& atoms,
    size_t chain_start,
    size_t chain_end
) {
    std::vector<std::pair<size_t, size_t>> output;
    if (chain_start >= chain_end || chain_end > atoms.size()) {
        return output;
    }
    size_t start = chain_start;
    int previousResidueIndex = atoms[chain_start].residue_index;

    for (size_t i = chain_start + 1; i < chain_end; i++) {
        if (!startsNewResidue(atoms[i], atoms[i - 1])) {
            continue;
        }

        int currentResidueIndex = atoms[i].residue_index;
        if (currentResidueIndex != previousResidueIndex + 1) {
            output.emplace_back(start, i);
            start = i;
        }
        previousResidueIndex = currentResidueIndex;
    }
    output.emplace_back(start, chain_end);
    return output;
}

std::vector<std::pair<size_t, size_t>> identifyCompleteBackboneRegions(
    const tcb::span<AtomCoordinate>& atoms
) {
    std::vector<std::pair<size_t, size_t>> output;
    if (atoms.empty()) {
        return output;
    }

    size_t residueStart = 0;
    size_t runStart = atoms.size();
    size_t completeResiduesInRun = 0;

    auto flushRun = [&](size_t residueEnd) {
        if (runStart < atoms.size() && completeResiduesInRun >= 2) {
            output.emplace_back(runStart, residueEnd);
        }
        runStart = atoms.size();
        completeResiduesInRun = 0;
    };

    for (size_t i = 1; i <= atoms.size(); i++) {
        bool endOfResidue = (i == atoms.size()) || startsNewResidue(atoms[i], atoms[i - 1]);
        if (!endOfResidue) {
            continue;
        }

        bool hasN = false;
        bool hasCA = false;
        bool hasC = false;
        for (size_t j = residueStart; j < i; j++) {
            if (atoms[j].atom == "N") {
                hasN = true;
            } else if (atoms[j].atom == "CA") {
                hasCA = true;
            } else if (atoms[j].atom == "C") {
                hasC = true;
            }
        }

        if (hasN && hasCA && hasC) {
            if (runStart == atoms.size()) {
                runStart = residueStart;
            }
            completeResiduesInRun++;
        } else {
            flushRun(residueStart);
        }
        residueStart = i;
    }

    flushRun(atoms.size());
    return output;
}

std::vector<BackboneRegion> identifyBackboneRegions(const tcb::span<AtomCoordinate>& atoms) {
    std::vector<BackboneRegion> output;
    if (atoms.empty()) {
        return output;
    }

    struct ResidueRegion {
        size_t start;
        size_t end;
        bool completeBackbone;
    };

    std::vector<ResidueRegion> residues;
    size_t residueStart = 0;
    for (size_t i = 1; i <= atoms.size(); i++) {
        bool endOfResidue = (i == atoms.size()) || startsNewResidue(atoms[i], atoms[i - 1]);
        if (!endOfResidue) {
            continue;
        }

        bool hasN = false;
        bool hasCA = false;
        bool hasC = false;
        for (size_t j = residueStart; j < i; j++) {
            if (atoms[j].atom == "N") {
                hasN = true;
            } else if (atoms[j].atom == "CA") {
                hasCA = true;
            } else if (atoms[j].atom == "C") {
                hasC = true;
            }
        }
        residues.push_back({residueStart, i, hasN && hasCA && hasC});
        residueStart = i;
    }

    size_t i = 0;
    while (i < residues.size()) {
        if (!residues[i].completeBackbone) {
            size_t start = residues[i].start;
            size_t end = residues[i].end;
            i++;
            while (i < residues.size() && !residues[i].completeBackbone) {
                end = residues[i].end;
                i++;
            }
            output.push_back({start, end, false});
            continue;
        }

        size_t runStartIndex = i;
        size_t start = residues[i].start;
        size_t end = residues[i].end;
        i++;
        while (i < residues.size() && residues[i].completeBackbone) {
            end = residues[i].end;
            i++;
        }
        if (i - runStartIndex >= 2) {
            output.push_back({start, end, true});
        } else {
            output.push_back({start, end, false});
        }
    }

    return output;
}

namespace {

template <typename T>
void appendBinaryValue(std::string& output, const T& value) {
    size_t offset = output.size();
    output.resize(offset + sizeof(T));
    memcpy(&output[offset], &value, sizeof(T));
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

} // namespace

bool serializeAtomCoordinates(
    const std::vector<AtomCoordinate>& atoms,
    std::string& output
) {
    output.clear();
    uint32_t atomCount = static_cast<uint32_t>(atoms.size());
    appendBinaryValue(output, atomCount);
    for (const auto& atom : atoms) {
        if (atom.atom.size() > 255 || atom.residue.size() > 255) {
            return false;
        }
        uint8_t atomNameLen = static_cast<uint8_t>(atom.atom.size());
        uint8_t residueNameLen = static_cast<uint8_t>(atom.residue.size());
        appendBinaryValue(output, atomNameLen);
        appendBinaryValue(output, residueNameLen);
        appendBinaryValue(output, atom.atom_index);
        appendBinaryValue(output, atom.residue_index);
        appendBinaryValue(output, atom.coordinate.x);
        appendBinaryValue(output, atom.coordinate.y);
        appendBinaryValue(output, atom.coordinate.z);
        appendBinaryValue(output, atom.occupancy);
        appendBinaryValue(output, atom.tempFactor);
        output.append(atom.atom);
        output.append(atom.residue);
    }
    return true;
}

bool deserializeAtomCoordinates(
    const char* data,
    size_t size,
    std::vector<AtomCoordinate>& atoms
) {
    atoms.clear();
    size_t offset = 0;
    uint32_t atomCount = 0;
    if (!readBinaryValue(data, size, offset, atomCount)) {
        return false;
    }
    atoms.reserve(atomCount);
    for (uint32_t i = 0; i < atomCount; i++) {
        uint8_t atomNameLen = 0;
        uint8_t residueNameLen = 0;
        AtomCoordinate atom;
        if (!readBinaryValue(data, size, offset, atomNameLen) ||
            !readBinaryValue(data, size, offset, residueNameLen) ||
            !readBinaryValue(data, size, offset, atom.atom_index) ||
            !readBinaryValue(data, size, offset, atom.residue_index) ||
            !readBinaryValue(data, size, offset, atom.coordinate.x) ||
            !readBinaryValue(data, size, offset, atom.coordinate.y) ||
            !readBinaryValue(data, size, offset, atom.coordinate.z) ||
            !readBinaryValue(data, size, offset, atom.occupancy) ||
            !readBinaryValue(data, size, offset, atom.tempFactor)) {
            return false;
        }
        size_t required = static_cast<size_t>(atomNameLen) + static_cast<size_t>(residueNameLen);
        if (offset + required > size) {
            return false;
        }
        atom.atom.assign(data + offset, data + offset + atomNameLen);
        offset += atomNameLen;
        atom.residue.assign(data + offset, data + offset + residueNameLen);
        offset += residueNameLen;
        atoms.push_back(std::move(atom));
    }
    return offset == size;
}
