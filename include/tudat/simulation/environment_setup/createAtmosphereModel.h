/*    Copyright (c) 2010-2019, Delft University of Technology
 *    All rigths reserved
 *
 *    This file is part of the Tudat. Redistribution and use in source and
 *    binary forms, with or without modification, are permitted exclusively
 *    under the terms of the Modified BSD license. You should have received
 *    a copy of the license with this file. If not, please or visit:
 *    http://tudat.tudelft.nl/LICENSE.
 */

#ifndef TUDAT_CREATEATMOSPHEREMODEL_H
#define TUDAT_CREATEATMOSPHEREMODEL_H

#include <string>
#include <map>

#include <memory>
#include <boost/date_time/posix_time/time_period.hpp>

#include "body.h"
#include "tudat/io/basicInputOutput.h"
#include "tudat/astro/aerodynamics/atmosphereModel.h"
#include "tudat/astro/aerodynamics/exponentialAtmosphere.h"
#include "tudat/astro/aerodynamics/customConstantTemperatureAtmosphere.h"
#include "tudat/astro/basic_astro/physicalConstants.h"
#include "tudat/math/interpolators/interpolator.h"
#include "tudat/basics/identityElements.h"
#include <boost/variant.hpp>


namespace tudat
{
namespace simulation_setup
{
using namespace aerodynamics;

// ---- coefficient storage: rows = coeffs, cols = {C,S} ----
using StokesBlock = Eigen::Matrix< double, Eigen::Dynamic, 2, Eigen::RowMajor >;

struct FileMeta
{
    // Epochs in seconds since J2000 (or your preferred convention)
    double start_epoch{};
    double end_epoch{};
    std::string source_tag; // filename / description
};

// Degree-major mapping (n = 0..nmax, m = 0..n)
inline std::size_t nm_to_index_deg_major( int n, int m )
{
    // index(n,m) = n(n+1)/2 + m
    if(n < 0 || m < 0 || m > n) throw std::out_of_range( "nm_to_index: invalid (n,m)" );
    return static_cast< std::size_t >(n * ( n + 1 ) / 2 + m);
}

inline int index_to_n_deg_major( std::size_t k )
{
    // invert k = n(n+1)/2 + m  => n = floor((sqrt(8k+1)-1)/2)
    const double nd = ( std::sqrt( 8.0 * static_cast< double >(k) + 1.0 ) - 1.0 ) * 0.5;
    return static_cast< int >(std::floor( nd + 1e-12 ));
}

inline std::pair< int, int > index_to_nm_deg_major( std::size_t k )
{
    const int n = index_to_n_deg_major( k );
    const std::size_t base = static_cast< std::size_t >(n * ( n + 1 ) / 2);
    const int m = static_cast< int >(k - base);
    return { n, m };
}

// ============= Core Data Models (Pure data storage) =============

// ---- Stokes Dataset ----
class ComaStokesDataset
{
public:
    struct FileMeta
    {
        double start_epoch{};
        double end_epoch{};
        std::string source_tag;
        double referenceRadius{};
    };

    // Factory method
    static ComaStokesDataset create(
        std::vector<FileMeta> files,
        std::vector<double> radii,
        std::vector<double> lons,
        int nmax)
    {
        if (files.empty() || radii.empty() || lons.empty() || nmax < 0)
            throw std::runtime_error("StokesDataset: invalid metadata.");

        // Check that no radius exceeds the reference radius from any file
        for (const auto& file : files)
        {
            for (const auto& radius : radii)
            {
                if (radius > file.referenceRadius)
                {
                    throw std::runtime_error("StokesDataset: The max radius in the vector (" +
                                           std::to_string(radius) +
                                           ") should not be larger than the reference radius (" +
                                           std::to_string(file.referenceRadius) + ").");
                }
            }
        }

        ComaStokesDataset g;
        g.files_ = std::move(files);
        g.radii_ = std::move(radii);
        g.lons_ = std::move(lons);
        g.nmax_ = nmax;

        g.n_files_ = g.files_.size();
        g.n_radii_ = g.radii_.size();
        g.n_lons_ = g.lons_.size();
        g.n_coeffs_ = static_cast<std::size_t>((nmax + 1) * (nmax + 2) / 2);

        const std::size_t totalRows = g.n_files_ * g.n_radii_ * g.n_lons_ * g.n_coeffs_;
        g.data_.setZero(static_cast<Eigen::Index>(totalRows), 2);
        return g;
    }

    // Pure data accessors
    std::size_t nFiles() const { return n_files_; }
    std::size_t nRadii() const { return n_radii_; }
    std::size_t nLongitudes() const { return n_lons_; }
    std::size_t nCoeffs() const { return n_coeffs_; }
    int nmax() const { return nmax_; }
    const std::vector<double>& radii() const { return radii_; }
    const std::vector<double>& lons() const { return lons_; }
    const std::vector<FileMeta>& files() const { return files_; }

    // Convenience accessor for reference radius
    double getReferenceRadius(std::size_t f = 0) const
    {
        if (f >= n_files_) throw std::out_of_range("File index out of range");
        return files_[f].referenceRadius;
    }

    // Block access
    auto block(std::size_t f, std::size_t r, std::size_t l)
    {
        const std::size_t start = startRow_(f, r, l);
        return data_.block(static_cast<Eigen::Index>(start), 0,
                          static_cast<Eigen::Index>(n_coeffs_), 2);
    }

    auto block(std::size_t f, std::size_t r, std::size_t l) const
    {
        const std::size_t start = startRow_(f, r, l);
        return data_.block(static_cast<Eigen::Index>(start), 0,
                          static_cast<Eigen::Index>(n_coeffs_), 2);
    }

    // Coefficient matrices
    std::pair<Eigen::MatrixXd, Eigen::MatrixXd>
    getCoefficientMatrices(std::size_t f, std::size_t r, std::size_t l) const
    {
        Eigen::MatrixXd cosine = Eigen::MatrixXd::Zero(nmax_ + 1, nmax_ + 1);
        Eigen::MatrixXd sine = Eigen::MatrixXd::Zero(nmax_ + 1, nmax_ + 1);

        const auto blk = block(f, r, l);
        for (std::size_t k = 0; k < n_coeffs_; ++k)
        {
            auto nm = index_to_nm_deg_major(k);
            int n = nm.first;
            int m = nm.second;
            cosine(n, m) = blk(static_cast<Eigen::Index>(k), 0);
            sine(n, m) = blk(static_cast<Eigen::Index>(k), 1);
        }
        return {cosine, sine};
    }

    // Single coefficient access
    void setCoeff(std::size_t f, std::size_t r, std::size_t l,
                  int n, int m, double C, double S)
    {
        const std::size_t k = nm_to_index_deg_major(n, m);
        if (k >= n_coeffs_) throw std::out_of_range("setCoeff: (n,m) exceeds nmax");
        const std::size_t row = startRow_(f, r, l) + k;
        data_(static_cast<Eigen::Index>(row), 0) = C;
        data_(static_cast<Eigen::Index>(row), 1) = S;
    }

    std::pair<double, double> getCoeff(std::size_t f, std::size_t r, std::size_t l,
                                       int n, int m) const
    {
        const std::size_t k = nm_to_index_deg_major(n, m);
        if (k >= n_coeffs_) throw std::out_of_range("getCoeff: (n,m) exceeds nmax");
        const std::size_t row = startRow_(f, r, l) + k;
        return {data_(static_cast<Eigen::Index>(row), 0),
                data_(static_cast<Eigen::Index>(row), 1)};
    }

    const StokesBlock& data() const { return data_; }
    StokesBlock& data() { return data_; }

private:
    std::size_t startRow_(std::size_t f, std::size_t r, std::size_t l) const
    {
        if (f >= n_files_ || r >= n_radii_ || l >= n_lons_)
            throw std::out_of_range("StokesDataset: index out of range.");
        const std::size_t cell = ((f * n_radii_) + r) * n_lons_ + l;
        return cell * n_coeffs_;
    }

    std::vector<FileMeta> files_;
    std::vector<double> radii_;
    std::vector<double> lons_;
    int nmax_{};
    std::size_t n_files_{}, n_radii_{}, n_lons_{}, n_coeffs_{};
    StokesBlock data_;
};

// ---- Poly Dataset ----
class ComaPolyDataset
{
public:
    struct FileMeta
    {
        double referenceRadius{};
        Eigen::VectorXd powersInvRadius;
        std::vector<std::pair<double, double>> timePeriods;
        int maxDegreeSH{};
        Eigen::Index numRadialTerms{};
        Eigen::Index numIntervals{};
        std::string sourcePath;
    };

    // Simple accessors
    std::size_t getNumFiles() const { return numPolyCoefFiles_; }
    const FileMeta& getFileMeta(std::size_t f) const
    {
        boundsCheck_(f);
        return fileMeta_[f];
    }

    const Eigen::MatrixXd& getPolyCoefficients(std::size_t f) const
    {
        boundsCheck_(f);
        return polyCoefficients_[f];
    }

    const Eigen::ArrayXXi& getSHDegreeAndOrderIndices(std::size_t f) const
    {
        boundsCheck_(f);
        return SHDegreeAndOrderIndices_[f];
    }

    // Convenience accessors
    double getReferenceRadius(std::size_t f) const
    {
        boundsCheck_(f);
        return fileMeta_[f].referenceRadius;
    }

    const Eigen::VectorXd& getPowersInvRadius(std::size_t f) const
    {
        boundsCheck_(f);
        return fileMeta_[f].powersInvRadius;
    }

    int getMaxDegreeSH(std::size_t f) const
    {
        boundsCheck_(f);
        return fileMeta_[f].maxDegreeSH;
    }

    // Column access methods
    Eigen::VectorXd columnForNM(std::size_t f, int n, int m) const
    {
        auto [ok, col] = findColumn_(f, n, m);
        if (!ok) throw std::out_of_range("columnForNM: (n,m) not found");
        return polyCoefficients_[f].col(col);
    }

    double value(std::size_t f, Eigen::Index termIndex, int n, int m) const
    {
        auto [ok, col] = findColumn_(f, n, m);
        if (!ok) throw std::out_of_range("value: (n,m) not found");
        if (termIndex < 0 || termIndex >= polyCoefficients_[f].rows())
            throw std::out_of_range("value: termIndex out of range");
        return polyCoefficients_[f](termIndex, col);
    }

    void clear()
    {
        numPolyCoefFiles_ = 0;
        polyCoefficients_.clear();
        SHDegreeAndOrderIndices_.clear();
        fileMeta_.clear();
        nmToColCache_.clear();
    }

protected:
    // Allow friend classes to populate data
    friend class ComaPolyDatasetReader;

    void setData(std::size_t numFiles,
                 std::vector<Eigen::MatrixXd> polyCoeffs,
                 std::vector<Eigen::ArrayXXi> shIndices,
                 std::vector<FileMeta> meta)
    {
        numPolyCoefFiles_ = numFiles;
        polyCoefficients_ = std::move(polyCoeffs);
        SHDegreeAndOrderIndices_ = std::move(shIndices);
        fileMeta_ = std::move(meta);

        // Build caches
        nmToColCache_.resize(numFiles);
        for (std::size_t f = 0; f < numFiles; ++f)
            buildNmMap_(f);
    }

private:
    struct PairHash
    {
        std::size_t operator()(const std::pair<int, int>& p) const noexcept
        {
            return (static_cast<std::size_t>(p.first) << 32) ^
                   static_cast<std::size_t>(p.second);
        }
    };

    void boundsCheck_(std::size_t f) const
    {
        if (f >= numPolyCoefFiles_)
            throw std::out_of_range("file index out of range");
    }

    void buildNmMap_(std::size_t f)
    {
        nmToColCache_[f].clear();
        const auto& sh = SHDegreeAndOrderIndices_[f];
        for (Eigen::Index c = 0; c < sh.cols(); ++c)
        {
            int n = sh(0, c);
            int m = sh(1, c);
            nmToColCache_[f][{n, m}] = static_cast<int>(c);
        }
    }

    std::pair<bool, int> findColumn_(std::size_t f, int n, int m) const
    {
        boundsCheck_(f);
        const auto& map = nmToColCache_[f];
        auto it = map.find({n, m});
        if (it == map.end()) return {false, -1};
        return {true, it->second};
    }

    std::size_t numPolyCoefFiles_{0};
    std::vector<Eigen::MatrixXd> polyCoefficients_;
    std::vector<Eigen::ArrayXXi> SHDegreeAndOrderIndices_;
    std::vector<FileMeta> fileMeta_;
    std::vector<std::unordered_map<std::pair<int, int>, int, PairHash>> nmToColCache_;
};

// ============= I/O Components (Separate from data) =============

// ---- Reader for Poly Coefficients ----
class ComaPolyDatasetReader
{
public:
    static ComaPolyDataset readFromFiles(const std::vector<std::string>& filePaths)
    {
        if (filePaths.empty())
            throw std::invalid_argument("ComaPolyDatasetReader: empty file list");

        const std::size_t n = filePaths.size();
        std::vector<Eigen::MatrixXd> polyCoefficients(n);
        std::vector<Eigen::ArrayXXi> SHDegreeAndOrderIndices(n);
        std::vector<ComaPolyDataset::FileMeta> fileMeta(n);

        for (std::size_t fileIdx = 0; fileIdx < n; ++fileIdx)
        {
            readSingleFile(filePaths[fileIdx], fileIdx,
                          polyCoefficients, SHDegreeAndOrderIndices, fileMeta);
        }

        ComaPolyDataset dataset;
        dataset.setData(n, std::move(polyCoefficients),
                       std::move(SHDegreeAndOrderIndices),
                       std::move(fileMeta));
        return dataset;
    }

private:
    static void readSingleFile(const std::string& filePath,
                               std::size_t fileIdx,
                               std::vector<Eigen::MatrixXd>& polyCoefficients,
                               std::vector<Eigen::ArrayXXi>& SHDegreeAndOrderIndices,
                               std::vector<ComaPolyDataset::FileMeta>& fileMeta)
    {
        // ===== Implementation adapted from old ComaPolyDataset::readInputFiles =====
        fileMeta[fileIdx].sourcePath = filePath;

        std::ifstream file(filePath);
        if (!file.is_open())
        {
            std::cerr << "[ERROR] Could not open file '" << filePath << "'.\n";
            std::exit(EXIT_FAILURE);
        }

        std::string line;
        std::vector<std::string> tokens;

        int maxDegreeSH = 0;
        Eigen::Index numTerms = 0, numCoefs = 0, numRadialTerms = 0, numIntervals = 0;
        Eigen::VectorXd powers;

        // ----- Parse header -----
        while (std::getline(file, line))
        {
            if (line.empty()) continue;
            if (line[0] != '#') break;

            std::string headerLine = line.substr(1);
            boost::trim(headerLine);
            boost::split(tokens, headerLine, boost::is_any_of(", \t"), boost::token_compress_on);
            if (tokens.empty()) continue;

            const std::string& key = tokens[0];

            if (boost::iequals(key, "N(SH)"))
            {
                maxDegreeSH = std::stoi(line.substr(line.find_last_of(" \t") + 1));
                numCoefs = (maxDegreeSH + 1) * (maxDegreeSH + 1);
            }
            else if (boost::icontains(key, "PWRS"))
            {
                std::string tail = line.substr(line.find("PWRS"));
                boost::trim(tail);
                std::vector<std::string> pwrtok;
                boost::split(pwrtok, tail, boost::is_any_of(", \t"), boost::token_compress_on);
                std::size_t start = (!pwrtok.empty() &&
                                     !std::all_of(pwrtok[0].begin(), pwrtok[0].end(), ::isdigit))
                                    ? 1
                                    : 0;
                auto count = static_cast<Eigen::Index>(pwrtok.size() - start);
                powers.resize(count);
                for (Eigen::Index j = 0; j < count; ++j)
                    powers[j] = std::stod(pwrtok[start + j]);
            }
            else if (boost::iequals(key, "R"))
            {
                double R_km = std::stod(line.substr(line.find_last_of(" \t") + 1));
                fileMeta[fileIdx].referenceRadius = R_km * 1000.0; // Convert km to meters
            }
            else if (line.find("N(r)") != std::string::npos && line.find("N(T)") != std::string::npos)
            {
                std::string content = line.substr(1); // strip '#'
                boost::trim(content);
                boost::split(tokens, content, boost::is_any_of(", \t"), boost::token_compress_on);

                if (tokens.size() >= 2)
                {
                    int a = std::stoi(tokens[tokens.size() - 2]);
                    int b = std::stoi(tokens[tokens.size() - 1]);

                    if (line.find("N(r)") < line.find("N(T)"))
                    {
                        numRadialTerms = a;
                        numIntervals = b;
                    }
                    else
                    {
                        numIntervals = a;
                        numRadialTerms = b;
                    }
                    numTerms = numRadialTerms * numIntervals;
                }
                else
                {
                    std::cerr << "[ERROR] Failed to parse N(T)/N(r) line: " << line << std::endl;
                    std::exit(EXIT_FAILURE);
                }
            }
            // OPTIONAL: parse time period lines and fill fileMeta[fileIdx].timePeriods
        }

        // ----- Validation -----
        if (numTerms <= 0 || numCoefs <= 0 || powers.size() == 0)
        {
            std::cerr << "[ERROR] Header parsing failed in file: " << filePath << std::endl;
            std::cerr << "  numTerms = " << numTerms
                      << "\n  numCoefs = " << numCoefs
                      << "\n  powersInvRadius.size() = " << powers.size() << std::endl;
            std::exit(EXIT_FAILURE);
        }

        fileMeta[fileIdx].maxDegreeSH = maxDegreeSH;
        fileMeta[fileIdx].numRadialTerms = numRadialTerms;
        fileMeta[fileIdx].numIntervals = numIntervals;
        fileMeta[fileIdx].powersInvRadius = powers;

        // --- Allocate containers ---
        Eigen::MatrixXd& currentPolyCoefficients = polyCoefficients[fileIdx];
        Eigen::ArrayXXi& currentShDegreeAndOrder = SHDegreeAndOrderIndices[fileIdx];

        currentPolyCoefficients.resize(numTerms, numCoefs);
        currentShDegreeAndOrder.resize(2, numCoefs);

        // ----- Read poly coefficients data block -----
        Eigen::Index coefIndex = -1;
        do
        {
            boost::trim(line);
            if (line.empty() || line[0] == '#') continue;

            boost::split(tokens, line, boost::is_any_of(", \t"), boost::token_compress_on);
            if (static_cast<Eigen::Index>(tokens.size()) == numTerms + 2)
            {
                ++coefIndex;
                currentShDegreeAndOrder(0, coefIndex) = std::stoi(tokens[0]); // n
                currentShDegreeAndOrder(1, coefIndex) = std::stoi(tokens[1]); // m
                for (Eigen::Index j = 0; j < numTerms; ++j)
                    currentPolyCoefficients(j, coefIndex) = std::stod(tokens[j + 2]);
            }
        } while (std::getline(file, line));

        file.close();
        // The (n,m)->col map is built by ComaPolyDataset::setData()
    }
};

// ---- Writer for Stokes Coefficients ----
class ComaStokesDatasetWriter
{
public:
    static void writeCsvForFile(const ComaStokesDataset& dataset,
                                std::size_t f,
                                const std::string& outputPath)
    {
        if (f >= dataset.nFiles())
            throw std::out_of_range("writeCsvForFile: file index out of range");

        auto csvEscape = [](std::ostream& os, const std::string& s)
        {
            bool needs = false;
            for (char c : s)
            {
                if (c == ',' || c == '"' || c == '\n' || c == '\r')
                { needs = true; break; }
            }
            if (!needs) { os << s; return; }
            os << '"';
            for (char c : s) os << (c == '"' ? "\"\"" : std::string(1, c));
            os << '"';
        };

        boost::filesystem::path p(outputPath);
        if (p.has_parent_path())
            boost::filesystem::create_directories(p.parent_path());

        std::ofstream os(outputPath, std::ios::binary);
        if (!os)
            throw std::runtime_error("writeCsvForFile: cannot open " + outputPath);
        os.imbue(std::locale::classic());

        const auto& fm = dataset.files()[f];

        auto set_sci = [&] {
            os.setf(std::ios::scientific, std::ios::floatfield);
            os << std::setprecision(17);
        };
        auto set_def = [&] {
            os.setf(std::ios::fmtflags(0), std::ios::floatfield);
            os << std::setprecision(17);
        };

        // Row 1: metadata
        os << "meta";
        set_sci();
        os << ",start_epoch=" << fm.start_epoch;
        os << ",end_epoch=" << fm.end_epoch;
        os << ",max_degree=" << dataset.nmax();
        os << ",max_order=" << dataset.nmax();
        os << ",n_radii=" << dataset.nRadii();
        os << ",n_lons=" << dataset.nLongitudes();
        os << ",n_coeffs=" << dataset.nCoeffs();
        os << ",source=";
        csvEscape(os, fm.source_tag);
        os << '\n';

        // Row 2: radii
        os << "radii [meter]";
        for (double r : dataset.radii())
        {
            set_def(); os << ',' << r;
        }
        os << '\n';

        // Row 3: longitudes
        os << "longitudes [degree]";
        for (double L : dataset.lons())
        {
            set_def(); os << ',' << L;
        }
        os << '\n';

        // Blocks
        for (std::size_t ri = 0; ri < dataset.nRadii(); ++ri)
        {
            for (std::size_t li = 0; li < dataset.nLongitudes(); ++li)
            {
                const std::size_t block_id = ri * dataset.nLongitudes() + li;

                os << "ID," << block_id << ',';
                set_def();
                os << "r_0=" << dataset.radii()[ri] << ','
                   << "l_0=" << dataset.lons()[li] << '\n';

                os << "n,m,C,S\n";

                const auto blk = dataset.block(f, ri, li);
                for (std::size_t k = 0; k < dataset.nCoeffs(); ++k)
                {
                    const auto nm = index_to_nm_deg_major(k);
                    const double C = blk(static_cast<Eigen::Index>(k), 0);
                    const double S = blk(static_cast<Eigen::Index>(k), 1);
                    os << nm.first << ',' << nm.second << ',';
                    set_sci();
                    os << C << ',' << S << '\n';
                }
            }
        }

        os.flush();
        if (!os) throw std::runtime_error("writeCsvForFile: write failed");
    }

    static void writeCsvAll(const ComaStokesDataset& dataset,
                           const std::string& outputDir,
                           const std::string& prefix = "stokes")
    {
        boost::filesystem::create_directories(outputDir);
        for (std::size_t f = 0; f < dataset.nFiles(); ++f)
        {
            boost::filesystem::path path = boost::filesystem::path(outputDir) /
                (prefix + "_file" + std::to_string(f) + ".csv");
            writeCsvForFile(dataset, f, path.string());
        }
    }
};

// ---- Reader for Stokes Coefficients (future) ----
class ComaStokesDatasetReader
{
public:
    static ComaStokesDataset readFromCsv(const std::string& csvPath)
    {
        std::ifstream ifs(csvPath);
        if (!ifs)
            throw std::runtime_error("readFromCsv: cannot open " + csvPath);

        // Parse metadata line
        std::string line;
        std::getline(ifs, line);
        if (line.empty() || line.substr(0, 4) != "meta")
            throw std::runtime_error("readFromCsv: invalid metadata line");

        auto parseMeta = [](const std::string& metaLine) -> std::tuple<double, double, int, int, int, int, int, std::string>
        {
            std::istringstream ss(metaLine);
            std::string token;
            double start_epoch = 0, end_epoch = 0;
            int max_degree = 0, max_order = 0, n_radii = 0, n_lons = 0, n_coeffs = 0;
            std::string source;

            while (std::getline(ss, token, ','))
            {
                if (token.find("start_epoch=") == 0)
                    start_epoch = std::stod(token.substr(12));
                else if (token.find("end_epoch=") == 0)
                    end_epoch = std::stod(token.substr(10));
                else if (token.find("max_degree=") == 0)
                    max_degree = std::stoi(token.substr(11));
                else if (token.find("max_order=") == 0)
                    max_order = std::stoi(token.substr(10));
                else if (token.find("n_radii=") == 0)
                    n_radii = std::stoi(token.substr(8));
                else if (token.find("n_lons=") == 0)
                    n_lons = std::stoi(token.substr(7));
                else if (token.find("n_coeffs=") == 0)
                    n_coeffs = std::stoi(token.substr(9));
                else if (token.find("source=") == 0)
                    source = token.substr(7);
            }
            return {start_epoch, end_epoch, max_degree, max_order, n_radii, n_lons, n_coeffs, source};
        };

        auto [start_epoch, end_epoch, max_degree, max_order, n_radii, n_lons, n_coeffs, source] = parseMeta(line);

        // Parse radii line
        std::getline(ifs, line);
        std::istringstream radiiStream(line);
        std::string token;
        std::getline(radiiStream, token, ','); // skip "radii [meter]"
        std::vector<double> radii;
        while (std::getline(radiiStream, token, ','))
            radii.push_back(std::stod(token));

        // Parse longitudes line
        std::getline(ifs, line);
        std::istringstream lonsStream(line);
        std::getline(lonsStream, token, ','); // skip "longitudes [degree]"
        std::vector<double> lons;
        while (std::getline(lonsStream, token, ','))
            lons.push_back(std::stod(token));

        // Create dataset with single file
        std::vector<ComaStokesDataset::FileMeta> files(1);
        files[0].start_epoch = start_epoch;
        files[0].end_epoch = end_epoch;
        files[0].source_tag = source;

        ComaStokesDataset dataset = ComaStokesDataset::create(std::move(files), radii, lons, max_degree);

        // Parse coefficient blocks
        while (std::getline(ifs, line))
        {
            if (line.find("ID,") == 0)
            {
                // Parse block ID line: "ID,block_id,r_0=...,l_0=..."
                std::istringstream idStream(line);
                std::string idToken;
                std::getline(idStream, idToken, ','); // skip "ID"
                std::getline(idStream, idToken, ','); // block_id
                int block_id = std::stoi(idToken);

                std::size_t ri = block_id / n_lons;
                std::size_t li = block_id % n_lons;

                // Skip header line "n,m,C,S"
                std::getline(ifs, line);

                // Read coefficients
                for (int k = 0; k < n_coeffs; ++k)
                {
                    std::getline(ifs, line);
                    std::istringstream coeffStream(line);
                    std::string nStr, mStr, cStr, sStr;
                    std::getline(coeffStream, nStr, ',');
                    std::getline(coeffStream, mStr, ',');
                    std::getline(coeffStream, cStr, ',');
                    std::getline(coeffStream, sStr, ',');

                    int n = std::stoi(nStr);
                    int m = std::stoi(mStr);
                    double C = std::stod(cStr);
                    double S = std::stod(sStr);

                    dataset.setCoeff(0, ri, li, n, m, C, S);
                }
            }
        }

        return dataset;
    }

    static ComaStokesDataset readFromCsvFolder(const std::string& dir,
                                               const std::string& prefix = "stokes")
    {
        namespace bf = boost::filesystem;

        if (!bf::exists(dir) || !bf::is_directory(dir))
            throw std::runtime_error("readFromCsvFolder: directory does not exist: " + dir);

        // Find all CSV files with the given prefix
        std::vector<std::string> csvFiles;
        for (bf::directory_iterator it(dir); it != bf::directory_iterator(); ++it)
        {
            const std::string filename = it->path().filename().string();
            if (filename.find(prefix + "_file") == 0 && filename.substr(filename.length() - 4) == ".csv")
                csvFiles.push_back(it->path().string());
        }

        if (csvFiles.empty())
            throw std::runtime_error("readFromCsvFolder: no CSV files found with prefix " + prefix);

        // Sort files by number
        std::sort(csvFiles.begin(), csvFiles.end(), [&prefix](const std::string& a, const std::string& b)
        {
            auto extractNum = [&prefix](const std::string& path) -> int
            {
                bf::path p(path);
                std::string name = p.stem().string();
                std::string numStr = name.substr(prefix.length() + 5); // +5 for "_file"
                return std::stoi(numStr);
            };
            return extractNum(a) < extractNum(b);
        });

        // Read first file to get structure
        ComaStokesDataset firstDataset = readFromCsv(csvFiles[0]);

        if (csvFiles.size() == 1)
            return firstDataset;

        // Create multi-file dataset
        std::vector<ComaStokesDataset::FileMeta> allFiles;
        for (const auto& csvFile : csvFiles)
        {
            ComaStokesDataset singleDataset = readFromCsv(csvFile);
            allFiles.push_back(singleDataset.files()[0]);
        }

        ComaStokesDataset multiDataset = ComaStokesDataset::create(
            std::move(allFiles), firstDataset.radii(), firstDataset.lons(), firstDataset.nmax());

        // Copy data from all files
        for (std::size_t f = 0; f < csvFiles.size(); ++f)
        {
            ComaStokesDataset singleDataset = readFromCsv(csvFiles[f]);
            for (std::size_t ri = 0; ri < multiDataset.nRadii(); ++ri)
            {
                for (std::size_t li = 0; li < multiDataset.nLongitudes(); ++li)
                {
                    auto srcBlock = singleDataset.block(0, ri, li);
                    auto destBlock = multiDataset.block(f, ri, li);
                    destBlock = srcBlock;
                }
            }
        }

        return multiDataset;
    }
};

// ============= Processing/Transformation Components =============

class StokesCoefficientsEvaluator
{
public:
private:
    static double radialPolyvalAndTemporalIFFT(const double r, const double alf,
                                               const Eigen::MatrixXd& P, const Eigen::ArrayXd& pw)
    {
        using namespace Eigen;

        int N = P.rows();

        if (N % 2 == 0)
        {
            return
            ((ArrayXd(N) << cos(ArrayXd::LinSpaced(N/2+1, 0.0, double(N/2)) * alf),
                           sin(ArrayXd::LinSpaced(N/2-1, 1.0, double(N/2-1)) * alf)).finished().matrix().transpose() *
             P * pow(r, pw).matrix()
            ).value() / double(N);
        }
        else
        {
            return
            ((ArrayXd(N) << cos(ArrayXd::LinSpaced(N/2+1, 0.0, double(N/2)) * alf),
                           sin(ArrayXd::LinSpaced(N/2, 1.0, double(N/2)) * alf)).finished().matrix().transpose() *
             P * pow(r, pw).matrix()
            ).value() / double(N);
        }
    }

    static double reducedToTemporalIFFT(const double alf, const Eigen::VectorXd& P)
    {
        using namespace Eigen;

        int N = P.rows();

        if (N % 2 == 0)
        {
            return
            ((ArrayXd(N) << cos(ArrayXd::LinSpaced(N/2+1, 0.0, double(N/2)) * alf),
                           sin(ArrayXd::LinSpaced(N/2-1, 1.0, double(N/2-1)) * alf)).finished().matrix().transpose() * P
            ).value() / double(N);
        }
        else
        {
            return
            ((ArrayXd(N) << cos(ArrayXd::LinSpaced(N/2+1, 0.0, double(N/2)) * alf),
                           sin(ArrayXd::LinSpaced(N/2, 1.0, double(N/2)) * alf)).finished().matrix().transpose() * P
            ).value() / double(N);
        }
    }

public:
    static void evaluate2D(
        const double radius_m,            // meter
        const double solarLongitude,    // radians
        const Eigen::ArrayXXd& polyCoefficients,
        const Eigen::ArrayXXi& atDegreeAndOrder,
        const Eigen::VectorXd& atPowersInvRadius,
        double refRadius_m,               // meter
        Eigen::MatrixXd& cosineCoefficients,
        Eigen::MatrixXd& sineCoefficients,
        int maxDegree,
        int maxOrder)
    {
        // --- Unit conversion ---
        const double radius_km = radius_m / 1000.0; // Conversion from m to km
        const double refRadius = refRadius_m / 1000.0; // Conversion from m to km

        const int maxDegAvailable = atDegreeAndOrder.row(0).maxCoeff();
        const int maxOrdAvailable = atDegreeAndOrder.row(1).abs().maxCoeff();

        if (maxDegree < 0) maxDegree = maxDegAvailable;
        if (maxOrder < 0) maxOrder = maxOrdAvailable;

        if (maxDegree > maxDegAvailable || maxOrder > maxOrdAvailable)
        {
            std::ostringstream err;
            err << "[FATAL] Requested maxDegree=" << maxDegree
                << ", maxOrder=" << maxOrder
                << " exceeds available (degree=" << maxDegAvailable
                << ", order=" << maxOrdAvailable << ")";
            throw std::runtime_error(err.str());
        }

        const int numRadialTerms = atPowersInvRadius.size();
        const int numIntervals = polyCoefficients.rows() / numRadialTerms;

        cosineCoefficients = Eigen::MatrixXd::Zero(maxDegree + 1, maxOrder + 1);
        sineCoefficients   = Eigen::MatrixXd::Zero(maxDegree + 1, maxOrder + 1);

        double VR = (refRadius < 1.0e-10) ? 0.0 : 1.0 / refRadius;

        if (radius_km <= refRadius || refRadius < 1.0e-10)
        {
            for (int i = 0; i < polyCoefficients.cols(); ++i)
            {
                const int l = atDegreeAndOrder(0, i);
                const int m = atDegreeAndOrder(1, i);
                const int absM = std::abs(m);

                if (l > maxDegree || absM > maxOrder)
                    continue;

                const Eigen::MatrixXd polyCoefs =
                    polyCoefficients.col(i).reshaped(numIntervals, numRadialTerms).matrix();

                double value = radialPolyvalAndTemporalIFFT(1.0/radius_km - VR, solarLongitude,
                                                          polyCoefs, atPowersInvRadius.array());

                if (m >= 0)
                    cosineCoefficients(l, absM) = value;
                else
                    sineCoefficients(l, absM) = value;
            }
        }
        else
        {
            for (int i = 0; i < polyCoefficients.cols(); ++i)
            {
                const int l = atDegreeAndOrder(0, i);
                const int m = atDegreeAndOrder(1, i);
                const int absM = std::abs(m);

                if (l > maxDegree || absM > maxOrder)
                    continue;

                const Eigen::VectorXd polyCoefs =
                    polyCoefficients.block(0, i, numIntervals, 1).matrix();

                double value = reducedToTemporalIFFT(solarLongitude, polyCoefs);

                if (m >= 0)
                    cosineCoefficients(l, absM) = value;
                else
                    sineCoefficients(l, absM) = value;
            }

            // Add logarithmic term for 1/r² decay
            applyDecayTerm(cosineCoefficients, radius_m, refRadius_m);
        }
    }

    /// Apply logarithmic decay term for 1/r² behavior when radius exceeds reference radius
    /// @param cosineCoefficients Matrix of cosine coefficients to modify
    /// @param radius_m Current radius in meters
    /// @param referenceRadius_m Reference radius in meters
    static void applyDecayTerm(Eigen::MatrixXd& cosineCoefficients, double radius_m, double referenceRadius_m)
    {
        // Convert to km for the logarithmic calculation (consistent with existing formula)
        const double radius_km = radius_m / 1000.0;
        const double referenceRadius_km = referenceRadius_m / 1000.0;

        // Apply logarithmic decay term to C(0,0) coefficient
        cosineCoefficients(0, 0) += 2.0 * std::log2((referenceRadius_km < 1.0e-10) ? 1.0/radius_km : referenceRadius_km/radius_km);
    }
};

class ComaDatasetTransformer
{
public:
    static ComaStokesDataset transformPolyToStokes(
        const ComaPolyDataset& polyDataset,
        const std::vector<double>& radii_m,
        const std::vector<double>& solLongitudes_deg,
        const int requestedMaxDegree = -1,
        const int requestedMaxOrder = -1)
    {
        // Validate inputs
        validateTransformInputs(polyDataset, radii_m, solLongitudes_deg,
                               requestedMaxDegree, requestedMaxOrder);

        // Determine effective maxima
        auto [effNmax, effMmax] = determineEffectiveMaxima(
            polyDataset, requestedMaxDegree, requestedMaxOrder);

        // Build file metadata
        auto files = buildFileMeta(polyDataset);

        // Create empty Stokes dataset
        ComaStokesDataset stokesDataset = ComaStokesDataset::create(
            std::move(files), radii_m, solLongitudes_deg, effNmax);

        // Fill dataset with computed coefficients
        fillStokesDataset(polyDataset, stokesDataset, effNmax, effMmax);

        return stokesDataset;
    }

private:
    static void validateTransformInputs(
        const ComaPolyDataset& polyDataset,
        const std::vector<double>& radii_m,
        const std::vector<double>& solLongitudes_deg,
        const int requestedMaxDegree,
        const int requestedMaxOrder)
    {
        if (radii_m.empty() || solLongitudes_deg.empty())
            throw std::invalid_argument("transformPolyToStokes: radii and longitudes must be non-empty.");
        if (requestedMaxDegree < -1 || requestedMaxOrder < -1)
            throw std::invalid_argument("transformPolyToStokes: requested maxima must be >= -1.");
        if (polyDataset.getNumFiles() == 0)
            throw std::invalid_argument("transformPolyToStokes: polyDataset has no files.");
    }

    static std::pair<int, int> determineEffectiveMaxima(
        const ComaPolyDataset& polyDataset,
        const int requestedMaxDegree,
        const int requestedMaxOrder)
    {
        int globalMaxDeg = 0;
        int globalMaxOrd = 0;

        const std::size_t F = polyDataset.getNumFiles();
        std::vector<int> perFileMaxDeg(F, 0);
        std::vector<int> perFileMaxOrd(F, 0);

        for (std::size_t f = 0; f < F; ++f)
        {
            const int fMaxDeg = polyDataset.getMaxDegreeSH(f);
            const int fMaxOrd = polyDataset.getSHDegreeAndOrderIndices(f).row(1).abs().maxCoeff();

            perFileMaxDeg[f] = fMaxDeg;
            perFileMaxOrd[f] = fMaxOrd;

            globalMaxDeg = std::max(globalMaxDeg, fMaxDeg);
            globalMaxOrd = std::max(globalMaxOrd, fMaxOrd);
        }

        const int effNmax = requestedMaxDegree < 0 ? globalMaxDeg : requestedMaxDegree;
        const int effMmax = requestedMaxOrder < 0 ? globalMaxOrd : requestedMaxOrder;

        if (effNmax < 0 || effMmax < 0)
            throw std::invalid_argument("determineEffectiveMaxima: resolved negative maxima.");
        // Ensure every file can support requested maxima
        for (std::size_t f = 0; f < F; ++f)
        {
            if (effNmax > perFileMaxDeg[f] || effMmax > perFileMaxOrd[f])
            {
                std::ostringstream oss;
                oss << "Requested (nmax=" << effNmax << ", mmax=" << effMmax
                    << ") exceeds availability in file #" << f
                    << " [" << polyDataset.getFileMeta(f).sourcePath << "]: "
                    << "(maxDegree=" << perFileMaxDeg[f]
                    << ", maxOrder=" << perFileMaxOrd[f] << ").";
                throw std::invalid_argument(oss.str());
            }
        }

        return {effNmax, effMmax};
    }

    static std::vector<ComaStokesDataset::FileMeta> buildFileMeta(
        const ComaPolyDataset& polyDataset)
    {
        std::vector<ComaStokesDataset::FileMeta> files;
        files.reserve(polyDataset.getNumFiles());
        for (std::size_t f = 0; f < polyDataset.getNumFiles(); ++f)
        {
            ComaStokesDataset::FileMeta fm{};
            fm.source_tag = polyDataset.getFileMeta(f).sourcePath;
            fm.referenceRadius = polyDataset.getFileMeta(f).referenceRadius;
            if (!polyDataset.getFileMeta(f).timePeriods.empty())
            {
                fm.start_epoch = polyDataset.getFileMeta(f).timePeriods.front().first;
                fm.end_epoch   = polyDataset.getFileMeta(f).timePeriods.front().second;
            }
            else
            {
                fm.start_epoch = 0.0;
                fm.end_epoch   = 0.0;
            }
            files.push_back(std::move(fm));
        }
        return files;
    }

    static void fillStokesDataset(
        const ComaPolyDataset& polyDataset,
        ComaStokesDataset& stokesDataset,
        const int effNmax,
        const int effMmax)
    {
        Eigen::MatrixXd C(effNmax + 1, effMmax + 1);
        Eigen::MatrixXd S(effNmax + 1, effMmax + 1);

        for (std::size_t f = 0; f < stokesDataset.nFiles(); ++f)
        {
            const auto& P   = polyDataset.getPolyCoefficients(f);
            const auto& nm  = polyDataset.getSHDegreeAndOrderIndices(f);
            const auto& pw  = polyDataset.getPowersInvRadius(f);
            const double R0 = polyDataset.getReferenceRadius(f); // (unit as in files; matches old code)

            for (std::size_t ri = 0; ri < stokesDataset.nRadii(); ++ri)
            {
                const double r_m = stokesDataset.radii()[ri]; // radius in meters

                for (std::size_t li = 0; li < stokesDataset.nLongitudes(); ++li)
                {
                    const double Lrad = stokesDataset.lons()[li] * M_PI / 180.0;

                    // Evaluate
                    StokesCoefficientsEvaluator::evaluate2D(
                        r_m, Lrad, P.array(), nm, pw, R0, C, S, effNmax, effMmax);

                    // Store back
                    for (int n = 0; n <= effNmax; ++n)
                    {
                        const int m_hi = std::min(n, effMmax);
                        for (int m = 0; m <= m_hi; ++m)
                        {
                            if (n >= C.rows() || m >= C.cols())
                                throw std::runtime_error("Evaluator returned C/S smaller than requested (nmax,mmax).");
                            stokesDataset.setCoeff(f, ri, li, n, m, C(n, m), S(n, m));
                        }
                    }
                }
            }
        }
    }
};

// ============= High-Level Processing Interface =============

class ComaModelFileProcessor
{
public:
    enum class FileType
    {
        PolyCoefficients,
        StokesCoefficients
    };

    // Constructor for polynomial coefficient files
    explicit ComaModelFileProcessor(std::vector<std::string> filePaths)
        : filePaths_(std::move(filePaths)), fileType_(FileType::PolyCoefficients)
    {
        if (filePaths_.empty())
            throw std::invalid_argument("ComaModelFileProcessor: empty file list");
    }

    // Constructor for Stokes coefficient files (SH files)
    explicit ComaModelFileProcessor(const std::string& inputDir, const std::string& prefix = "stokes")
        : fileType_(FileType::StokesCoefficients), shInputDir_(inputDir), shPrefix_(prefix)
    {
        // Load the SH dataset immediately
        preloadedSHDataset_ = ComaStokesDatasetReader::readFromCsvFolder(inputDir, prefix);
    }

    // Create poly dataset from files (only available for poly coefficient files)
    ComaPolyDataset createPolyCoefDataset() const
    {
        if (fileType_ != FileType::PolyCoefficients)
        {
            throw std::runtime_error("createPolyCoefDataset: not available when processor is constructed from SH files. "
                                   "Use a processor constructed from polynomial coefficient files instead.");
        }
        return ComaPolyDatasetReader::readFromFiles(filePaths_);
    }

    // Create Stokes dataset
    ComaStokesDataset createSHDataset(
        const std::vector<double>& radii_m,
        const std::vector<double>& solLongitudes_deg,
        const int requestedMaxDegree = -1,
        const int requestedMaxOrder = -1) const
    {
        if (fileType_ == FileType::StokesCoefficients)
        {
            // When constructed from SH files, return the preloaded dataset
            // Note: radii_m, solLongitudes_deg, and degree/order parameters are ignored
            // as they are already determined by the loaded SH files
            return preloadedSHDataset_;
        }
        else
        {
            // When constructed from poly files, transform as before
            const ComaPolyDataset polyDataset = createPolyCoefDataset();
            return ComaDatasetTransformer::transformPolyToStokes(
                polyDataset, radii_m, solLongitudes_deg,
                requestedMaxDegree, requestedMaxOrder);
        }
    }

    // Create SH files (combines dataset creation and writing)
    void createSHFiles(
        const std::string& outputDir,
        const std::vector<double>& radii_m,
        const std::vector<double>& solLongitudes_deg,
        const int requestedMaxDegree = -1,
        const int requestedMaxOrder = -1) const
    {
        const ComaStokesDataset stokesDataset = createSHDataset(
            radii_m, solLongitudes_deg,
            requestedMaxDegree, requestedMaxOrder);

        ComaStokesDatasetWriter::writeCsvAll(stokesDataset, outputDir);
    }

    // Get the file type of this processor
    FileType getFileType() const
    {
        return fileType_;
    }

private:
    std::vector<std::string> filePaths_;
    FileType fileType_;

    // For SH file processor
    std::string shInputDir_;
    std::string shPrefix_;
    ComaStokesDataset preloadedSHDataset_;
};

// Stream operator for FileType enum (needed for Boost Test)
inline std::ostream& operator<<(std::ostream& os, const ComaModelFileProcessor::FileType& type)
{
    switch (type)
    {
        case ComaModelFileProcessor::FileType::PolyCoefficients:
            return os << "PolyCoefficients";
        case ComaModelFileProcessor::FileType::StokesCoefficients:
            return os << "StokesCoefficients";
        default:
            return os << "Unknown";
    }
}






//  List of wind models available in simulations
/*
 *  List of wind models available in simulations. Wind models not defined by this
 *  given enum cannot be used for automatic model setup.
 */
enum WindModelTypes { constant_wind_model, custom_wind_model, coma_wind_model };

//  Class for providing settings for wind model.
/*
 *  Class for providing settings for automatic wind model creation. This class is a
 *  functional (base) class for settings of wind models that require no information in
 *  addition to their type. Wind model classes defining requiring additional information
 *  must be created using an object derived from this class.
 */
//! @get_docstring(WindModelSettings.__docstring__)
class WindModelSettings
{
public:
    //  Constructor
    /*
     * Constructor
     * \param windModelType Type of wind model that is to be created
     */
    WindModelSettings( const WindModelTypes windModelType,
                       const reference_frames::AerodynamicsReferenceFrames associatedFrame = reference_frames::vertical_frame ):
        windModelType_( windModelType ), associatedFrame_( associatedFrame )
    {
    }

    //  Destructor
    virtual ~WindModelSettings( )
    {
    }

    //  Function to retrieve type of wind model that is to be created
    /*
     * Function to retrieve type of wind model that is to be created
     * \return Type of wind model that is to be created
     */
    WindModelTypes getWindModelType( )
    {
        return windModelType_;
    }

    reference_frames::AerodynamicsReferenceFrames getAssociatedFrame( )
    {
        return associatedFrame_;
    }

protected:
    //  Type of wind model that is to be created
    WindModelTypes windModelType_;

    reference_frames::AerodynamicsReferenceFrames associatedFrame_;
};

class ConstantWindModelSettings : public WindModelSettings
{
public:
    ConstantWindModelSettings( const Eigen::Vector3d constantWindVelocity,
                               const reference_frames::AerodynamicsReferenceFrames associatedFrame = reference_frames::vertical_frame ):
        WindModelSettings( constant_wind_model, associatedFrame ), constantWindVelocity_( constantWindVelocity )
    {
    }

    Eigen::Vector3d getConstantWindVelocity( )
    {
        return constantWindVelocity_;
    }

private:
    Eigen::Vector3d constantWindVelocity_;
};

//  Class to define settings for a custom, user-defined, wind model
class CustomWindModelSettings : public WindModelSettings
{
public:
    //  Constructor
    /*
     * Constructor
     * \param windFunction Function that returns wind vector as a function of altitude, longitude, latitude and time (in that
     * order).
     */
    CustomWindModelSettings( const std::function< Eigen::Vector3d( const double, const double, const double, const double ) > windFunction,
                             const reference_frames::AerodynamicsReferenceFrames associatedFrame = reference_frames::vertical_frame ):
        WindModelSettings( custom_wind_model, associatedFrame ), windFunction_( windFunction )
    {
    }

    //  Destructor
    ~CustomWindModelSettings( )
    {
    }

    //  Function to retrieve function that returns wind vector as a function of altitude, longitude, latitude and time
    /*
     * Function to retrieve function that returns wind vector as a function of altitude, longitude, latitude and time
     * \return Function that returns wind vector as a function of altitude, longitude, latitude and time
     */
    std::function< Eigen::Vector3d( const double, const double, const double, const double ) > getWindFunction( )
    {
        return windFunction_;
    }

    //  Function to reset function that returns wind vector as a function of altitude, longitude, latitude and time
    /*
     * Function to reset function that returns wind vector as a function of altitude, longitude, latitude and time
     * \param windFunction New function that returns wind vector as a function of altitude, longitude, latitude and time
     */
    void setWindFunction( const std::function< Eigen::Vector3d( const double, const double, const double, const double ) > windFunction )
    {
        windFunction_ = windFunction;
    }

protected:
    //  Function that returns wind vector as a function of altitude, longitude, latitude and time (in that order).
    std::function< Eigen::Vector3d( const double, const double, const double, const double ) > windFunction_;
};

//  List of atmosphere models available in simulations
/*
 *  List of atmosphere models available in simulations. Atmosphere models not defined by this
 *  given enum cannot be used for automatic model setup.
 */
enum AtmosphereTypes
{
    exponential_atmosphere,
    custom_constant_temperature_atmosphere,
    tabulated_atmosphere,
    nrlmsise00,
    mars_dtm_atmosphere,
    scaled_atmosphere,
    coma_model
};

//  Class for providing settings for atmosphere model.
/*
 *  Class for providing settings for automatic atmosphere model creation. This class is a
 *  functional (base) class for settings of atmosphere models that require no information in
 *  addition to their type. Atmosphere model classes defining requiring additional information
 *  must be created using an object derived from this class.
 */
//! @get_docstring(AtmosphereSettings.__docstring__)
class AtmosphereSettings
{
public:
    //  Constructor, sets type of atmosphere model.
    /*
     *  Constructor, sets type of atmosphere model. Settings for atmosphere models requiring
     *  additional information should be defined in a derived class.
     *  \param atmosphereType Type of atmosphere model that is to be created.
     */
    AtmosphereSettings( const AtmosphereTypes atmosphereType ):
        atmosphereType_( atmosphereType )
    {
    }

    //  Destructor
    virtual ~AtmosphereSettings( )
    {
    }

    //  Function to return type of atmosphere model that is to be created.
    /*
     *  Function to return type of atmosphere model that is to be created.
     *  \return Type of atmosphere model that is to be created.
     */
    AtmosphereTypes getAtmosphereType( )
    {
        return atmosphereType_;
    }

    //  Function to return settings for the atmosphere's wind model.
    /*
     *  Function to return settings for the atmosphere's wind model.
     *  \return Settings for the atmosphere's wind model.
     */
    std::shared_ptr< WindModelSettings > getWindSettings( )
    {
        return windSettings_;
    }

    //  Function to (re)set settings for the atmosphere's wind model.
    /*
     *  Function to (re)set settings for the atmosphere's wind model.
     *  \param windSettings Settings for the atmosphere's wind model.
     */
    void setWindSettings( const std::shared_ptr< WindModelSettings > windSettings )
    {
        windSettings_ = windSettings;
    }

private:
    //   Type of atmosphere model that is to be created.
    AtmosphereTypes atmosphereType_;

    //  Settings for the atmosphere's wind model.
    std::shared_ptr< WindModelSettings > windSettings_;
};

//  AtmosphereSettings for defining an exponential atmosphere.
//! @get_docstring(ExponentialAtmosphereSettings.__docstring__)
class ExponentialAtmosphereSettings : public AtmosphereSettings
{
public:
    //  Default constructor.
    /*
     *  Default constructor, taking full atmosphere model parameters.
     *  \param densityScaleHeight Scale height for density profile of atmosphere.
     *  \param constantTemperature Constant atmospheric temperature.
     *  \param densityAtZeroAltitude Atmospheric density at ground level.
     *  \param specificGasConstant Specific gas constant for (constant) atmospheric chemical
     *  composition.
     *  \param ratioOfSpecificHeats Ratio of specific heats for (constant) atmospheric chemical
     *  composition.
     */
    ExponentialAtmosphereSettings( const double densityScaleHeight,
                                   const double constantTemperature,
                                   const double densityAtZeroAltitude,
                                   const double specificGasConstant = physical_constants::SPECIFIC_GAS_CONSTANT_AIR,
                                   const double ratioOfSpecificHeats = 1.4 ):
        AtmosphereSettings( exponential_atmosphere ), densityScaleHeight_( densityScaleHeight ),
        constantTemperature_( constantTemperature ), densityAtZeroAltitude_( densityAtZeroAltitude ),
        specificGasConstant_( specificGasConstant ), ratioOfSpecificHeats_( ratioOfSpecificHeats ),
        bodyWithPredefinedExponentialAtmosphere_( undefined_body )
    {
    }

    //  Default constructor.
    /*
     *  Default constructor, taking only the name of the body for which to load the predefined
     *  exponential atmosphere model parameters.
     *  \param bodyWithPredefinedExponentialAtmosphere Enumeration denoting the name of the body for which the
     *  predefined atmosphere model is to be loaded.
     */
    ExponentialAtmosphereSettings( const BodiesWithPredefinedExponentialAtmospheres bodyWithPredefinedExponentialAtmosphere ):
        AtmosphereSettings( exponential_atmosphere ), bodyWithPredefinedExponentialAtmosphere_( bodyWithPredefinedExponentialAtmosphere )
    {
        // Check that the body name inserted is available
        switch(bodyWithPredefinedExponentialAtmosphere)
        {
            case earth:
            case mars:
                // all is good
                break;
            default:
                throw std::runtime_error(
                        "Error while creating exponential atmosphere. The body name provided "
                        "does not match any predefined atmosphere model. Available models for: "
                        "Earth, Mars." );
        }
    }

    //  Function to return scale heigh for density profile of atmosphere.
    /*
     *  Function to return scale heigh for density profile of atmosphere.
     *  \return Scale heigh for density profile of atmosphere.
     */
    double getDensityScaleHeight( )
    {
        return densityScaleHeight_;
    }

    //  Function to return constant atmospheric temperature.
    /*
     *  Function to return constant atmospheric temperature.
     *  \return Constant atmospheric temperature.
     */
    double getConstantTemperature( )
    {
        return constantTemperature_;
    }

    //  Function to return atmospheric density at ground level.
    /*
     *  Function to return atmospheric density at ground level.
     *  \return Atmospheric density at ground level.
     */
    double getDensityAtZeroAltitude( )
    {
        return densityAtZeroAltitude_;
    }

    //  Function to return specific gas constant for (constant) atmospheric chemical
    /*
     *  Function to return specific gas constant for (constant) atmospheric chemical
     *  \return Specific gas constant for (constant) atmospheric chemical
     */
    double getSpecificGasConstant( )
    {
        return specificGasConstant_;
    }

    //  Function to return ratio of specific heats for (constant) atmospheric chemical
    /*
     *  Function to return ratio of specific heats for (constant) atmospheric chemical
     *  \return Specific gas constant for (constant) atmospheric chemical
     */
    double getRatioOfSpecificHeats( )
    {
        return ratioOfSpecificHeats_;
    }

    //  Function to return the name of the body for which to load the predefined
    //  atmosphere model parameters.
    BodiesWithPredefinedExponentialAtmospheres getBodyName( )
    {
        return bodyWithPredefinedExponentialAtmosphere_;
    }

private:
    //  Scale heigh for density profile of atmosphere.
    double densityScaleHeight_;

    //  Constant atmospheric temperature.
    double constantTemperature_;

    //  Atmospheric density at ground level.
    double densityAtZeroAltitude_;

    //  Specific gas constant for (constant) atmospheric chemical
    double specificGasConstant_;

    //  Ratio of specific heats for (constant) atmospheric chemical
    double ratioOfSpecificHeats_;

    //  Enumeration denoting the name of the body for which to load the predefined
    //  atmosphere model.
    BodiesWithPredefinedExponentialAtmospheres bodyWithPredefinedExponentialAtmosphere_;
};

//  AtmosphereSettings for defining custom constant temperature atmosphere.
class CustomConstantTemperatureAtmosphereSettings : public AtmosphereSettings
{
public:
    //  Typedef for density function.
    typedef std::function< double( const double, const double, const double, const double ) > DensityFunction;

    //  Default constructor.
    /*
     *  Default constructor setting all parameters manually.
     *  \param densityFunction Function to retireve the density at the current altitude,
     *      longitude, latitude and time.
     *  \param constantTemperature Constant atmospheric temperature.
     *  \param specificGasConstant The constant specific gas constant of the atmosphere.
     *  \param ratioOfSpecificHeats The constant ratio of specific heats of the atmosphere.
     */
    CustomConstantTemperatureAtmosphereSettings( const DensityFunction& densityFunction,
                                                 const double constantTemperature,
                                                 const double specificGasConstant = physical_constants::SPECIFIC_GAS_CONSTANT_AIR,
                                                 const double ratioOfSpecificHeats = 1.4 ):
        AtmosphereSettings( custom_constant_temperature_atmosphere ), densityFunction_( densityFunction ),
        constantTemperature_( constantTemperature ), specificGasConstant_( specificGasConstant ),
        ratioOfSpecificHeats_( ratioOfSpecificHeats )
    {
    }

    //  Constructor.
    /*
     *  Constructor which uses one of the built-in density functions as input.
     *  \param densityFunctionType Enumeration denoting which density function to implement.
     *  \param constantTemperature Constant atmospheric temperature.
     *  \param specificGasConstant The constant specific gas constant of the atmosphere.
     *  \param ratioOfSpecificHeats The constant ratio of specific heats of the atmosphere.
     *  \param modelSpecificParameters Vector of parameters to be used to set up the density
     *      function. Both meaning and number of parameters depends on the model.
     */
    CustomConstantTemperatureAtmosphereSettings( const AvailableConstantTemperatureAtmosphereModels densityFunctionType,
                                                 const double constantTemperature,
                                                 const double specificGasConstant,
                                                 const double ratioOfSpecificHeats,
                                                 const std::vector< double >& modelSpecificParameters ):
        AtmosphereSettings( custom_constant_temperature_atmosphere ), densityFunctionType_( densityFunctionType ),
        constantTemperature_( constantTemperature ), specificGasConstant_( specificGasConstant ),
        ratioOfSpecificHeats_( ratioOfSpecificHeats ), modelSpecificParameters_( modelSpecificParameters )
    {
    }

    //  Get the function to compute the density at the current conditions.
    /*
     *  Get the function to compute the density at the current conditions.
     *  \return Function to compute the density at the current conditions.
     */
    DensityFunction getDensityFunction( )
    {
        return densityFunction_;
    }

    //  Get the type of function to compute the density at the current conditions.
    /*
     *  Get the type of function to compute the density at the current conditions.
     *  \return Type of function to compute the density at the current conditions.
     */
    AvailableConstantTemperatureAtmosphereModels getDensityFunctionType( )
    {
        return densityFunctionType_;
    }

    //  Get constant temperature.
    /*
     *  Returns the atmospheric temperature (constant, property of exponential atmosphere) in
     *  Kelvin.
     *  \return constantTemperature Constant atmospheric temperature in exponential atmosphere.
     */
    double getConstantTemperature( )
    {
        return constantTemperature_;
    }

    //  Get specific gas constant.
    /*
     *  Returns the specific gas constant of the atmosphere in J/(kg K), its value is assumed constant,
     *  due to the assumption of constant atmospheric composition.
     *  \return Specific gas constant in exponential atmosphere.
     */
    double getSpecificGasConstant( )
    {
        return specificGasConstant_;
    }

    //  Get ratio of specific heats.
    /*
     *  Returns the ratio of specific hears of the atmosphere, its value is assumed constant,
     *  due to the assumption of constant atmospheric composition.
     *  \return Ratio of specific heats exponential atmosphere.
     */
    double getRatioOfSpecificHeats( )
    {
        return ratioOfSpecificHeats_;
    }

    //  Get model specific parameters.
    /*
     *  Get model specific parameters.
     *  \return Vector of parameters to be used to set up the density function.
     */
    std::vector< double > getModelSpecificParameters( )
    {
        return modelSpecificParameters_;
    }

private:
    //  Function to compute the density at the current conditions.
    /*
     *  Function to compute the density at the current conditions. Note that the independent
     *  variables need to be altitude, longitude, latitude and time, in this precise order.
     */
    DensityFunction densityFunction_;

    //  Enumeration denoting which density function to implement.
    /*
     *  Enumeration denoting which density function to implement
     */
    AvailableConstantTemperatureAtmosphereModels densityFunctionType_;

    //  Constant temperature.
    /*
     *  The atmospheric temperature (constant, property of exponential atmosphere) in Kelvin.
     */
    const double constantTemperature_;

    //  Specific gas constant.
    /*
     *  Specific gas constant of the atmosphere, its value is assumed constant, due to the
     *  assumption of constant atmospheric composition.
     */
    const double specificGasConstant_;

    //  Ratio of specific heats at constant pressure and constant volume.
    /*
     *  Ratio of specific heats of the atmosphere at constant pressure and constant volume.
     *  This value is set to a constant, implying constant atmospheric composition.
     */
    const double ratioOfSpecificHeats_;

    //  Vector of parameters to be used to set up the density function.
    /*
     *  Vector of parameters to be used to set up the density function. Both meaning and number of parameters depends on the model.
     */
    std::vector< double > modelSpecificParameters_;
};

//  AtmosphereSettings for defining an NRLMSISE00 atmosphere reading space weather data from a text file.
class NRLMSISE00AtmosphereSettings : public AtmosphereSettings
{
public:
    //  Constructor.
    /*
     *  Constructor.
     *  \param spaceWeatherFile File containing space weather data, as in
     *  https://celestrak.com/SpaceData/sw19571001.txt
     */
    NRLMSISE00AtmosphereSettings( const std::string& spaceWeatherFile,
                                  const bool useStormConditions = false,
                                  const bool useAnomalousOxygen = true ):
        AtmosphereSettings( nrlmsise00 ), spaceWeatherFile_( spaceWeatherFile ), useStormConditions_( useStormConditions ),
        useAnomalousOxygen_( useAnomalousOxygen )
    {
    }

    //  Function to return file containing space weather data.
    /*
     *  Function to return file containing space weather data.
     *  \return Filename containing space weather data.
     */
    std::string getSpaceWeatherFile( )
    {
        return spaceWeatherFile_;
    }

    //  Function to return geomagnetic activity setting.
    /*
     *  Function to return geomagnetic activity setting.
     *  \return Geomagnetic activity value.
     */
    bool getUseStormConditions( )
    {
        return useStormConditions_;
    }

    bool getUseAnomalousOxygen( )
    {
        return useAnomalousOxygen_;
    }

private:
    //  File containing space weather data.
    /*
     *  File containing space weather data, as in https://celestrak.com/SpaceData/sw19571001.txt
     */
    std::string spaceWeatherFile_;

    //  Boolean denoting how to deal with geomagnetic activity setting.
    /*
     *  Controls the geomagnetic activity behavior.
     *  If true, it uses full vector of Ap values (switch for geomagnetic activity in nrlmsise to -1).
     *  If false, it only uses daily value of Ap (switch for geomagnetic activity in nrlmsise to 0).
     */
    bool useStormConditions_;

    bool useAnomalousOxygen_;
};

class MarsDtmAtmosphereSettings : public AtmosphereSettings
{
public:
    MarsDtmAtmosphereSettings( const std::string& spaceWeatherFile = "" ):
        AtmosphereSettings( mars_dtm_atmosphere ), spaceWeatherFile_( spaceWeatherFile )
    {
    }

    std::string getSpaceWeatherFile( )
    {
        return spaceWeatherFile_;
    }

private:
    //  File containing space weather data.
    /*
     *  File containing space weather data, as in https://celestrak.com/SpaceData/sw19571001.txt
     */
    std::string spaceWeatherFile_;
};


//  AtmosphereSettings for defining an atmosphere with tabulated data from file.
// //! @get_docstring(TabulatedAtmosphereSettings.__docstring__)
class TabulatedAtmosphereSettings : public AtmosphereSettings
{
public:
    //  Default constructor.
    /*
     *  Default constructor.
     *  \param atmosphereTableFile Map of files containing information on the atmosphere. The order of both
     *      independent and dependent parameters needs to be specified in the independentVariablesNames and
     *      dependentVariablesNames vectors, respectively. Note that specific gas constant and specific heat ratio
     *      will be given the default constant values for Earth, unless they are included in the file map.
     *  \param independentVariablesNames List of independent parameters describing the atmosphere.
     *  \param dependentVariablesNames List of dependent parameters output by the atmosphere.
     *  \param specificGasConstant The constant specific gas constant of the atmosphere.
     *  \param ratioOfSpecificHeats The constant ratio of specific heats of the atmosphere.
     *  \param boundaryHandling List of methods for interpolation behavior when independent variable is out of range.
     *  \param defaultExtrapolationValue List of default values to be used for extrapolation, in case of
     *      use_default_value or use_default_value_with_warning as methods for boundaryHandling.
     */
    TabulatedAtmosphereSettings(
            const std::map< int, std::string >& atmosphereTableFile,
            const std::vector< AtmosphereIndependentVariables >& independentVariablesNames = { altitude_dependent_atmosphere },
            const std::vector< AtmosphereDependentVariables >& dependentVariablesNames = { density_dependent_atmosphere,
                                                                                           pressure_dependent_atmosphere,
                                                                                           temperature_dependent_atmosphere },
            const double specificGasConstant = physical_constants::SPECIFIC_GAS_CONSTANT_AIR,
            const double ratioOfSpecificHeats = 1.4,
            const std::vector< interpolators::BoundaryInterpolationType >& boundaryHandling = {},
            const std::vector< std::vector< std::pair< double, double > > >& defaultExtrapolationValue = {} ):
        AtmosphereSettings( tabulated_atmosphere ), atmosphereFile_( atmosphereTableFile ),
        independentVariables_( independentVariablesNames ), dependentVariables_( dependentVariablesNames ),
        specificGasConstant_( specificGasConstant ), ratioOfSpecificHeats_( ratioOfSpecificHeats ), boundaryHandling_( boundaryHandling ),
        defaultExtrapolationValue_( defaultExtrapolationValue )
    {
    }

    //  Constructor with single boundary handling parameters.
    /*
     *  Constructor with single boundary handling parameters. The specifier is assumed to be the same for
     *  each (in)dependent variable.
     *  \param atmosphereTableFile Map of files containing information on the atmosphere. The order of both
     *      independent and dependent parameters needs to be specified in the independentVariablesNames and
     *      dependentVariablesNames vectors, respectively. Note that specific gas constant and specific heat ratio
     *      will be given the default constant values for Earth, unless they are included in the file map.
     *  \param independentVariablesNames List of independent parameters describing the atmosphere.
     *  \param dependentVariablesNames List of dependent parameters output by the atmosphere.
     *  \param specificGasConstant The constant specific gas constant of the atmosphere.
     *  \param ratioOfSpecificHeats The constant ratio of specific heats of the atmosphere.
     *  \param boundaryHandling Method for interpolation behavior when independent variable is out of range.
     *  \param defaultExtrapolationValue Default value to be used for extrapolation, in case of use_default_value or
     *      use_default_value_with_warning as methods for boundaryHandling.
     */
    TabulatedAtmosphereSettings( const std::map< int, std::string >& atmosphereTableFile,
                                 const std::vector< AtmosphereIndependentVariables >& independentVariablesNames,
                                 const std::vector< AtmosphereDependentVariables >& dependentVariablesNames,
                                 const double specificGasConstant,
                                 const double ratioOfSpecificHeats,
                                 const interpolators::BoundaryInterpolationType boundaryHandling,
                                 const double defaultExtrapolationValue = IdentityElement::getAdditionIdentity< double >( ) ):
        TabulatedAtmosphereSettings(
                atmosphereTableFile,
                independentVariablesNames,
                dependentVariablesNames,
                specificGasConstant,
                ratioOfSpecificHeats,
                std::vector< interpolators::BoundaryInterpolationType >( independentVariablesNames.size( ), boundaryHandling ),
                std::vector< std::vector< std::pair< double, double > > >(
                        dependentVariablesNames.size( ),
                        std::vector< std::pair< double, double > >(
                                independentVariablesNames.size( ),
                                std::make_pair( defaultExtrapolationValue, defaultExtrapolationValue ) ) ) )
    {
    }

    //  Constructor compatible with old version.
    /*
     *  Constructor compatible with old version.
     *  \param atmosphereTableFile File containing atmospheric properties. The file name of the atmosphere table. The file
     *      should contain four columns of data, containing altitude (first column), and the associated density, pressure and
     *      density values in the second, third and fourth columns.
     *  \param dependentVariablesNames List of dependent parameters output by the atmosphere.
     *  \param specificGasConstant The constant specific gas constant of the atmosphere.
     *  \param ratioOfSpecificHeats The constant ratio of specific heats of the atmosphere.
     *  \param boundaryHandling Method for interpolation behavior when independent variable is out of range.
     *  \param defaultExtrapolationValue Default value to be used for extrapolation, in case of use_default_value or
     *      use_default_value_with_warning as methods for boundaryHandling.
     */
    TabulatedAtmosphereSettings(
            const std::string& atmosphereTableFile,
            const std::vector< AtmosphereDependentVariables >& dependentVariablesNames = { density_dependent_atmosphere,
                                                                                           pressure_dependent_atmosphere,
                                                                                           temperature_dependent_atmosphere },
            const double specificGasConstant = physical_constants::SPECIFIC_GAS_CONSTANT_AIR,
            const double ratioOfSpecificHeats = 1.4,
            const interpolators::BoundaryInterpolationType boundaryHandling = interpolators::use_boundary_value,
            const double defaultExtrapolationValue = IdentityElement::getAdditionIdentity< double >( ) ):
        TabulatedAtmosphereSettings( { { 0, atmosphereTableFile } },
                                     { altitude_dependent_atmosphere },
                                     dependentVariablesNames,
                                     specificGasConstant,
                                     ratioOfSpecificHeats,
                                     { boundaryHandling },
                                     std::vector< std::vector< std::pair< double, double > > >(
                                             dependentVariablesNames.size( ),
                                             std::vector< std::pair< double, double > >(
                                                     1,
                                                     std::make_pair( defaultExtrapolationValue, defaultExtrapolationValue ) ) ) )
    {
    }

    //  Constructor with no specific gas constant nor ratio of specific heats.
    /*
     *  Constructor with no specific gas constant nor ratio of specific heats.
     *  \param atmosphereTableFile Map of files containing information on the atmosphere. The order of both
     *      independent and dependent parameters needs to be specified in the independentVariablesNames and
     *      dependentVariablesNames vectors, respectively. Note that specific gas constant and specific heat ratio
     *      will be given the default constant values for Earth, unless they are included in the file map.
     *  \param independentVariablesNames List of independent parameters describing the atmosphere.
     *  \param dependentVariablesNames List of dependent parameters output by the atmosphere.
     *  \param boundaryHandling List of methods for interpolation behavior when independent variable is out of range.
     *  \param defaultExtrapolationValue List of default values to be used for extrapolation, in case of
     *      use_default_value or use_default_value_with_warning as methods for boundaryHandling.
     */
    TabulatedAtmosphereSettings( const std::map< int, std::string >& atmosphereTableFile,
                                 const std::vector< AtmosphereIndependentVariables >& independentVariablesNames,
                                 const std::vector< AtmosphereDependentVariables >& dependentVariablesNames,
                                 const std::vector< interpolators::BoundaryInterpolationType >& boundaryHandling,
                                 const std::vector< std::vector< std::pair< double, double > > >& defaultExtrapolationValue = {} ):
        AtmosphereSettings( tabulated_atmosphere ), atmosphereFile_( atmosphereTableFile ),
        independentVariables_( independentVariablesNames ), dependentVariables_( dependentVariablesNames ),
        specificGasConstant_( physical_constants::SPECIFIC_GAS_CONSTANT_AIR ), ratioOfSpecificHeats_( 1.4 ),
        boundaryHandling_( boundaryHandling ), defaultExtrapolationValue_( defaultExtrapolationValue )
    {
    }

    //  Constructor with no specific gas constant nor ratio of specific heats.
    /*
     *  Constructor with no specific gas constant nor ratio of specific heats. These two values will be given
     *  the default Earth value, or are specified inside the atmosphere table file (and thus, inside the
     *  dependent variables vector).
     *  \param atmosphereTableFile Map of files containing information on the atmosphere. The order of both
     *      independent and dependent parameters needs to be specified in the independentVariablesNames and
     *      dependentVariablesNames vectors, respectively. Note that specific gas constant and specific heat ratio
     *      will be given the default constant values for Earth, unless they are included in the file map.
     *  \param independentVariablesNames List of independent parameters describing the atmosphere.
     *  \param dependentVariablesNames List of dependent parameters output by the atmosphere.
     *  \param boundaryHandling List of methods for interpolation behavior when independent variable is out of range.
     *  \param defaultExtrapolationValue List of default values to be used for extrapolation, in case of
     *      use_default_value or use_default_value_with_warning as methods for boundaryHandling.
     */
    TabulatedAtmosphereSettings( const std::map< int, std::string >& atmosphereTableFile,
                                 const std::vector< AtmosphereIndependentVariables >& independentVariablesNames,
                                 const std::vector< AtmosphereDependentVariables >& dependentVariablesNames,
                                 const std::vector< interpolators::BoundaryInterpolationType >& boundaryHandling,
                                 const std::vector< double >& defaultExtrapolationValue ):
        AtmosphereSettings( tabulated_atmosphere ), atmosphereFile_( atmosphereTableFile ),
        independentVariables_( independentVariablesNames ), dependentVariables_( dependentVariablesNames ),
        specificGasConstant_( physical_constants::SPECIFIC_GAS_CONSTANT_AIR ), ratioOfSpecificHeats_( 1.4 ),
        boundaryHandling_( boundaryHandling )
    {
        // Assign default values
        defaultExtrapolationValue_.resize( dependentVariablesNames.size( ) );
        for(unsigned int i = 0; i < dependentVariablesNames.size( ); i++)
        {
            for(unsigned int j = 0; j < independentVariablesNames.size( ); j++)
            {
                if(boundaryHandling_.at( j ) == interpolators::use_default_value ||
                    boundaryHandling_.at( j ) == interpolators::use_default_value_with_warning)
                {
                    defaultExtrapolationValue_.at( i ).push_back(
                            std::make_pair( defaultExtrapolationValue.at( i ), defaultExtrapolationValue.at( i ) ) );
                }
                else
                {
                    defaultExtrapolationValue_.at( i ).push_back( std::make_pair( IdentityElement::getAdditionIdentity< double >( ),
                                                                                  IdentityElement::getAdditionIdentity< double >( ) ) );
                }
            }
        }
    }

    //  Constructor with no specific gas constant nor ratio of specific heats, and with
    //  single boundary handling parameters.
    /*
     *  Constructor with no specific gas constant nor ratio of specific heats. These two values will be given
     *  the default Earth value, or are specified inside the atmosphere table file (and thus, inside the
     *  dependent variables vector). Only one boundary handling parameter is specified, which is then repeated for
     *  dimension.
     *  \param atmosphereTableFile Map of files containing information on the atmosphere. The order of both
     *      independent and dependent parameters needs to be specified in the independentVariablesNames and
     *      dependentVariablesNames vectors, respectively. Note that specific gas constant and specific heat ratio
     *      will be given the default constant values for Earth, unless they are included in the file map.
     *  \param independentVariablesNames List of independent parameters describing the atmosphere.
     *  \param dependentVariablesNames List of dependent parameters output by the atmosphere.
     *  \param boundaryHandling Method for interpolation behavior when independent variable is out of range.
     *  \param defaultExtrapolationValue Default value to be used for extrapolation, in case of use_default_value or
     *      use_default_value_with_warning as methods for boundaryHandling.
     */
    TabulatedAtmosphereSettings( const std::map< int, std::string >& atmosphereTableFile,
                                 const std::vector< AtmosphereIndependentVariables >& independentVariablesNames,
                                 const std::vector< AtmosphereDependentVariables >& dependentVariablesNames,
                                 const interpolators::BoundaryInterpolationType boundaryHandling,
                                 const double defaultExtrapolationValue = IdentityElement::getAdditionIdentity< double >( ) ):
        TabulatedAtmosphereSettings(
                atmosphereTableFile,
                independentVariablesNames,
                dependentVariablesNames,
                physical_constants::SPECIFIC_GAS_CONSTANT_AIR,
                1.4,
                std::vector< interpolators::BoundaryInterpolationType >( independentVariablesNames.size( ), boundaryHandling ),
                std::vector< std::vector< std::pair< double, double > > >(
                        dependentVariablesNames.size( ),
                        std::vector< std::pair< double, double > >(
                                independentVariablesNames.size( ),
                                std::make_pair( defaultExtrapolationValue, defaultExtrapolationValue ) ) ) )
    {
    }

    //  Function to return file containing atmospheric properties.
    /*
     *  Function to return file containing atmospheric properties.
     *  \return Map of filenames containing atmospheric properties.
     */
    std::map< int, std::string > getAtmosphereFile( )
    {
        return atmosphereFile_;
    }

    //  Function to return file containing atmospheric properties.
    /*
     *  Function to return file containing atmospheric properties.
     *  \return Filename containing atmospheric properties.
     */
    std::string getAtmosphereFile( const unsigned int fileIndex )
    {
        return atmosphereFile_.at( fileIndex );
    }

    //  Function to return independent variables names.
    /*
     *  Function to return independent variables names.
     *  \return Independent variables.
     */
    std::vector< AtmosphereIndependentVariables > getIndependentVariables( )
    {
        return independentVariables_;
    }

    //  Function to return dependent variables names.
    /*
     *  Function to return dependent variables names.
     *  \return Dependent variables.
     */
    std::vector< AtmosphereDependentVariables > getDependentVariables( )
    {
        return dependentVariables_;
    }

    //  Function to return specific gas constant of the atmosphere.
    /*
     *  Function to return specific gas constant of the atmosphere.
     *  \return Specific gas constant of the atmosphere.
     */
    double getSpecificGasConstant( )
    {
        return specificGasConstant_;
    }

    //  Function to return ratio of specific heats of the atmosphere.
    /*
     *  Function to return ratio of specific heats of the atmosphere.
     *  \return Ratio of specific heats of the atmosphere at constant pressure and constant volume.
     */
    double getRatioOfSpecificHeats( )
    {
        return ratioOfSpecificHeats_;
    }

    //  Function to return boundary handling method.
    /*
     *  Function to return boundary handling method.
     *  \return Boundary handling method for when independent variables are outside specified range.
     */
    std::vector< interpolators::BoundaryInterpolationType > getBoundaryHandling( )
    {
        return boundaryHandling_;
    }

    //  Function to return default extrapolation value.
    /*
     *  Function to return boundary handling method.
     *  \return Boundary handling method for when independent variables are outside specified range.
     */
    std::vector< std::vector< std::pair< double, double > > > getDefaultExtrapolationValue( )
    {
        return defaultExtrapolationValue_;
    }

private:
    //  File containing atmospheric properties.
    /*
     *  File containing atmospheric properties, file should contain
     *  columns of atmospheric data with at least density, pressure and temperature,
     *  (whose order is specified in dependentVariables), and with at least one
     *  indendent variables.
     */
    std::map< int, std::string > atmosphereFile_;

    //  A vector of strings containing the names of the independent variables contained in the atmosphere file
    /*
     * A vector of strings containing the names of the independent variables contained in the atmosphere file,
     * in the correct order (from left, being the first entry in the vector, to the right).
     */
    std::vector< AtmosphereIndependentVariables > independentVariables_;

    //  A vector of strings containing the names of the variables contained in the atmosphere file
    /*
     * A vector of strings containing the names of the variables contained in the atmosphere file,
     * in the correct order (from left, being the first entry in the vector, to the right).
     */
    std::vector< AtmosphereDependentVariables > dependentVariables_;

    //  Specific gas constant of the atmosphere.
    /*
     * Specific gas constant of the atmosphere.
     */
    double specificGasConstant_;

    //  Ratio of specific heats of the atmosphere at constant pressure and constant volume.
    /*
     *  Ratio of specific heats of the atmosphere at constant pressure and constant volume.
     */
    double ratioOfSpecificHeats_;

    //  Behavior of interpolator when independent variable is outside range.
    /*
     *  Behavior of interpolator when independent variable is outside range.
     */
    std::vector< interpolators::BoundaryInterpolationType > boundaryHandling_;

    //  Default value to be used for extrapolation.
    /*
     *  Default value to be used for extrapolation.
     */
    std::vector< std::vector< std::pair< double, double > > > defaultExtrapolationValue_;
};

class ScaledAtmosphereSettings : public AtmosphereSettings
{
public:
    ScaledAtmosphereSettings( const std::shared_ptr< AtmosphereSettings > baseSettings,
                              const double scaling,
                              const bool isScalingAbsolute ):
        AtmosphereSettings( scaled_atmosphere ), baseSettings_( baseSettings ), scaling_( [ = ]( const double ) {
            return scaling;
        } ),
        isScalingAbsolute_( isScalingAbsolute )
    {
    }

    ScaledAtmosphereSettings( const std::shared_ptr< AtmosphereSettings >& baseSettings,
                              const std::function< double( const double ) > scaling,
                              const bool isScalingAbsolute ):
        AtmosphereSettings( scaled_atmosphere ), baseSettings_( baseSettings ), scaling_( scaling ), isScalingAbsolute_( isScalingAbsolute )
    {
    }

    std::shared_ptr< AtmosphereSettings > getBaseSettings( )
    {
        return baseSettings_;
    }

    std::function< double( const double ) > getScaling( )
    {
        return scaling_;
    }

    bool getIsScalingAbsolute( ) const
    {
        return isScalingAbsolute_;
    }

protected:
    std::shared_ptr< AtmosphereSettings > baseSettings_;

    std::function< double( const double ) > scaling_;

    bool isScalingAbsolute_;
};


/**
 * \class ComaSettings
 * \brief Configuration settings for coma atmosphere models
 *
 * This class can be initialized with either polynomial coefficients or
 * pre-computed Stokes coefficients. It provides a unified interface for
 * passing data to the ComaModel while maintaining flexibility in input types.
 */
class ComaSettings final : public AtmosphereSettings
{
public:
    // Type alias for cleaner code
    using DataVariant = boost::variant<ComaPolyDataset, ComaStokesDataset>;

    /**
     * \brief Constructor with polynomial coefficient data
     * \param polyData Pre-loaded polynomial coefficient dataset
     * \param requestedDegree Maximum spherical harmonic degree (-1 for auto)
     * \param requestedOrder Maximum spherical harmonic order (-1 for auto)
     */
    explicit ComaSettings(const ComaPolyDataset& polyData,
                const int requestedDegree = -1,
                const int requestedOrder = -1)
        : AtmosphereSettings(coma_model),
          data_(polyData),
          requestedDegree_(requestedDegree),
          requestedOrder_(requestedOrder)
    {
        validateAndSetDefaults();
    }

    /**
     * \brief Constructor with Stokes coefficient data
     * \param stokesData Pre-computed Stokes coefficient dataset
     * \param requestedDegree Maximum spherical harmonic degree (-1 for auto)
     * \param requestedOrder Maximum spherical harmonic order (-1 for auto)
     */
    explicit ComaSettings(const ComaStokesDataset& stokesData,
                          const int requestedDegree = -1,
                          const int requestedOrder = -1)
        : AtmosphereSettings(coma_model),
          data_(stokesData),
          requestedDegree_(requestedDegree),
          requestedOrder_(requestedOrder)
    {
        validateAndSetDefaults();
    }

    /**
     * \brief Get the underlying data (poly or Stokes coefficients)
     * \return Variant containing either ComaPolyDataset or ComaStokesDataset
     */
    const DataVariant& getData() const
    {
        return data_;
    }

    /**
     * \brief Check if settings contain polynomial coefficient data
     */
    bool hasPolyData() const
    {
        return data_.type() == typeid(ComaPolyDataset);
    }

    /**
     * \brief Check if settings contain Stokes coefficient data
     */
    bool hasStokesData() const
    {
        return data_.type() == typeid(ComaStokesDataset);
    }

    /**
     * \brief Get polynomial dataset if available
     * \throws std::runtime_error if data is not polynomial type
     */
    const ComaPolyDataset& getPolyDataset() const
    {
        if (auto* p = boost::get<ComaPolyDataset>(&data_)) return *p;
        throw std::runtime_error("ComaSettings does not contain polynomial data");
    }

    /**
     * \brief Get Stokes dataset if available
     * \throws std::runtime_error if data is not Stokes type
     */
    const ComaStokesDataset& getStokesDataset() const
    {
        if (auto* p = boost::get<ComaStokesDataset>(&data_)) return *p;
        throw std::runtime_error("ComaSettings does not contain Stokes data");
    }

    /**
     * \brief Get requested maximum degree
     */
    int getRequestedDegree() const
    {
        return requestedDegree_;
    }

    /**
     * \brief Get requested maximum order
     */
    int getRequestedOrder() const
    {
        return requestedOrder_;
    }

    /**
     * \brief Get the effective maximum degree available in the data
     */
    int getAvailableMaxDegree() const
    {
        return availableMaxDegree_;
    }

    /**
     * \brief Get the effective maximum order available in the data
     */
    int getAvailableMaxOrder() const
    {
        return availableMaxOrder_;
    }

private:
    /**
     * \brief Validate settings and set defaults for degree/order
     */
    void validateAndSetDefaults()
    {
        // Determine available maxima from data
        if (hasPolyData())
        {
            const auto& poly = getPolyDataset();
            availableMaxDegree_ = determineMaxDegreeFromPoly(poly);
            availableMaxOrder_ = determineMaxOrderFromPoly(poly);
        }
        else if (hasStokesData())
        {
            const auto& stokes = getStokesDataset();
            availableMaxDegree_ = stokes.nmax();
            // For Stokes data, order equals degree in the triangular storage
            availableMaxOrder_ = stokes.nmax();
        }

        // Set defaults if -1
        if (requestedDegree_ < 0)
        {
            requestedDegree_ = availableMaxDegree_;
        }
        if (requestedOrder_ < 0)
        {
            requestedOrder_ = availableMaxOrder_;
        }

        // Validate requested values don't exceed available
        if (requestedDegree_ > availableMaxDegree_)
        {
            throw std::invalid_argument(
                "Requested degree " + std::to_string(requestedDegree_) +
                " exceeds available maximum " + std::to_string(availableMaxDegree_));
        }
        if (requestedOrder_ > availableMaxOrder_)
        {
            throw std::invalid_argument(
                "Requested order " + std::to_string(requestedOrder_) +
                " exceeds available maximum " + std::to_string(availableMaxOrder_));
        }
    }

    /**
     * \brief Determine maximum degree from polynomial dataset
     */
    static int determineMaxDegreeFromPoly(const ComaPolyDataset& poly)
    {
        int maxDeg = 0;
        for (std::size_t f = 0; f < poly.getNumFiles(); ++f)
        {
            maxDeg = std::max(maxDeg, poly.getMaxDegreeSH(f));
        }
        return maxDeg;
    }

    /**
     * \brief Determine maximum order from polynomial dataset
     */
    static int determineMaxOrderFromPoly(const ComaPolyDataset& poly)
    {
        int maxOrd = 0;
        for (std::size_t f = 0; f < poly.getNumFiles(); ++f)
        {
            const auto& indices = poly.getSHDegreeAndOrderIndices(f);
            int fileMaxOrd = indices.row(1).abs().maxCoeff();
            maxOrd = std::max(maxOrd, fileMaxOrd);
        }
        return maxOrd;
    }

    // Data members
    DataVariant data_;                // Holds either poly or Stokes data
    int requestedDegree_;             // User-requested max degree
    int requestedOrder_;              // User-requested max order
    int availableMaxDegree_{0};      // Maximum available in data
    int availableMaxOrder_{0};       // Maximum available in data
};


//! @get_docstring(exponentialAtmosphereSettings,2)
inline std::shared_ptr< AtmosphereSettings > exponentialAtmosphereSettings(
        const double densityScaleHeight,
        const double densityAtZeroAltitude,
        const double constantTemperature,
        const double specificGasConstant = physical_constants::SPECIFIC_GAS_CONSTANT_AIR,
        const double ratioOfSpecificHeats = 1.4 )
{
    return std::make_shared< ExponentialAtmosphereSettings >(
            densityScaleHeight,
            constantTemperature,
            densityAtZeroAltitude,
            specificGasConstant,
            ratioOfSpecificHeats );
}

//! @get_docstring(exponentialAtmosphereSettings,1)
inline std::shared_ptr< AtmosphereSettings > exponentialAtmosphereSettings( const double densityScaleHeight,
                                                                            const double densityAtZeroAltitude )
{
    return std::make_shared< ExponentialAtmosphereSettings >( densityScaleHeight, TUDAT_NAN, densityAtZeroAltitude, TUDAT_NAN, TUDAT_NAN );
}

//! @get_docstring(exponentialAtmosphereSettings,0)
inline std::shared_ptr< AtmosphereSettings > exponentialAtmosphereSettings( const std::string& bodyName )
{
    BodiesWithPredefinedExponentialAtmospheres bodyId;
    if(bodyName == "Earth")
    {
        bodyId = BodiesWithPredefinedExponentialAtmospheres::earth;
    }
    else if(bodyName == "Mars")
    {
        bodyId = BodiesWithPredefinedExponentialAtmospheres::mars;
    }
    else
    {
        throw std::runtime_error(
                "Error while creating exponential atmosphere. The body name provided "
                "does not match any predefined atmosphere model. Available models for: "
                "Earth, Mars." );
    }
    return std::make_shared< ExponentialAtmosphereSettings >( bodyId );
}

//! @get_docstring(nrlmsise00AtmosphereSettings)
inline std::shared_ptr< AtmosphereSettings > nrlmsise00AtmosphereSettings( const std::string dataFile = paths::getSpaceWeatherDataPath( ) +
                                                                                   "/sw19571001.txt",
                                                                           const bool useStormConditions = true,
                                                                           const bool useAnomalousOxygen = true )
{
    return std::make_shared< NRLMSISE00AtmosphereSettings >( dataFile, useStormConditions, useAnomalousOxygen );
}

inline std::shared_ptr< AtmosphereSettings > marsDtmAtmosphereSettings( )
{
    return std::make_shared< MarsDtmAtmosphereSettings >( );
}


typedef std::function< double( const double, const double, const double, const double ) > DensityFunction;
//! @get_docstring(customConstantTemperatureAtmosphereSettings,0)
inline std::shared_ptr< AtmosphereSettings > customConstantTemperatureAtmosphereSettings(
        const std::function< double( const double ) > densityFunction,
        const double constantTemperature,
        const double specificGasConstant = physical_constants::SPECIFIC_GAS_CONSTANT_AIR,
        const double ratioOfSpecificHeats = 1.4 )
{
    DensityFunction fullDensityFunction = [ = ]( const double altitude, const double, const double, const double ) {
        return densityFunction( altitude );
    };
    return std::make_shared< CustomConstantTemperatureAtmosphereSettings >(
            fullDensityFunction,
            constantTemperature,
            specificGasConstant,
            ratioOfSpecificHeats );
}

//! @get_docstring(customConstantTemperatureAtmosphereSettings,1)
inline std::shared_ptr< AtmosphereSettings > customConstantTemperatureAtmosphereSettings(
        const DensityFunction densityFunction,
        const double constantTemperature,
        const double specificGasConstant = physical_constants::SPECIFIC_GAS_CONSTANT_AIR,
        const double ratioOfSpecificHeats = 1.4 )
{
    return std::make_shared< CustomConstantTemperatureAtmosphereSettings >(
            densityFunction,
            constantTemperature,
            specificGasConstant,
            ratioOfSpecificHeats );
}

//! @get_docstring(scaledAtmosphereSettings,0)
inline std::shared_ptr< AtmosphereSettings > scaledAtmosphereSettings( const std::shared_ptr< AtmosphereSettings > baseSettings,
                                                                       const std::function< double( const double ) > scaling,
                                                                       const bool isScalingAbsolute )
{
    return std::make_shared< ScaledAtmosphereSettings >( baseSettings, scaling, isScalingAbsolute );
}

//! @get_docstring(scaledAtmosphereSettings,1)
inline std::shared_ptr< AtmosphereSettings > scaledAtmosphereSettings( const std::shared_ptr< AtmosphereSettings > baseSettings,
                                                                       const double scaling,
                                                                       const bool isScalingAbsolute )
{
    return std::make_shared< ScaledAtmosphereSettings >( baseSettings, scaling, isScalingAbsolute );
}

inline std::shared_ptr< AtmosphereSettings > tabulatedAtmosphereSettings(
        const std::string& atmosphereTableFile,
        const std::vector< AtmosphereDependentVariables >& dependentVariablesNames = { density_dependent_atmosphere,
                                                                                       pressure_dependent_atmosphere,
                                                                                       temperature_dependent_atmosphere },
        const double specificGasConstant = physical_constants::SPECIFIC_GAS_CONSTANT_AIR,
        const double ratioOfSpecificHeats = 1.4 )
{
    return std::make_shared< TabulatedAtmosphereSettings >( atmosphereTableFile,
                                                            dependentVariablesNames,
                                                            specificGasConstant,
                                                            ratioOfSpecificHeats,
                                                            interpolators::throw_exception_at_boundary );
}


//@get_docstring(ComaSettings,0)
inline std::shared_ptr< AtmosphereSettings > comaSettings(
        const ComaPolyDataset& polyData,
        const int requestedDegree = -1,
        const int requestedOrder = -1 )
{
    return std::make_shared< ComaSettings >( polyData, requestedDegree, requestedOrder );
}

//@get_docstring(ComaSettings,1)
inline std::shared_ptr< AtmosphereSettings > comaSettings(
        const ComaStokesDataset& stokesData,
        const int requestedDegree = -1,
        const int requestedOrder = -1 )
{
    return std::make_shared< ComaSettings >( stokesData, requestedDegree, requestedOrder );
}



// === Coma processing: factory-style helpers (header-inline) ===
inline std::shared_ptr<ComaModelFileProcessor> comaModelFileProcessorFromPolyFiles(
    const std::vector<std::string>& filePaths)
{
    return std::make_shared<ComaModelFileProcessor>(filePaths);
}

inline std::shared_ptr<ComaModelFileProcessor> comaModelFileProcessorFromSHFiles(
    const std::string& inputDir,
    const std::string& prefix = "stokes")
{
    return std::make_shared<ComaModelFileProcessor>(inputDir, prefix);
}





//! @get_docstring(customWindModelSettings)
inline std::shared_ptr< WindModelSettings > customWindModelSettings(
        const std::function< Eigen::Vector3d( const double, const double, const double, const double ) > windFunction,
        const reference_frames::AerodynamicsReferenceFrames associatedFrame = reference_frames::vertical_frame )
{
    return std::make_shared< CustomWindModelSettings >( windFunction, associatedFrame );
}

//! @get_docstring(constantWindModelSettings)
inline std::shared_ptr< WindModelSettings > constantWindModelSettings(
        const Eigen::Vector3d constantWindVelocity,
        const reference_frames::AerodynamicsReferenceFrames associatedFrame = reference_frames::vertical_frame )
{
    return std::make_shared< ConstantWindModelSettings >( constantWindVelocity, associatedFrame );
}

//  Function to create a wind model.
/*
 *  Function to create a wind model based on model-specific settings for the wind model.
 *  \param windSettings Settings for the wind model that is to be created, defined
 *  a pointer to an object of class (derived from) WindModelSettings.
 *  \param body Name of the body for which the wind model is to be created.
 *  \return Wind model created according to settings in windSettings.
 */
std::shared_ptr< aerodynamics::WindModel > createWindModel( const std::shared_ptr< WindModelSettings > windSettings,
                                                            const std::string& body );

//  Function to create an atmosphere model.
/*
 *  Function to create an atmosphere model based on model-specific settings for the atmosphere.
 *  \param atmosphereSettings Settings for the atmosphere model that is to be created, defined
 *  a pointer to an object of class (derived from) AtmosphereSettings.
 *  \param body Name of the body for which the atmosphere model is to be created.
 *  \return Atmosphere model created according to settings in atmosphereSettings.
 */
std::shared_ptr< aerodynamics::AtmosphereModel > createAtmosphereModel( const std::shared_ptr< AtmosphereSettings > atmosphereSettings,
                                                                        const std::string& body,
                                                                        const SystemOfBodies& bodies );
} // namespace simulation_setup
} // namespace tudat

#endif  // TUDAT_CREATEATMOSPHEREMODEL_H