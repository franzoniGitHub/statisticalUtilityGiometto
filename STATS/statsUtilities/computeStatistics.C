// computeStatistics.C
//
// Horizontal- and time-averaged statistics (first order, second order, autocorrelations)
// for the STATS case. All paths are relative to the directory this program is run from
// (../stats, ../output, ../z.dat, ../Umean), exactly as in the previous version.
//
// Build (C++17; OpenMP is optional but recommended):
//   g++ -std=c++17 -O3 -march=native -fopenmp ComputeStatistics.C -o computeStatistics
// Threads: export OMP_NUM_THREADS=8. The job is largely I/O-bound, so using more threads than
// the disk can feed will not help.
//
// Overview of the algorithm
//   * Every time step is read ONCE. b, U and nut are parsed cell by cell with a buffered reader
//     (no big in-memory fields), and all the sums needed for means, variances, covariances and
//     autocorrelations are accumulated in the same pass.
//   * One-pass variances are made numerically safe by accumulating shifted values
//     (x - s_k), where s_k is the horizontal average of step 0 at height k. Then
//     var = <(x-s)^2> - <x-s>^2 with negligible cancellation.
//   * Time steps are independent: they are distributed over OpenMP threads, each with its own
//     accumulators, merged at the end.
//   * DIS depends on ../Umean (via forSecondOrder.sh), so it is read in a second, scalar-only pass.
//   * Autocorrelations are computed only at the 5 selected heights, averaging over ALL reference
//     points in the periodic direction (not only the first one).
//
// Cell ordering assumption: x fastest, then y, then z (single-block blockMesh, no renumberMesh).

#include <algorithm>
#include <array>
#include <atomic>
#include <cctype>
#include <chrono>
#include <cmath>
#include <cstdio>
#include <cstdlib>
#include <cstring>
#include <fstream>
#include <iomanip>
#include <iostream>
#include <limits>
#include <memory>
#include <optional>
#include <sstream>
#include <stdexcept>
#include <string>
#include <vector>

#include "../SETTINGS.h"

namespace {

// ===========================================================================
// Grid, paths and small utilities
// ===========================================================================

const long cellsPerLevel = long(Nx) * long(Ny);
const long nCells        = cellsPerLevel * long(Nz);
constexpr int nSelected  = 5;
const std::array<int, nSelected> selectedLevels = {k1, k2, k3, k4, k5};

std::string statsFile(const std::string& field, int step)
{
    return "../stats/" + field + std::to_string(step);
}

// Output file with full double precision; throws if it cannot be opened.
std::ofstream openOut(const std::string& path)
{
    std::ofstream f(path);
    if (!f) throw std::runtime_error("Cannot open output file " + path + " (does the directory exist?)");
    f << std::setprecision(std::numeric_limits<double>::max_digits10);
    return f;
}

void runScript(const std::string& command)
{
    const int rc = std::system(command.c_str());
    if (rc != 0) throw std::runtime_error("'" + command + "' failed with return code " + std::to_string(rc));
}

enum class Answer { Yes, No, Quit };

Answer ask(const std::string& question)
{
    std::cout << question;
    char answer = '\0';
    if (!(std::cin >> answer)) return Answer::Quit;
    if (answer == 'Y') return Answer::Yes;
    if (answer == 'n') return Answer::No;
    return Answer::Quit;
}

class Timer {
public:
    double seconds() const
    {
        return std::chrono::duration<double>(std::chrono::steady_clock::now() - start_).count();
    }
private:
    std::chrono::steady_clock::time_point start_ = std::chrono::steady_clock::now();
};

// d/dz on a non-uniform grid: second-order three-point formulas, one-sided at the boundaries.
std::vector<double> ddz(const std::vector<double>& z, const std::vector<double>& f)
{
    const int n = int(z.size());
    std::vector<double> d(n, 0.0);
    {   // bottom
        const double h1 = z[1] - z[0], h2 = z[2] - z[1];
        d[0] = -(2*h1 + h2) / (h1*(h1 + h2)) * f[0] + (h1 + h2) / (h1*h2) * f[1] - h1 / (h2*(h1 + h2)) * f[2];
    }
    for (int i = 1; i < n - 1; ++i) {
        const double h1 = z[i] - z[i-1], h2 = z[i+1] - z[i];
        d[i] = -h2 / (h1*(h1 + h2)) * f[i-1] + (h2 - h1) / (h1*h2) * f[i] + h1 / (h2*(h1 + h2)) * f[i+1];
    }
    {   // top
        const double a = z[n-1] - z[n-2], b = z[n-2] - z[n-3];
        d[n-1] = (2*a + b) / (a*(a + b)) * f[n-1] - (a + b) / (a*b) * f[n-2] + a / (b*(a + b)) * f[n-3];
    }
    return d;
}

// Plain difference f[hi]-f[lo]: one-sided at the boundaries, central inside.
// Only used for the two diagnostic columns of the stationary-balance files.
double rawDifference(const std::vector<double>& f, int i)
{
    const int n = int(f.size());
    const int lo = (i == 0) ? 0 : i - 1;
    const int hi = (i == n - 1) ? n - 1 : i + 1;
    return f[hi] - f[lo];
}

// ===========================================================================
// Buffered reader for OpenFOAM ASCII volScalarField / volVectorField files
// ===========================================================================
// Parses the header instead of skipping a fixed number of lines, checks that the list has
// the expected number of entries, and reads the values with strtod from a 1 MB buffer.

class FoamFieldReader {
public:
    FoamFieldReader(const std::string& path, long expectedCount, bool vectorField)
        : path_(path), file_(std::fopen(path.c_str(), "rb"), &std::fclose), buf_(kBufSize + 1, '\0')
    {
        if (!file_) throw std::runtime_error("Cannot open " + path);
        parseHeader(expectedCount, vectorField);
    }

    double scalar()
    {
        skipSpace();
        ensure(kMaxToken);
        const char* start = buf_.data() + pos_;
        char* end = nullptr;
        const double v = std::strtod(start, &end);
        if (end == start) fail("expected a number (truncated file, or fewer entries than declared?)");
        pos_ = size_t(end - buf_.data());
        return v;
    }

    void vector(double& x, double& y, double& z)
    {
        expect('(');
        x = scalar(); y = scalar(); z = scalar();
        expect(')');
    }

    // Closing parenthesis of the list: detects files with more entries than declared.
    void finish() { expect(')'); }

private:
    static constexpr size_t kBufSize  = size_t(1) << 20;
    static constexpr size_t kMaxToken = 64;

    void refill()
    {
        const size_t remaining = len_ - pos_;
        std::memmove(buf_.data(), buf_.data() + pos_, remaining);
        len_ = remaining;
        pos_ = 0;
        const size_t got = std::fread(buf_.data() + len_, 1, kBufSize - len_, file_.get());
        if (got == 0) eof_ = true;
        len_ += got;
        buf_[len_] = '\0';
    }

    void ensure(size_t n) { while (!eof_ && len_ - pos_ < n) refill(); }

    void skipSpace()
    {
        for (;;) {
            if (pos_ >= len_) { if (eof_) return; refill(); continue; }
            if (!std::isspace(static_cast<unsigned char>(buf_[pos_]))) return;
            ++pos_;
        }
    }

    void expect(char c)
    {
        skipSpace();
        if (pos_ >= len_ || buf_[pos_] != c) fail(std::string("expected '") + c + "'");
        ++pos_;
    }

    std::string word()
    {
        skipSpace();
        std::string w;
        for (;;) {
            if (pos_ >= len_) { if (eof_) break; refill(); continue; }
            const char ch = buf_[pos_];
            if (std::isspace(static_cast<unsigned char>(ch))) break;
            w += ch;
            ++pos_;
        }
        return w;
    }

    long integer()
    {
        skipSpace();
        ensure(kMaxToken);
        const char* start = buf_.data() + pos_;
        char* end = nullptr;
        const long v = std::strtol(start, &end, 10);
        if (end == start) fail("expected the list size");
        pos_ = size_t(end - buf_.data());
        return v;
    }

    void parseHeader(long expectedCount, bool vectorField)
    {
        for (;;) {
            const std::string w = word();
            if (w.empty()) fail("no 'internalField' entry found");
            if (w == "format") {
                if (word().rfind("binary", 0) == 0) fail("binary format is not supported (use writeFormat ascii)");
            } else if (w == "internalField") {
                break;
            }
        }
        const std::string kind = word();
        if (kind != "nonuniform") fail("internalField is '" + kind + "': only nonuniform lists are supported");
        const std::string type = word();
        const std::string wanted = vectorField ? "List<vector>" : "List<scalar>";
        if (type != wanted) fail("expected " + wanted + ", found '" + type + "'");
        const long count = integer();
        if (count != expectedCount)
            fail("list has " + std::to_string(count) + " entries, expected Nx*Ny*Nz = " + std::to_string(expectedCount));
        expect('(');
    }

    [[noreturn]] void fail(const std::string& message) const
    {
        throw std::runtime_error(path_ + ": " + message);
    }

    std::string path_;
    std::unique_ptr<FILE, int (*)(FILE*)> file_;
    std::vector<char> buf_;
    size_t pos_ = 0, len_ = 0;
    bool eof_ = false;
};

// ===========================================================================
// Accumulators
// ===========================================================================

// Sums kept per height k. B/UX/UY/UZ are shifted values (x - s_k); products use shifted values too.
enum Sum : int {
    S_B, S_UX, S_UY, S_UZ,            // first moments (contiguous: S_UX + c)
    S_NUT,                            // nut (not shifted)
    S_BB, S_UXUX, S_UYUY, S_UZUZ,     // second moments (contiguous: S_UXUX + c)
    S_UXUZ,
    S_UXB, S_UYB, S_UZB,              // (contiguous: S_UXB + c)
    S_DIS,                            // DIS (not shifted)
    N_SUMS
};

struct Accumulators {
    std::vector<std::array<double, N_SUMS>> level;  // [k][sum], value-initialized to zero
    std::vector<double> corrX;                      // [(s*3 + c)*Nx + lag]
    std::vector<double> corrY;                      // [(s*3 + c)*Ny + lag]

    explicit Accumulators(bool withCorrelations) : level(Nz)
    {
        if (withCorrelations) {
            corrX.assign(size_t(nSelected) * 3 * Nx, 0.0);
            corrY.assign(size_t(nSelected) * 3 * Ny, 0.0);
        }
    }

    double*       cx(int s, int c)       { return corrX.data() + (size_t(s)*3 + c) * Nx; }
    double*       cy(int s, int c)       { return corrY.data() + (size_t(s)*3 + c) * Ny; }
    const double* cx(int s, int c) const { return corrX.data() + (size_t(s)*3 + c) * Nx; }
    const double* cy(int s, int c) const { return corrY.data() + (size_t(s)*3 + c) * Ny; }

    void merge(const Accumulators& other)
    {
        for (int k = 0; k < Nz; ++k)
            for (int m = 0; m < N_SUMS; ++m) level[k][m] += other.level[k][m];
        for (size_t i = 0; i < other.corrX.size(); ++i) corrX[i] += other.corrX[i];
        for (size_t i = 0; i < other.corrY.size(); ++i) corrY[i] += other.corrY[i];
    }
};

struct Needs {
    bool b = false, nut = false, u = false, correlations = false;
};

struct Shifts {
    std::vector<double> b;
    std::array<std::vector<double>, 3> u;
    Shifts() : b(Nz, 0.0) { for (auto& v : u) v.assign(Nz, 0.0); }
};

// ===========================================================================
// Per-time-step work
// ===========================================================================

// Periodic two-point correlations of one horizontal plane a(i,j) (i fastest), summed over all
// reference points: cx[r] += sum_ij a(i,j) a(i+r,j),  cy[r] += sum_ij a(i,j) a(i,j+r).
void accumulateCorrelations(const double* a, double* cx, double* cy, std::vector<double>& row)
{
    for (int j = 0; j < Ny; ++j) {
        const double* line = a + size_t(j) * Nx;
        std::copy(line, line + Nx, row.begin());        // row = line repeated twice,
        std::copy(line, line + Nx, row.begin() + Nx);   // so row[i + r] wraps around
        for (int r = 0; r < Nx; ++r) {
            const double* shifted = row.data() + r;
            double s = 0.0;
            for (int i = 0; i < Nx; ++i) s += line[i] * shifted[i];
            cx[r] += s;
        }
    }
    for (int r = 0; r < Ny; ++r) {
        double s = 0.0;
        for (int j = 0; j < Ny; ++j) {
            const double* p = a + size_t(j) * Nx;
            const double* q = a + size_t((j + r) % Ny) * Nx;
            for (int i = 0; i < Nx; ++i) s += p[i] * q[i];
        }
        cy[r] += s;
    }
}

// Reads b, nut and U of one time step (as requested by 'need') and adds their sums to 'acc'.
void accumulateStep(int step, const Needs& need, const Shifts& shift, Accumulators& acc)
{
    std::optional<FoamFieldReader> bReader, nutReader, uReader;
    if (need.b)   bReader.emplace(statsFile("b", step), nCells, false);
    if (need.nut) nutReader.emplace(statsFile("nut", step), nCells, false);
    if (need.u)   uReader.emplace(statsFile("U", step), nCells, true);

    std::vector<double> plane, row;
    if (need.correlations) {
        plane.resize(3 * size_t(cellsPerLevel));
        row.resize(2 * size_t(Nx));
    }

    for (int k = 0; k < Nz; ++k) {
        const bool storePlane = need.correlations &&
            std::find(selectedLevels.begin(), selectedLevels.end(), k) != selectedLevels.end();
        const double b0 = shift.b[k];
        const double u0 = shift.u[0][k], v0 = shift.u[1][k], w0 = shift.u[2][k];

        std::array<double, N_SUMS> s{};  // sums over this plane (better rounding than adding cell by cell)
        for (long p = 0; p < cellsPerLevel; ++p) {
            double db = 0.0;
            if (bReader) {
                db = bReader->scalar() - b0;
                s[S_B]  += db;
                s[S_BB] += db * db;
            }
            if (nutReader) s[S_NUT] += nutReader->scalar();
            if (uReader) {
                double ux, uy, uz;
                uReader->vector(ux, uy, uz);
                ux -= u0; uy -= v0; uz -= w0;
                s[S_UX] += ux;         s[S_UY] += uy;         s[S_UZ] += uz;
                s[S_UXUX] += ux * ux;  s[S_UYUY] += uy * uy;  s[S_UZUZ] += uz * uz;
                s[S_UXUZ] += ux * uz;
                s[S_UXB] += ux * db;   s[S_UYB] += uy * db;   s[S_UZB] += uz * db;
                if (storePlane) {
                    plane[p] = ux;
                    plane[cellsPerLevel + p] = uy;
                    plane[2 * cellsPerLevel + p] = uz;
                }
            }
        }
        for (int m = 0; m < N_SUMS; ++m) acc.level[k][m] += s[m];

        if (storePlane)
            for (int sel = 0; sel < nSelected; ++sel)
                if (selectedLevels[sel] == k)
                    for (int c = 0; c < 3; ++c)
                        accumulateCorrelations(plane.data() + c * cellsPerLevel, acc.cx(sel, c), acc.cy(sel, c), row);
    }

    if (bReader)   bReader->finish();
    if (nutReader) nutReader->finish();
    if (uReader)   uReader->finish();
}

void accumulateDissipationStep(int step, Accumulators& acc)
{
    FoamFieldReader reader(statsFile("DIS", step), nCells, false);
    for (int k = 0; k < Nz; ++k) {
        double s = 0.0;
        for (long p = 0; p < cellsPerLevel; ++p) s += reader.scalar();
        acc.level[k][S_DIS] += s;
    }
    reader.finish();
}

// Runs stepFn(step, localAccumulators) for all steps, in parallel over time steps.
template <class StepFn>
Accumulators runOverSteps(int nsteps, bool withCorrelations, const std::string& label, StepFn stepFn)
{
    Accumulators total(withCorrelations);
    std::atomic<int> done{0};
    std::atomic<bool> failed{false};
    std::string errorMessage;
    Timer timer;

    #pragma omp parallel
    {
        Accumulators local(withCorrelations);

        #pragma omp for schedule(dynamic, 1)
        for (int step = 0; step < nsteps; ++step) {
            if (failed) continue;  // exceptions cannot leave an OpenMP region: skip the remaining work
            try {
                stepFn(step, local);
            } catch (const std::exception& e) {
                #pragma omp critical(computeStatistics_error)
                {
                    if (!failed) { errorMessage = e.what(); failed = true; }
                }
            }
            const int d = ++done;
            if (d % 20 == 0 || d == nsteps) {
                #pragma omp critical(computeStatistics_print)
                std::cout << "   " << label << ": " << d << "/" << nsteps << " steps done ("
                          << std::fixed << std::setprecision(1) << timer.seconds() << " s)"
                          << std::defaultfloat << std::endl;
            }
        }

        #pragma omp critical(computeStatistics_merge)
        total.merge(local);
    }

    if (failed) throw std::runtime_error(errorMessage);
    return total;
}

// Horizontal averages of step 0, used as shifts for the numerically safe one-pass moments.
Shifts computeShifts(const Needs& need)
{
    Needs firstStep = need;
    firstStep.nut = false;
    firstStep.correlations = false;
    Accumulators acc(false);
    accumulateStep(0, firstStep, Shifts{}, acc);

    Shifts shift;
    for (int k = 0; k < Nz; ++k) {
        shift.b[k] = acc.level[k][S_B] / double(cellsPerLevel);
        for (int c = 0; c < 3; ++c) shift.u[c][k] = acc.level[k][S_UX + c] / double(cellsPerLevel);
    }
    return shift;
}

// ===========================================================================
// Results
// ===========================================================================

struct Profiles {
    std::vector<double> z, bMean, bRms, covUxUz, tke, nutNu, dUxdz, dis, dbdz;
    std::array<std::vector<double>, 3> uMean, uRms, covUb;

    Profiles()
    {
        for (auto* v : {&z, &bMean, &bRms, &covUxUz, &tke, &nutNu, &dUxdz, &dis, &dbdz}) v->assign(Nz, 0.0);
        for (int c = 0; c < 3; ++c) { uMean[c].assign(Nz, 0.0); uRms[c].assign(Nz, 0.0); covUb[c].assign(Nz, 0.0); }
    }
};

// Everything is normalized as in the previous version: b/b_norm, U/U_norm, z/Z_norm.
void finalizeMoments(const Accumulators& acc, const Shifts& shift, int nsteps, Profiles& p)
{
    const double n = double(cellsPerLevel) * double(nsteps);
    for (int k = 0; k < Nz; ++k) {
        const auto& S = acc.level[k];
        const double mb = S[S_B] / n;
        const std::array<double, 3> mu = {S[S_UX] / n, S[S_UY] / n, S[S_UZ] / n};

        p.bMean[k] = (shift.b[k] + mb) / b_norm;
        p.bRms[k]  = std::sqrt(std::max(0.0, S[S_BB] / n - mb * mb)) / b_norm;
        for (int c = 0; c < 3; ++c) {
            p.uMean[c][k] = (shift.u[c][k] + mu[c]) / U_norm;
            p.uRms[c][k]  = std::sqrt(std::max(0.0, S[S_UXUX + c] / n - mu[c] * mu[c])) / U_norm;
            p.covUb[c][k] = (S[S_UXB + c] / n - mu[c] * mb) / (U_norm * b_norm);
        }
        p.covUxUz[k] = (S[S_UXUZ] / n - mu[0] * mu[2]) / (U_norm * U_norm);
        p.tke[k]     = 0.5 * (p.uRms[0][k] * p.uRms[0][k] + p.uRms[1][k] * p.uRms[1][k] + p.uRms[2][k] * p.uRms[2][k]);
        p.nutNu[k]   = S[S_NUT] / (n * nu);
    }
    p.dUxdz = ddz(p.z, p.uMean[0]);
    p.dbdz  = ddz(p.z, p.bMean);
}

void finalizeDissipation(const Accumulators& acc, int nsteps, Profiles& p)
{
    const double n = double(cellsPerLevel) * double(nsteps);
    const double scale = Z_norm / std::pow(U_norm, 3.0);
    for (int k = 0; k < Nz; ++k) p.dis[k] = acc.level[k][S_DIS] / n * scale;
}

// ===========================================================================
// Input
// ===========================================================================

void validateSettings()
{
    if (Nx < 1 || Ny < 1) throw std::runtime_error("Nx and Ny must be positive");
    if (Nz < 3) throw std::runtime_error("Nz must be at least 3 (needed by the vertical derivatives)");
    for (int k : selectedLevels)
        if (k < 0 || k >= Nz) throw std::runtime_error("Selected height index " + std::to_string(k) + " is outside [0, Nz)");
}

int readNSteps(const std::string& path)
{
    std::ifstream f(path);
    int n = 0;
    if (!f || !(f >> n) || n < 1) throw std::runtime_error("Unable to read a valid number of steps from " + path);
    return n;
}

// z.dat holds the Nz+1 face heights; the profiles use the cell centres.
void readZ(Profiles& p)
{
    std::ifstream f("../z.dat");
    if (!f) throw std::runtime_error("Unable to open z.dat in the STATS directory");
    double previous = 0.0, current = 0.0;
    if (!(f >> previous)) throw std::runtime_error("z.dat is empty");
    for (int k = 0; k < Nz; ++k) {
        if (!(f >> current)) throw std::runtime_error("z.dat has fewer than Nz+1 values");
        p.z[k] = (current + previous) / (2.0 * Z_norm);
        previous = current;
    }
}

// Used when the second-order statistics run without the first-order option (column 14 of ALLmatrix).
void readNutNu(Profiles& p)
{
    std::ifstream f("../output/firstOrder/nut_nu.dat");
    double z = 0.0;
    for (int k = 0; k < Nz; ++k) {
        if (!(f >> z >> p.nutNu[k])) {
            std::cout << "Warning: could not read ../output/firstOrder/nut_nu.dat, the nut/nu column of ALLmatrix.dat will be zero.\n";
            std::fill(p.nutNu.begin(), p.nutNu.end(), 0.0);
            return;
        }
    }
}

// ===========================================================================
// Output
// ===========================================================================

void writeFirstOrder(const Profiles& p)
{
    auto bmean = openOut("../output/firstOrder/bmean.dat");
    auto umean = openOut("../output/firstOrder/Umean.dat");
    auto nutnu = openOut("../output/firstOrder/nut_nu.dat");
    for (int k = 0; k < Nz; ++k) {
        bmean << p.z[k] << "  " << p.bMean[k] << "\n";
        umean << p.z[k] << "  " << p.uMean[0][k] << "  " << p.uMean[1][k] << "  " << p.uMean[2][k] << "\n";
        nutnu << p.z[k] << "  " << p.nutNu[k] << "\n";
    }
}

// Mean velocity as an OpenFOAM volVectorField, for the utility that computes the dissipation.
void writeOpenFoamUmean(const Profiles& p)
{
    auto f = openOut("../Umean");
    f << "FoamFile{version     2.0;    format      ascii;    class       volVectorField;    object      Umean;}\n";
    f << "dimensions      [0 1 -1 0 0 0 0]; internalField   nonuniform List<vector>\n" << nCells << "\n";
    f << "(\n";
    for (int k = 0; k < Nz; ++k) {
        // The value is the same for every cell of a plane: format the line once, write it Nx*Ny times.
        std::ostringstream line;
        line << std::setprecision(std::numeric_limits<double>::max_digits10)
             << "(" << p.uMean[0][k] * U_norm << "  " << p.uMean[1][k] * U_norm << "  " << p.uMean[2][k] * U_norm << ")\n";
        const std::string s = line.str();
        for (long c = 0; c < cellsPerLevel; ++c) f.write(s.data(), std::streamsize(s.size()));
    }
    f << ");\n";
    f << "boundaryField{    front    {        type            cyclic;    }\n";
    f << "    back    {        type            cyclic;    }\n";
    f << "    left    {        type            cyclic;    }\n";
    f << "    right    {        type            cyclic;    }\n";
    f << "    floor    {        type            zeroGradient;    }\n";
    f << "    ceiling    {        type            zeroGradient;    }\n";
    f << "}";
    if (!f) throw std::runtime_error("Error while writing ../Umean");
}

// Stationary momentum and buoyancy balances. Columns:
//   z   forcing   balancing flux divergence   raw difference of the flux   raw difference of z
void writeStationaryBalances(const Profiles& p)
{
    std::vector<double> tauXZ(Nz), taubZ(Nz);
    for (int k = 0; k < Nz; ++k) {
        tauXZ[k] = p.dUxdz[k] / std::sqrt(Gr) - p.covUxUz[k];
        taubZ[k] = p.dbdz[k] / (std::sqrt(Gr) * Pr) - p.covUb[2][k];
    }
    const std::vector<double> dtauXZ = ddz(p.z, tauXZ);
    const std::vector<double> dtaubZ = ddz(p.z, taubZ);

    auto momentum = openOut("../output/secondOrder/momentumStationary.dat");
    auto energy   = openOut("../output/secondOrder/energyStationary.dat");
    for (int k = 0; k < Nz; ++k) {
        const double dz = rawDifference(p.z, k);
        momentum << p.z[k] << "   " << p.bMean[k] * std::sin(alpha) << "   " << -dtauXZ[k] << "   "
                 << rawDifference(tauXZ, k) << "   " << dz << "\n";
        energy   << p.z[k] << "   " << p.uMean[0][k] * std::sin(alpha) << "   " << dtaubZ[k] << "   "
                 << rawDifference(taubZ, k) << "   " << dz << "\n";
    }
}

void writeSecondOrder(const Profiles& p)
{
    {
        auto all = openOut("../output/secondOrder/ALLmatrix.dat");
        all << "z\tmeanb\tmeanUx\tmeanUy\tmeanUz\trmsb\trmsUx\trmsUy\trmsUz\tcovarUxUz\tcovarUxb\tcovarUyb\tcovarUzb\tTKE\tnut/nu\tdUx/dz\tDIS\tdb/dz\n";
        const std::array<const std::vector<double>*, 18> columns = {
            &p.z, &p.bMean, &p.uMean[0], &p.uMean[1], &p.uMean[2],
            &p.bRms, &p.uRms[0], &p.uRms[1], &p.uRms[2],
            &p.covUxUz, &p.covUb[0], &p.covUb[1], &p.covUb[2],
            &p.tke, &p.nutNu, &p.dUxdz, &p.dis, &p.dbdz};
        for (int k = 0; k < Nz; ++k) {
            for (const auto* col : columns) all << (*col)[k] << "\t";
            all << "\n";
        }
    }

    auto brms       = openOut("../output/secondOrder/brms.dat");
    auto urms       = openOut("../output/secondOrder/Urms.dat");
    auto uxuzCovar  = openOut("../output/secondOrder/UxUzcovar.dat");
    auto ubCovar    = openOut("../output/secondOrder/Ubcovar.dat");
    auto tke        = openOut("../output/secondOrder/TKE.dat");
    auto tkeBalance = openOut("../output/secondOrder/TKEbalance.dat");
    for (int k = 0; k < Nz; ++k) {
        brms      << p.z[k] << "  " << p.bRms[k] << "\n";
        urms      << p.z[k] << "  " << p.uRms[0][k] << "  " << p.uRms[1][k] << "  " << p.uRms[2][k] << "\n";
        uxuzCovar << p.z[k] << "  " << p.covUxUz[k] << "\n";
        ubCovar   << p.z[k] << "  " << p.covUb[0][k] << "  " << p.covUb[1][k] << "  " << p.covUb[2][k] << "\n";
        tke       << p.z[k] << "  " << p.tke[k] << "\n";
        // z   shear production   along-slope buoyancy flux   slope-normal buoyancy flux   dissipation
        tkeBalance << p.z[k] << "  " << -p.covUxUz[k] * p.dUxdz[k] << "  " << std::sin(alpha) * p.covUb[0][k]
                   << "  " << std::cos(alpha) * p.covUb[2][k] << "  " << p.dis[k] << "\n";
    }

    writeStationaryBalances(p);
}

// Two-point correlation coefficient rho(r) = C(r)/C(0), with C(r) = <u'(x) u'(x+r)>.
// First column: separation r / L_domain = j/N (j = 0 gives rho = 1). Rows beyond N/2 mirror the
// first half because the domain is periodic.
void writeAutocorrelation(const Accumulators& acc, int nsteps)
{
    const double n = double(cellsPerLevel) * double(nsteps);
    for (int s = 0; s < nSelected; ++s) {
        const int k = selectedLevels[s];
        std::array<std::vector<double>, 3> rhoX, rhoY;
        for (int c = 0; c < 3; ++c) {
            const double m = acc.level[k][S_UX + c] / n;  // mean of the shifted velocity
            const double* cx = acc.cx(s, c);
            const double* cy = acc.cy(s, c);
            const double varX = cx[0] / n - m * m;
            const double varY = cy[0] / n - m * m;
            rhoX[c].resize(Nx);
            rhoY[c].resize(Ny);
            for (int r = 0; r < Nx; ++r) rhoX[c][r] = (cx[r] / n - m * m) / varX;
            for (int r = 0; r < Ny; ++r) rhoY[c][r] = (cy[r] / n - m * m) / varY;
        }
        const std::string id = std::to_string(s + 1);
        auto fy = openOut("../output/autocorrelation/alongY/selection" + id + "Y.dat");
        for (int j = 0; j < Ny; ++j)
            fy << double(j) / Ny << "  " << rhoY[0][j] << "  " << rhoY[1][j] << "  " << rhoY[2][j] << "\n";
        auto fx = openOut("../output/autocorrelation/alongX/selection" + id + "X.dat");
        for (int i = 0; i < Nx; ++i)
            fx << double(i) / Nx << "  " << rhoX[0][i] << "  " << rhoX[1][i] << "  " << rhoX[2][i] << "\n";
    }
}

}  // namespace

// ===========================================================================
// Main
// ===========================================================================

int main()
{
    try {
        std::cout << "The program \"computeStatistics.C\" is running\n";

        if (clear_bool) {
            switch (ask("You chose to clear the STATS case and the environment (also the output directories will be deleted): please, confirm (Y/n). Give any other answer to quit:\n")) {
                case Answer::Yes:  runScript("bash clearAll.sh"); break;
                case Answer::No:   break;
                case Answer::Quit: std::cout << "Closing the program.\n"; return 0;
            }
        }
        if (environment_bool) {
            switch (ask("You chose to create the environment: please, confirm (Y/n). Give any other answer to quit:\n")) {
                case Answer::Yes:  runScript("bash environment.sh"); break;
                case Answer::No:   break;
                case Answer::Quit: std::cout << "Closing the program.\n"; return 0;
            }
        }

        if (!firstOrder_bool && !secondOrder_bool && !autocorrelation_bool) {
            std::cout << "No statistics requested in SETTINGS.h. Bye!\n";
            return 0;
        }

        validateSettings();
        const int nsteps = readNSteps("../stats/nSteps.dat");
        std::cout << "The number of time steps is " << nsteps << "\n";

        Profiles profiles;
        std::cout << "Reading the z values from z.dat\n";
        readZ(profiles);

        Needs need;
        need.u            = true;
        need.b            = firstOrder_bool || secondOrder_bool;
        need.nut          = firstOrder_bool;
        need.correlations = autocorrelation_bool;

        std::cout << "Computing the horizontal averages of step 0 (reference values for the one-pass moments)\n";
        const Shifts shift = computeShifts(need);

        std::cout << "Reading all time steps (means, second-order moments"
                  << (need.correlations ? ", autocorrelations" : "") << ")\n";
        const Accumulators acc = runOverSteps(nsteps, need.correlations, "main pass",
            [&](int step, Accumulators& local) { accumulateStep(step, need, shift, local); });
        finalizeMoments(acc, shift, nsteps, profiles);

        if (firstOrder_bool) {
            std::cout << "Writing the first order output files.\n";
            writeFirstOrder(profiles);
            writeOpenFoamUmean(profiles);
        } else if (secondOrder_bool) {
            readNutNu(profiles);
        }

        if (secondOrder_bool) {
            std::cout << "Computing the second order statistics\n";
            runScript("bash forSecondOrder.sh");
            if (readNSteps("../stats/nSteps2.dat") != nsteps)
                throw std::runtime_error("The general number of time steps is different from that from forSecondOrder.sh");

            const Accumulators disAcc = runOverSteps(nsteps, false, "dissipation pass",
                [](int step, Accumulators& local) { accumulateDissipationStep(step, local); });
            finalizeDissipation(disAcc, nsteps, profiles);

            std::cout << "Writing the second order output files.\n";
            writeSecondOrder(profiles);
        }

        if (autocorrelation_bool) {
            std::cout << "Writing the autocorrelation output files.\n";
            writeAutocorrelation(acc, nsteps);
        }

        std::cout << "END of this program. Bye!\n";
        return 0;
    } catch (const std::exception& e) {
        std::cerr << "Error: " << e.what() << "\nClosing the program.\n";
        return 1;
    } catch (...) {
        std::cerr << "Unknown exception caught. Closing the program.\n";
        return 1;
    }
}
