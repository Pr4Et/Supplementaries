// eer2mrc.cpp
// Transforms eer file to mrc with super-resolution, motion corrected.
// Designed by Shahar Seifer, Elbaum lab, Weizmann Institute of Science, 2026
// Based on https://github.com/fei-company/EerReaderLib/tree/master
//
//Either run the exe file from a folder that contains a file "input.eer" or use the full structure:
// eer2mrc [0] [input_file.eer] [output_file.mrc] [upscaleFactor] [timestep_S] [pixelSize_A]



#include <algorithm>
#include <array>
#include <chrono>
#include <cmath>
#include <cstdint>
#include <filesystem>
#include <fstream>
#include <iostream>
#include <string>
#include <vector>
#include <cstdio>
#include <tiffio.h>

//for motioncorrection
#include <fftw3.h>
#include <cstring>
#include <stdexcept>
#include <sstream>   // <-- required for std::ostringstream

//for multithread
#include <unordered_map>
#include <mutex>
#include <omp.h>


#include "../ElectronCountedFramesDecompressor.h"
namespace fs = std::filesystem;

// ------------------------------- MRC writer (single 2D slice, 8-bit) -------------------------------
#pragma pack(push, 1)
struct MRCHeader {
    int32_t nx, ny, nz, mode;
    int32_t nxstart, nystart, nzstart;
    int32_t mx, my, mz;
    float   xlen, ylen, zlen;
    float   alpha, beta, gamma;
    int32_t mapc, mapr, maps;
    float   dmin, dmax, dmean;
    int32_t ispg, nsymbt;
    int32_t extra[25];
    float   originX, originY, originZ;
    char    map[4];
    uint8_t machst[4];
    float   rms;
    int32_t nlabels;
    char    labels[10][80];
};
#pragma pack(pop)
static_assert(sizeof(MRCHeader) == 1024, "MRC header must be 1024 bytes");

struct ImageStats { float dmin, dmax, dmean, rms; };

template <typename T>
ImageStats ComputeStats(const std::vector<T>& data) {
    if (data.empty()) return { 0,0,0,0 };
    auto [minIt, maxIt] = std::minmax_element(data.begin(), data.end());
    double sum = 0.0; for (auto v : data) sum += (double)v;
    double mean = sum / (double)data.size();
    double var = 0.0; for (auto v : data) { double dv = (double)v - mean; var += dv * dv; }
    double rms = std::sqrt(var / (double)data.size());
    return { (float)*minIt, (float)*maxIt, (float)mean, (float)rms };
}



void SaveMRC(const fs::path& outPath,
    uint32_t width, uint32_t height,
    const std::vector<uint8_t>& image,
    double pixelSizeAngstrom /* <=0 if unknown */)
{
    if (image.size() != (size_t)width * (size_t)height)
        throw std::runtime_error("SaveMRC: image size does not match width*height.");

    MRCHeader hdr{};
    hdr.nx = (int32_t)width;  hdr.ny = (int32_t)height; hdr.nz = 1;
    hdr.mode = 0; // uint8
    hdr.nxstart = hdr.nystart = hdr.nzstart = 0;
    hdr.mx = hdr.nx; hdr.my = hdr.ny; hdr.mz = hdr.nz;

    float pix = (pixelSizeAngstrom > 0.0) ? (float)pixelSizeAngstrom : 1.0f;
    hdr.xlen = pix * (float)hdr.mx;  hdr.ylen = pix * (float)hdr.my;  hdr.zlen = pix * (float)hdr.mz;

    hdr.alpha = hdr.beta = hdr.gamma = 90.0f;
    hdr.mapc = 1; hdr.mapr = 2; hdr.maps = 3;

    const auto st = ComputeStats(image);
    hdr.dmin = st.dmin; hdr.dmax = st.dmax; hdr.dmean = st.dmean;

    hdr.ispg = 0; hdr.nsymbt = 0;
    std::fill(std::begin(hdr.extra), std::end(hdr.extra), 0);
    hdr.originX = hdr.originY = hdr.originZ = 0.0f;

    hdr.map[0] = 'M'; hdr.map[1] = 'A'; hdr.map[2] = 'P'; hdr.map[3] = ' ';
    hdr.machst[0] = 0x44; hdr.machst[1] = 0x41; hdr.machst[2] = 0x00; hdr.machst[3] = 0x00;

    hdr.rms = st.rms;
    hdr.nlabels = 1;
    for (auto& lab : hdr.labels) std::memset(lab, 0, sizeof(lab));
    {
        std::string label = "EER -> MRC (uint8)";
        if (pixelSizeAngstrom > 0.0) {
            char buf[64]; std::snprintf(buf, sizeof(buf), " | pixel=%.6f A", pixelSizeAngstrom);
            label += buf;
        }
        std::strncpy(hdr.labels[0], label.c_str(), sizeof(hdr.labels[0]) - 1);
    }

    std::ofstream out(outPath, std::ios::binary);
    if (!out) throw std::runtime_error("SaveMRC: cannot open output file.");
    out.write(reinterpret_cast<const char*>(&hdr), sizeof(hdr));
    out.write(reinterpret_cast<const char*>(image.data()), (std::streamsize)image.size());
    if (!out) throw std::runtime_error("SaveMRC: failed writing MRC file.");
}


#ifndef TIFFTAG_ACQ_METADATA
#define TIFFTAG_ACQ_METADATA 65001 // Acquisition metadata (XML, ISO-8601 timestamp)
#endif
#ifndef TIFFTAG_FRAME_METADATA
#define TIFFTAG_FRAME_METADATA 65002 // Frame metadata (XML, ISO-8601 timestamp)
#endif

// --------- Helper: safe substring extraction for <item name="...">value</item> ----------
static bool XmlFindItemValue(const std::string& xml, const std::string& name, std::string& outVal) {
    // parser: looks for <item name="name">...</item> (unit attribute may exist)
    // Works for both acquisition (65001) and frame (65002) metadata layouts.
    const std::string key = "<item name=\"" + name + "\"";
    size_t p = xml.find(key);
    if (p == std::string::npos) return false;
    p = xml.find('>', p);
    if (p == std::string::npos) return false;
    size_t q = xml.find("</item>", p + 1);
    if (q == std::string::npos) return false;
    outVal = xml.substr(p + 1, q - (p + 1));
    // Trim whitespace
    while (!outVal.empty() && (outVal.back() == '\n' || outVal.back() == '\r' || outVal.back() == ' ' || outVal.back() == '\t')) outVal.pop_back();
    size_t s = 0; while (s < outVal.size() && (outVal[s] == '\n' || outVal[s] == '\r' || outVal[s] == ' ' || outVal[s] == '\t')) s++;
    if (s > 0) outVal.erase(0, s);
    return true;
}

// --------- Helper: read private TIFF string tag as std::string (handles non-terminated) ----------
static bool TiffReadPrivateString(TIFF* tif, uint16_t tag, std::string& out) {
    // Strategy: get raw field as pointer+count when possible; fall back to char* if libtiff returns that.
    // Many libtiff builds will return char* for ASCII, but spec says "non-zero-terminated". We handle both.
    const TIFFField* fld = TIFFFieldWithTag(tif, tag);
    if (!fld) return false;
    // Try the generic interface that returns count + data:
    uint64_t count = 0;
    void* data = nullptr;
    /*if (TIFFGetField(tif, tag, &data) == 1 && data) {
        // Frequently returns char*; try to detect length by the field info if available.
        const char* c = reinterpret_cast<const char*>(data);
        // Attempt to get count if libtiff knows it; if not, copy up to first '\0'
        // Since spec says not terminated, use string constructor with strlen as a safe fallback.
        // However, many implementations *do* provide a terminating zero. We’ll be robust:
        out.assign(c, std::strlen(c));
        return true;
    }*/
    // Some libtiff builds expose a (count, data) pair for undefined/private types:
    if (TIFFGetField(tif, tag, &count, &data) == 1 && data && count > 0) {
        out.assign(reinterpret_cast<const char*>(data), reinterpret_cast<const char*>(data) + static_cast<size_t>(count));
        return true;
    }
    return false;
}

// --------- Helper: parse ISO-8601 timestamp with timezone to epoch milliseconds ----------
static bool ParseIso8601ToMillis(const std::string& s, int64_t& out_ms) {
    // Accepts formats like: YYYY-MM-DDTHH:MM:SS[.ffffff][Z|(+/-)HH:MM]
    // We implement a minimal parser sufficient for EER timestamps.
    int Y, m, d, H, M, tz_h = 0, tz_m = 0; double S;
    char c1, c2, c3, c4, c5;
    // First, split date/time and timezone
    size_t tpos = s.find('T');
    if (tpos == std::string::npos) return false;
    std::string date = s.substr(0, tpos);
    std::string time = s.substr(tpos + 1);
    // Parse date: YYYY-MM-DD
    if (sscanf(date.c_str(), "%d-%d-%d", &Y, &m, &d) != 3) return false;
    // Detect timezone delimiter
    int tz_sign = 0; size_t zpos = time.find('Z');
    size_t plus = time.find('+'); size_t minus = time.rfind('-');
    size_t tzpos = std::string::npos;
    if (zpos != std::string::npos) tzpos = zpos;
    else if (plus != std::string::npos) { tzpos = plus; tz_sign = +1; }
    else if (minus != std::string::npos && minus > 1) { tzpos = minus; tz_sign = -1; } // avoid the leading '-' for seconds fraction
    std::string tmain = (tzpos == std::string::npos) ? time : time.substr(0, tzpos);
    std::string tzone = (tzpos == std::string::npos) ? "" : time.substr(tzpos);

    // Parse time: HH:MM:SS[.fraction]
    // We’ll tolerate micro/nanosecond fraction; sscanf to seconds then scan fraction.
    if (sscanf(tmain.c_str(), "%d:%d:%lf", &H, &M, &S) != 3) return false;
    // Parse timezone offset
    int tz_offset_minutes = 0;
    if (!tzone.empty() && tzone != "Z") {
        if (sscanf(tzone.c_str(), "%c%d:%d", &c1, &tz_h, &tz_m) == 3) {
            if (c1 == '+') tz_sign = +1; else if (c1 == '-') tz_sign = -1;
            tz_offset_minutes = tz_sign * (tz_h * 60 + tz_m);
        }
        else return false;
    }

    // Convert to time_t (UTC) by interpreting the time as local-of-timezone and then subtracting offset:
    struct tm tm_utc {};
    tm_utc.tm_year = Y - 1900;
    tm_utc.tm_mon = m - 1;
    tm_utc.tm_mday = d;
    tm_utc.tm_hour = H;
    tm_utc.tm_min = M;
    tm_utc.tm_sec = int(S);
    // timegm (non-standard) alternative: use time_t as if UTC
#if defined(_WIN32)
    time_t t = _mkgmtime(&tm_utc);
#else
    time_t t = timegm(&tm_utc);
#endif
    if (t == (time_t)-1) return false;
    // Seconds within the second (fraction)
    double sec_frac = S - double(int(S));
    // Apply timezone offset (minutes)
    int64_t ms = int64_t(t) * 1000 + int64_t(std::llround(sec_frac * 1000.0));
    ms -= int64_t(tz_offset_minutes) * 60 * 1000; // convert to true UTC epoch ms
    out_ms = ms;
    return true;
}

// --------- Extract acquisition/frames timestamps ---------
struct EerTiming {
    std::vector<int64_t> frame_ms; // per-frame timestamps (epoch ms)
    int64_t t0_ms = -1;            // acquisition start (from tag 65001)
    double exposure_s = -1.0;
    int numberOfFrames = -1;
};

static EerTiming ExtractEerTimestamps(const char* filepath) {
    EerTiming T;
    TIFF* tif = TIFFOpen(filepath, "r");
    if (!tif) return T;

    // IFD 0: acquisition metadata (65001)
    TIFFSetDirectory(tif, 0);
    std::string acqXml;
    if (TiffReadPrivateString(tif, TIFFTAG_ACQ_METADATA, acqXml)) {
        std::string v;
        if (XmlFindItemValue(acqXml, "timestamp", v)) {
            ParseIso8601ToMillis(v, T.t0_ms); // best effort
        }
        if (XmlFindItemValue(acqXml, "exposureTime", v)) {
            try { T.exposure_s = std::stod(v); }
            catch (...) {}
        }
        if (XmlFindItemValue(acqXml, "numberOfFrames", v)) {
            try { T.numberOfFrames = std::stoi(v); }
            catch (...) {}
        }
    }

    // Iterate all frame IFDs (skip final-image IFD if present; those have Compression=1)
    uint16_t comp = 0;
    do {
        TIFFGetField(tif, TIFFTAG_COMPRESSION, &comp);
        if (comp == 65000 || comp == 65001 || comp == 65002) {
            std::string frameXml;
            if (TiffReadPrivateString(tif, TIFFTAG_FRAME_METADATA, frameXml)) {
                std::string v;
                if (XmlFindItemValue(frameXml, "timestamp", v)) {
                    int64_t ms = -1;
                    if (ParseIso8601ToMillis(v, ms)) T.frame_ms.push_back(ms);
                }
            }
        }
    } while (TIFFReadDirectory(tif)); // advance to next IFD

    TIFFClose(tif);

    // Fallback if no per-frame timestamps were present:
    if (T.frame_ms.empty() && T.t0_ms >= 0 && T.exposure_s > 0 && T.numberOfFrames > 0) {
        double dt_s = T.exposure_s / double(T.numberOfFrames);
        int64_t dt_ms = int64_t(std::llround(dt_s * 1000.0));
        T.frame_ms.resize(size_t(T.numberOfFrames));
        for (int i = 0; i < T.numberOfFrames; ++i) T.frame_ms[size_t(i)] = T.t0_ms + int64_t(i) * dt_ms;
    }
    return T;
}

// --------- Write an MRC stack (float32, mode=2) ----------
void SaveMRCStackFloat(const fs::path& outPath,
    uint32_t width, uint32_t height, uint32_t nz,
    const std::vector<float>& volume,
    double pixelSizeAngstrom,
    const std::string& label0)
{
    if (volume.size() != size_t(width) * size_t(height) * size_t(nz))
        throw std::runtime_error("SaveMRCStackFloat: volume size mismatch.");

    MRCHeader hdr{};
    hdr.nx = (int32_t)width;
    hdr.ny = (int32_t)height;
    hdr.nz = (int32_t)nz;
    hdr.mode = 2; // float32

    hdr.nxstart = hdr.nystart = hdr.nzstart = 0;
    hdr.mx = hdr.nx; hdr.my = hdr.ny; hdr.mz = hdr.nz;
    float pix = (pixelSizeAngstrom > 0.0) ? (float)pixelSizeAngstrom : 1.0f;
    hdr.xlen = pix * (float)hdr.mx;
    hdr.ylen = pix * (float)hdr.my;
    hdr.zlen = pix * (float)hdr.mz;

    hdr.alpha = hdr.beta = hdr.gamma = 90.0f;
    hdr.mapc = 1; hdr.mapr = 2; hdr.maps = 3;

    // Stats across whole volume
    double vmin = +1e300, vmax = -1e300, vsum = 0.0, vsq = 0.0;
    for (float v : volume) { vmin = std::min(vmin, (double)v); vmax = std::max(vmax, (double)v); vsum += v; vsq += (double)v * (double)v; }
    double n = (double)volume.size();
    double mean = (n > 0 ? vsum / n : 0.0);
    double rms = (n > 0 ? std::sqrt(std::max(0.0, vsq / n - mean * mean)) : 0.0);
    hdr.dmin = (float)vmin; hdr.dmax = (float)vmax; hdr.dmean = (float)mean;

    hdr.ispg = 0; hdr.nsymbt = 0;
    std::fill(std::begin(hdr.extra), std::end(hdr.extra), 0);
    hdr.originX = hdr.originY = hdr.originZ = 0.0f;

    hdr.map[0] = 'M'; hdr.map[1] = 'A'; hdr.map[2] = 'P'; hdr.map[3] = ' ';
    hdr.machst[0] = 0x44; hdr.machst[1] = 0x41; hdr.machst[2] = 0x00; hdr.machst[3] = 0x00;
    hdr.rms = (float)rms;

    hdr.nlabels = 1;
    for (auto& lab : hdr.labels) std::memset(lab, 0, sizeof(lab));
    std::string label = label0.empty() ? "EER -> MRC stack (float32)" : label0;
    std::strncpy(hdr.labels[0], label.c_str(), sizeof(hdr.labels[0]) - 1);

    std::ofstream out(outPath, std::ios::binary);
    if (!out) throw std::runtime_error("SaveMRCStackFloat: cannot open output file.");
    out.write(reinterpret_cast<const char*>(&hdr), sizeof(hdr));
    out.write(reinterpret_cast<const char*>(volume.data()), (std::streamsize)(volume.size() * sizeof(float)));
    if (!out) throw std::runtime_error("SaveMRCStackFloat: failed writing MRC stack.");
}

// --------- Build time-binned stack (no motion correction yet) ----------
struct TimeBinningResult {
    uint32_t w = 0, h = 0, nz = 0;
    std::vector<float> volume;           // size = w*h*nz
    std::vector<int64_t> binStart_ms;    // per-slice absolute start time
    double step_s = 1.0;
};




static int InitOpenMPAndReport(int preferred_threads = 0)
{
#ifdef _OPENMP
    omp_set_dynamic(0); // do not let the runtime change thread count
    if (preferred_threads > 0)
        omp_set_num_threads(preferred_threads);

    int maxT = omp_get_max_threads();
    int procs = omp_get_num_procs();
    std::cout << "[OMP] _OPENMP=" << _OPENMP
        << "  procs=" << procs
        << "  max_threads=" << maxT << std::endl;
    return maxT;
#else
    std::cout << "[OMP] DISABLED at compile time (_OPENMP not defined)" << std::endl;
    return 1;
#endif
}


// Build time-binned stack using "one thread per bin"
//   a) Preload all compressed frames once (serial)
//   b) Build frame->bin mapping (serial)
//   c) Parallel over bins; each thread decodes only frames that belong to its bin
static TimeBinningResult BuildTimeBinnedStack_ByBins(
    const char* eerPath,
    double step_s,
    int upscaleFactor,
    double /*pixelSize_A*/)   // kept for signature compatibility
{
    if (step_s <= 0.0)
        throw std::invalid_argument("BuildTimeBinnedStack_ByBins: step_s must be > 0.");
    if (upscaleFactor <= 0)
        throw std::invalid_argument("BuildTimeBinnedStack_ByBins: upscaleFactor must be > 0.");

    TimeBinningResult R;
    R.step_s = step_s;

    // ---- 0) Timestamps (serial) ----
    EerTiming T = ExtractEerTimestamps(eerPath);
    if (T.frame_ms.empty())
        throw std::runtime_error("No frame timestamps available; cannot time-bin.");

    // ---- 1) Preload all frames once (serial) ----
    // USE THE ONE-ARG FILENAME CTOR HERE (preload/index phase)
    ElectronCountedFramesDecompressor loader(eerPath);

    unsigned width = 0, height = 0, nFrames = 0;
    loader.getSize(width, height, nFrames);

    // Truncate to timestamps length if mismatch (conservative)
    if (nFrames != T.frame_ms.size()) {
        const uint32_t n2 = std::min<uint32_t>(nFrames, (uint32_t)T.frame_ms.size());
        std::cerr << "Warning: decoder frames=" << nFrames
            << " timestamps=" << T.frame_ms.size()
            << " -> truncating to " << n2 << "\n";
        nFrames = n2;
    }

    // Upscaled dimensions
    const uint32_t w = width * (uint32_t)upscaleFactor;
    const uint32_t h = height * (uint32_t)upscaleFactor;

    // Capture frame settings (usually constant across the stack)
    const EerFrameSettings settings = loader.frameSettings();

    // Copy compressed bytes for each frame for later parallel decode
    std::vector<std::vector<unsigned char>> frames(nFrames);
    for (uint32_t i = 0; i < nFrames; ++i) {
        frames[i] = loader.copyFrameData(i);
        if (frames[i].empty())
            throw std::runtime_error("Preload error: empty frame after copyFrameData().");
    }

    // ---- 2) Build bins (serial) ----
    const int64_t t0 = *std::min_element(T.frame_ms.begin(), T.frame_ms.begin() + nFrames);
    const int64_t t1 = *std::max_element(T.frame_ms.begin(), T.frame_ms.begin() + nFrames);
    int64_t step_ms = (int64_t)std::llround(step_s * 1000.0);
    if (step_ms <= 0) step_ms = 1;
    const uint32_t nbins = (uint32_t)((t1 - t0) / step_ms) + 1;

   

    int maxThreads = InitOpenMPAndReport(); // prints info once

    // Frame -> bin map
    std::vector<uint32_t> frameBin(nFrames);
    for (uint32_t fi = 0; fi < nFrames; ++fi) {
        const int64_t dt = T.frame_ms[fi] - t0;
        int bin = (int)(dt / step_ms);
        if (bin < 0) bin = 0;
        if (bin >= (int)nbins) bin = (int)nbins - 1;
        frameBin[fi] = (uint32_t)bin;
    }

    // Bin -> list of frames
    std::vector<std::vector<uint32_t>> binFrames(nbins);
    for (uint32_t fi = 0; fi < nFrames; ++fi) {
        binFrames[frameBin[fi]].push_back(fi);
    }

    std::cout << "==== DEBUG TIMING INFO ====\n";
    std::cout << "nFrames = " << nFrames << "\n";
    std::cout << "t0 (ms) = " << t0 << "\n";
    std::cout << "t1 (ms) = " << t1 << "\n";
    std::cout << "Delta t (ms) = " << (t1 - t0) << "\n";
    std::cout << "step_ms = " << step_ms << "\n";
    std::cout << "nbins   = " << nbins << "\n";
    for (size_t b = 0; b < nbins; b++) {
        std::cout << "  bin " << b << " has " << binFrames[b].size() << " frames\n";
    }
    std::cout << "==========================\n";




    // ---- 3) Allocate output volume ----
    R.w = w; R.h = h; R.nz = nbins;
    R.volume.assign((size_t)w * (size_t)h * (size_t)nbins, 0.0f);
    R.binStart_ms.resize(nbins);
    for (uint32_t b = 0; b < nbins; ++b)
        R.binStart_ms[b] = t0 + (int64_t)b * step_ms;

    // ---- 4) Parallel over bins ----
#pragma omp parallel for schedule(dynamic)
    for (int b = 0; b < (int)nbins; ++b)
    {
        if (binFrames[b].empty())
            continue;

        float* slice = R.volume.data() + (size_t)b * (size_t)w * (size_t)h;

        std::vector<ElectronPos> hits;   // <-- requires the struct definition
        hits.reserve(200000);            // optional, just avoids reallocations

        for (uint32_t idx = 0; idx < binFrames[b].size(); ++idx)
        {
            uint32_t fi = binFrames[b][idx];
            const auto& fbytes = frames[fi];

            ElectronCountedFramesDecompressor dec_local(settings, fbytes);

            unsigned cap = dec_local.nElectronFractionUpperLimit(0, 1);
            hits.resize(cap);                           // allocate storage

            unsigned nElect = dec_local.decompressCoordinateList(hits.data(), 0);

            hits.resize(nElect);                        // shrink to actual valid hits

            for (unsigned e = 0; e < nElect; ++e)
            {
                unsigned x = hits[e].x;
                unsigned y = hits[e].y;
                if (x < w && y < h)
                    slice[(size_t)y * w + x] += 1.0f;
            }
        }
    }
    return R;
}


void PrintAllTiffTags(const char* path) {
    TIFF* tif = TIFFOpen(path, "r");
    if (!tif) { std::fprintf(stderr, "Cannot open %s\n", path); return; }
    do { TIFFPrintDirectory(tif, stderr, TIFFPRINT_NONE); } while (TIFFReadDirectory(tif));
    TIFFClose(tif);
}

// ---------------- EER subpixel info (safe) ----------------
#ifndef TIFFTAG_POSSKIPBITS
#define TIFFTAG_POSSKIPBITS 65007
#endif
#ifndef TIFFTAG_HORZSUBBITS
#define TIFFTAG_HORZSUBBITS 65008
#endif
#ifndef TIFFTAG_VERTSUBBITS
#define TIFFTAG_VERTSUBBITS 65009
#endif

struct EerDecoderSettings { uint16_t posSkipBits = 7, horzSubBits = 2, vertSubBits = 2; };

bool GetEerDecoderSettings(TIFF* tif, EerDecoderSettings& s) {
    uint16_t compression = 0;
    if (TIFFGetField(tif, TIFFTAG_COMPRESSION, &compression) != 1) return false;

    if (compression == 65000) { s.posSkipBits = 8; s.horzSubBits = 2; s.vertSubBits = 2; return true; } // 8/2/2
    if (compression == 65001) { s.posSkipBits = 7; s.horzSubBits = 2; s.vertSubBits = 2; return true; } // 7/2/2
    if (compression == 65002) {
        uint16_t v = 0;
        if (TIFFGetField(tif, TIFFTAG_POSSKIPBITS, &v) == 1) s.posSkipBits = v;
        if (TIFFGetField(tif, TIFFTAG_HORZSUBBITS, &v) == 1) s.horzSubBits = v;
        if (TIFFGetField(tif, TIFFTAG_VERTSUBBITS, &v) == 1) s.vertSubBits = v;
        return true;
    }
    // Unknown compression -> treat as non-EER
    return false;
}

// ---------------- Brief summary (NO metadata read) ----------------
bool PrintBriefEerSummary(const char* filepath,
    int& suggestedUpscale,   // OUTPUT
    double& pxWidth_A,       // OUTPUT (Å)
    double& pxHeight_A)      // OUTPUT (Å)
{
    suggestedUpscale = 1;
    pxWidth_A = pxHeight_A = -1.0;  // “unknown”

    TIFF* tif = TIFFOpen(filepath, "r");
    if (!tif) {
        std::cerr << "ERROR: Cannot open file: " << filepath << "\n";
        return false;
    }

    TIFFSetDirectory(tif, 0);

    uint32_t w = 0, h = 0; uint16_t comp = 0;
    TIFFGetField(tif, TIFFTAG_IMAGEWIDTH, &w);
    TIFFGetField(tif, TIFFTAG_IMAGELENGTH, &h);
    TIFFGetField(tif, TIFFTAG_COMPRESSION, &comp);

    std::cout << "---- EER / TIFF Summary ----\n";
    std::cout << "Image size: " << w << " x " << h << "\n";
    std::cout << "Compression: " << comp << " (65000/65001/65002 => EER)\n";

    // We intentionally skip reading AcquisitionMetadata (65001) to avoid ABI issues.
    std::cout << "Pixel size extraction from metadata: SKIPPED (use CLI [pixelSize_A] if needed)\n";

    // Upscaling suggestion from compression / sub-pixel bits
    if (comp == 65000 || comp == 65001) {
        suggestedUpscale = 4; // 2 subpixel bits per axis in both modes
        std::cout << "Mode " << comp << ": fixed 2 subpixel bits -> suggest 4×\n";
    }
    else if (comp == 65002) {
        uint16_t posSkip = 7, hsub = 2, vsub = 2;
        TIFFGetField(tif, TIFFTAG_POSSKIPBITS, &posSkip);
        TIFFGetField(tif, TIFFTAG_HORZSUBBITS, &hsub);
        TIFFGetField(tif, TIFFTAG_VERTSUBBITS, &vsub);
        int upscaleX = 1 << hsub, upscaleY = 1 << vsub;
        suggestedUpscale = std::min(upscaleX, upscaleY);
        std::cout << "Mode 65002: PosSkipBits=" << posSkip
            << " HorzSubBits=" << hsub
            << " VertSubBits=" << vsub
            << "-> suggest " << suggestedUpscale << "×\n";
    }
    else {
        std::cout << "Not an EER compression; upscale = 1×\n";
    }

    std::cout << "--------------------------------\n";
    TIFFClose(tif);
    return true;
}


// ======================================================================
// Downsample by N (box-mean) from float (full-res) -> double (low-res)
// ======================================================================
// General N×N (box-mean) downsample: float (full-res) -> double (low-res)
static void DownsampleN_toDouble(const float* src, int w, int h,
    double* dst, int w2, int h2, int f)
{
    // w2 == w/f, h2 == h/f
    for (int y2 = 0; y2 < h2; y2++) {
        const int y0 = y2 * f;
        for (int x2 = 0; x2 < w2; x2++) {
            const int x0 = x2 * f;
            double sum = 0.0;
            for (int yy = 0; yy < f; yy++) {
                const float* row = src + (size_t)(y0 + yy) * w + x0;
                for (int xx = 0; xx < f; xx++) sum += (double)row[xx];
            }
            dst[(size_t)y2 * w2 + x2] = sum / (double)(f * f);
        }
    }
}

// ======================================================================
// FFT-based Phase Correlation (double-precision) to estimate shift
// Input: ref, img are double (downsampled). Output shifts in downsampled pixels.
// ======================================================================
static void PhaseCorrelationFFTW_double_TS(const double* ref,
    const double* img,
    int w, int h,
    double& dx, double& dy)
{
    const int N = w * h;
    const int nFreq = w * (h / 2 + 1);

    // Allocate per-call (thread-local) FFTW buffers
    fftw_complex* F_ref = (fftw_complex*)fftw_malloc(sizeof(fftw_complex) * nFreq);
    fftw_complex* F_img = (fftw_complex*)fftw_malloc(sizeof(fftw_complex) * nFreq);
    fftw_complex* F_prod = (fftw_complex*)fftw_malloc(sizeof(fftw_complex) * nFreq);
    double* corr = (double*)fftw_malloc(sizeof(double) * N);
    if (!F_ref || !F_img || !F_prod || !corr) {
        if (F_ref)  fftw_free(F_ref);
        if (F_img)  fftw_free(F_img);
        if (F_prod) fftw_free(F_prod);
        if (corr)   fftw_free(corr);
        throw std::runtime_error("FFTW alloc failed");
    }

    // Plans (guard plan creation just in case this FFTW build isn't thread-safe for planning)
    fftw_plan pRef, pImg, pCorr;
#pragma omp critical(fftw_plan_creation)
    {
        pRef = fftw_plan_dft_r2c_2d(h, w, const_cast<double*>(ref), F_ref, FFTW_ESTIMATE);
        pImg = fftw_plan_dft_r2c_2d(h, w, const_cast<double*>(img), F_img, FFTW_ESTIMATE);
        pCorr = fftw_plan_dft_c2r_2d(h, w, F_prod, corr, FFTW_ESTIMATE);
    }
    if (!pRef || !pImg || !pCorr)
    {
        if (pRef)  fftw_destroy_plan(pRef);
        if (pImg)  fftw_destroy_plan(pImg);
        if (pCorr) fftw_destroy_plan(pCorr);
        fftw_free(F_ref); fftw_free(F_img); fftw_free(F_prod); fftw_free(corr);
        throw std::runtime_error("FFTW plan creation failed");
    }

    // Execute
    fftw_execute(pRef);
    fftw_execute(pImg);

    // Cross-power spectrum (normalize)
    for (int i = 0; i < nFreq; i++) {
        const double a = F_ref[i][0], b = F_ref[i][1];
        const double c = F_img[i][0], d = F_img[i][1];
        const double re = a * c + b * d;
        const double im = b * c - a * d;
        const double mag = std::sqrt(re * re + im * im) + 1e-20;
        F_prod[i][0] = re / mag;
        F_prod[i][1] = im / mag;
    }

    fftw_execute(pCorr);

    // Peak locate + 1D parabolic subpixel
    int peakX = 0, peakY = 0; double peakVal = -1e300;
    for (int y = 0; y < h; y++) {
        const double* row = corr + (size_t)y * w;
        for (int x = 0; x < w; x++) {
            const double v = row[x];
            if (v > peakVal) { peakVal = v; peakX = x; peakY = y; }
        }
    }
    auto parabolic = [](double L, double C, double R)->double {
        const double den = (L - 2.0 * C + R);
        if (fabs(den) < 1e-12) return 0.0;
        return 0.5 * (L - R) / den;
        };
    double subX = 0.0, subY = 0.0;
    if (peakX > 0 && peakX < w - 1) {
        const double L = corr[(size_t)peakY * w + (peakX - 1)];
        const double C = corr[(size_t)peakY * w + (peakX)];
        const double R = corr[(size_t)peakY * w + (peakX + 1)];
        subX = parabolic(L, C, R);
    }
    if (peakY > 0 && peakY < h - 1) {
        const double U = corr[(size_t)(peakY - 1) * w + peakX];
        const double C = corr[(size_t)(peakY)*w + peakX];
        const double D = corr[(size_t)(peakY + 1) * w + peakX];
        subY = parabolic(U, C, D);
    }

    double cx = (double)peakX + subX;
    double cy = (double)peakY + subY;
    if (cx > w / 2.0) cx -= w;
    if (cy > h / 2.0) cy -= h;
    dx = -cx; dy = -cy;

    // Cleanup
    fftw_destroy_plan(pRef);
    fftw_destroy_plan(pImg);
    fftw_destroy_plan(pCorr);
    fftw_free(F_ref);
    fftw_free(F_img);
    fftw_free(F_prod);
    fftw_free(corr);
}

// ======================================================================
// Savitzky–Golay smoothing (5-point, quadratic) on motion curves (double)
// ======================================================================
static void SmoothTrajectory(std::vector<double>& v)
{
    const int n = (int)v.size();
    if (n < 5) return;

    std::vector<double> out(n);
    for (int i = 2; i < n - 2; i++) {
        out[i] = (-3.0 * v[i - 2] + 12.0 * v[i - 1] + 17.0 * v[i]
            + 12.0 * v[i + 1] - 3.0 * v[i + 2]) / 35.0;
    }
    out[0] = v[0];
    out[1] = v[1];
    out[n - 2] = v[n - 2];
    out[n - 1] = v[n - 1];
    v.swap(out);
}

// ======================================================================
// Bilinear shift on full-resolution float slices (dx,dy can be fractional)
// ======================================================================
static void BilinearShift(const float* src, int w, int h,
    float* dst, double dx, double dy)
{
    // dst(x,y) <-- src(x - dx, y - dy)
    for (int y = 0; y < h; y++) {
        const double yy = (double)y - dy;
        const int y0 = (int)std::floor(yy);
        const double fy = yy - (double)y0;

        if (y0 < 0 || y0 >= h - 1) {
            std::memset(dst + (size_t)y * w, 0, sizeof(float) * (size_t)w);
            continue;
        }
        const float* row0 = src + (size_t)y0 * w;
        const float* row1 = row0 + w;

        for (int x = 0; x < w; x++) {
            const double xx = (double)x - dx;
            const int x0 = (int)std::floor(xx);
            const double fx = xx - (double)x0;

            if (x0 < 0 || x0 >= w - 1) {
                dst[(size_t)y * w + x] = 0.0f;
                continue;
            }

            const float v00 = row0[x0];
            const float v01 = row0[x0 + 1];
            const float v10 = row1[x0];
            const float v11 = row1[x0 + 1];

            const float v0 = (float)((1.0 - fx) * v00 + fx * v01);
            const float v1 = (float)((1.0 - fx) * v10 + fx * v11);
            dst[(size_t)y * w + x] = (float)((1.0 - fy) * v0 + fy * v1);
        }
    }
}

// ======================================================================
// MotionCor-style correction (CPU, FFTW double, bilinear shift on float)
// slices: contiguous (w*h*nz) float volume
// ======================================================================
// Parallel MotionCor-style correction with configurable downsample (ds = 8 or 16)
// - slices: contiguous [w*h*nz] float volume (in-place)
// - ds: 8 (original) or 16 (faster, coarser); must divide both w and h
// - nthreads: 0 => let OMP decide; >0 => force that many threads
void MotionCorStyleCorrectOMP(float* slices, int w, int h, int nz,
    int ds /*=8*/, int nthreads /*=0*/)
{
    if (nz <= 1) return;
    if (ds <= 0 || (w % ds) || (h % ds))
        throw std::invalid_argument("MotionCorStyleCorrectOMP: ds must divide w and h");

#ifdef _OPENMP
    omp_set_dynamic(0);
    if (nthreads > 0) omp_set_num_threads(nthreads);
#endif

    const int w2 = w / ds;
    const int h2 = h / ds;

    // 1) Build downsampled stack (double) -- parallel over slices
    std::vector<double> dsStack((size_t)w2 * h2 * nz);
#ifdef _OPENMP
#pragma omp parallel for schedule(static)
#endif
    for (int k = 0; k < nz; k++) {
        const float* src = slices + (size_t)k * w * h;
        double* dst = dsStack.data() + (size_t)k * w2 * h2;
        DownsampleN_toDouble(src, w, h, dst, w2, h2, ds);
    }

    // 2) Reference = mean of downsampled slices (parallel reduce)
    std::vector<double> ref((size_t)w2 * h2, 0.0);
#ifdef _OPENMP
#pragma omp parallel
    {
        std::vector<double> local((size_t)w2 * h2, 0.0);
#pragma omp for nowait
        for (int k = 0; k < nz; k++) {
            const double* src = dsStack.data() + (size_t)k * w2 * h2;
            for (size_t i = 0; i < (size_t)w2 * h2; i++) local[i] += src[i];
        }
#pragma omp critical
        {
            for (size_t i = 0; i < (size_t)w2 * h2; i++) ref[i] += local[i];
        }
    }
#else
    for (int k = 0; k < nz; k++) {
        const double* src = dsStack.data() + (size_t)k * w2 * h2;
        for (size_t i = 0; i < (size_t)w2 * h2; i++) ref[i] += src[i];
    }
#endif
    const double inv_nz = 1.0 / (double)nz;
    for (size_t i = 0; i < (size_t)w2 * h2; i++) ref[i] *= inv_nz;

    // 3) Measure shifts (parallel over slices)
    std::vector<double> dx(nz, 0.0), dy(nz, 0.0);
#ifdef _OPENMP
#pragma omp parallel for schedule(static)
#endif
    for (int k = 0; k < nz; k++) {
        const double* img = dsStack.data() + (size_t)k * w2 * h2;
        double sx = 0.0, sy = 0.0;
        PhaseCorrelationFFTW_double_TS(ref.data(), img, w2, h2, sx, sy);
        dx[k] = sx * (double)ds;
        dy[k] = sy * (double)ds;
    }

    // 4) Smooth global trajectory (cheap; keep serial or parallel if nz is large)
    SmoothTrajectory(dx);
    SmoothTrajectory(dy);

    // 5) Apply shifts to full-res slices (parallel over slices)
    std::vector<float> tmp((size_t)w * h);
#ifdef _OPENMP
#pragma omp parallel for schedule(static)
#endif
    for (int k = 0; k < nz; k++) {
        const float* src = slices + (size_t)k * w * h;
        float* dst = slices + (size_t)k * w * h;
        // Use thread-local temp: allocate per-thread to avoid contention
        std::vector<float> localTmp((size_t)w * h);
        BilinearShift(src, w, h, localTmp.data(), dx[k], dy[k]);
        std::memcpy(dst, localTmp.data(), sizeof(float) * (size_t)w * h);
    }
}

// ------------------------------------ Program entry ------------------------------------
int main(int argc, char* argv[]) {
    std::cout << "By default I look for file input.eer in CWD: "
        << std::filesystem::current_path().string() << "\n";
	std::cout << "Usage: " << argv[0] << "[save motion? 1/0] [input.eer] [output.mrc] [upscaleFactor] [timestep_S] [pixelSize_A]\n";

    
    try {
        // Defaults
        fs::path inputEER = "input.eer";
        fs::path outputMRC = "output.mrc";
        int    upscaleFactor = 0;      // 0 => auto-suggest
        double pixelSizeAngstrom = -1.0;   // unknown
        double timeStepSec = 1.0; // <0 => no stack, write single summed image
        int save_motion = 0; // 1 = save motion stacks for validation

        bool doMotionCor = true;
        int  mcDownsample = 16; //  MotionCorStyleCorrect uses fixed 8× internally



        // Usage: Eer2Mrc.exe <input.eer> [output.mrc] [upscaleFactor] [pixelSize_A]
        if (argc >= 2) save_motion = std::stoi(argv[1]);
        if (argc >= 3) inputEER = fs::path(argv[2]);
        if (argc >= 4) outputMRC = fs::path(argv[3]);
        if (argc >= 5) upscaleFactor = std::stoi(argv[4]);
        if (argc >= 6) timeStepSec = std::stod(argv[5]);
        if (argc >= 7) pixelSizeAngstrom = std::stod(argv[6]);

        int suggested = 1; double pxW = -1, pxH = -1;
        PrintBriefEerSummary(inputEER.string().c_str(), suggested, pxW, pxH);

        if (upscaleFactor <= 0)           // auto
            upscaleFactor = suggested;

        if (pixelSizeAngstrom <= 0 && pxW > 0) // we skipped metadata; this will remain -1.0 unless you passed CLI
            pixelSizeAngstrom = pxW;

        if (upscaleFactor <= 0) {
            std::cerr << "Error: upscaleFactor must be a positive integer.\n";
            return EXIT_FAILURE;
        }


        if (timeStepSec > 0.0) {

            std::cout << "Building time-binned stack with Dt = "
                << timeStepSec << " s\n";

            //auto R = BuildTimeBinnedStack(inputEER.string().c_str(),
            //    timeStepSec,
            //    upscaleFactor,
            //    pixelSizeAngstrom);

            auto R = BuildTimeBinnedStack_ByBins(inputEER.string().c_str(),
                timeStepSec,
                upscaleFactor,
                pixelSizeAngstrom);



            // --- Save ORIGINAL (unshifted) stack ---
            fs::path outStack = outputMRC;
            if (outStack.extension() != ".mrc") {
                outStack = outStack.replace_extension(".mrc");
            }

            // Name: original stack
            fs::path stackOriginal = outStack.parent_path() /
                (outStack.stem().string() + "_stack_original.mrc");

            // Label for original stack
            std::ostringstream lbl;
            lbl.setf(std::ios::fixed);
            lbl.precision(3);
            lbl << "EER -> MRC stack (original); Dt=" << timeStepSec << " s; nz=" << R.nz;
            if (pixelSizeAngstrom > 0.0)  lbl << " ; pixel=" << pixelSizeAngstrom << " A";
            if (save_motion) {
                SaveMRCStackFloat(stackOriginal, R.w, R.h, R.nz, R.volume, pixelSizeAngstrom, lbl.str());
                std::cout << "Wrote ORIGINAL MRC stack: "
                    << stackOriginal.string() << " ("
                    << R.w << "x" << R.h << "x" << R.nz
                    << ", float32)\n";
            }

            // ======================================================
            // OPTIONAL: MOTION-CORRECTED STACK
            // ======================================================
            if (doMotionCor) {

                // Save a copy of original volume in case you need both
                std::vector<float> originalVolume = R.volume;

                std::cout << "Running MotionCor-style correction on stack: "
                    << R.w << "x" << R.h << "x" << R.nz << " ...\n";

                const auto t0 = std::chrono::steady_clock::now();

                // In-place correction
                MotionCorStyleCorrectOMP(R.volume.data(), (int)R.w, (int)R.h, (int)R.nz, /*ds=*/mcDownsample, /*nthreads=*/0);


                const auto t1 = std::chrono::steady_clock::now();
                const double sec =
                    std::chrono::duration<double>(t1 - t0).count();

                std::cout << "MotionCor-style correction finished in "
                    << sec << " s\n";

                // Name: corrected stack
                fs::path stackShifted = outStack.parent_path() /
                    (outStack.stem().string() + "_stack_shifted.mrc");

                // Label for shifted
                std::ostringstream lbl2;
                lbl2.setf(std::ios::fixed);
                lbl2.precision(3);
                lbl2 << "EER -> MRC stack (MotionCor-corrected); Dt="
                    << timeStepSec << " s; nz=" << R.nz;
                lbl2 << " ; MotionCor";
                if (pixelSizeAngstrom > 0.0)
                    lbl2 << " ; pixel=" << pixelSizeAngstrom << " A";

                if (save_motion) {
                    SaveMRCStackFloat(stackShifted, R.w, R.h, R.nz, R.volume, pixelSizeAngstrom, lbl2.str());
                    std::cout << "Wrote SHIFTED (motion-corrected) MRC stack: "
                        << stackShifted.string() << " ("
                        << R.w << "x" << R.h << "x" << R.nz
                        << ", float32)\n";
                }
            }


            // == = Compute SUM of all motion - corrected slices == =
            std::vector<float> sumImage((size_t)R.w * R.h, 0.0f);
            for (int k = 0; k < R.nz; k++) {
                const float* slice = R.volume.data() + (size_t)k * R.w * R.h;
                for (size_t p = 0; p < (size_t)R.w * R.h; p++) {
                    sumImage[p] += slice[p];
                }
            }
            std::cout << "Writing MRC: " << outputMRC.string() << std::endl;
            std::vector<float> volume = sumImage; // 2D -> treat as a single z-slice
            SaveMRCStackFloat(outputMRC, R.w, R.h, 1, volume, pixelSizeAngstrom, "float2D");
            std::cout << "Done.\n";

            return EXIT_SUCCESS;
        }
        else
        {
			// Decode EER single frame
            unsigned width = 0, height = 0, nrOfFrames = 0;
            ElectronCountedFramesDecompressor decompressor(inputEER.string());
            decompressor.getSize(width, height, nrOfFrames);
            std::cout << "EER size: " << width << "x" << height
                << " | frames: " << nrOfFrames << '\n';

            const uint32_t w = width * (uint32_t)upscaleFactor;
            const uint32_t h = height * (uint32_t)upscaleFactor;
            std::vector<uint8_t> image((size_t)w * (size_t)h, 0);

            std::cout << "Decoding EER with upscaleFactor=" << upscaleFactor << " ..." << std::endl;
            decompressor.decompressImage(image.data(), upscaleFactor);
            for (unsigned f = 1; f < nrOfFrames; ++f)
                decompressor.decompressImage_AddTo(image.data(), upscaleFactor);
            std::cout << "Electrons counted: " << decompressor.getNElectronsCounted() << std::endl;

            // If your decompressor can provide pixel size in Å, prefer that (pseudo):
            // double pxMeta = -1.0;
            // if (decompressor.tryGetPixelSizeAngstrom(pxMeta)) pixelSizeAngstrom = pxMeta;

            std::cout << "Writing MRC: " << outputMRC.string() << std::endl;
            SaveMRC(outputMRC, w, h, image, pixelSizeAngstrom);
            std::cout << "Done." << std::endl;
            return EXIT_SUCCESS;
        }
    }
    catch (const std::exception& e) {
        std::cerr << "Fatal: " << e.what() << std::endl;
        return EXIT_FAILURE;
    }
}