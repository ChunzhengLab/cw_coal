# cw_coal

## Project Overview

`cw_coal` is a modern C++20 library and command-line tool suite for simulating parton-to-hadron coalescence in high-energy nuclear physics. It implements multiple coalescence algorithms (brute-force, greedy, KD-tree accelerated), event I/O (AMPT/ROOT), particle ID inference, QA and correlation analysis, and Python-based visualization. Batch execution via Condor is also supported.

## Key Features

- **Coalescence Algorithms**:
  - **Brute-Force**: `BruteForceGlobal`, `BruteForceGreedy`, `BruteForceDualGreedy`
  - **KD-Tree Accelerated**: `KDTreeGlobal`, `KDTreeGreedy`, `KDTreeDualGreedy` (powered by [nanoflann])
- **Event Model**:
  - `Event`, `Parton`, `Hadron` classes with random generation, momentum/position handling, and constituent tracking
- **PID Inference & Assignment**:
  - Flavor-based meson/diquark identification and heavy-quarkonium inference
- **I/O Support**:
  - `EventReaderAMPT`, `EventRandomGen`, `EventWriter` for ROOT file formats
- **Analysis Tools**:
  - `analysis` CLI for QA histograms and two-particle correlations
- **Python Scripts**:
  - `scripts/visualize_event.py`, `investigate_baryon-pref_ratio.py`, `draw.py`
- **Batch Scheduling**:
  - Condor scripts under `condor/` with customizable job lists

## Command-Line Interface: `cwcoal`

The primary CLI tool, `cwcoal`, orchestrates event reading/generation, parton shuffling, coalescence, PID assignment, and output writing—all via simple options.

### Synopsis
```bash
cwcoal [OPTIONS]
```

### Options

| Option                         | Description                                                                                       |
|--------------------------------|---------------------------------------------------------------------------------------------------|
| `-h, --help`                   | Show help message and exit                                                                       |
| `-i, --data-input <file>`      | Input ROOT file with events (AMPT/previous output). If omitted, use `--toymode` to generate events |
| `-o, --data-output <file>`     | Output ROOT file path for hadronized events                                                      |
| `-a, --algorithm <string>`     | Coalescence algorithm: `KDTreeGlobal`, `KDTreeGreedy`, `BruteForceGlobal`, `BruteForceGreedy` (default: `KDTreeGlobal`) |
| `-n, --events <N>`             | Number of events to process or generate (default: 10)                                            |
| `-p, --partons <N>`            | Partons per event (default: sample from histogram if unset)                                      |
| `-b, --bn <int>`               | Net baryon number sum (default: 0)                                                               |
| `-s, --savedir <dir>`          | Directory for all output files (default: `.`)                                                    |
| `-r, --baryon-preference <dbl>`| Scale factor for baryon vs. meson preference (default: 1.0)                                       |
| `-F, --shuffle-fraction <dbl>` | Fraction of partons to shuffle across events (0–1)                                               |
| `-T, --toymode`                | Enable toy-mode: generate random events instead of reading input                                  |

### Workflow
1. **Event Source**: read from `--data-input` or generate with `--toymode`.
2. **Shuffle**: optionally mix partons across events via `--shuffle-fraction`.
3. **Coalescence**: choose an algorithm (`--algorithm`) to form hadrons.
4. **PID Assignment**: infer PDG codes and assign to hadrons.
5. **Output**: write hadronized events to `--data-output` and produce QA/CVE histograms in `--savedir`.

### Examples

Generate 100 toy events with 50 partons each using KDTreeGreedy:
```bash
cwcoal --toymode --events 100 --partons 50 \
       --algorithm KDTreeGreedy \
       --data-output toy_hadrons.root \
       --savedir results
```

Process 500 AMPT events with brute-force global coalescence and 20% shuffling:
```bash
cwcoal --data-input ampt_events.root \
       --events 500 \
       --algorithm BruteForceGlobal \
       --shuffle-fraction 0.2 \
       --data-output ampt_hadrons.root \
       --savedir results
```

Afterwards, run QA and CVE analysis:
```bash
analysis --input results/ampt_hadrons.root \
         --output results/qa_cve.root
```

## Repository Structure
```
cw_coal/
├── CMakeLists.txt        # Build configuration
├── include/              # Public headers
│   ├── core/             # Event, Parton, Hadron, PID, CombinerBase
│   ├── io/               # I/O interfaces
│   └── ana/              # Analyzer interfaces (QA, CVE)
├── src/                  # Source implementation
│   ├── core/             # Core classes
│   ├── combiner/         # Algorithm implementations
│   ├── io/               # Readers, writers, random generator
│   └── app/              # CLI tools (`cwcoal`, `analysis`)
├── test/                 # Unit tests (Google Test)
├── scripts/              # Python visualization and analysis
├── refdata/              # Reference data (ZPC lists)
└── condor/               # Condor batch scripts and lists
```

## Installation & Build

### Prerequisites
- CMake ≥ 3.10
- C++20 compiler (e.g. GCC 9+, Clang 10+)
- [nanoflann] for KD-tree accelerated combiners
- [ROOT] for I/O and data structures
- Python 3 (optional, for scripts)

### Build Steps
```bash
git clone <repo_url>
cd cw_coal
mkdir build && cd build
cmake ..
make -j$(nproc)
make install  # optional: install into system or custom prefix
```

## Testing
```bash
cd build
ctest --output-on-failure
```

## License
MIT License

[nanoflann]: https://github.com/jlblancoc/nanoflann
[ROOT]: https://root.cern/

