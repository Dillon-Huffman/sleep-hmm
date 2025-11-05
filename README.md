# sleep-hmm

A real-time sleep scoring framework using Hidden Markov Models for closed-loop sleep manipulation in mice.

## 🔬 Overview

This repository contains the implementation for real-time sleep state classification using Hidden Markov Models (HMMs). The system analyzes EEG and EMG data to automatically classify sleep states, enabling closed-loop experimental paradigms for sleep research.

## 📁 Repository Structure

```
sleep-hmm/
├── matlab-build-hmm/           # HMM model training and construction
├── matlab-simulate-model/      # Model validation and simulation
└── labview-realtime-hmm/       # Real-time implementation
    ├── Types/                  # Custom LabVIEW data types
    └── Utils/                  # Utility VIs and helper functions
```

### `matlab-build-hmm/`

MATLAB code for building and training HMM models.
#### Main File: BuildRealTimeModel.m
- Train HMM models from EEG/EMG sleep recordings
- Reads LabVIEW generated features from CSV file
- Saves model parameters to CSV for real-time implementation

### `matlab-simulate-model/` (not yet implemented)

~~MATLAB code for model simulation and validation.~~
~~#### Main Script: simulateModel.m~~
~~- Test HMM model performance on historical datasets~~
~~- Validate classifier accuracy against manual expert scoring~~
~~- Analyze sleep architecture and state transitions~~
~~- Generate performance metrics (Cohen's Kappa, confusion matrices)~~
  
### `labview-realtime-hmm/`
LabVIEW implementation for real-time sleep scoring.
#### Main File: Acquisition_BassShaker.vi

- Real-time EEG/EMG signal acquisition and processing
- Live HMM model execution and sleep state classification (Wake/NREM/REM)
- Closed Loop Digital triggering + Analog Waveform Generation

#### Subdirectories
- **Types/**: Custom LabVIEW data types and type definitions
- **Utils/**: Utility VIs and helper functions for signal processing


## 🚀 Getting Started

### Prerequisites

- MATLAB (for model building and validation)
- LabVIEW (for real-time implementation)
- EEG/EMG data acquisition hardware
- NI USB DAQs

### Workflow
1. **Collect baseline recording** → Use `labview-realtime-hmm/Acquisition_BassShaker.vi` with "Modeling" set to "Building" to collect 6-8 hour baseline data set
2. **Build Models** → Use `matlab-build-hmm/buildRealTimeModel.m` to train HMM models on LabVIEW generated features
3. **Deploy in Real-time** → Use `labview-realtime-hmm/Acquisition_BassShaker.vi`, importing MATLAB-generated model from step 2 and setting "modeling" to "implementing"

## 📊 Performance
Based on 24-hour score sets from C57BL/6 mice:
- **Cohen's Kappa >75%** agreement with manual expert scoring
- Long-term stability up to **3 weeks** after model construction
- Robust performance during sleep disruption (tested in selective REM sleep restriction trials)

## 📄 Citation

This software implements methods described in:

```bibtex
@article{huffman2021realtime,
  title={A real-time sleep scoring framework for closed-loop sleep manipulation in mice},
  author={Huffman, Dillon and Ajwad, Asma'a and Yaghouby, Farid and O'Hara, Bruce F and Sunderam, Sridhar},
  journal={Journal of Sleep Research},
  volume={30},
  number={4},
  pages={e13262},
  year={2021},
  doi={10.1111/jsr.13262}
}
```

**Paper**: [https://doi.org/10.1111/jsr.13262](https://doi.org/10.1111/jsr.13262)

## 📝 License

Please refer to the [LICENSE](LICENSE) file for usage terms and conditions.

## 📧 Contact

For questions or collaboration inquiries related to the methods or data, reach out to Dillon Huffman (dillon.huffman@sigsoln.com)
