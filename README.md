# Online Heart Rate Estimation Framework (WFPV-OVD-PISC)

[![License: MIT](https://img.shields.io/badge/License-MIT-yellow.svg)](https://opensource.org/licenses/MIT)

This repository contains the source code for the paper:  
**"A Robust and Efficient Framework for Online Heart Rate Estimation from PPG Signals During Intensive Physical Exercise"**  
*Submitted to Biomedical Signal Processing and Control*.

## 🚀 Key Features
*   **High Accuracy:** 1.50 BPM Average Absolute Error on IEEE SPC Dataset.
*   **Ultra-Low Latency:** Fixed-lag Online Viterbi Decoding (6s latency).
*   **Efficiency:** Processes 115 mins of data in **2.46s** (MATLAB).
*   **Real-World Ready:** Includes a **Java implementation** tested on Android Smartwatches.

## 📂 Project Structure
*   `/Matlab_Source`: The core algorithm implementation and evaluation scripts (reproduces Table 2 & 3 in the paper).
*   `/Android_Java`: The ported Java class files suitable for Android WearOS deployment.
*   `/Data`: (Optional) Links to the IEEE SPC dataset or sample data.

## 🛠️ Usage
### MATLAB
1. Run `main.m` to process the dataset.
2. Adjust parameters in `config.m` if needed.

### Java (Android)
The `HeartRateEstimator.java` class is a self-contained module.
```java
// Example usage
float[] ppgWindow = ...; // 8s PPG data
float heartRate = estimator.estimate(ppgWindow);

## 🔗 Citation
If you use this code in your research, please cite our paper:
