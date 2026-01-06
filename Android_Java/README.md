# Android/Java Implementation

This directory contains the lightweight Java implementation of the heart rate estimation framework, specifically ported and optimized for **Android WearOS** and embedded Java environments. 

It replicates the functionality of the MATLAB research code (WFPV + Online Viterbi + PISC) but is engineered for **zero-dependency** deployment and minimal memory footprint.

## 📂 File Structure

The implementation is decoupled into three core files for easy integration:

*   **`OnlineHeartRateEstimator.java`**  
    The main entry point. It manages the algorithm state (buffers, Viterbi path history) and executes the core pipeline (Signal Enhancement $\to$ OVD $\to$ PISC).

*   **`DSPUtils.java`**  
    A static utility class containing low-level signal processing primitives, including:
    *   Fast Fourier Transform (FFT)
    *   Butterworth Band-pass Filter implementation
    *   Downsampling and Moving Average routines
    *   Wiener Filter logic

*   **`HeartRateResult.java`**  
    A data model class that encapsulates the estimation output, including:
    *   Estimated Heart Rate (BPM)
    *   Confidence score (optional, internal metric)
    *   Status flags (e.g., if correction was triggered)

## 🚀 Integration Guide

### 1. Copy Files
Simply copy the three `.java` files into your Android project's source path (e.g., `src/main/java/com/yourcompany/algo/`). No external libraries (like NumPy or OpenCV) are required.

### 2. Initialization
Initialize the estimator once when your measurement service starts.

```java
import com.yourcompany.algo.OnlineHeartRateEstimator;
import com.yourcompany.algo.HeartRateResult;

// Initialize with sampling rate (e.g., 25Hz after downsampling or 125Hz raw)
OnlineHeartRateEstimator estimator = new OnlineHeartRateEstimator();
```
### 3. Real-Time Processing
Feed the algorithm with PPG and Accelerometer data windows.
```java
// Example: Processing a window of data
// Input arrays should match the window size defined in the paper (e.g., 8 seconds)

double[] ppgSignal = ...; // Raw PPG data array
double[] accX = ...;      // Accelerometer X-axis
double[] accY = ...;      // Accelerometer Y-axis
double[] accZ = ...;      // Accelerometer Z-axis

// Execute estimation
HeartRateResult result = estimator.execute(ppgSignal, accX, accY, accZ);

if (result.bpm != null) {
    Log.d("HR_Algo", "Current Heart Rate: " + result.bpm);
}
```
