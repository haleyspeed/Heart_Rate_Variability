# Heart_Rate_Variability
Takes 2-13 channels of EEG and EMG activity and determines covariance of brain, respiratory, and cardiac function. 

### Simple Peak/Trough Detection
- Makes use of predefined templates to determine the threshold (using LabChart or Clampfit) for peak/trough detection on each channel
- Baseline can be determined by inter-peak average, overall average, a user-defined window, or from pre-defined template
- Suitable for short recordings (<2-4 hours)
- Currently limited to 2-channel recordings

### Built-In Functions
- Peak amplitude
- Peak to peak latency
- Annotated plots for all analyses
- Weighted averaging for long-term recordings
- Capability for unlimited channels of data

### Future updates
- Spike coherence across channels
- Peak to trough latency
- Trough amplitude

