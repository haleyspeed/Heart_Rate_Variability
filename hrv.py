""" Heartrate Variability (Cardio-Respiratory Physiologyical Synchrony)
for multi-channel recordings

Written by Haley Speed (haleygeek@gmail.com)
Github: https://github.com/haleyspeed/Heart_Rate_Variability
"""

#------------------------- User Input ---------------------------------------#
file_in = 'Synchrony_EEG_ECG_RESP_K6527_HmzNdufsCrep_10_06_16_TISel.edf'

# Channel configuration
unfiltered_breath_channel = 8
filtered_breath_channel = 9
ecg_channel = 12
left_eeg_channel = 10
right_eeg_channel = 11
emg_channel = 13

# Set up analysis window
start_mins = 0.0001 # The time (in min) to start analyzing the recording
length_mins = "all"  # Length of Recording you wannt to analyze, or 'all'

sample_rate = 1000 # In Hz

# From template designer "Refine Event Search"
unfiltered_breath_channel_baseline = -0.026
unfiltered_breath_channel_std = 0.022
filtered_breath_channel_baseline = -0.023
filtered_breath_channel_std = 0.021
ecg_channel_baseline = -0.001
ecg_channel_std = 0.003

# Thresholds for event detection (* SD)
thresh_ufbreath = 2    # Default = 2 * SD
thresh_breath = 2    # Default = 2 * SD
thresh_ecg = 4   # Default = 4 * SD Avoids contamination of Q-R peak by P or T

#----------------------------- Import Packages -------------------------------#

import pyedflib
import os
import numpy as np
import matplotlib.pyplot as plt
import pandas as pd

#--------------------------- Function Definitions ----------------------------#

def dp2min(datapoints):
    """ Takes either a list or an int and converts
    to minutes based on sample rate"""
    if type(datapoints) == list:
        minutes = []
        for dp in datapoints:
            mins = dp / (60 * sample_rate)
            minutes.append(mins)
    else:
        minutes = datapoints / (60 * sample_rate)
    return minutes

def dp2s(datapoints):
    """ Takes either a list or an int and converts
    to milliseconds based on sample rate"""
    if type(datapoints) == list:
        s = []
        for dp in datapoints:
            s = dp / sample_rate
            s.append(int(s))
    else:
        s = datapoints/ sample_rate
    return s

def dp2ms(datapoints):
    """ Takes either a list or an int and converts
    to milliseconds based on sample rate"""
    if type(datapoints) == list:
        milliseconds = []
        for dp in datapoints:
            ms = dp / sample_rate * 1000
            milliseconds.append(int(ms))
    else:
        milliseconds = datapoints/ sample_rate *1000
    return milliseconds

def min2dp (minutes):
    """ Takes either a list or a float and converts
    to datapoints (int) based on sample rate"""
    if type(minutes) == list:
        datapoints = []
        for mins in minutes:
            dp = mins * (60 * sample_rate)
            datapoints.append(int(dp))
    else:
        datapoints = int(minutes * 60 * sample_rate)
    return datapoints

def s2dp (seconds):
    """ Takes either a list or a float and converts
    to datapoints (int) based on sample rate"""
    if type(seconds) == list:
        datapoints = []
        for sec in seconds:
            dp = sec * sample_rate
            datapoints.append(int(dp))
    else:
        datapoints = int(seconds * sample_rate)
    return datapoints

def ms2dp (milliseconds):
    """ Takes either a list or a float and converts
    to datapoints (int) based on sample rate"""
    if type(seconds) == list:
        datapoints = []
        for millisec in milliseconds:
            dp = millisec * sample_rate * 1000
            datapoints.append(int(dp))
    else:
        datapoints = int(milliseconds * sample_rate * 1000)
    return datapoints

def get_edf (dir_in, file_in, file_open):
    """ Opens a Labchart EDF file from the specified file and directory"""
    if file_open == 0:
        file_name = os.path.join(dir_in, file_in)
        f = pyedflib.EdfReader(file_name)
        file_open = 1
    else:
        pass
    n = f.signals_in_file
    signal_labels = f.getSignalLabels()
    sigbufs = np.zeros((n, f.getNSamples()[0]))
    for i in np.arange(n):
        sigbufs[i,:] = f.readSignal(i)
    return signal_labels, sigbufs

def get_plot_length (plot_length, sigbufs):
    """ Takes in a plot length in minutes (whole trace or subsection of a trace)
        and returns it in datapoints"""
    if plot_length == 'all':
        plot_end_time = len(sigbufs[0,:])
    elif plot_length > len(sigbufs[0,:]):
        plot_end_time = len(sigbufs[0,:])
    elif plot_length < 0:
        print ("USER ERROR: Cannot plot a negative plot_length.")
    else:
        plot_end_time = plot_length * (60*sample_rate)
    return int(plot_end_time)


def get_events(x, y, threshold, baseline):
    """ Searches the trace for the specified length (in datapoints)
    based on the user-specified threshold and baseline, with option
    to plot the results (True/False)"""
    events = []
    if threshold == 0: # first pass assumes 0 (ground) in template search
        threshold = ((max(y) - min(y))*.1)
    baseline_x = []
    baseline_y = []

    i = 0
    while i < len(x)-1:
        if y[i] < threshold and y[i] > -1*threshold:
            while y[i] < threshold and i < len(x)-1:
                baseline_x.append (x[i])
                baseline_y.append (y[i])
                i = i + 1
        elif y[i] >= threshold:
            event_x = []
            event_y = []
            peak_x = 0
            peak_y = 0
            while y[i] > threshold and i < len(x)-1:
                event_x.append (x[i])
                event_y.append (y[i])
                i = i + 1
            event = {'signal': signal, 'x': event_x, 'y': event_y}
            events.append (event)
        else:
            i = i + 1
    baseline_mean = np.mean(baseline_y)
    baseline_std = np.std(baseline_y)
    return [events, baseline_mean, baseline_std]


def get_event_stats (events, baseline):
    df = pd.DataFrame({'peaks_y': [], 'peaks_x' : [], 'peaks_amp': [],
    'widths' : [], 'half_width_x':[], 'half_width_y':[]})
    for event in events:
        peak_y = max (event['y'])
        peak_amp = max (event['y']) - baseline
        width = max(event['x']) - min(event['x'])
        half_peak = peak_y/2
        half_x = 0
        half_y = 0
        i = 0
        while i < len(event['x'])-1:
            if event['y'][i] < event['y'][i + 1]:
                peak_x = event['x'][i]
                i = i + 1
            else:
                i = len(event['x'])
        row = {'peaks_y': peak_y, 'peaks_x' : peak_x, 'peaks_amp': peak_amp,
             'widths' : width, 'half_width_x': half_x, 'half_width_y' : half_y}
        df = df.append(row, ignore_index = True)
    return df


#--------------------------- Default variables -------------------------------#

# Assumes that the python script (hrv.py) and edf file are in the same dir
dir_in = str(os.getcwd())

# Channel refers to Labchart numbering, which starts at 1
channel_list = {'Unfiltered_Breath' : unfiltered_breath_channel,
                'Filtered_Breath' : filtered_breath_channel,
                'ECG': ecg_channel,
                'Left_EEG' : left_eeg_channel,
                'Right_EEG' : right_eeg_channel,
                'EMG': emg_channel}
amp_channels =  {'Unfiltered_Breath' : unfiltered_breath_channel,
                'Filtered_Breath' : filtered_breath_channel,
                'ECG': ecg_channel}
freq_channels = {'Left_EEG' : left_eeg_channel,
                'Right_EEG' : right_eeg_channel,
                'EMG': emg_channel}

# Signal refers to Python numbering, which starts at 0
signal_list = channel_list.copy()
label_list = []
for keys, values in signal_list.items():
    signal_list[keys] = values - 1
    label_list.append(keys)

amp_signals = amp_channels.copy()
amp_labels = []
for keys, values in amp_signals.items():
    amp_signals[keys] = values - 1
    amp_labels.append(keys)

freq_signals = freq_channels.copy()
freq_labels = []
for keys, values in freq_signals.items():
    freq_signals[keys] = values - 1
    freq_labels.append(keys)

file_open = 0   # 0 = EDF file hasn't been opened, 1 = it has been opened

# Read in the raw traces for each channel
signal_labels, sigbufs = get_edf (dir_in, file_in, file_open)

# Traces per channel
x_dp = list(range (0,len(sigbufs[0]))) # x values for all channels in datapoints
x_mins = [dp2min(x) for x in x_dp] # x values for all channels in minutes
x_s = [dp2s(x) for x in x_dp] # x values for all channels in seconds
y_breath_unfiltered = sigbufs[signal_list['Unfiltered_Breath'],:]
y_breath = sigbufs[signal_list['Filtered_Breath'],:]
y_ecg = sigbufs[signal_list['ECG'],:]
y_eeg_l = sigbufs[signal_list['Left_EEG'],:]
y_eeg_r = sigbufs[signal_list['Right_EEG'],:]
y_emg = sigbufs[signal_list['EMG'],:]
y_amp = [y_breath_unfiltered, y_breath, y_ecg]
y_freq = [y_eeg_l, y_eeg_r, y_emg]
y_list = y_amp + y_freq

# Window variables
start_dp = min2dp(start_mins)
start_s = dp2s(start_dp)
start_ms = dp2ms(start_dp)

# Get plot length in datapoints
end_mins = start_mins + get_plot_length (length_mins, sigbufs)
end_dp = min2dp (end_mins)
end_s = dp2s (end_dp)
end_ms = dp2ms (end_dp)

# Event detection
thresholds = {'Unfiltered_Breath': thresh_ufbreath,
              'Filtered_Breath': thresh_breath,
              'ECG': thresh_ecg}

baseline = [unfiltered_breath_channel_baseline,
            filtered_breath_channel_baseline,
            ecg_channel_baseline]
std = [unfiltered_breath_channel_std,
        filtered_breath_channel_std,
        ecg_channel_std]


#------------------------------- Event Search -------------------------------#

baselines = []
stds = []
events_list = []
desc_list = []
i = 0
for signal in amp_signals:
    result = get_events(x_dp,
                        y_amp[i],
                        thresholds[amp_labels[i]] * std[i],
                        baseline[i])
    """ result = events, baseline_mean, baseline_std """
    events_list.append(result[0])
    baselines.append(result[1])
    stds.append(result[2])
    
    # Get peaks for each event
    df_events = get_event_stats (result[0], 0)
    desc_list.append(df_events)
    df_events.to_csv(file_in.replace('.edf', '_' + str(amp_labels[i])
        + '_events.csv'))
    i = i + 1

# Compare event features to template and eliminate as needed
# From here move back into notebooks
# Finish out analysis as written in notebooks
