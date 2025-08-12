# -*- coding: utf-8 -*-
"""
Created on Sun Oct 22 17:16:32 2023

@author: Camille Gontier
"""

# Packages ####################################################################

# Importation of all the required packages
import h5py
import matplotlib.pyplot as plt
import spikeinterface.extractors as se
from probeinterface import Probe, ProbeGroup
from probeinterface.plotting import plot_probe, plot_probe_group
from probeinterface import generate_multi_columns_probe
import numpy as np
import spikeinterface as si  # import core only
import spikeinterface.preprocessing as spre
from spikeinterface.preprocessing import bandpass_filter
import spikeinterface.sorters as ss
import spikeinterface.postprocessing as spost
import spikeinterface.qualitymetrics as sqm
import spikeinterface.comparison as sc
import spikeinterface.exporters as sexp
import spikeinterface.curation as scur
import spikeinterface.widgets as sw
from pprint import pprint
import spikeinterface
from spikeinterface.postprocessing import compute_principal_components
# from spikeinterface.qualitymetrics import (compute_snrs, compute_firing_rates,
#     compute_isi_violations, calculate_pc_metrics, compute_quality_metrics)
# %matplotlib qt

# Parameters ##################################################################

# Path and name of the file to be read
file_name = "data/MP12_Baseline3.h5"

# My data were recorded from 2 chips, each having 60 channels. Here, I am 
# specifying which chip I will study
chip_id = 1
num_channels = 60 # One chip at a time

sampling_frequency = 25000.  # in Hz

# Correspondance between the lines in the data sructure and the electrodes
channel_indices = np.array([10,12,15,16,19,21,
                    7,9,11,14,17,20,22,24,
                    5,6,8,13,18,23,25,26,
                    2,1,3,4,27,28,30,29,
                    59,60,58,57,34,33,31,32,
                    56,55,53,48,43,38,36,35,
                    54,52,50,47,44,41,39,37,
                    51,49,46,45,42,40])-1

# Parallelization, might not be required
global_job_kwargs = dict(n_jobs=4, chunk_duration="1s")
si.set_global_job_kwargs(**global_job_kwargs)

# Data reading ################################################################

# Transforms the file to be read into a Python object. To adapt depending 
# on the format of the recording
f = h5py.File(file_name)
recording_np = f['Data']['Recording_0']['AnalogStream']['Stream_1']['ChannelData'][()]
num_timepoints = recording_np.shape[1]

# Recording object ############################################################

# The recordings need to be transformed into a NumpyRecording object,
# then saved as a .raw object, and then read again using 
# BinaryRecordingExtractor. No idea why, that's the only workaround 
# I could get to work.

recording = se.NumpyRecording(recording_np[(num_channels*(chip_id-1)):(num_channels*chip_id),:].T, sampling_frequency=sampling_frequency)
file_paths = ['traces0.raw']
se.BinaryRecordingExtractor.write_recording(recording, file_paths)
recording = se.BinaryRecordingExtractor(file_paths=file_paths, sampling_frequency=sampling_frequency, num_channels=num_channels, dtype=recording.dtype)

# Plots the first 5 seconds of recordings 
# w_ts = sw.plot_traces(recording, time_range=(0, 5))

# Probe creation ##############################################################

# Specifies the geometry of the probe. In my case, chips were square, with
# 8x8 electrodes.
square_probe = generate_multi_columns_probe(num_columns=8,
                                            num_contact_per_column=[6,8,8,8,8,8,8,6],
                                            xpitch=200, ypitch=200,
                                            y_shift_per_column=[200,0,0,0,0,0,0,200],
                                            contact_shapes='circle', contact_shape_params={'radius': 15})

square_probe.create_auto_shape('rect')

square_probe.set_device_channel_indices(channel_indices)
recording = recording.set_probe(square_probe)

# Displays the electrodes
# plot_probe(square_probe, with_contact_id=True, with_device_index=True)
# plt.show()
# plt.figure
plot_probe(square_probe)
plt.title('')
plt.savefig("probe.svg", dpi=300) 

# Data preprocessing ##########################################################

# Basic preprocessing, i.e. filtering, etc.
recording_cmr = recording
recording_f = bandpass_filter(recording, freq_min=300, freq_max=6000)
recording_cmr = spre.common_reference(recording_f, reference='global', operator='median')
recording_preprocessed = recording_cmr.save(format='binary')

# Sorting algorithm ###########################################################

# Prints all the installed and available sorting algorithms
# print('Installed sorters', ss.installed_sorters())

# I am using mountainsort5, but a different one can be used
sorting = ss.run_sorter(sorter_name="mountainsort5", recording=recording_preprocessed)

print(sorting)
w_rs = sw.plot_rasters(sorting, time_range=(0, 40))
w_rs.ax.set_ylabel("Units")
w_rs.ax.set_xlabel("Time [s]")
w_rs.ax.set_yticks([])
w_rs.figure.savefig("microg1_1.svg",dpi=300)

###############################################################################

print('Unit ids = {}'.format(sorting.get_unit_ids()))
st = sorting.get_unit_spike_train(unit_id=16, segment_index=0)
print('Num. events for unit 1seg0 = {}'.format(len(st)/40))

# Quality control #############################################################

# we = si.extract_waveforms(recording=recording_preprocessed,
#                           sorting=sorting,
#                           folder='waveforms_mearec',
#                           sparse=False,
#                           ms_before=1,
#                           ms_after=2.,
#                           max_spikes_per_unit=500,
#                           n_jobs=1,
#                           chunk_durations='1s')
# print(we)

# firing_rates = compute_firing_rates(we)
# print(firing_rates)
# isi_violation_ratio, isi_violations_count = compute_isi_violations(we)
# print(isi_violation_ratio)
# snrs = compute_snrs(we)
# print(snrs)

# pc = compute_principal_components(waveform_extractor=we, load_if_exists=True,
#                                      n_components=3, mode='by_channel_local')
# print(pc)

# pc_metrics = calculate_pc_metrics(pc, metric_names=['nearest_neighbor'])
# print(pc_metrics)
