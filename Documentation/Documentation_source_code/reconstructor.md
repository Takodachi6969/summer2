# Reconstructor Class Documentation
# Overview
The rpcCoincidence and Reconstructor classes are designed to process and analyze data from proAnubis setup.
It uses heavly the [RTools](RTools documentation.html)

# Table of Contents
- [Reconstructor Class Documentation](#reconstructor-class-documentation)
- [Overview](#overview)
- [Table of Contents](#table-of-contents)
- [`rpcCoincidence` Class](#rpccoincidence-class)
  - [Initialization](#initialization)
  - [Parameters](#parameters)
  - [Methods](#methods)
  - [Description](#description)
- [Reconstructor Class](#reconstructor-class)
  - [Initialization](#initialization-1)
  - [Parameters](#parameters-1)
  - [Attributes](#attributes)
  - [Methods](#methods-1)
    - [self.`populate_hits`](#selfpopulate_hits)
    - [self.`update_event`](#selfupdate_event)
      - [Description](#description-1)
      - [Parameters](#parameters-2)
    - [self.`make_cluster`](#selfmake_cluster)
    - [Description](#description-2)
    - [Returns](#returns)
    - [Related Functions](#related-functions)
    - [self.`reconstruct_and_extrapolate`](#selfreconstruct_and_extrapolate)
    - [Description](#description-3)
    - [Parameters](#parameters-3)
    - [Related Functions](#related-functions-1)
    - [self.`reconstruct_and_findtof`](#selfreconstruct_and_findtof)
    - [Description](#description-4)
    - [Parameters](#parameters-4)
    - [self.`extract_angles_phi_eta_timed_DZ_modified`](#selfextract_angles_phi_eta_timed_dz_modified)
    - [Description](#description-5)
    - [Parameters](#parameters-5)
    - [self.`apply_systematic_correction`](#selfapply_systematic_correction)
    - [Description](#description-6)
    - [Parameters](#parameters-6)
  - [Usage Example](#usage-example)


# `rpcCoincidence` Class
## Initialization
```python
def __init__(self, event_num, time_bin, hits):
```

## Parameters
- **event_num (int)**: The event number.
- **time_bin (int)**: The time bin associated with the event.
- **hits (list)**: The list of hits recorded in the event.

## Methods
`__str__`
```python
def __str__(self):
```

## Description
Returns a string representation of the rpcCoincidence object.

# Reconstructor Class
## Initialization
```python
def __init__(self, event_chunk, processsed_event, tolerance=None, coincidence_window=15, tof_correction=True):
```

## Parameters
- **event_chunk (list)**: List of events to be processed.
- **processsed_event (int)**: The number of processed events.
- **tolerance (list, optional)**: List of tolerance values for reconstruction. Defaults to a range of 20 values if not provided.
- **coincidence_window (int, optional)**: The window size for coincidence detection. Default is 15.
- **tof_correction (bool, optional)**: Boolean to apply time-of-flight correction. Default is True.

## Attributes

- **event_chunk (list)**: List of events to be processed.
- **etaHits (list)**: List of eta hits for each RPC.
- **phiHits (list)**: List of phi hits for each RPC.
- **window_size (int)**: Window size for coincidence detection.
- **tof_correction (bool)**: Boolean to apply time-of-flight correction.
- **processedEvents (int)**: Number of processed events.
- **tol (list)**: List of tolerance values.
- **dT (list)**: List of time differences.
- **recon (list)**: List of reconstructions.
- **possible_reconstructions (list)**: List of possible reconstructions per RPC.
- **successful_reconstructions (list)**: List of successful reconstructions per RPC for each tolerance value.
- **successful_reconstructed_coords (dict)**: Dictionary of successful reconstructed coordinates for each RPC.
- **failed_reconstructed_coords (dict)**: Dictionary of failed reconstructed coordinates for each RPC.
- **possible_reconstructions_coords (dict)**: Dictionary of possible reconstructed coordinates for each RPC.
- **eta_histogram (numpy array)**: Histogram for eta angles.
- **phi_histogram (numpy array)**: Histogram for phi angles.
- **solid_theta_histogram (numpy array)**: Histogram for solid theta angles.
- **solid_phi_histogram (numpy array)**: Histogram for solid phi angles.
- **tdcstatus (list)**: Status of the TDCs.
  
## Methods
### self.`populate_hits`

Populates the eta and phi hits from the event chunk.

```python
def populate_hits(self):
        self.etaHits = [[] for rpc in range(6)]
        self.phiHits = [[] for rpc in range(6)]
        skip_event = False
        for idx, event in enumerate(self.event_chunk):
            for tdc in range(5):
                if event.tdcEvents[tdc].qual != 0:
                    skip_event = True
                    break 

            if skip_event:
                continue 
            for tdc in range(5):
                for word in event.tdcEvents[tdc].words:
                    rpc, thisHit = ATools.tdcChanToRPCHit(word,tdc, self.processedEvents + idx)
                    if thisHit.channel == [0]:
                        continue
                    if thisHit.eta:
                        self.etaHits[rpc].append(thisHit)

                    else:
                        self.phiHits[rpc].append(thisHit)
                        
```
---

### self.`update_event`
```python
def update_event(self, event_chunk, processed_event):
        self.event_chunk = event_chunk
        self.processedEvents = processed_event
        self.etaHits = [[] for rpc in range(6)]
        self.phiHits = [[] for rpc in range(6)]
```
#### Description
Updates the event data to a new chunk without changing other stored instances.

#### Parameters
- **event_chunk (list)**: List of new events to be processed.
- **processedEvents (int)**: Number of new processed events.

---

### self.`make_cluster`
```python
def make_cluster(self):
        coincident_hits = RTools.FindCoincidentHits(self.etaHits, self.phiHits, self.window_size, tof_correction=self.tof_correction)
        clustered = RTools.cluster(coincident_hits)
        return clustered
```

### Description
Creates clusters of coincident hits from the eta and phi hits.

### Returns
- clustered (list): List of clustered hits.

### Related Functions
[FindCoincidentHits](Rtools documentation.html#rtoolsfindcoincidenthits)
[Cluster](Rtools documentation.html#rtoolscluster)

---

### self.`reconstruct_and_extrapolate`
```python
def reconstruct_and_extrapolate(self, dataset, chi2_region = [0, 100]):
        # Ensure RPC is a list, even if it's a single integer
        if self.tdcstatus[3] == True:
            for rpc in range(6):
                for i, data in enumerate(dataset):
                    if RTools.count_entries(data) < 100:
                        E_recon = RTools.reconstruct_timed_Chi2_ByRPC(data, 3, rpc)
                        if E_recon:
                            if len(E_recon[2]) >= 5:
                                if E_recon[4] > chi2_region[0] and E_recon[4] < chi2_region[1]:


                                    muon_coords = RTools.does_muon_hit_RPC(E_recon[0], E_recon[1], rpc)
                                    if muon_coords:
                                        self.possible_reconstructions[rpc] += 1
                                        self.possible_reconstructions_coords[rpc][int(muon_coords[0] / 2.7625)][int(muon_coords[1] / 2.9844)] += 1
                                        for idx, t in enumerate(self.tol):
                                            check = RTools.does_RPC_detect_muon(muon_coords, E_recon[7], t)
                                            if check:
                                                self.successful_reconstructions[rpc][idx] += 1
                                                self.successful_reconstructed_coords[rpc][int(muon_coords[0] / 2.7625)][int(muon_coords[1] / 2.9844)] += 1
                                            else:
                                                self.failed_reconstructed_coords[rpc][int(muon_coords[0] / 2.7625)][int(muon_coords[1] / 2.9844)] += 1

```
### Description
Performs reconstruction and extrapolation of the dataset.

### Parameters
- **dataset (list)**: List of data to be reconstructed.
- **chi2_region (list, optional)**: Chi-squared region for filtering. Default is [0, 100].

### Related Functions
[RTools](RTools documentation.html)

---

### self.`reconstruct_and_findtof`
```python
def reconstruct_and_findtof(self, dataset, rpc_comparisons):
        for i, data in enumerate(dataset):
            if RTools.count_entries(data) < 100:
                E_recon = RTools.reconstruct_timed_Chi2_ByRPC(data, 3, -1, rpc_indicies=rpc_comparisons)
                if E_recon:
                    if len(E_recon[2]) >= 6:
                        self.dT.append(E_recon[5])
                        self.recon.append(E_recon)
```

### Description
Performs reconstruction and finds the time-of-flight differences.
### Parameters
- **dataset (list)**: List of data to be reconstructed.
- **rpc_comparisons (list)**: List of RPC comparisons for the reconstruction.

### self.`extract_angles_phi_eta_timed_DZ_modified`

### Description
Extracts and analyzes angles and time differences from filtered events. The structure is very similar to `reconstruct_and_extrapolate`, with extra steps to turn the reconstructed line into angles

### Parameters
- **filtered_events (list)**: List of filtered events.
- **max_length (int, optional)**: Maximum length for reconstruction. Default is None.
- **exact_length (bool, optional)**: Flag to enforce exact length in reconstruction. Default is False.


---

### self.`apply_systematic_correction`
```python
def apply_systematic_correction(self, residEta, residPhi):
        for rpc in range(6):
            for i, etahit in enumerate(self.etaHits[rpc]):
                self.etaHits[rpc][i].time += residEta[rpc][etahit.channel]
            for j , phihit in enumerate(self.phiHits[rpc]):
                self.phiHits[rpc][j].time += residPhi[rpc][phihit.channel]
```

### Description
Applies systematic corrections to the eta and phi hits.

### Parameters
- **residEta (list)**: List of eta residuals for correction.
- **residPhi (list)**: List of phi residuals for correction.

## Usage Example
```python
# Initialize an event chunk and processed event number
event_chunk = [...]  # List of events
processed_event = 100

# Create a Reconstructor object
reconstructor = Reconstructor(event_chunk, processed_event)

# Populate hits from the event chunk
reconstructor.populate_hits()

# Update event data
new_event_chunk = [...]  # New list of events
new_processed_event = 200
reconstructor.update_event(new_event_chunk, new_processed_event)

# Create clusters of coincident hits
clusters = reconstructor.make_cluster()

# Perform reconstruction and extrapolation
reconstructor.reconstruct_and_extrapolate(clusters)

# Perform reconstruction and find time-of-flight differences
rpc_comparisons = [...]  # List of RPC comparisons
reconstructor.reconstruct_and_findtof(clusters, rpc_comparisons)

# Extract angles and time differences
filtered_clusters = [...]  # List of filtered events
reconstructor.extract_angles_phi_eta_timed_DZ_modified(filtered_clusters)

# Apply systematic corrections
residEta = [...]  # List of eta residuals
residPhi = [...]  # List of phi residuals
reconstructor.apply_systematic_correction(residEta, residPhi)

# Plot time-of-flight offsets
reconstructor.plot_tof_offset(rpc_comparisons)
```