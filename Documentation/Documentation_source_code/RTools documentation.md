# Reconstruction Tools Documentation
This package contains everything needed to do recontruction algorithm in the proAnubis setup. This includes but not limited to , temporal and spatial coincidence, track fitting, efficiency and angular distribution calculations

# Table of content
- [Reconstruction Tools Documentation](#reconstruction-tools-documentation)
- [Table of content](#table-of-content)
  - [RTools.`find_tdc_event_count`](#rtoolsfind_tdc_event_count)
    - [**Description**](#description)
    - [**Parameters**](#parameters)
    - [**Returns**](#returns)
    - [**Note**](#note)
    - [**related function**](#related-function)
  - [RTools.`find_tof_time`](#rtoolsfind_tof_time)
    - [**Definition**](#definition)
  - [Rtools.`FindCoincidentHits`](#rtoolsfindcoincidenthits)
    - [**Description**](#description-1)
    - [**Parameters**](#parameters-1)
    - [**Returns**](#returns-1)
    - [**related functions**](#related-functions)
  - [find\_tof\_time](#find_tof_time)
  - [Rtools.`cluster`](#rtoolscluster)
    - [**Description**](#description-2)
    - [**Parameters**](#parameters-2)
    - [**Returns**](#returns-2)
    - [**Note**](#note-1)
  - [The weird data structure was designed to adapt to the reconstruction fit, which was written by someone else...](#the-weird-data-structure-was-designed-to-adapt-to-the-reconstruction-fit-which-was-written-by-someone-else)
  - [Rtools.`extract_coords_timed_Chi2`](#rtoolsextract_coords_timed_chi2)
    - [**Description**](#description-3)
    - [**Parameters**](#parameters-3)
    - [**Returns**](#returns-3)
  - [Rtools.`generate_hit_coords_combo_Chi2`](#rtoolsgenerate_hit_coords_combo_chi2)
    - [**Description**](#description-4)
    - [**Definition**](#definition-1)
    - [**Parameters**](#parameters-4)
    - [Returns](#returns-4)
  - [Rtools.`extract_DT_DZ_Chi2`](#rtoolsextract_dt_dz_chi2)
    - [**Description**](#description-5)
    - [**Parameters**](#parameters-5)
    - [**Returns**](#returns-5)
    - [**Note**](#note-2)
  - [This function contains magic number for the speed of signal propagation, which needs to be changed if the fit is different](#this-function-contains-magic-number-for-the-speed-of-signal-propagation-which-needs-to-be-changed-if-the-fit-is-different)
  - [Rtools.`fit_event_chi2`](#rtoolsfit_event_chi2)
    - [**Description**](#description-6)
    - [**Parameters**](#parameters-6)
  - [**Returns**](#returns-6)
    - [External links](#external-links)
  - [Rtools.`reconstruct_timed_Chi2_ByRPC`](#rtoolsreconstruct_timed_chi2_byrpc)
    - [**Description**](#description-7)
    - [**Parameters**](#parameters-7)
    - [**Returns**](#returns-7)
    - [**Explainations and Pseudocode**](#explainations-and-pseudocode)
  - [Rtools.`does_muon_hit_RPC`](#rtoolsdoes_muon_hit_rpc)
    - [**Description**](#description-8)
    - [**Parameters**](#parameters-8)
    - [**Returns**](#returns-8)
  - [Rtools.`does_RPC_detect_muon`](#rtoolsdoes_rpc_detect_muon)
    - [**Description**](#description-9)
    - [**Parameters**](#parameters-9)
    - [**Returns**](#returns-9)
    - [**Note**](#note-3)
  - [Rtools.`check_event_attributes`](#rtoolscheck_event_attributes)
    - [**Description**](#description-10)
    - [**Parameters**](#parameters-10)
    - [**Returns**](#returns-10)
  - [Rtools.`filter_events`](#rtoolsfilter_events)
    - [**Description**](#description-11)
    - [**Parameters**](#parameters-11)
    - [**Returns**](#returns-11)
  - [Rtools.`count_entries`](#rtoolscount_entries)
    - [***Description***](#description-12)
    - [**Parameters**](#parameters-12)
    - [**Returns**](#returns-12)
  - [Rtools.`find_tdc_event_count`](#rtoolsfind_tdc_event_count-1)
    - [**Description**](#description-13)
    - [**Parameters**](#parameters-13)
    - [**Returns**](#returns-13)
  - [Rtools.`compile_and_plot_tof`](#rtoolscompile_and_plot_tof)
    - [**Description**](#description-14)
  - [**Parameters**](#parameters-14)
  - [**Returns**](#returns-14)
  - [Rtools.`compile_and_plot_tof_chunk`](#rtoolscompile_and_plot_tof_chunk)
    - [**Description**](#description-15)
    - [**Parameters**](#parameters-15)
    - [**Returns**](#returns-15)
  - [Example Usage](#example-usage)

## RTools.`find_tdc_event_count`


### **Description**
calculates the total number of words for each TDC across a chunk of events. It processes the events and sums the lengths of the `words` lists for each TDC.

### **Parameters**
- **event_chunk** (`list`): A list of proanubis event objects. Each event object contains an attribute `tdcEvents`, which is a list of TDC objects. Each TDC object has an attribute `words`, which is a list of words.

### **Returns**
- **event_number** (`list`): A list of lists. Each inner list corresponds to a TDC and contains a single integer representing the total count of words for that TDC across all events in the given chunk.

### **Note**
- Function assumes exactly **5** TDCs
- The returned `event_number` list contains 5 lists, each with a single integer value representing the total count of words for the corresponding TDC.

### **related function**
[Realigner](Realignment documentation.html)

## RTools.`find_tof_time`
### **Definition**
```python
def find_tof_time(eta, phi, slope = 0.05426554612593516, offSet = 15.8797407836404):
    if (len(set([eta.eta, phi.eta])) == 1):
        return 0
    else:
        return slope*(phi.channel-eta.channel)-offSet

```

---

## Rtools.`FindCoincidentHits`

### **Description**
This function processes eta and phi hits from different RPC arrays to find coincident hits based on time windows and optional Time-of-Flight (TOF) correction.

### **Parameters**
- **etaHits (list)**: A list of lists containing eta RPChit objects for each RPC. Each hit object must have attributes time, channel, event_num, and eta.
- **phiHits (list)**: A list of lists containing phi RPChit objects for each RPC. Each hit object must have attributes time, channel, event_num, and eta.
- **time_window (float)**: The time window in which to consider hits as coincident.
- **tof_correction (bool, optional)**: Whether to apply TOF correction. Default is True.
- **slope (float, optional)**: Slope parameter for TOF correction. Default is 0.05426554612593516.
- **offSet (float, optional)**: Offset parameter for TOF correction. Default is 15.8797407836404.

### **Returns**
- **coincident_hits (list)**: A list of lists where each sublist represents an event with coincident hits. Each sublist contains the event number, time bin, and details of the hits.
``` python
[event_num, time_bin,[[hit.rpc, hit.channel, hit.time, hit.eta] for hit in unique_hits]]
```

### **related functions**
(#populatehits)
[rpcHits](ATools documentation.html#class-rpchit)
[find_tof_time](#rtoolsfind_tof_time)
---

## Rtools.`cluster`
### **Description**
This function clusters the hits within each RPC into spatial clusters in both phi and eta directions. Putting all adjacent strips together

### **Parameters**
- **coincident_hits (list)**: A list of coincident hit events. Each event is a list containing the event number, time bin, and hit details.

### **Returns**
- **coincident_hits_clustered (list)**: A list of events with clustered hits. Each event contains the event number, time bin, and clusters of hits for each RPC.
```python
#rpc_phi_clusters is a list of sublists [hit.rpc, hit.channel, hit.time, hit.eta]
[event_num, time_bin, [rpc_phi_clusters, rpc_eta_clusters]]
```

### **Note**
The weird data structure was designed to adapt to the reconstruction fit, which was written by someone else...
---

## Rtools.`extract_coords_timed_Chi2`

### **Description**
Converts spatial clusters in RPCs into x and y coordinates, with z given by the RPC number.

### **Parameters**
- **clusters (list)**: A list representing a cluster, the output of function [cluster](#rtoolscluster). It contains the event number, time bin, and clusters of hits.
- **max_cluster_size (int)**: Maximum allowed size for a cluster.
  
### **Returns**
- **coords (list)**: A list of coordinates with errors. Each entry represents an RPC and contains lists of x and y coordinates, with variances and times.
```python
    [x_coords] = [[x,err_x,x_time],...]
    
    RPC_coords = [x_coords,y_coords]

    coords = [[RPC1_coords],[RPC2_coords],[RPC3_coords],...]
```
---


## Rtools.`generate_hit_coords_combo_Chi2`
### **Description**
Generates combinations of hit coordinates across RPCs for chi-squared fitting.

### **Definition**
```python
def generate_hit_coords_combo_Chi2(coords, RPC_heights, max_length=None, exact_length=False, combinations=None, hit_coords=None, depth=0):
    if combinations is None:
        combinations = []
    if hit_coords is None:
        hit_coords = []
    if max_length is None:
        max_length = len(coords)

    # Base case: If we've reached the end of the coords or the length condition is met
    if depth == len(coords) or len(hit_coords) == max_length:
        if not exact_length or len(hit_coords) == max_length:
            combinations.append(hit_coords.copy())
        return combinations

    # Extract x and y values for the current depth
    x_values = coords[depth][0]
    y_values = coords[depth][1]

    # If there are no x or y values at this depth, move to the next depth
    if not x_values or not y_values:
        return generate_hit_coords_combo_Chi2(coords, RPC_heights, max_length, exact_length, combinations, hit_coords, depth + 1)

    # Iterate over all combinations of x and y values
    for x in x_values:
        for y in y_values:
            if x is not None and y is not None and isinstance(x[0], (int, float)) and isinstance(y[0], (int, float)):
                hit_coords.append([x, y, RPC_heights[depth]])
                generate_hit_coords_combo_Chi2(coords, RPC_heights, max_length, exact_length, combinations, hit_coords, depth + 1)
                hit_coords.pop()

    return combinations
```

### **Parameters**
- **coords (list)**: A list of coordinates with errors.
- **RPC_heights (list)**: A list of heights for each RPC.
- **max_length (int, optional)**: Maximum length of the combination.
- **exact_length (bool, optional)**: Whether to generate combinations of exact length.
  
### Returns
- **combinations (list)**: A list of hit coordinate combinations.

---

## Rtools.`extract_DT_DZ_Chi2`
### **Description**
Extracts the time difference (dT) and height difference (dZ) between specified pairs of RPCs.

### **Parameters**
- **coords (list)**: A list of coordinates with time information for each RPC.
- **rpc_indices (list of lists)**: A list of pairs of indices specifying the RPCs to compare.

### **Returns**
- **dT_all (list)**: A list of time differences for each pair of RPC indices.
- **dZ_all (list)**: A list of height differences for each pair of RPC indices.

### **Note**
This function contains magic number for the speed of signal propagation, which needs to be changed if the fit is different
---

## Rtools.`fit_event_chi2`
### **Description**
Fits a line to the event coordinates using Singular Value Decomposition (SVD) and calculates the chi-squared value.

### **Parameters**
- **combo (list)**: A list of coordinates with errors for each RPC, output from [generate_hit_coords_combo_chi2](#rtoolsgenerate_hit_coords_combo_chi2).
- **rpc_indicies (list of lists, optional)**: A list of pairs of indices specifying the RPCs to compare.

## **Returns**
- **p0 (list)**: Centroid of the fitted line.
- **d (list)**: Direction vector of the fitted line.
- **Chi2 (float)**: Chi-squared value of the fit.
- **coordinates (list)**: List of coordinates used in the fit.
- **dT (list)**: List of time differences.
- **dZ (list)**: List of height differences.

### External links
[SVD fits](https://en.wikipedia.org/wiki/Singular_value_decomposition)

---

## Rtools.`reconstruct_timed_Chi2_ByRPC`
### **Description**
Reconstructs an event by excluding one RPC and fitting the remaining coordinates.

### **Parameters**
- **event (list)**: A list representing an event.
- **max_cluster_size (int)**: Maximum allowed size for a cluster.
- **RPC_excluded (int)**: The RPC to exclude from the reconstruction.
- **rpc_indicies (list of lists, optional)**: A list of pairs of indices specifying the RPCs to compare.

### **Returns**
- **reconstruction_result (tuple or None)**: If successful, returns a tuple containing the **optimized centroid**, **direction vector**, **coordinates**, **combinations**, **chi-squared value**, **time differences**, **height differences**, and **test coordinates**. Otherwise, returns None.

### **Explainations and Pseudocode**

1. **Initialize Parameters and Variables**:
    - Set `max_Chi2` to 10.
    - Define `RPC_heights` as the heights of the middle point of each RPC.

2. **Extract Coordinates**:
    - Call `extract_coords_timed_Chi2` with `event` and `max_cluster_size`.
    - Store the returned coordinates in `coords`.

3. **Exclude Specified RPC**:
    - Initialize `test_coords` to -1.
    - If `RPC_excluded` is not -1:
        - Set `test_coords` to the coordinates of the excluded RPC (`coords[RPC_excluded]`).
        - Set the coordinates of the excluded RPC to `[[], [], "N"]` in `coords`.

4. **Check for Sufficient Data**:
    - Count the number of empty RPCs in `coords`.
    - If more than 4 RPCs are empty, return `None`.

5. **Check Cross-Chamber Requirement**:
    - Initialize `cross_chamberness` to 0.
    - Check if there are sufficient clusters in the top, middle, and bottom chambers.
    - If `cross_chamberness` is less than 2, return `None`.

6. **Generate Coordinate Combinations**:
    - Call `generate_hit_coords_combo_Chi2` with `coords` and `RPC_heights`.
    - Store the returned combinations in `combinations`.
    - If the number of combinations exceeds 20, return `None`.

7. **Fit Each Combination**:
    - Initialize `Chi2_current` to infinity.
    - Initialize `optimised_coords`, `optimised_d`, and `optimised_centroid` to `None`.
    - Initialize `dT` to a list containing `None`.
    - Loop through each combination in `combinations`:
        - If the length of the combination is less than 5, skip to the next combination.
        - Call `fit_event_chi2` with the combination and `rpc_indicies`.
        - If the returned `Chi2` is less than `Chi2_current`:
            - Update `dZ`, `dT`, `Chi2_current`, `optimised_centroid`, `optimised_d`, and `optimised_coords` with the new values.

8. **Adjust Direction Vector Based on Time Differences**:
    - If `dT[-1]` is not `None`:
        - If `dT[-1]` is greater than 0 and `optimised_d[2]` is less than 0, multiply `optimised_d` by -1.
        - If `dT[-1]` is less than 0 and `optimised_d[2]` is greater than 0, multiply `optimised_d` by -1.

9. **Check Chi-Squared Value**:
    - If `Chi2_current` is less than `max_Chi2`, return the optimized results.
    - Otherwise, return `None`.

---


## Rtools.`does_muon_hit_RPC`
### **Description**
Checks if a reconstructed muon hits a specified RPC based on its trajectory.

### **Parameters**
- **optimised_centroid (list)**: Centroid of the reconstructed muon trajectory.
- **optimised_d (list)**: Direction vector of the reconstructed muon trajectory.
- **RPC (int)**: The RPC to check.

### **Returns**
- **hit_coords (list or None)**: If the muon hits the RPC, returns the **coordinates** of the hit. Otherwise, returns **None**.

---

## Rtools.`does_RPC_detect_muon`
### **Description**
Checks if an RPC detects a muon based on the reconstructed hit coordinates and a given tolerance.

### **Parameters**
- **muon_coords (list)**: Expected coordinates of the muon hit.
- **test_coords (list)**: Actual coordinates reconstructed by the RPC.
- **tol (float)**: Tolerance for the detection check.

### **Returns**
- detection_result (bool or float): If the RPC detects the muon, returns the **separation distance**. Otherwise, returns **False**.

### **Note**

Tolerance is defines in cm as the maximum distance the algorithm look out for detected muons on the RPC. 


---

## Rtools.`check_event_attributes`
### **Description**
Checks if an event meets the minimum required number of RPCs and chambers with hits.

### **Parameters**
- **event (list)**: A list representing an event.
- **min_chamber_number (int)**: Minimum required number of chambers with hits.
- **min_RPC_number (int)**: Minimum required number of RPCs with hits.

### **Returns**
- **valid_event (bool)**: **True** if the event meets the requirements, otherwise **False**.

---

## Rtools.`filter_events`
### **Description**
Filters a list of events based on the minimum required number of RPCs and chambers with hits.

### **Parameters**
- **events (list)**: A list of events.
- **min_chamber_number (int)**: Minimum required number of chambers with hits.
- **min_RPC_number (int)**: Minimum required number of RPCs with hits.

### **Returns**
- **filtered_events (list)**: A list of events that meet the requirements.

---

## Rtools.`count_entries`
### ***Description***
Recursively counts the total number of entries in a nested list.

### **Parameters**
- **lst (list)**: A nested list.

### **Returns**
- **total_count (int)**: Total number of entries in the list.

---

## Rtools.`find_tdc_event_count`
### **Description**
Calculates the total number of words for each TDC (Time-to-Digital Converter) across a chunk of events.

### **Parameters**
- **event_chunk (list)**: A list of event objects. Each event object contains a list of TDC objects with a words attribute.

### **Returns**
- **event_number (list)**: A list of lists. Each inner list corresponds to a TDC and contains a single integer representing the total count of words for that TDC across all events in the chunk.

## Rtools.`compile_and_plot_tof`
### **Description**
Compiles time-of-flight (TOF) data and generates Gaussian fits for specified RPC pairs, saving the results to a PDF.

## **Parameters**
- dTs (list): A list of time differences for events.
- rpc_indicies (list of lists): A list of pairs of RPC indices for TOF analysis.
- pdf_filename (str, optional): The filename for the output PDF. Default is "Data_output/compiled_tof_plots.pdf".

## **Returns**
- **pdf_filename (str)**: The filename of the generated PDF.
- **plots**: plots for the tof analysis
---


## Rtools.`compile_and_plot_tof_chunk`
### **Description**
Compiles TOF data, splits it into chunks, and generates Gaussian fits for each chunk, saving the results to a PDF.

### **Parameters**
- **dTs (list)**: A list of time differences for events.
- **rpc_indicies (list of lists)**: A list of pairs of RPC indices for TOF analysis.
- **num_chunks (int, optional)**: Number of chunks to split the data into. Default is 10.
- **pdf_filename (str, optional)**: The filename for the output PDF. Default is "Data_output/tof_chunks.pdf".

### **Returns**
- pdf_filename (str): The filename of the generated PDF.
- **plots**: plots for the tof analysis




## Example Usage
To find coincident hits in eta and phi data:

```python
coincident_hits = FindCoincidentHits(etaHits, phiHits, time_window=5)
```

To cluster the coincident hits:

```python
clustered_hits = cluster(coincident_hits)
```

To reconstruct the event using a chi-squared fitting method:

```python
result = reconstruct_timed_Chi2_modified(clustered_hits[0], max_cluster_size=3)
if result:
    centroid, direction, coordinates, combinations, chi2, dT, dZ = result
```

To compile and plot TOF data:

```python
compile_and_plot_tof(reconstructor.dT,rpc_comparison, pdf_filename='Data_output/tof.pdf')
```

To writein to a function

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