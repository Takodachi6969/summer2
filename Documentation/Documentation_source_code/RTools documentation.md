# Reconstruction Tools Documentation

## RTools.`find_tdc_event_count`


### Description
calculates the total number of words for each TDC (Time-to-Digital Converter) across a chunk of events. It processes the events and sums the lengths of the `words` lists for each TDC.

### Parameters
- **event_chunk** (`list`): A list of proanubis event objects. Each event object contains an attribute `tdcEvents`, which is a list of TDC objects. Each TDC object has an attribute `words`, which is a list of words.

### Returns
- **event_number** (`list`): A list of lists. Each inner list corresponds to a TDC and contains a single integer representing the total count of words for that TDC across all events in the given chunk.

### Note
- Function assumes exactly **5** TDCs
- The returned `event_number` list contains 5 lists, each with a single integer value representing the total count of words for the corresponding TDC.


## Rtools.`FindCoincidentHits`

### Description
This function processes eta and phi hits from different RPC arrays to find coincident hits based on time windows and optional Time-of-Flight (TOF) correction.

### **Parameters**
- **etaHits (list)**: A list of lists containing eta hit objects for each RPC. Each hit object must have attributes time, channel, event_num, and eta.
- **phiHits (list)**: A list of lists containing phi hit objects for each RPC. Each hit object must have attributes time, channel, event_num, and eta.
- **time_window (float)**: The time window in which to consider hits as coincident.
- **tof_correction (bool, optional)**: Whether to apply TOF correction. Default is True.
- **slope (float, optional)**: Slope parameter for TOF correction. Default is 0.05426554612593516.
- **offSet (float, optional)**: Offset parameter for TOF correction. Default is 15.8797407836404.

### **Returns**
- **coincident_hits (list)**: A list of lists where each sublist represents an event with coincident hits. Each sublist contains the event number, time bin, and details of the hits.
``` python
[event_num, time_bin,[[hit.rpc, hit.channel, hit.time, hit.eta] for hit in unique_hits]]
```

---

## Rtools.`cluster`
### **Description**
This function clusters the hits within each RPC into spatial clusters in both phi and eta directions.

### **Parameters**
- **coincident_hits (list)**: A list of coincident hit events. Each event is a list containing the event number, time bin, and hit details.

### **Returns**
- **coincident_hits_clustered (list)**: A list of events with clustered hits. Each event contains the event number, time bin, and clusters of hits for each RPC.
```python
[event_num, time_bin, [rpc_phi_clusters, rpc_eta_clusters]]
```

---

## Rtools.`extract_coords_timed_Chi2`

### **Description**
Converts spatial clusters in RPCs into x and y coordinates, with z given by the RPC number.

### **Parameters**
- **event (list)**: A list representing an event. It contains the event number, time bin, and clusters of hits.
- **max_cluster_size (int)**: Maximum allowed size for a cluster.
  
### **Returns**
- **coords (list)**: A list of coordinates with errors. Each entry represents an RPC and contains lists of x and y coordinates, with variances and times.
---


## Rtools.`generate_hit_coords_combo_Chi2`
### **Description**
Generates combinations of hit coordinates across RPCs for chi-squared fitting.

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

---

## Rtools.`fit_event_chi2`
### **Description**
Fits a line to the event coordinates using Singular Value Decomposition (SVD) and calculates the chi-squared value.

### **Parameters**
- **coordinates_with_error (list)**: A list of coordinates with errors for each RPC.
- **rpc_indicies (list of lists, optional)**: A list of pairs of indices specifying the RPCs to compare.

## **Returns**
- **p0 (list)**: Centroid of the fitted line.
- **d (list)**: Direction vector of the fitted line.
- **Chi2 (float)**: Chi-squared value of the fit.
- **coordinates (list)**: List of coordinates used in the fit.
- **dT (list)**: List of time differences.
- **dZ (list)**: List of height differences.

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

## **check_event_attributes**
### **Description**
Checks if an event meets the minimum required number of RPCs and chambers with hits.

### **Parameters**
- **event (list)**: A list representing an event.
- **min_chamber_number (int)**: Minimum required number of chambers with hits.
- **min_RPC_number (int)**: Minimum required number of RPCs with hits.

### **Returns**
- **valid_event (bool)**: **True** if the event meets the requirements, otherwise **False**.

---

## **filter_events**
### **Description**
Filters a list of events based on the minimum required number of RPCs and chambers with hits.

### **Parameters**
- **events (list)**: A list of events.
- **min_chamber_number (int)**: Minimum required number of chambers with hits.
- **min_RPC_number (int)**: Minimum required number of RPCs with hits.

### **Returns**
- **filtered_events (list)**: A list of events that meet the requirements.

---

## **count_entries**
### ***Description***
Recursively counts the total number of entries in a nested list.

### **Parameters**
- **lst (list)**: A nested list.

### **Returns**
- **total_count (int)**: Total number of entries in the list.

---

## **find_tdc_event_count**
### **Description**
Calculates the total number of words for each TDC (Time-to-Digital Converter) across a chunk of events.

### **Parameters**
- **event_chunk (list)**: A list of event objects. Each event object contains a list of TDC objects with a words attribute.

### **Returns**
- **event_number (list)**: A list of lists. Each inner list corresponds to a TDC and contains a single integer representing the total count of words for that TDC across all events in the chunk.

## **compile_and_plot_tof**
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


## **compile_and_plot_tof_chunk**
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