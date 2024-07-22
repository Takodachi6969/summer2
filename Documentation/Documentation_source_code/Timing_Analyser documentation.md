# `Timing_Analyser` Class Documentation
## Overview
The `Timing_Analyser` class is designed to process and analyze timing data from event chunks. This class helps in calculating and visualizing time of flight analysis, residuals, and other related metrics for the pro_anubis detectors. It provides functionalities to update events, read TDC (Time-to-Digital Converter) time differences, calculate residuals, check eta trigger, and plot various data for analysis.

## Class Definition

### Initialization
```python
def __init__(self, event_chunk, processsed_event, diffHists=None, scDiffs=None, normDiffs=None):
```
- ### **Parameters**
- `event_chunk`: List of events to be processed.
- `processsed_event`: processed_event_number.
- `diffHists`: (Optional) Dictionary of histograms for storing time differences.
- `scDiffs`: (Optional) List for storing scaled differences.
- `normDiffs`: (Optional) List for storing normalized differences.
- ### **Attributes**:
- `totDiffs`: Dictionary storing total differences for each RPC.
- `nDiffs`: Dictionary storing the number of differences for each RPC.
- `diffHists`: Dictionary of histograms for storing time differences, initialized if not provided.
- `scDiffs`: List for storing scaled differences, initialized if not provided.
- `normDiffs`: List for storing normalized differences, initialized if not provided.
- `residEtaLatest`: List for storing latest eta residuals.
- `residPhiLatest`: List for storing latest phi residuals.
- `count`: List for storing event counts.
- `rpc_involvement`: Dictionary for tracking RPC involvement in events.


## **Methods**
### TAnalyser.`update_event`
```python
def update_event(self, event_chunk, processed_event):
```
### **Description**
Updates the event data to a new chunk, without changing other stored instances

### **Parameters**:
- `event_chunk`: List of new events to be processed.
- `processed_event`: List of new processed events.
---

### TAnalyser.`readTDCTimeDiffs`
### **Description**
This method processes the event chunk to calculate the time differences between eta and phi hits for each RPC's stirips and updates the respective histograms and difference lists.

### **Related_functions**

[RPC_hit](ATools documentation.html#class-rpchit)    
[TTools.`tdcEventToRPCData`](Ttools documentation.html#tdceventtorpcdata)

---

### TAnalyser.`Calculate_Residual_and_plot_TDC_Time_Diffs`
```python
def Calculate_Residual_and_plot_TDC_Time_Diffs(self, outDict, pdf_filename='plots.pdf', max_itr=1):
```
### **Description**
This function calculates residuals of each eta and phi crossing after time walk correction and generates the intra time of flight plots for each RPC, before and after applying the correction. The plots are saved in a specified PDF file. The function accounts for bad channels in the dataset and iteratively refines the residual calculations.

### **Parameters**
- **self**: An instance of the class containing the function. This instance should have the following attributes:
- **scDiffs**: A 2D array to store scalar differences for each channel pairs.
- **normDiffs**: A 2D array to store normalized differences for each channel.
- **residEtaLatest**: A list to store the latest eta residuals.
- **residPhiLatest**: A list to store the latest phi residuals.
- **outDict (dict)**: A dictionary containing the output data. Expected structure:
- **diffHists**: A 3D list or array containing histograms for each RPC, phi, and eta channel.
- **pdf_filename (str, optional)**: The name of the PDF file to save the plots. Default is 'plots.pdf'.
- **max_itr (int, optional)**: The maximum number of iterations for residual calculations. Default is 1.
- 
### **Returns**
- **residEtaLatest (list)**: The list of the latest eta residuals.
- **residPhiLatest (list)**: The list of the latest phi residuals.

### **Further information**
The residual is calculated by averaging the eta and phi offset along one strip direction, and projected into a 2D array. When iteration equals 1, eta were projected using phi, and any more iteration reprojects that to phi, eta and so on

### **Note**
The time walk correction is defined as 
```python
correction = slope * (eta - phi) + offset
```
where the slope and offset were defined using "magic number". Although they should not change, if more events used, a better fit could be found for these corrections, and might require adjustment in the future

---

### TAnalyser.`check_eta_trigger`
```python
def check_eta_trigger(self):
```
### **Description**
The `check_eta_trigger` function is designed to process a chunk of event data and determine if any events meet specific trigger criteria based on the number of hits recorded in RPCs per event. It identifies events where hits occur in at least four different RPCs and classifies them accordingly. The function also identifies and reports events that fail to meet the minimum criteria. The data of all triggers are also recorded through self.`count` and self.`rpc_involvement`

### **Parameters**
- **self**: An instance of the class that contains this method. The class is expected to have the following attributes:
- **event_chunk**: A list of event objects to be processed.
processedEvents: A counter for the number of processed events.
- **count**: A dictionary that stores lists of event numbers, keyed by the number of RPCs involved.
- **rpc_involvement**: A dictionary that tracks the number of times each RPC is involved in events, keyed by the number of RPCs involved.

### **Returns**
- **tuple**: A tuple containing:
- **bool**: Indicates whether all events met the criteria (i.e., True if all events had hits in at least 4 RPCs, False otherwise).
- **list**: A list of tuples, each containing an event number and the corresponding count of RPCs involved, for events that did not meet the criteria (i.e., had hits in fewer than 4 RPCs).

### **Example Usage**
```python
result, failed_events = self.check_eta_trigger()

if result:
    print("All events met the criteria.")
else:
    print("Some events failed the criteria:")
    for event_num, count in failed_events:
        print(f"Event {event_num} involved {count} RPC(s).")

```

### **Note**
Only good quality data are tested in the quality. This can be disabled manually or another check condition can be added. 


### **Related function**
[ATools.`tdcChanToRPCHit`](ATools documentation#atoolsfind_tdc_alignment_metric)

---

### TAnalyser.`plot_rpc_involvement_histogram`

### Description
The `plot_rpc_involvement_histogram` function generates a histogram that visualizes the normalized involvement of each RPC in events, categorized by the number of RPCs involved in each event. The histogram helps in understanding the distribution and proportion of events involving each RPC.

### Parameters
- **self**: An instance of the class that contains this method. The class is expected to have an attribute `rpc_involvement` which is a dictionary that tracks the number of times each RPC is involved in events, keyed by the number of RPCs involved.

### Example Usage
```python
self.plot_rpc_involvement_histogram()

```