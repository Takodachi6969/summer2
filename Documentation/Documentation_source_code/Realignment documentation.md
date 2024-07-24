# Realignent Documentation
# Overview
This is the documentation of the realignment

# Table of content

- [Realignent Documentation](#realignent-documentation)
- [Overview](#overview)
- [Table of content](#table-of-content)
    - [fReader.`get_aligned_events`](#freaderget_aligned_events)
      - [Description:](#description)
      - [Parameters:](#parameters)
      - [Returns:](#returns)
      - [Usage Example:](#usage-example)
    - [self.`doRealign`](#selfdorealign)
      - [Description:](#description-1)
      - [Parameters:](#parameters-1)
      - [Returns:](#returns-1)
    - [self.`check_alignment_status`](#selfcheck_alignment_status)
      - [Description:](#description-2)
      - [Parameters:](#parameters-2)
    - [self.`update_adjustment_window`](#selfupdate_adjustment_window)
      - [Description:](#description-3)
      - [Parameters:](#parameters-3)
    - [self.`monitor_tdc3_state`](#selfmonitor_tdc3_state)
      - [Description:](#description-4)
      - [Parameters:](#parameters-4)
      - [Returns:](#returns-2)
    - [self.`InsertFakeEvents`](#selfinsertfakeevents)
      - [Description:](#description-5)
      - [Parameters:](#parameters-5)
    - [Code Description: `get_aligned_events`](#code-description-get_aligned_events)
      - [1. Initialization](#1-initialization)
      - [2. Reading and Processing Events](#2-reading-and-processing-events)
      - [3. Aligning Events](#3-aligning-events)
      - [4. Updating Alignment Status](#4-updating-alignment-status)
      - [5. Updating Adjustment Window](#5-updating-adjustment-window)
      - [6. Monitoring TDC State](#6-monitoring-tdc-state)
      - [7. Returning Results](#7-returning-results)

### fReader.`get_aligned_events`
**fReader.`get_aligned_events`(order = [(0,1), (1,2), (2,3), (3,4)], interval = 100, extract_tdc_mets = False, recordtimes = False)**

#### Description:
output chunks of interval length events

#### Parameters:

1. **order**: 
   - **Type**: `list of Tuples`
   - **Description**: Specifies your order.
   - **Default Value**: `[(0,1), (1,2), (2,3), (3,4)]`

2. **interval**: 
   - **Type**: `int`
   - **Description**: The number of events to process in one chunk, larger chunk leads to better alignment.
   - **Default Value**: `100`

3. **extract_tdc_mets**: 
   - **Type**: `bool`
   - **Description**: Flag that determines whether to extract TDC metrics, where the TDC states will be monitored.
   - **Default Value**: `False`
4. **recordtimes**: 
   - **Type**: `bool`
   - **Description**: Flag that determines whether to record the minimum time and word for each event
   - **Default Value**: `False`

#### Returns:

1. **When `extract_tdc_mets` is `True`**:
    - **evts_chunk**: 
      - **Type**: `list`
      - **Description**: A list of processed events that are not necessarily globally aligned.
    - **tdc_mets**: 
      - **Type**: `list of lists`
      - **Description**: A list of TDC metrics for each tdc, only outputted every 2500 events. .
      - **Example Value**: `[[TDC0 list], [TDC1 list], [TDC2 list], [TDC3 lsit], [TDC4 list]]`
    - **TDC_error_time**: 
      - **Type**: `list of lists`
      - **Description**: A list for eah TDC recording the minimum event time and event word in each event
      - **Example Value**: `[[(min_time, min_word), event number], [TDC1 list], [TDC2 list], [TDC3 lsit], [TDC4 list]]`

2. **When `extract_tdc_mets` is `False` and `self.global_alignment` is `True`**:
    - **evts_chunk**: 
      - **Type**: `list`
      - **Description**: A list of processed events.

3. **When `extract_tdc_mets` is `False` and `self.global_alignment` is `False`**:
    - **None**: 
      - **Description**: The function returns `None`.

#### Usage Example:

```python
import importlib
importlib.reload(rawFileReader) # Reload fReader

# Set your monitoring chunk size
interval = 100 

# Reload the class object
fReader = rawFileReader.fileReader(file_path)

# Define the alignment order
order = [(0,1), (1,2), (2,3), (3,4)]

# Set the maximum number of events to process
max_process_event_chunk = 100 

# Initialize the counter for processed events
processedEvents = 0 

# Process events using a progress bar

while processedEvents < max_process_event_chunk:
   processedEvents += 1
   event_chunk = fReader.get_aligned_events(order=order, interval=interval)
```

 see also [AnalysisTools](ATools documentation.html) for other complimentary functions used

---

### self.`doRealign`

#### Description:
Aligns the events based on the given order and updates the alignment metrics.

#### Parameters:
- **event_chunk**: 
  - **Type**: `list`
  - **Description**: The chunk of events to be realigned.

- **order**: 
  - **Type**: `list of Tuples`
  - **Description**: Specifies the order for alignment.

- **skipChans**: 
  - **Type**: `list`
  - **Description**: Channels to be skipped during alignment.
  - **Default Value**: `[0]`

#### Returns:
- **aligned**: 
  - **Type**: `bool`
  - **Description**: Whether the events are aligned.

- **realigned**: 
  - **Type**: `bool`
  - **Description**: Whether the events were realigned.

---

### self.`check_alignment_status`

#### Description:
Checks the status of the alignment and updates the global alignment and last bad status.

#### Parameters:
- **aligned**: 
  - **Type**: `bool`
  - **Description**: Whether the events are aligned.

- **realigned**: 
  - **Type**: `bool`
  - **Description**: Whether the events were realigned.

---

### self.`update_adjustment_window`

#### Description:
Updates the adjustment window based on the alignment status.

#### Parameters:
- **realigned**: 
  - **Type**: `bool`
  - **Description**: Whether the events were realigned.

---

### self.`monitor_tdc3_state`

#### Description:
Monitors the state of the TDCs and updates TDC metrics.

#### Parameters:
- **recordtimes**: 
  - **Type**: `bool`
  - **Description**: Flag that determines whether to record the minimum time and word for each event.
  - **Default Value**: `False`

#### Returns:
- **TDC_error_time**: 
  - **Type**: `list of lists`
  - **Description**: A list for each TDC recording the minimum event time and event word.

- **tdc_mets**: 
  - **Type**: `list of lists`
  - **Description**: A list of TDC metrics for each TDC.

---

### self.`InsertFakeEvents`

#### Description:
Inserts fake events into the event builder based on the insertion list.

#### Parameters:
- **insertion_list**: 
  - **Type**: `list`
  - **Description**: List of the number of fake events to insert for each TDC.



### Code Description: `get_aligned_events`

The `get_aligned_events` function is designed to process a chunk of events, align them according to the alignment metric, monitor TDC (Time-to-Digital Converter) states, and optionally extract TDC metrics.

#### 1. Initialization
- The function initializes three main variables: 
  - `evts_chunk` to store the processed events.
  - `tdc_mets` to store TDC metrics (initialized to zero).
  - `TDC_error_time` to record error times for each TDC (initialized as empty lists).
- A counter `i` is set to 0 to keep track of the number of events processed.

#### 2. Reading and Processing Events
- The function does the following for each pair ordered
- The function enters a while loop, which runs until `i` reaches the specified `interval` value.
- Within the loop, `self.readBlock()` is called to read a block of data. If the read fails, the loop breaks.
- If the block contains events (`self.hasEvents()`), each event is appended to `evts_chunk`, and counters `i` and `self.tdc_monitoring_counter` are incremented.

#### 3. Aligning Events
- The function calls `self.doRealign(evts_chunk, order)` to align the events in `evts_chunk` according to the specified `order`.
- `doRealign` performs the alignment by iterating through the specified `order`, calculating alignment metrics using `ATools.calcAvgAlign`, and attempting to realign events if the metrics indicate misalignment.
- `ATools.calcAvgAlign` finds the alignment metric by choosing two different RPCs operated under different TDC. Then the minimum hit time hits are used in each event, the average separation in eta and phi channels are calculated for a length of interval events. If the events are aligned, the averaged distance would be small because they represent real particles, otherwise the average distance is large since it is purely stochastic. 
- If realignment is found, the pairwise difference is recorded, where Atools.`ConstructEventInsertionList` is used to calculate how many fakes events needed to be inserted. fake events are inserted using `self.InsertFakeEvents` where the variable is processed internally. 


#### 4. Updating Alignment Status
- After alignment, the function calls `self.check_alignment_status(aligned, realigned)` to update the global alignment status and the status of the last processed events.

#### 5. Updating Adjustment Window
- The function then calls `self.update_adjustment_window(realigned)` to adjust the window size used for alignment based on the recent alignment status. If the events are not aligned and no alignment were found, the window increases

#### 6. Monitoring TDC State
- The function extends `self.tdc_monitoring_event_buffer` with the newly processed events.
- If `self.tdc_monitoring_counter` exceeds 2500, the function calls `self.monitor_tdc3_state(recordtimes=True)` to monitor the state of the TDCs, update TDC metrics, and reset the counter and buffer.

#### 7. Returning Results
- If `extract_tdc_mets` is `True`, the function returns `evts_chunk`, `tdc_mets`, and `TDC_error_time`.
- If `extract_tdc_mets` is `False` and `self.global_alignment` is `True`, the function returns `evts_chunk`.
- Otherwise, the function returns `None`.

---