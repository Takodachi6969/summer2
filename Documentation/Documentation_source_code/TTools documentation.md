# TTools Functions

This documentation provides an overview of various functions used for TDC time of flight analysis. These functions help in converting TDC events to RPC data, plotting channels, and alignment times, and displaying strip time information. The functions are designed for handling and visualizing data specific to proAnubis systems.

## Table of Contents
- [tdcEventToRPCData](#tdceventtorpcdata)
- [rpcHitToTdcChan](#rpchittotdcchan)
- [plot_tdc_error_times](#plot_tdc_error_times)
- [plot_tdc_error_times_custom_ranges](#plot_tdc_error_times_custom_ranges)
- [plot_tdc_error_channels](#plot_tdc_error_channels)
- [plot_tdc_error_channels_custom_ranges](#plot_tdc_error_channels_custom_ranges)
- [plot_tdc_alignment_channels_custom_ranges](#plot_tdc_alignment_channels_custom_ranges)
- [plot_tdc_alignment_times_custom_ranges](#plot_tdc_alignment_times_custom_ranges)
- [show_strip_time_info](#show_strip_time_info)

## TTools.`tdcEventToRPCData`

### **Description**
This function maps RPC hits to TDC channels.

### **Parameters**
- **rpc (int)**: RPC index.
- **rpcChan (int)**: RPC channel.
- **eta (bool)**: Eta parameter to distinguish between channels.
- 
### **Returns**
- **tdc (int)** and **tdcChannel (int)**: Corresponding TDC index and channel, or None if invalid combination.

### **Example usage**
```python
tdcEvenToRPCData(0, 12, True)
```
---

## TTools.`plot_tdc_error_times`

### **Description**
This function plots histograms of min_time for all TDCs.

### **Parameters**
**TDC_error_time (list of lists)**: TDC error times. Each sublist is the data for tdc 0 to 5 respectively

### **Usage Example**

```python
[(event_time, event_word), event_num] = [(253, 605028605), 1]
TDC_error_time = [[[(event_time, event_word), event_num]] for _ in range(5)]
plot_tdc_error_times(TDC_error_time)
```

### **Note**

This is horrible Naming, should be the min_time 

---

## TTools.`plot_tdc_error_times_custom_ranges`
### **Description**
This function plots histograms of min_time for different TDCs over custom process event ranges and saves the plots to a PDF.

### **Parameters**
- **TDC_error_time (list of list of tuples)**: Data for TDC error times.
- **ranges (list of tuples)**: Custom process time ranges.
- **output_pdf (str, optional)**: Output PDF file name. Default is TDC_error_times.pdf.

### **Usage Example**

```python
plot_tdc_error_times_custom_ranges(TDC_error_time, [(0, 15000), (25000, 40000)])
```

### **Note**

This is horrible Naming, should be the min_time 

---

## TTools.`plot_tdc_error_channels`
### **Description**
This function plots histograms of good and bad time channels for selected TDCs.

### **Parameters**
- **TDC_error_time (list of lists)**: TDC error times.
- **tdcs_to_plot (list of int, optional)**: TDC indices to plot. Default is [0, 1, 2, 3, 4].
  
### **Usage Example**

```python
plot_tdc_error_channels(TDC_error_time, [0, 1, 2])
```
### **Note**

This is horrible Naming, should be the min_time 

---

## TTools.`plot_tdc_error_channels_custom_ranges`
### **Description**
This function plots histograms of good and bad time channels for different TDCs over custom process time ranges and saves the plots to a PDF.

### **Parameters**
- **TDC_error_time (list of list of tuples)**: Data for TDC error times.
- **ranges (list of tuples)**: Custom process time ranges.
- **tdcs_to_plot (list of int, optional)**: TDC indices to plot. Default is plotting all TDCs (0-4).
- **output_pdf (str, optional)**: Output PDF file name. Default is TDC_error_channels.pdf.
- 
### **Usage Example**

```python
plot_tdc_error_channels_custom_ranges(TDC_error_time, [(0, 15000), (25000, 40000)], [0, 1])
```

### **Note**

This is horrible Naming, should be the min_time 

---

## TTools.`plot_tdc_alignment_channels_custom_ranges`
### **Description**
This function plots histograms of alignment channels for different TDCs over custom process time ranges and saves the plots to a PDF.

### **Parameters**
- **TDC_alignment_time (list of list of tuples)**: Data for TDC alignment times.
- **ranges (list of tuples)**: Custom process time ranges.
- **tdcs_to_plot (list of int, optional)**: TDC indices to plot. Default is plotting all TDCs (0-4).
- **output_pdf (str, optional)**: Output PDF file name. Default is TDC_alignment_channels.pdf.

### **Usage Example**
```python
plot_tdc_alignment_channels_custom_ranges(TDC_alignment_time, [(0, 15000), (25000, 40000)], [0, 1])

```

---

## TTools.`plot_tdc_alignment_times_custom_ranges`
### **Description**
This function plots histograms of min_time for different TDCs over custom process time ranges and saves the plots to a PDF.

### **Parameters**
- **TDC_alignment_time (list of list of tuples)**: Data for TDC alignment times.
- **ranges (list of tuples)**: Custom process time ranges.
- **output_pdf (str, optional)**: Output PDF file name. Default is TDC_alignment_times.pdf.
- 
### **Usage Example**
```python
plot_tdc_alignment_times_custom_ranges(TDC_alignment_time, [(0, 15000), (25000, 40000)])
```

---

## TTools.`show_strip_time_info`
### **Description**
This function displays the strip time information by plotting the difference between eta and phi hit times. 

### **Parameters**
- **outDict (dict)**: Output dictionary containing difference histograms, this is populated by TAnalyser.
- **ph (int)**: Phi index.
- **et (int)**: Eta index.
- **rpc (int)**: RPC index.

### **Usage Example**
```python
outDict = {'totDiffs':TAnalyser.totDiffs,
                    'nDiffs':TAnalyser.nDiffs,
                    'diffHists':TAnalyser.diffHists} 
show_strip_time_info(outDict, 12, 22, 0)

```
### **Related Functions**
[TAnlayser](Timing_Analyser documentation.html)

