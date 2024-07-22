# Analysis Tools Documentation

## Table of Contents

- [Analysis Tools Documentation](#analysis-tools-documentation)
  - [Table of Contents](#table-of-contents)
    - [Class `rpcHit`](#class-rpchit)
    - [ATools.`calcAvgAlign`](#atoolscalcavgalign)
      - [Description:](#description)
      - [Parameters:](#parameters)
      - [Returns:](#returns)
    - [ATools.`find_tdc_alignment_metric`](#atoolsfind_tdc_alignment_metric)
      - [Description:](#description-1)
      - [Parameters:](#parameters-1)
      - [Returns:](#returns-1)
    - [ATools.`testAlign`](#atoolstestalign)
      - [Description:](#description-2)
      - [Parameters:](#parameters-2)
      - [Returns:](#returns-2)
    - [ATools.`ConstructEventInsertionList`](#atoolsconstructeventinsertionlist)
      - [Description:](#description-3)
      - [Parameters:](#parameters-3)
      - [Returns:](#returns-3)


### Class `rpcHit`

The `rpcHit` class is designed to represent a hit in a RPC. It stores the channel, time, eta, event number, and RPC identifier.

```python
class rpcHit():
    def __init__(self, channel, time, eta, event_num, rpc):
        self.rpc = rpc
        self.time = time
        self.channel = channel
        self.eta = eta
        self.event_num = event_num

    def __str__(self):
        return f"rpcHit(channel={self.channel}, time={self.time}, eta={self.eta}, event_num={self.event_num}, rpc={self.rpc})"
```
### ATools.`calcAvgAlign`

#### Description:
Calculates the average alignment metric for a given chunk of events by evaluating the average distance between hits across 2 RPCs operated by different TDC.

#### Parameters:
- **event_chunk**: 
  - **Type**: `list`
  - **Description**: The chunk of events to be aligned.

- **offSet**: 
  - **Type**: `int`
  - **Description**: Offset used for alignment calculation.
  - **Default Value**: `0`

- **i**: 
  - **Type**: `int`
  - **Description**: Index for the first eta hit.
  - **Default Value**: `1`

- **j**: 
  - **Type**: `int`
  - **Description**: Index for the second eta hit.
  - **Default Value**: `2`

- **k**: 
  - **Type**: `int`
  - **Description**: Index for the first phi hit.
  - **Default Value**: `0`

- **l**: 
  - **Type**: `int`
  - **Description**: Index for the second phi hit.
  - **Default Value**: `2`

- **tdc1**: 
  - **Type**: `int`
  - **Description**: TDC index for the first set of hits.
  - **Default Value**: `0`

- **tdc0**: 
  - **Type**: `int`
  - **Description**: TDC index for the second set of hits.
  - **Default Value**: `1`

- **processedEvents**: 
  - **Type**: `int`
  - **Description**: Number of processed events.
  - **Default Value**: `0`

- **skipChans**: 
  - **Type**: `list`
  - **Description**: Channels to be skipped during alignment calculation.
  - **Default Value**: `[]`

#### Returns:
- **alignment_metric**: 
  - **Type**: `float`
  - **Description**: The average alignment metric for the event chunk. If no valid metrics are found, returns `-1`.

---

### ATools.`find_tdc_alignment_metric`

#### Description:
Determines the alignment metric indices for two given TDCs, to ensure at least 2 RPCs are crossed with the same TDC.

#### Parameters:
- **tdc0**: 
  - **Type**: `int`
  - **Description**: Index of the first TDC.

- **tdc1**: 
  - **Type**: `int`
  - **Description**: Index of the second TDC.

#### Returns:
- **alignment_indices**: 
  - **Type**: `Tuple[int, int, int, int]`
  - **Description**: Indices for eta and phi hits used in alignment. Returns `i, j, k, l` where `i` and `j` are for eta hits, and `k` and `l` are for phi hits.

---

### ATools.`testAlign`

#### Description:
Tests the alignment of hits between two RPCs by calculating the distance between two channels in eta and phi direction using minimum time.

#### Parameters:
- **rpc1Hits**: 
  - **Type**: `list`
  - **Description**: List of hits from the first RPC.

- **rpc2Hits**: 
  - **Type**: `list`
  - **Description**: List of hits from the second RPC.

- **skipChans**: 
  - **Type**: `list`
  - **Description**: Channels to be skipped during alignment calculation.
  - **Default Value**: `[]`

#### Returns:
- **time_difference**: 
  - **Type**: `int`
  - **Description**: The absolute difference in minimum channels between the two RPCs. Returns `-1` if alignment cannot be determined.

---

### ATools.`ConstructEventInsertionList`

#### Description:
Constructs a list indicating the number of fake events to be inserted for each TDC based on alignment updates.

#### Parameters:
- **updates**: 
  - **Type**: `list of int`
  - **Description**: List of update values for alignment, where each value corresponds to the adjustment needed for the respective TDC.

- **order**: 
  - **Type**: `list of Tuples`
  - **Description**: Specifies the order of TDC alignment pairs. Each tuple contains two TDC indices that are being aligned.

#### Returns:
- **insertion_list**: 
  - **Type**: `list of int`
  - **Description**: A list indicating the number of fake events to be inserted for each TDC. Ensures no negative values by adjusting the list with the minimum value if necessary.
---