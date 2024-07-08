from PIL import Image
# import anubisPlotUtils as anPlot
import json
import numpy as np
import os
import matplotlib
import matplotlib.pyplot as plt
import matplotlib.colors as colors
matplotlib.use('TkAgg')  # or 'Qt5Agg', 'GTK3Agg', etc.
import sys
# import ANUBIS_triggered_functions as ANT
import pandas as pd
import matplotlib.backends.backend_pdf
from matplotlib.ticker import MultipleLocator
from functools import reduce
import matplotlib.pyplot as plt
import numpy as np
# from scipy.stats import normpip install pillow
sys.path.insert(1, 'Osiris Temp\processing\python')
import rawFileReader
import csv
from itertools import groupby
import math
from scipy.optimize import curve_fit

import importlib
# %matplotlib inline
from collections import defaultdict
import cProfile
import pstats
import io


class rpcHit():
    def __init__(self, channel, time, eta, event_num, rpc):
        self.rpc = rpc
        self.time = time
        self.channel = channel
        self.eta = eta
        self.event_num = event_num


def tdcChanToRPCHit(word, tdc, event_num):
        tdcChannel = (word >> 24) & 0x7f
        tdcHitTime = word & 0xfffff
        eta = False
        rpcChan = -1
        if tdc == 0:
            if tdcChannel < 32:
                rpcChan = tdcChannel
                eta = True
                rpc = 0
            elif tdcChannel < 96:
                rpcChan = tdcChannel - 32
                eta = False
                rpc = 0
            else:
                rpcChan = tdcChannel - 96
                eta = True
                rpc = 1
        elif tdc == 1:
            if tdcChannel < 64:
                rpcChan = tdcChannel
                eta = False
                rpc = 1
            elif tdcChannel < 96:
                rpcChan = tdcChannel - 64
                eta = True
                rpc = 2
            else:
                rpcChan = tdcChannel - 96
                eta = False
                rpc = 2
        elif tdc == 2:
            if tdcChannel < 32:
                rpcChan = tdcChannel + 32
                eta = False
                rpc = 2
            elif tdcChannel < 64:
                rpcChan = tdcChannel - 32
                eta = True
                rpc = 3
            elif tdcChannel < 128:
                rpcChan = tdcChannel - 64
                eta = False
                rpc = 3
        elif tdc == 3:
            if tdcChannel < 32:
                rpcChan = tdcChannel
                eta = True
                rpc = 4
            elif tdcChannel < 96:
                rpcChan = tdcChannel - 32
                eta = False
                rpc = 4
            else:
                rpcChan = tdcChannel - 96
                eta = True
                rpc = 5
        elif tdc == 4:
            rpcChan = tdcChannel
            eta = False
            rpc = 5
        return rpc, rpcHit(rpcChan, tdcHitTime * 0.8, eta, event_num, rpc)
def Unpack_event(eventList, processedEvents):
    allEtaHits = [[] for rpc in range(6)]
    allPhiHits = [[] for rpc in range(6)]
    for event_num, event in enumerate(eventList):
        etaHits = [[] for rpc in range(6)]
        phiHits = [[] for rpc in range(6)]
        for i in range(5):
            for word in event.tdcEvents[i].words:
                rpc, thisHit = tdcChanToRPCHit(word, i, processedEvents - len(eventList) + event_num + 1)
                if thisHit.eta:
                    etaHits[thisHit.rpc].append(thisHit)
                else:
                    phiHits[thisHit.rpc].append(thisHit)
        for rpc in range(6):
            allEtaHits[rpc].extend(etaHits[rpc])
            allPhiHits[rpc].extend(phiHits[rpc])
    return allEtaHits, allPhiHits
def FindCoincidentHits(etaHits, phiHits, time_window, tof_correction = False):
    channels = []
    
    # Combine etaHits and phiHits from all RPCs
    for RPC in range(6):
        channels += [hit for hit in etaHits[RPC] if 150 <= hit.time <= 300]
        channels += [hit for hit in phiHits[RPC] if 150 <= hit.time <= 300]

    # Sort events by event number and then by time within each event
    event_sorted = sorted(channels, key=lambda rpcHit: (rpcHit.event_num, rpcHit.time))
    
    # Group by event number
    grouped_and_sorted = {key: list(group) 
                          for key, group in groupby(event_sorted, lambda rpcHit: rpcHit.event_num)}
    
    coincident_hits = []

    for event_num, hits in grouped_and_sorted.items():
        temp_hits = []

        for i in range(len(hits) - 1):
            if tof_correction:
                if hits[i].eta:
                    correction = find_tof_time(hits[i], hits[i+1])
                else:
                    correction = find_tof_time(hits[i+1], hits[i])
            else:
                correction = 0
            
            if abs(hits[i+1].time - hits[i].time + correction) <= time_window:
                temp_hits.append(hits[i])
                temp_hits.append(hits[i+1])

        if temp_hits:
            # Remove duplicates and sort by time
            unique_hits = { (hit.channel, hit.time, hit.eta, hit.event_num, hit.rpc): hit for hit in temp_hits }.values()
            time_bin = min(hit.time for hit in unique_hits)
            
            coincident_hits.append([
                event_num,
                time_bin,
                [[hit.rpc, hit.channel, hit.time, hit.eta] for hit in unique_hits]
            ])


    return coincident_hits
def cluster(coincident_hits):
    coincident_hits_clustered = []

    for coincidence_event in coincident_hits:

        coincident_event_clustered = [coincidence_event[0], coincidence_event[1], []]

        hit_locations = coincidence_event[2]
        phi_locations = [x for x in hit_locations if x[3] == False]
        eta_locations = [x for x in hit_locations if x[3] == True]

        phi_locations = sorted(phi_locations, key=lambda x: x[1])
        eta_locations = sorted(eta_locations, key=lambda x: x[1])

        for RPC in range(6):
            rpc_phi_clusters = []
            rpc_eta_clusters = []

            i = 0
            for index, hit in enumerate([x for x in phi_locations if x[0] == RPC]):
                if index == 0:
                    previous_element = hit[1]
                    rpc_phi_clusters.append([hit])
                else:
                    if abs(hit[1] - previous_element) > 1:
                        rpc_phi_clusters.append([hit])
                        i += 1
                    else:
                        rpc_phi_clusters[i].append(hit)
                    previous_element = hit[1]

            j = 0
            for index, hit in enumerate([x for x in eta_locations if x[0] == RPC]):
                if index == 0:
                    previous_element = hit[1]
                    rpc_eta_clusters.append([hit])
                else:
                    if abs(hit[1] - previous_element) > 1:
                        rpc_eta_clusters.append([hit])
                        j += 1
                    else:
                        rpc_eta_clusters[j].append(hit)
                    previous_element = hit[1]

            rpc_combined = [rpc_phi_clusters, rpc_eta_clusters]

            coincident_event_clustered[2].append(rpc_combined)

        coincident_hits_clustered.append(coincident_event_clustered)

    return coincident_hits_clustered

def find_tof_time(eta, phi):
    slope = 0.173
    offSet = 13.2
    if (len(set([eta.eta, phi.eta])) == 1):
        return 0
    else:
        return slope*(phi.channel-eta.channel)-offSet

def check_event_attributes(event,min_chamber_number,min_RPC_number):
    #Used in filter_events() function to decide whether or not to save an event as described by the user's inputs.

    # event = ['Event x',TIMEBIN, [[[RPC1_PHI_CLUSTERS],[RPC1_ETA_CLUSTERS]],[[...],[...]],...]

    RPC_counter = 0
    chamber_counter = 0
    condition_1 = False
    condition_2 = False
    condition_3 = False

    for RPC in range(6):
        if RPC<3:
            #Checking triplet layer.
            if event[2][RPC][0] and event[2][RPC][1]:
                #Reqiure       phi^              and               eta^       strips to go off
                RPC_counter+=1 
                #If RPC has two eta and phi strips going off then consider it "hit"
                if not condition_1:
                    #Count triplet chamber being hit.
                    chamber_counter+=1
                    condition_1 = True
        elif RPC == 3:
            #Singlet layer
            if event[2][RPC][0] and event[2][RPC][1]:
                RPC_counter+=1
                if not condition_2:
                    chamber_counter+=1
                    condition_2 = True
        else:
            #Doublet layer
            if event[2][RPC][0] and event[2][RPC][1]:
                RPC_counter+=1
                if not condition_3:
                    chamber_counter+=1
                    condition_3 = True

    return RPC_counter >= min_RPC_number and chamber_counter >= min_chamber_number
def filter_events(events,min_chamber_number,min_RPC_number):
    #Initiliase array of filtered events
    filtered_events = []

    for event in events:
        if check_event_attributes(event,min_chamber_number,min_RPC_number):
            filtered_events.append(event)

    # print(f"Number of events in filter = {len(filtered_events)}")
    
    return filtered_events
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
def extract_DT_DZ_Chi2(coords):

    #coords = [[[x0,var,time],[y0,var],z0],[[x1,var,time],[y1,var],z1],...,[[x5,var,time],[y5,var],z5]]

    times = [[RPC,x[0][2]] for RPC, x in enumerate(coords) if isinstance(x[2], (float, int))]

    #Should already be sorted, but just in case.
    #Sort times by RPC, with RPC at lowest height at first entry.

    if len(times) > 1:

        times_sorted = sorted(times, key=lambda x: x[0])

        #print(times_sorted)

        dT = times_sorted[-1][1]-times_sorted[0][1]
        #if dT>0 this implies the particles hit the higher RPC after the lower one, so the particle is travelling upwards here.
        #Vice-versa for dT < 0 

        RPC_heights = [0.6,1.8,3.0,61.8,121.8,123] #Heights of middle point of each RPC, measured from the bottom of the Triplet Low RPC. Units are cm.

        first_RPC = times_sorted[0][0]
        last_RPC = times_sorted[-1][0]

        dZ = RPC_heights[last_RPC] - RPC_heights[first_RPC]
    
        return dT, dZ
    else:
        pass
def fit_event_chi2(coordinates_with_error):
    #Coordinates = [[[x0,var,time],[y0,var],z0],[[x1,var,time],[y1,var],z1],...,[[x5,var,time],[y5,var],z5]]
    #Z coordinate given by height of relevant RPC.
    #Using SVD

    # Calculate dT for event, in ns
    dT, dZ = extract_DT_DZ_Chi2(coordinates_with_error)
    
    coordinates = []

    for coords in coordinates_with_error:
        coordinates.append([coords[0][0],coords[1][0],coords[2]])

    centroid = np.mean(coordinates, axis=0)
    subtracted = coordinates-centroid

    # performing SVD
    _, _, V = np.linalg.svd(subtracted)
    
    # find the direction vector (which is the right singular vector corresponding to the largest singular value)
    direction = V[0, :]

    # A line is defined by the average and its direction
    p0 = centroid
    d = direction

    #Work out Chi2. Minimise this to find best fit (from possible combos)

    Chi2 = 0

    i = 0 

    for point in coordinates_with_error:
        
        i+=2
        
        z = point[2]
        x = point[0][0]
        y = point[1][0]
        x_var = point[0][1]
        y_var = point[1][1]

        z_0 = centroid[2]

        # t = (z-z_0)/d_z

        t = (z-z_0)/d[2]

        # Find expected (x,y) coordinates at that height.

        x_traj = centroid[0] + t*d[0]
        y_traj = centroid[1] + t*d[1]

        Chi2_x = (x-x_traj)**2 / x_var
        Chi2_y = (y-y_traj)**2 / y_var

        Chi2+= Chi2_x
        Chi2+= Chi2_y

    # i is number of fitted points. There are 4 fitted paramters, 2 for each x and y. 
    doF = i - 4

    Chi2 = Chi2/ doF

    return p0, d, Chi2, coordinates, dT, dZ
def extract_coords_timed_Chi2(event,max_cluster_size):

    #This function converts spatially clusters in RPCs into x and y coordinates (z given by RPC number)
    # event = ['Event x',TIMEBIN, [[[RPC1_PHI_CLUSTERS],[RPC1_ETA_CLUSTERS]],[[...],[...]],...]

    #Extract x and y coords of cluster in event
    distance_per_phi_channel = 2.7625 #cm
    distance_per_eta_channel = 2.9844 #cm
    
    coords = []

    for RPC in range(6):
        
        x_clusters = [x for x in event[2][RPC][0] if len(x)<=max_cluster_size] #phi direction
        y_clusters = [y for y in event[2][RPC][1] if len(y)<=max_cluster_size] #eta direction

        #Finding size of largest cluster, consider coordinates bad if largest cluster is larger than 6.
        x_clusters_lengths = [len(x) for x in event[2][RPC][0]]
        y_clusters_lengths = [len(y) for y in event[2][RPC][1]]

        max_length = max(max(x_clusters_lengths, default=0), max(y_clusters_lengths, default=0))

        x_coords = []
        y_coords = []
        
        for x_cluster in x_clusters:
            # Extract phi channels and times from the cluster
            phi_channels = [x[1] for x in x_cluster]
            phi_times = [t[2] for t in x_cluster]

            # Convert the channel number into a measurement along the RPC
            x_values = [(phi_channel + 0.5) * distance_per_phi_channel for phi_channel in phi_channels]

            # Variance in x coord
            x_var = (1 * distance_per_phi_channel) ** 2 / 12

            # Find the index of the minimum time
            min_time_index = phi_times.index(min(phi_times))

            # Append the x value corresponding to the minimum time
            x_coords.append([x_values[min_time_index], x_var, min(phi_times)])

        for y_cluster in y_clusters:
            #y_cluster = [[RPC,CHANNEL,TIME,'eta'],...]
            eta_channels_corrected = [31-y[1] for y in y_cluster] #corrected for labelling from 0 to 31.
            eta_times = [t[2] for t in y_cluster]
            y_values = [(channel_num+0.5)*distance_per_eta_channel for channel_num in eta_channels_corrected]
            
            y_var = (1*distance_per_eta_channel)**2 /12
            
            # Find the index of the minimum time
            min_time_index = eta_times.index(min(eta_times))
            
            y_coords.append([y_values[min_time_index],y_var,min(eta_times)])

        if x_coords and y_coords:

            coords.append([x_coords, y_coords])

        else:
            coords.append([[],[],"N"])

    #[x_coords] = [[x,err_x,x_time],...]
    
    #RPC_coords = [x_coords,y_coords]

    #coords = [[RPC1_coords],[RPC2_coords],[RPC3_coords],...]
    return coords
def reconstruct_timed_Chi2_ByRPC(event,max_cluster_size, RPC_excluded):

    #timed tag indicates that timing information from RPC is used to determine direction of vertical transversal of "particle" in the event.

    max_Chi2 = 10

    # event = ['Event x',TIMEBIN, [[[RPC1_PHI_CLUSTERS],[RPC1_ETA_CLUSTERS]],[[...],[...]],...]
    RPC_heights = [0.6,1.8,3.0,61.8,121.8,123] #Heights of middle point of each RPC, measured from the bottom of the Triplet Low RPC. Units are cm.


    #Extract x and y coords of cluster in event
    coords = extract_coords_timed_Chi2(event,max_cluster_size)

     # Filter out coords of RPC under test 

    test_coords = coords[RPC_excluded]

    coords[RPC_excluded] = [[],[],"N"] 
    
    # if RPC_excluded[1] < 6:
    #     coords[RPC_excluded[1]] = [[],[],"N"] 

    # Count the number of empty RPCs
    empty_RPC_count = sum(1 for item in coords if item == [[], [],'N'])

    # If less than 3 elements of coords are occupied, exit the function
    if empty_RPC_count > 3:
        #print("Failed to reconstruct, not enough coords")
        return None  # Exit the function
    
    cross_chamberness = 0

    if coords[0] != [[], [], 'N'] or coords[1] != [[], [], 'N'] or coords[2] != [[], [], 'N']:
        cross_chamberness += 1

    if coords[3] != [[], [], 'N']:
        cross_chamberness += 1

    if coords[4] != [[], [], 'N'] or coords[5] != [[], [], 'N']:
        cross_chamberness += 1

    if cross_chamberness < 2:
        #print("Failed to reconstruct, too few chambers")
        return None

    #ITERATING OVER EVERY POSSIBLE COMBINATION OF x,y,z over all 3 RPCs (limited to one x,y per RPC).
    #Doesn't look particularly nice, but there are not many coordinates to loop over usually....

    combinations = generate_hit_coords_combo_Chi2(coords,RPC_heights)

    #Now for each combo in combinations, attempt to reconstruct a path. See which one gives the best trajectory.

    #If success, print parameters of fitting function.
    #If fail, print reconstruction failed.

    Chi2_current = np.inf
    optimised_coords = None
    optimised_d= None
    optimised_centroid= None
    dT = np.inf

    for ind,combo in enumerate(combinations):

        centroid, d, Chi2, coordinates, delta_T, delta_Z= fit_event_chi2(combo)
        if Chi2 < Chi2_current:

            # If new fit is better than old then replace old fit properties.
            dZ = delta_Z 
            dT = delta_T
            Chi2_current = Chi2
            optimised_centroid = centroid
            optimised_d = d
            optimised_coords = coordinates

    #if dT>0 this implies the particles hit the higher RPC after the lower one, so the particle is travelling upwards here.
    #Vice-versa for dT < 0.

    #dT = 0 case?

    if dT != np.inf:

        if dT > 0:
            if optimised_d[2] < 0:
                optimised_d = np.multiply(optimised_d,-1)
        else:
            if optimised_d[2] > 0:
                optimised_d = np.multiply(optimised_d,-1)

        if Chi2_current<max_Chi2:
            return optimised_centroid, optimised_d, optimised_coords, combinations, Chi2_current, dT, dZ, test_coords

    else:
        #print("Failed to reconstruct, Chi2 too large")
        #return optimised_centroid, optimised_d, optimised_coords, combinations, residuals_current
        return None
def check_event_attributes_by_RPC(event, min_chamber_number, min_RPC_number, RPC_excluded):
    # Used in filter_events() function to decide whether or not to save an event as described by the user's inputs.
    # The user selects an RPC to exclude from the filter.
    # e.g. say we want to exclude RPC 4 and the user selects a min_RPC_number of 4.
    # The function will check if the event has at least 4 RPCs hit.

    # USING ONLY ETA FILTER HERE SINCE THIS IS WHAT TRIGGERS THE CHANNEL!

    # event = ['Event x', TIMEBIN, [[[RPC1_PHI_CLUSTERS],[RPC1_ETA_CLUSTERS]],[[...],[...]],...]

    RPC_counter = 0
    chamber_counter = 0
    condition_1 = False
    condition_2 = False
    condition_3 = False

    # Check if any of the excluded RPCs have ETA clusters
    # for excluded_RPC in RPC_excluded:
    if event[2][RPC_excluded][1] != []:
        return False

    for RPC in range(6):
        if RPC != RPC_excluded:
            # If the RPC is not excluded, we check for ETA clusters.
            if event[2][RPC][1]:
                RPC_counter += 1

                if RPC < 3:
                    # Checking triplet layer.
                    if not condition_1:
                        chamber_counter += 1
                        condition_1 = True
                elif RPC == 3:
                    # Checking singlet layer.
                    if not condition_2:
                        chamber_counter += 1
                        condition_2 = True
                else:
                    # Checking doublet layer.
                    if not condition_3:
                        chamber_counter += 1
                        condition_3 = True

    return RPC_counter >= min_RPC_number and chamber_counter >= min_chamber_number
def filter_events_by_RPC(events,min_chamber_number,min_RPC_number,RPC_excluded):
    #Initiliase array of filtered events
    filtered_events = []

    for event in events:
        if check_event_attributes_by_RPC(event,min_chamber_number,min_RPC_number,RPC_excluded):
            filtered_events.append(event)
        else :
            event[2][RPC_excluded][1] = []
            event[2][RPC_excluded][0] = []
            filtered_events.append(event)
                

   # print(f"Number of events in filter = {len(filtered_events)}")
    
    return filtered_events
def reconstruct_timed_Chi2(event,max_cluster_size, max_length=None, exact_length=False):

    #timed tag indicates that timing information from RPC is used to determine direction of vertical transversal of "particle" in the event.

    max_Chi2 = 10

    # event = ['Event x',TIMEBIN, [[[RPC1_PHI_CLUSTERS],[RPC1_ETA_CLUSTERS]],[[...],[...]],...]
    RPC_heights = [0.6,1.8,3.0,61.8,121.8,123] #Heights of middle point of each RPC, measured from the bottom of the Triplet Low RPC. Units are cm.

    #Extract x and y coords of cluster in event

    coords = extract_coords_timed_Chi2(event,max_cluster_size)

    # Count the number of empty RPCs
    empty_RPC_count = sum(1 for item in coords if item == [[], [],'N'])

    # If less than 3 elements of coords are occupied, exit the function
    if empty_RPC_count > 3:
        #print("Failed to reconstruct, not enough coords")
        return None  # Exit the function
    
    #NEED TO CHECK IF STILL CROSS CHAMBER! 

    cross_chamberness = 0

    if coords[0] != [[], [], 'N'] or coords[1] != [[], [], 'N'] or coords[2] != [[], [], 'N']:
        cross_chamberness += 1

    if coords[3] != [[], [], 'N']:
        cross_chamberness += 1

    if coords[4] != [[], [], 'N'] or coords[5] != [[], [], 'N']:
        cross_chamberness += 1

    # print(coords[0])
    # print(coords[1])
    # print(coords[2])
    # print(coords[3])
    # print(coords[4])
    # print(coords[5])
    # print(cross_chamberness)

    if cross_chamberness < 2:
        #print("Failed to reconstruct, too few chambers")
        return None

    #ITERATING OVER EVERY POSSIBLE COMBINATION OF x,y,z over all 3 RPCs (limited to one x,y per RPC).
    #Doesn't look particularly nice, but there are not many coordinates to loop over usually....

    combinations = generate_hit_coords_combo_Chi2(coords,RPC_heights, max_length=max_length, exact_length=exact_length)

    #Now for each combo in combinations, attempt to reconstruct a path. See which one gives the best trajectory.

    #If success, print parameters of fitting function.
    #If fail, print reconstruction failed.

    Chi2_current = np.inf
    optimised_coords = None
    optimised_d= None
    optimised_centroid= None
    dT = np.inf

    for ind,combo in enumerate(combinations):

        centroid, d, Chi2, coordinates, delta_T, delta_Z= fit_event_chi2(combo)
        if Chi2 < Chi2_current:

            # If new fit is better than old then replace old fit properties.
            dZ = delta_Z 
            dT = delta_T
            Chi2_current = Chi2
            optimised_centroid = centroid
            optimised_d = d
            optimised_coords = coordinates

    #if dT>0 this implies the particles hit the higher RPC after the lower one, so the particle is travelling upwards here.
    #Vice-versa for dT < 0.

    #dT = 0 case?
    
    if dT != np.inf:

        if dT > 0:
            if optimised_d[2] < 0:
                optimised_d = np.multiply(optimised_d,-1)
        else:
            if optimised_d[2] > 0:
                optimised_d = np.multiply(optimised_d,-1)

        if Chi2_current<max_Chi2:
            return optimised_centroid, optimised_d, optimised_coords, combinations, Chi2_current, dT, dZ

        else:
            #print("Failed to reconstruct, Chi2 too large")
            #return optimised_centroid, optimised_d, optimised_coords, combinations, residuals_current
            return None
def extract_angles_phi_eta_timed_DZ(filtered_events, max_length=None, exact_length=False):

    #Input is filtered_events, output of ANT.filter_events() function

    angles_eta = []
    angles_phi = []
    delta_times = []
    dZ = []
    chi2_values = []

    for i,filtered_event in enumerate(filtered_events):

        # print(f"Index= {i}") 
        
        result = reconstruct_timed_Chi2(filtered_event,3, max_length=max_length, exact_length=exact_length)

        if result is not None:

            delta_times.append(result[5])

            chi2_values.append(result[4])

            # if abs(result[5])>40:
            #     print(f"Index of abnormal event= {i}")
            # #Only save angles that actually were reconstructed well

            dZ.append(result[6])
            
            # a.b = |a||b|cos(x)

            #eta angle. 
            #work out the projection of the direction vector in the plane.
            
            v_parr_eta = np.array([0,result[1][1],result[1][2]])

            theta_eta = np.arccos(np.dot(v_parr_eta,[0,0,1]) / np.linalg.norm(v_parr_eta))

            if theta_eta > np.pi / 2:
                theta_eta= np.pi - theta_eta
            
            if v_parr_eta[1] > 0:
                theta_eta*=-1

            angles_eta.append(theta_eta)

            # Phi angles
            #work out the projection of the direction vector in the plane.
            
            v_parr_phi = np.array([result[1][0],0,result[1][2]])

            theta_phi = np.arccos(np.dot(v_parr_phi,[0,0,1]) / np.linalg.norm(v_parr_phi))

            if theta_phi > np.pi / 2:
                theta_phi= np.pi - theta_phi
            
            if v_parr_phi[0] < 0:
                theta_phi*=-1

            angles_phi.append(theta_phi)

    return angles_eta, angles_phi, delta_times, dZ, chi2_values
def interactive_muon_plot(centroid,d,event_coords):
    #Coefficients = [a,b,c]
    #event_coords = [[x0,y0,z0],[x1,y1,z1],...,[x5,y5,z5]]

    # Generate line coordinates
    t_values = np.linspace(-150, 150, 100)
    line_coordinates = centroid.reshape(-1, 1) + d.reshape(-1, 1) * t_values.reshape(1, -1)

    x_coords = [x[0] for x in event_coords]
    y_coords = [y[1] for y in event_coords]
    z_coords = [z[2] for z in event_coords]


    RPC_origins = [[0,0,0],[0,0,1.2],[0,0,2.4],[0,0,61.2],[0,0,121.2],[0,0,122.4]]
    RPC_dimensions = [180,99,1.2]

    #Calculate vertices for RPCs
    rpc_vertices = []
    for i in range(6):
        rpc_vertices.append(calculate_cuboid_vertices(RPC_origins[i],RPC_dimensions))

    # Configure Plotly to be rendered inline in the notebook.
    plotly.offline.init_notebook_mode()

    # Configure the trace for the RPCs
    rpc_0 = go.Mesh3d(
        x=[vertex[0] for vertex in rpc_vertices[0]],
        y=[vertex[1] for vertex in rpc_vertices[0]],
        z=[vertex[2] for vertex in rpc_vertices[0]],
        i=[0, 1, 2, 3, 0, 4, 5, 6, 7, 4, 5, 1, 2, 6, 7, 3],  # Indices for faces
        j=[1, 2, 3, 0, 4, 5, 6, 7, 5, 6, 2, 3, 7, 4, 0, 1],  # Indices for faces
        k=[2, 3, 0, 1, 6, 7, 4, 5, 1, 5, 6, 2, 3, 7, 4, 0],  # Indices for faces
        opacity=0.1,  # Set the opacity as needed
        color='green'  # Set the color of the cuboid
    )

    rpc_1 = go.Mesh3d(
        x=[vertex[0] for vertex in rpc_vertices[1]],
        y=[vertex[1] for vertex in rpc_vertices[1]],
        z=[vertex[2] for vertex in rpc_vertices[1]],
        i=[0, 1, 2, 3, 0, 4, 5, 6, 7, 4, 5, 1, 2, 6, 7, 3],  # Indices for faces
        j=[1, 2, 3, 0, 4, 5, 6, 7, 5, 6, 2, 3, 7, 4, 0, 1],  # Indices for faces
        k=[2, 3, 0, 1, 6, 7, 4, 5, 1, 5, 6, 2, 3, 7, 4, 0],  # Indices for faces
        opacity=0.1,  # Set the opacity as needed
        color='green'  # Set the color of the cuboid
    )

    rpc_2 = go.Mesh3d(
        x=[vertex[0] for vertex in rpc_vertices[2]],
        y=[vertex[1] for vertex in rpc_vertices[2]],
        z=[vertex[2] for vertex in rpc_vertices[2]],
        i=[0, 1, 2, 3, 0, 4, 5, 6, 7, 4, 5, 1, 2, 6, 7, 3],  # Indices for faces
        j=[1, 2, 3, 0, 4, 5, 6, 7, 5, 6, 2, 3, 7, 4, 0, 1],  # Indices for faces
        k=[2, 3, 0, 1, 6, 7, 4, 5, 1, 5, 6, 2, 3, 7, 4, 0],  # Indices for faces
        opacity=0.1,  # Set the opacity as needed
        color='green'  # Set the color of the cuboid
    )

    rpc_3 = go.Mesh3d(
        x=[vertex[0] for vertex in rpc_vertices[3]],
        y=[vertex[1] for vertex in rpc_vertices[3]],
        z=[vertex[2] for vertex in rpc_vertices[3]],
        i=[0, 1, 2, 3, 0, 4, 5, 6, 7, 4, 5, 1, 2, 6, 7, 3],  # Indices for faces
        j=[1, 2, 3, 0, 4, 5, 6, 7, 5, 6, 2, 3, 7, 4, 0, 1],  # Indices for faces
        k=[2, 3, 0, 1, 6, 7, 4, 5, 1, 5, 6, 2, 3, 7, 4, 0],  # Indices for faces
        opacity=0.1,  # Set the opacity as needed
        color='green'  # Set the color of the cuboid
    )

    rpc_4 = go.Mesh3d(
        x=[vertex[0] for vertex in rpc_vertices[4]],
        y=[vertex[1] for vertex in rpc_vertices[4]],
        z=[vertex[2] for vertex in rpc_vertices[4]],
        i=[0, 1, 2, 3, 0, 4, 5, 6, 7, 4, 5, 1, 2, 6, 7, 3],  # Indices for faces
        j=[1, 2, 3, 0, 4, 5, 6, 7, 5, 6, 2, 3, 7, 4, 0, 1],  # Indices for faces
        k=[2, 3, 0, 1, 6, 7, 4, 5, 1, 5, 6, 2, 3, 7, 4, 0],  # Indices for faces
        opacity=0.1,  # Set the opacity as needed
        color='green'  # Set the color of the cuboid
    )

    rpc_5 = go.Mesh3d(
        x=[vertex[0] for vertex in rpc_vertices[5]],
        y=[vertex[1] for vertex in rpc_vertices[5]],
        z=[vertex[2] for vertex in rpc_vertices[5]],
        i=[0, 1, 2, 3, 0, 4, 5, 6, 7, 4, 5, 1, 2, 6, 7, 3],  # Indices for faces
        j=[1, 2, 3, 0, 4, 5, 6, 7, 5, 6, 2, 3, 7, 4, 0, 1],  # Indices for faces
        k=[2, 3, 0, 1, 6, 7, 4, 5, 1, 5, 6, 2, 3, 7, 4, 0],  # Indices for faces
        opacity=0.1,  # Set the opacity as needed
        color='green'  # Set the color of the cuboid
    )

    # Configure the trace.
    trace = go.Scatter3d(
        x=x_coords,  # <-- Put your data instead
        y=y_coords,  # <-- Put your data instead
        z=z_coords, # <-- Put your data instead
        mode='markers',
        marker={
            'size': 10,
            'opacity': 1.0,
        }
    )

    #Plot zenith angle
    # zenith = go.Vector(...)

    # Extract x, y, z coordinates from line coordinates
    x_line = line_coordinates[0]
    y_line = line_coordinates[1]
    z_line = line_coordinates[2]
    
    # Configure the trace for the line
    trace_line = go.Scatter3d(
        x=x_line,
        y=y_line,
        z=z_line,
        mode='lines',
        line={
            'color': 'red',
            'width': 2,
        },
        name='Line'
    )

    trace_centroid = go.Scatter3d(
        x=[centroid[0]],
        y=[centroid[1]],
        z=[centroid[2]],
        mode='markers',
        marker={
            'size': 5,
            'color': 'green',
        },
        name='Centroid'
    )

    # Configure the layout.
    layout = go.Layout(
        margin={'l': 0, 'r': 0, 'b': 0, 't': 0},
        scene=dict(
            xaxis=dict(showgrid=True, gridcolor='rgb(211,211,211)', gridwidth=2,range=[-5,185], title = "x/cm"),
            yaxis=dict(showgrid=True, gridcolor='rgb(211,211,211)', gridwidth=2,range=[-5, 105],title = "y/cm"),
            zaxis=dict(showgrid=True, gridcolor='rgb(211,211,211)', gridwidth=2,range=[-5, 125],title="z/cm")
        )
    )

    #Include trace_centroid in data= [] to plot centroid on the plot.

    data = [trace,trace_line,rpc_0,rpc_1,rpc_2,rpc_3,rpc_4,rpc_5]

    plot_figure = go.Figure(data=data, layout=layout)

    # Render the plot.
    plotly.offline.iplot(plot_figure)
def calculate_cuboid_vertices(origin, dimensions):
    # Function to calculate vertices of cuboid from origin and dimensions
    x_min, y_min, z_min = origin
    x_max = x_min + dimensions[0]
    y_max = y_min + dimensions[1]
    z_max = z_min + dimensions[2]

    vertices = [
        [x_min, y_min, z_min], [x_max, y_min, z_min],
        [x_max, y_max, z_min], [x_min, y_max, z_min],
        [x_min, y_min, z_max], [x_max, y_min, z_max],
        [x_max, y_max, z_max], [x_min, y_max, z_max]
    ]

    return vertices
def find_tdc_alignment_metric(tdc0, tdc1):
    if tdc0 > tdc1:
        tdc0, tdc1 = tdc1, tdc0
    i, j, k, l = None, None, None, None
    if tdc0 == 0:
        if tdc1 == 1:
            i, j, k, l = 1, 2, 0, 1
        if tdc1 == 2:
            i, j, k, l = 3, 0, 3, 0
        if tdc1 == 3:
            i, j, k, l = 4, 0, 4, 0
        if tdc1 == 4:
            i, j, k, l = -1, -1, 5, 0
    if tdc0 == 1:
        if tdc1 == 2:
            i, j, k, l = 3, 2, 3, 1
        if tdc1 == 3:
            i, j, k, l = 5, 2, 4, 1
        if tdc1 == 4:
            i, j, k, l = -1, -1, 5, 1
    if tdc0 == 2:
        if tdc1 == 3:
            i, j, k, l = 4, 3, 4, 3
        if tdc1 == 4:
            i, j, k, l = -1, -1, 5, 3
    if tdc0 == 3:
        if tdc1 == 4:
            i, j, k, l = -1, -1, 5, 4
        
    return i, j, k, l
        
def get_element(i, j, data):
    if i == j:
        return None
    if i < j:
        return data.get((i, j), None)
    else:
        return data.get((j, i), None)    
#Reconstruction, with timing of RPC hit taken into account in trajectory.
def testAlign(rpc1Hits, rpc2Hits, skipChans = []):
    minTimes = [300,300]
    minChans = [-1,-1]
    if len(rpc1Hits)<1 or len(rpc2Hits)<1:
        return -1
    for hit in rpc1Hits:
        if hit.time<minTimes[0] and hit.channel not in skipChans:
            minTimes[0]=hit.time
            minChans[0]=hit.channel
    for hit in rpc2Hits:
        if hit.time<minTimes[1] and hit.channel not in skipChans:
            minTimes[1]=hit.time
            minChans[1]=hit.channel
    if -1 in minChans:
        return -1
    return abs(minChans[1]-minChans[0])
def count_entries(matrix):
    total_entries = 0
    for row in matrix:
        for cell in row:
            total_entries += len(cell)
    return total_entries

def reconstruct_timed_Chi2_modified(event,max_cluster_size, max_length=None, exact_length=False):

    #timed tag indicates that timing information from RPC is used to determine direction of vertical transversal of "particle" in the event.

    max_Chi2 = 10

    # event = ['Event x',TIMEBIN, [[[RPC1_PHI_CLUSTERS],[RPC1_ETA_CLUSTERS]],[[...],[...]],...]
    RPC_heights = [0.6,1.8,3.0,61.8,121.8,123] #Heights of middle point of each RPC, measured from the bottom of the Triplet Low RPC. Units are cm.

    #Extract x and y coords of cluster in event

    coords = extract_coords_timed_Chi2(event,max_cluster_size)

    # Count the number of empty RPCs
    empty_RPC_count = sum(1 for item in coords if item == [[], [],'N'])

    # If less than 3 elements of coords are occupied, exit the function
    if empty_RPC_count > 3:
        #print("Failed to reconstruct, not enough coords")
        return None  # Exit the function
    
    #NEED TO CHECK IF STILL CROSS CHAMBER! 

    cross_chamberness = 0

    if coords[0] != [[], [], 'N'] or coords[1] != [[], [], 'N'] or coords[2] != [[], [], 'N']:
        cross_chamberness += 1

    if coords[3] != [[], [], 'N']:
        cross_chamberness += 1

    if coords[4] != [[], [], 'N'] or coords[5] != [[], [], 'N']:
        cross_chamberness += 1

    # print(coords[0])
    # print(coords[1])
    # print(coords[2])
    # print(coords[3])
    # print(coords[4])
    # print(coords[5])
    # print(cross_chamberness)

    # if cross_chamberness < 2:
    #     #print("Failed to reconstruct, too few chambers")
    #     return None

    #ITERATING OVER EVERY POSSIBLE COMBINATION OF x,y,z over all 3 RPCs (limited to one x,y per RPC).
    #Doesn't look particularly nice, but there are not many coordinates to loop over usually....

    combinations = generate_hit_coords_combo_Chi2(coords,RPC_heights, max_length=max_length, exact_length=exact_length)

    #Now for each combo in combinations, attempt to reconstruct a path. See which one gives the best trajectory.

    #If success, print parameters of fitting function.
    #If fail, print reconstruction failed.

    Chi2_current = np.inf
    optimised_coords = None
    optimised_d= None
    optimised_centroid= None
    dT = np.inf

    for ind,combo in enumerate(combinations):

        centroid, d, Chi2, coordinates, delta_T, delta_Z= fit_event_chi2(combo)
        if Chi2 < Chi2_current:

            # If new fit is better than old then replace old fit properties.
            dZ = delta_Z 
            dT = delta_T
            Chi2_current = Chi2
            optimised_centroid = centroid
            optimised_d = d
            optimised_coords = coordinates

    #if dT>0 this implies the particles hit the higher RPC after the lower one, so the particle is travelling upwards here.
    #Vice-versa for dT < 0.

    #dT = 0 case?
    
    if dT != np.inf:
        
        if optimised_d[2] > 0:
            optimised_d = np.multiply(optimised_d,-1)

        if Chi2_current<max_Chi2:
            return optimised_centroid, optimised_d, optimised_coords, combinations, Chi2_current, dT, dZ

        else:
            #print("Failed to reconstruct, Chi2 too large")
            #return optimised_centroid, optimised_d, optimised_coords, combinations, residuals_current
            return None

def does_muon_hit_RPC(optimised_centroid, optimised_d, RPC):

    RPC_heights = [0.6,1.8,3.0,61.8,121.8,123] 
    #Heights of middle point of each RPC, measured from the bottom of the Triplet Low RPC. Units are cm.

    # x_bar = x_centroid + d_vector * t
    # Find value of paramter t when the muon trajectory passes through the RPC height.
    
    z_0 = optimised_centroid[2]
    z = RPC_heights[RPC]

    # t = (z-z_0)/d_z

    t = (z-z_0)/optimised_d[2]

    # Find expected (x,y) coordinates at that height.

    x = optimised_centroid[0] + t*optimised_d[0]
    y = optimised_centroid[1] + t*optimised_d[1]

    # Check if these (x,y) coordinates lie within the RPC. 

    #Extract x and y coords of cluster in event
    distance_per_phi_channel = 2.7625 #cm
    distance_per_eta_channel = 2.9844 #cm

    # Max y (eta side) is 31.5 * distance_per_eta_channel
    # Max x (phi side) is 63.5 * distance_per_phi_channel

    if 0 < x < 63.5*distance_per_phi_channel and 0 < y < 31.5*distance_per_eta_channel:
        #Return coordinates where you expect the muon to hit this RPC from the reconstructed event.
        return [x,y]
    else:
        #print("Muon does not hit RPC")
        return None    

def does_RPC_detect_muon(muon_coords,test_coords,tol):
    #Tolerance in units of cm. 

    #Could experiment with tolerance.
    if test_coords != [[],[],"N"]: 

        x_coords = test_coords[0]
        y_coords = test_coords[1]

        for x_set in x_coords:
            for y_set in y_coords:

                x = x_set[0]
                y = y_set[0]
    
                #If statement ensures only calculate the coords if the test_coords actually exist.

                #Offset is 2D vector that represents difference 
                offset = np.subtract(np.array([x,y]),muon_coords)

                separation = np.linalg.norm(offset)

                #print(separation)

                if separation <= tol:
                    #Say the RPC only successfully reconstructs an event 
                    #if the distance between expected hit and reconstructed hit is less than tolerance.

                    #print("RPC successfully detects hit!")
                    return separation
        
        #print("No RPC coordinates constructed pass near the expected point!")
        return False

    else:
        #print("No coordinates reconstructed by RPC")
        return False   
 
def calc_efficiency_RPC(dataset, RPC, tol):
    # Ensure RPC is a list, even if it's a single integer
    RPC_excluded = RPC

    # events = filter_events_by_RPC(dataset, 2, 5, RPC_excluded)
    datas = dataset

    possible_reconstructions = 0
    successful_reconstructions = [0 for i in range(len(tol))]

    for i, data in enumerate(datas):

        # print(f"Event index {i}")

        E_recon = reconstruct_timed_Chi2_ByRPC(data, 3, RPC_excluded)

        if E_recon:
            if len(E_recon[2]) >= 5:
                # Adding this check to see if other 5 RPCs are in reconstructed event.
                # This is necessary to ensure the reconstructed path is accurate.

                muon_coords = does_muon_hit_RPC(E_recon[0], E_recon[1], RPC_excluded)

                if muon_coords:
                    possible_reconstructions += 1
                    for idx, t in enumerate(tol):
                        check = does_RPC_detect_muon(muon_coords, E_recon[7], t)
                        if check:
                            successful_reconstructions[idx] += 1

    # print(possible_reconstructions)
    # print(successful_reconstructions)
    # if possible_reconstructions == 0:
    #     return -1
    # print(possible_reconstructions)
    # print(successful_reconstructions)

    if possible_reconstructions == 0:
        return -1
    return [i / possible_reconstructions for i in successful_reconstructions]
def tdcEventToRPCData(event,activeTDCs=[0,1,2,3,4], event_num = 0):
    rpcHits = [[] for rpc in range(6)]
    for tdc in activeTDCs:
        for word in event.tdcEvents[tdc].words:
            rpc, thisHit = tdcChanToRPCHit(word,tdc, event_num=event_num)
            rpcHits[rpc].append(thisHit)
    return rpcHits

def calcAvgAlign(event_chunk, offSet=0, i = 1, j = 2, k = 0, l = 2, tdc1 =0, tdc0 = 1, processedEvents = 0, skipChans = []):
    mets = []
    for idx, event in enumerate(event_chunk):
        etaHits = [[] for rpc in range(6)]
        phiHits = [[] for rpc in range(6)]
        if (idx+abs(offSet))<len(event_chunk):
            if offSet<=0:
                oneIdx = idx+abs(offSet)
                twoIdx = idx
            else:
                oneIdx = idx
                twoIdx = idx+offSet
            for word in event_chunk[oneIdx].tdcEvents[tdc1].words:
                rpc, thisHit = tdcChanToRPCHit(word,tdc1, processedEvents + idx)
                if thisHit.eta:
                    etaHits[rpc].append(thisHit)

                else:
                    phiHits[rpc].append(thisHit)
            for word in event_chunk[twoIdx].tdcEvents[tdc0].words:
                rpc, thisHit = tdcChanToRPCHit(word,tdc0, processedEvents + idx)
                if thisHit.eta:
                    etaHits[rpc].append(thisHit)

                else:
                    phiHits[rpc].append(thisHit)     
            if i != -1:  
                etOff = testAlign(etaHits[i],etaHits[j], skipChans = skipChans)
                phOff = testAlign(phiHits[k],phiHits[l], skipChans = skipChans)
                if etOff>=0 and phOff>=0:
                    mets.append(math.sqrt(etOff*etOff+phOff*phOff))
            else:
                phOff = testAlign(phiHits[k],phiHits[l], skipChans = skipChans)
                if phOff>=0:
                    mets.append(math.sqrt(phOff*phOff))
    if len(mets)>0:
        return sum(mets)/len(mets)
    else:
        return 100


def ConstructEventInsertionList(updates, order):
    insertion_list = [0 for _ in range(5)]
    for idx, update in enumerate(updates):
        i, j = order[idx]
        insertion_list[j] += (insertion_list[i] - update)
    min_value = min(insertion_list)
    if min_value < 0:
        addition = abs(min_value)
    else:
        addition = 0
    return [x + addition for x in insertion_list]
    
def plot_cluster_size_histograms(cluster_sizes):
    phi_cluster_sizes = []
    eta_cluster_sizes = []

    for event_cluster_sizes in cluster_sizes:
        for rpc_cluster_sizes in event_cluster_sizes:
            phi_cluster_sizes.extend(rpc_cluster_sizes[0])
            eta_cluster_sizes.extend(rpc_cluster_sizes[1])

    # Calculate bin edges
    phi_bin_edges = np.arange(0.5, max(phi_cluster_sizes) + 1.5, 1)
    eta_bin_edges = np.arange(0.5, max(eta_cluster_sizes) + 1.5, 1)

    # Plot the phi cluster size histogram
    plt.figure(figsize=(10, 6))
    plt.hist(phi_cluster_sizes, bins=phi_bin_edges, edgecolor='black', align='mid', density=True)
    plt.xlabel('Cluster Size')
    plt.ylabel('Proportion of Events')
    plt.title('Phi Cluster Size Distribution')
    plt.xticks(np.arange(1, max(phi_cluster_sizes) + 1))
    plt.xlim(0, 5)
    plt.grid(axis='y')
    plt.show()

    # Plot the eta cluster size histogram
    plt.figure(figsize=(10, 6))
    plt.hist(eta_cluster_sizes, bins=eta_bin_edges, edgecolor='black', align='mid', density=True)
    plt.xlabel('Cluster Size')
    plt.ylabel('Proportion of Events')
    plt.title('Eta Cluster Size Distribution')
    plt.xticks(np.arange(1, max(eta_cluster_sizes) + 1))
    plt.xlim(0, 10)
    plt.grid(axis='y')
    plt.show()