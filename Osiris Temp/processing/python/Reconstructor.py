import mplhep as hep
import numpy as np
hep.style.use([hep.style.ATLAS])
import sys
import importlib
sys.path.insert(1, 'C://Users//Peter//OneDrive - University of Cambridge//Desktop//summer2//Osiris Temp//processing//python')
import Analysis_tools as ATools
import Reconstruction_Tools as RTools
importlib.reload(RTools)


class rpcCoincidence():
    def __init__(self, event_num, time_bin, hits):
        self.event_num = event_num
        self.time_bin = time_bin
        self.hits = hits

    def __str__(self):
        return f"rpcCoincidence(event_num={self.event_num}, time_bin={self.time_bin}, hits={self.hits})"

class Reconstructor():
    def __init__(self, event_chunk, processsed_event, tolerance = None, coincidence_window = 15, tof_correction = True):
        self.event_chunk = event_chunk
        self.etaHits = [[] for rpc in range(6)]
        self.phiHits = [[] for rpc in range(6)]
        self.window_size = coincidence_window
        self.tof_correction = tof_correction
        self.processedEvents = processsed_event
        self.tol = [i for i in range(20)] if tolerance is None else tolerance
        self.dT = []
        self.recon = []
        
        
        
        self.possible_reconstructions = [0 for rppc in range(6)]
        self.successful_reconstructions = [[0 for i in range(len(self.tol))] for rpc in range(6)]
        self.successful_reconstructed_coords = {0:[[0 for etchan in range(32)] for phchan in range(64)],
                1:[[0 for etchan in range(32)] for phchan in range(64)], 
                2:[[0 for etchan in range(32)] for phchan in range(64)], 
                3:[[0 for etchan in range(32)] for phchan in range(64)],
                4:[[0 for etchan in range(32)] for phchan in range(64)],
                5:[[0 for etchan in range(32)] for phchan in range(64)]}
        self.failed_reconstructed_coords = {0:[[0 for etchan in range(32)] for phchan in range(64)],
                1:[[0 for etchan in range(32)] for phchan in range(64)], 
                2:[[0 for etchan in range(32)] for phchan in range(64)], 
                3:[[0 for etchan in range(32)] for phchan in range(64)],
                4:[[0 for etchan in range(32)] for phchan in range(64)],
                5:[[0 for etchan in range(32)] for phchan in range(64)]}
        self.possible_reconstructions_coords = {0:[[0 for etchan in range(32)] for phchan in range(64)],
                1:[[0 for etchan in range(32)] for phchan in range(64)], 
                2:[[0 for etchan in range(32)] for phchan in range(64)], 
                3:[[0 for etchan in range(32)] for phchan in range(64)],
                4:[[0 for etchan in range(32)] for phchan in range(64)],
                5:[[0 for etchan in range(32)] for phchan in range(64)]}
        
        
        self.eta_histogram = np.zeros(len(np.arange(-90.5, 91.5, 1)) - 1)
        self.phi_histogram = np.zeros(len(np.arange(-90.5, 91.5, 1)) - 1)
        self.solid_theta_histogram = np.zeros(len(np.arange(-180.5, 181.5, 1)) - 1)
        self.solid_phi_histogram = np.zeros(len(np.arange(-180.5, 181.5, 1)) - 1)
        self.tdcstatus = [True for tdc in range(5)]
        

    def populate_hits(self):
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
                        
    
    def update_event(self, event_chunk, processed_event):
        self.event_chunk = event_chunk
        self.processedEvents = processed_event
        self.etaHits = [[] for rpc in range(6)]
        self.phiHits = [[] for rpc in range(6)]
        
    
    def make_cluster(self):
        coincident_hits = RTools.FindCoincidentHits(self.etaHits, self.phiHits, self.window_size, tof_correction=self.tof_correction)
        clustered = RTools.cluster(coincident_hits)
        return clustered
    
    
    def reconstruct_and_extrapolate(self, dataset, chi2_region = [0, 100]):
        # Ensure RPC is a list, even if it's a single integer
        if self.tdcstatus[3] == True:
            for rpc in range(6):
                for i, data in enumerate(dataset):
                    E_recon = RTools.reconstruct_timed_Chi2_ByRPC(data, 3, rpc)
                    print(E_recon)
                    if E_recon:
                        if len(E_recon[2]) >= 5:
                            if E_recon[4] > chi2_region[0] and E_recon[4] < chi2_region[1]:
                                # self.chi2.append(E_recon[4])
                                # self.event_of_interest.append(E_recon)
                                # Adding this check to see if other 5 RPCs are in reconstructed event.
                                # This is necessary to ensure the reconstructed path is accurate.

                                muon_coords = RTools.does_muon_hit_RPC(E_recon[0], E_recon[1], rpc)
                                print(muon_coords)
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


    def reconstruct_and_findtof(self, dataset):
        for i, data in enumerate(dataset):
                E_recon = RTools.reconstruct_timed_Chi2_ByRPC(data, 3, -1)
                if E_recon:
                    if len(E_recon[2]) >= 6:
                        self.dT.append(E_recon[5])
                        self.recon.append(E_recon)
                        
    def extract_angles_phi_eta_timed_DZ_modified(self, filtered_events, max_length=None, exact_length=False):
        angles_eta = []
        angles_phi = []
        delta_times = []
        dZ = []
        chi2_values = []
        solid_angle = []
        
        bin_size = 1
        eta_phi_bin_edges = np.arange(-90.5, 91.5, bin_size)
        solid_phi_bin_edges = np.arange(-180.5, 181.5, bin_size)

        for i, filtered_event in enumerate(filtered_events):
            result = RTools.reconstruct_timed_Chi2_modified(filtered_event, 3, max_length=max_length, exact_length=exact_length)

            if result is not None:
                delta_times.append(result[5])
                chi2_values.append(result[4])
                dZ.append(result[6])

                v_parr_eta = np.array([0, result[1][1], result[1][2]])
                theta_eta = np.arccos(np.dot(v_parr_eta, [0, 0, 1]) / np.linalg.norm(v_parr_eta))

                if theta_eta > np.pi / 2:
                    theta_eta = np.pi - theta_eta
                if v_parr_eta[1] > 0:
                    theta_eta *= -1

                angles_eta.append(theta_eta)

                v_parr_phi = np.array([result[1][0], 0, result[1][2]])
                theta_phi = np.arccos(np.dot(v_parr_phi, [0, 0, 1]) / np.linalg.norm(v_parr_phi))

                if theta_phi > np.pi / 2:
                    theta_phi = np.pi - theta_phi
                if v_parr_phi[0] < 0:
                    theta_phi *= -1

                angles_phi.append(theta_phi)

                magnitude = np.linalg.norm(result[1])
                theta = np.arccos(result[1][2] / -magnitude)  # Polar angle
                phi = np.arctan2(result[1][1], result[1][0])  # Azimuthal angle

                # Normalize phi to the range [-π, π]
                if phi > np.pi:
                    phi -= 2 * np.pi
                elif phi < -np.pi:
                    phi += 2 * np.pi

                # Convert angles to degrees
                theta_deg = np.degrees(theta)
                
                phi_deg = np.degrees(phi)

                # Append the solid angle (theta, phi) to self.angles_solid
                solid_angle.append((theta_deg, phi_deg))

        angles_eta_degrees = [x * (180 / np.pi) for x in angles_eta]
        angles_phi_degrees = [x * (180 / np.pi) for x in angles_phi]
        angles_solid_theta_degrees = [x for x, y in solid_angle]
        angles_solid_phi_degrees = [y for x, y in solid_angle]

        self.eta_histogram += np.histogram(angles_eta_degrees, bins=eta_phi_bin_edges)[0]
        self.phi_histogram += np.histogram(angles_phi_degrees, bins=eta_phi_bin_edges)[0]
        self.solid_theta_histogram += np.histogram(angles_solid_theta_degrees, bins=solid_phi_bin_edges)[0]
        self.solid_phi_histogram += np.histogram(angles_solid_phi_degrees, bins=solid_phi_bin_edges)[0]