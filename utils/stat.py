import numpy as np
import matplotlib.pyplot as plt
from utils.logger import setup_logger
from utils.utils import printr, printg, printyellow, green, red, reset
import ruptures as rpt

def weighted_profile(q, x, w, bins=512, range=None):
    with np.errstate(divide='ignore', invalid='ignore'):
        commArgs = dict(bins=bins, range=[np.min(x), np.max(x)] if range is None else range)
        qmed, edges = np.histogram(x, weights=q*w, **commArgs)
        ncelle, edges = np.histogram(x, weights=w, **commArgs)
        cell_center = 0.5 * (edges[1:] + edges[:-1])
        qmed /= ncelle
    return qmed, cell_center

def power_law(x, a, b):
    return a * x**b

# Color strings for console printing
green = "\033[92m"
red = "\033[91m"
yellow = "\033[93m"
reset = "\033[0m"

class PowerLawSegmenter:
    class segment:
        def __init__(self, nstart, nend, slope, intercept, data=None, index=None):
            self.nstart = nstart
            self.nend = nend
            self.slope = slope
            self.intercept = intercept
            self.data = data
            self.index = index
        @property 
        def size(self):
            # Return the number of points in the segment.
            return self.nend - self.nstart
        def set_data(self, data):
            self.data = data
        def __getitem__(self, index):
            return (self.nstart, self.nend, self.slope, self.intercept)[index]

    def __init__(self, max_delta=0.1, allowed_duplicates=3, min_segment_length=3, ntests=3, verbose=False):
        """
        Initialize the segmenter with specified parameters.
        
        Parameters:
            max_delta (float): Maximum allowed change in slope (in log–log space) before starting a new segment.
            allowed_duplicates (int): If the number of repeated (duplicate) points is below this, keep only one.
            min_segment_length (int): Minimum number of points per segment.
            ntests (int): Number of subsequent points to test when a candidate outlier is detected.
            verbose (bool): If True, extra output is printed.
        """
        self.verbose = verbose
        self.logger = setup_logger("PowerLawSegmenter")
        self.logger.setLevel(self.logger.loglevels.TRACE if self.verbose else self.logger.loglevels.INFO)
        self.max_delta = max_delta
        self.allowed_duplicates = allowed_duplicates
        self.min_segment_length = min_segment_length
        self.ntests = ntests
        self.segments = None
        self.coefficients = None   # Initialize coefficients attribute
        self.processed_time = None
        self.processed_radius = None
        self.nsegments = 0  # initialize number of segments

    def test_slope(self, icand, time, radius, start, slope_old):
        """
        (Legacy method from iterative segmentation; not used in the automated approach below.)
        Tests if the slope deviation persists.
        """
        condition_exceeded = True  # assume the deviation persists
        n = len(time)
        passed = 0
        for jtest in range(icand + 1, min(icand + 1 + self.ntests, n)):
            seg_log_time_test = np.concatenate((time[start:icand], time[icand+1:jtest+1]))
            seg_log_radius_test = np.concatenate((radius[start:icand], radius[icand+1:jtest+1]))
            slope_test, intercept_test = np.polyfit(seg_log_time_test, seg_log_radius_test, 1)
            DeltaSlopeTest = np.abs(slope_test - slope_old)
            if DeltaSlopeTest / self.max_delta <= 1.0:
                self.print(f"      |_____  test {jtest-icand}/{self.ntests} (i = {jtest}): (s_test,s_old) = ({slope_test:.3e},{slope_old:.3e}),  Δs_test/Δs_max  = {DeltaSlopeTest/self.max_delta:.4e}" +
                      " | " + green + " Passed! " + reset)
                passed += 1
                if passed > self.ntests//3:
                    self.print(yellow + f"            |_____ Outlier. Including tested points in the segment. Passed: {passed}/{self.ntests}" + reset)
                    condition_exceeded = False                
                    break  # deviation did not persist; candidate is an outlier
            else:
                self.print(f"      |_____  test {jtest-icand}/{self.ntests} (i = {jtest}): (s_test,s_old) = ({slope_test:.3e},{slope_old:.3e}),  Δs_test/Δs_max  = {DeltaSlopeTest/self.max_delta:.4e}" +
                      " | " + red + " Failed! " + reset)
        return condition_exceeded, jtest, slope_test, intercept_test

    def print(self, *args, **kwargs):
        if self.verbose:
            return print(*args, **kwargs)
        else:
            return

    def add_segment(self, *args, **kwargs):
        self.segments.append(self.segment(*args, **kwargs))
        self.nsegments += 1 

    def merge_similar_segments(self, slope_threshold=0.05, verbose=True, live_plot=False, fig=None, ax=None):
        """
        Iteratively merge adjacent segments if their slopes are similar (ignoring intercept).
        
        Parameters:
            slope_threshold (float): Relative difference threshold for slopes.
            verbose (bool): If True, print merge information.
            live_plot (bool): If True, update the provided live plot with merging progress.
            fig, ax: matplotlib figure and axis to update if live_plot is True.
            
        Returns:
            Merged list of segments.
        """
        if self.segments is None:
            raise ValueError("No segments to merge. Run segment_data() first.")
            
        merged = True
        while merged:
            merged = False
            new_segments = []
            i = 0
            while i < len(self.segments):
                if i < len(self.segments) - 1:
                    seg1 = self.segments[i]
                    seg2 = self.segments[i+1]
                    rel_diff_slope = abs(seg1.slope - seg2.slope) / max(abs(seg1.slope), abs(seg2.slope), 1e-8)
                    if rel_diff_slope < slope_threshold:
                        start_idx = seg1.nstart
                        end_idx = seg2.nend
                        new_data_time = self.processed_time[start_idx:end_idx]
                        new_data_radius = self.processed_radius[start_idx:end_idx]
                        coeffs = np.polyfit(np.log10(new_data_time), np.log10(new_data_radius), 1)
                        new_slope = coeffs[0]
                        new_intercept = coeffs[1]
                        merged_seg = self.segment(start_idx, end_idx, new_slope, new_intercept,
                                                  data=(new_data_time, new_data_radius), index=seg1.index)
                        new_segments.append(merged_seg)
                        if verbose:
                            self.logger.info(f"Merging segments {i} and {i+1}: new slope = {new_slope:.4e}")
                        i += 2
                        merged = True
                    else:
                        new_segments.append(self.segments[i])
                        i += 1
                else:
                    new_segments.append(self.segments[i])
                    i += 1
            self.segments = new_segments

            if live_plot and ax is not None:
                ax.cla()
                ax.loglog(self.processed_time, self.processed_radius, 'kx', ms=2.5, label="Data")
                colors = ['r', 'b', 'g', 'm', 'c']
                for j, seg in enumerate(self.segments):
                    seg_time = self.processed_time[seg.nstart:seg.nend]
                    fit_log_radius = seg.slope * np.log10(seg_time) + seg.intercept
                    fit_radius = 10 ** fit_log_radius
                    ax.plot(seg_time, fit_radius, color=colors[j % len(colors)],
                            lw=2, label=f"Seg {j+1}: s={seg.slope:.4e}")
                ax.set_xlabel("Time")
                ax.set_ylabel("Radius")
                ax.set_title("Merging Segments")
                ax.legend()
                plt.pause(0.5)
        self.nsegments = len(self.segments)
        return self.segments

    def segment_data(self, time, radius, tol=1e-1, max_candidates=4, live_plot=False,
                     merge_segments=True, merge_slope_threshold=0.05):
        """
        Segments the data into pieces that follow a power law (linear in log-log space)
        using change-point detection and model selection (AIC). Optionally plots the ongoing analysis
        and, if requested, merges segments with similar slopes.
        
        Parameters:
            time (array): 1D array of time values.
            radius (array): 1D array of corresponding radius values.
            tol (float): Tolerance for duplicate preprocessing.
            max_candidates (int): Maximum number of candidate change points to try.
            live_plot (bool): If True, updates an interactive plot during analysis.
            merge_segments (bool): If True, merge adjacent segments with similar slopes.
            merge_slope_threshold (float): Relative slope difference threshold for merging.
            
        Returns:
            segments (list of segment objects): Each segment contains (start_index, end_index, slope, intercept).
        """
        self.logger.info("Starting segmentation analysis.")
        # Preprocess duplicates.
        new_time, new_radius = self.preprocess_duplicates(time, radius, tol)
        self.processed_time = new_time
        self.processed_radius = new_radius
        self.logger.info(f"Preprocessed duplicates: Original length = {len(time)}, Processed length = {len(new_time)}")
        
        log_time = np.log10(new_time)
        log_radius = np.log10(new_radius)
        signal = np.column_stack((log_time, log_radius))
        n = len(new_time)
        
        candidate_models = []
        aic_values = []
        cps_list = []
        
        # Initialize live plotting if enabled.
        if live_plot:
            plt.ion()
            fig, ax = plt.subplots()
        
        self.logger.info("Examining candidate models for segmentation.")
        # Loop over candidate models with different numbers of change points.
        for num_cps in range(0, max_candidates+1):
            self.logger.trace(f"Trying candidate model with {num_cps} change point(s).")
            algo = rpt.Binseg(model="l2").fit(signal)
            cps = algo.predict(n_bkps=num_cps+1)
            self.logger.trace(f"Candidate change points: {cps}")
            cps_list.append(cps)
            
            total_rss = 0.0
            segments_params = []
            start_idx = 0
            for seg_idx, cp in enumerate(cps):
                if cp - start_idx < self.min_segment_length:
                    self.logger.trace(yellow + f"Segment {seg_idx+1} too short: start index {start_idx}, end index {cp} (length {cp-start_idx}).{red} Applying heavy penalty.{reset}" + reset)
                    total_rss += 1e6
                    segments_params.append((start_idx, cp, np.nan, np.nan))
                else:
                    seg_x = new_time[start_idx:cp]
                    seg_y = new_radius[start_idx:cp]
                    coeffs = np.polyfit(np.log10(seg_x), np.log10(seg_y), 1)
                    slope = coeffs[0]
                    intercept = coeffs[1]
                    pred_log = slope * np.log10(seg_x) + intercept
                    rss = np.sum((np.log10(seg_y) - pred_log)**2)
                    total_rss += rss
                    segments_params.append((start_idx, cp, slope, intercept))
                    self.print(f"                                |_____ Segment {seg_idx+1}: indices {start_idx} to {cp}, slope = {slope:.4e}, intercept = {intercept:.4e}, RSS = {rss:.4e}")
                start_idx = cp
            
            num_segments = len(cps)
            # Total free parameters: 2 per segment plus (num_segments - 1) breakpoints.
            k = num_segments * 2 + (num_segments - 1)
            aic = n * np.log(total_rss / n) + 2 * k
            aic_values.append(aic)
            self.logger.trace(f"Candidate model with {num_segments} segment(s): Total RSS = {total_rss:.4e}, k = {k}, AIC = {aic:.4e}")
            candidate_models.append(segments_params)
            
            if live_plot:
                ax.cla()
                ax.loglog(new_time, new_radius, 'kx', ms=2.5, label="Data")
                colors = ['r', 'b', 'g', 'm', 'c']
                for j, seg in enumerate(segments_params):
                    seg_start, seg_end, slope, intercept = seg
                    if np.isnan(slope):
                        continue
                    seg_time = new_time[seg_start:seg_end]
                    fit_log_radius = slope * np.log10(seg_time) + intercept
                    fit_radius = 10 ** fit_log_radius
                    ax.plot(seg_time, fit_radius, color=colors[j % len(colors)],
                            lw=2, label=f"Seg {j+1}: s={slope:.4e}")
                ax.set_xlabel("Time")
                ax.set_ylabel("Radius")
                ax.set_title(f"Candidate model with {num_segments} segment(s), AIC = {aic:.4e}")
                ax.legend()
                plt.pause(0.5)
        
        best_index = np.argmin(aic_values)
        best_segments_params = candidate_models[best_index]
        self.logger.info(green + f"Selected candidate model with {len(cps_list[best_index])} segments and AIC = {aic_values[best_index]:.4e}" + reset)
        
        self.nsegments = len(best_segments_params)
        self.segments = []
        self.coefficients = []  # Create the list to store coefficients
        for seg in best_segments_params:
            start, end, slope, intercept = seg
            self.segments.append(
                self.segment(start, end, slope, intercept, data=(new_time[start:end], new_radius[start:end]),
                             index=len(self.segments))
            )
            self.coefficients.append((slope, intercept))
        
        for seg in self.segments:
            self.logger.trace(f"Final segment: start index = {seg.nstart}, end index = {seg.nend}, slope = {seg.slope:.4e}, intercept = {seg.intercept:.4e}")
        
        # If merging is requested, perform merging and update the live plot if enabled.
        if merge_segments:
            self.logger.info("Merging similar segments based on slopes only...")
            self.segments = self.merge_similar_segments(slope_threshold=merge_slope_threshold,
                                                        verbose=True,
                                                        live_plot=live_plot,
                                                        fig=fig if live_plot else None,
                                                        ax=ax if live_plot else None)
            self.logger.info(f"After merging, total segments: {len(self.segments)}")
            if live_plot and ax is not None:
                ax.cla()
                ax.loglog(new_time, new_radius, 'kx', ms=2.5, label="Data")
                colors = ['r', 'b', 'g', 'm', 'c']
                for j, seg in enumerate(self.segments):
                    seg_time = new_time[seg.nstart:seg.nend]
                    fit_log_radius = seg.slope * np.log10(seg_time) + seg.intercept
                    fit_radius = 10 ** fit_log_radius
                    ax.plot(seg_time, fit_radius, alpha=0.8, lw=2,
                            label=f"Seg {j+1}: s={seg.slope:.4e}")
                ax.set_xlabel("Time")
                ax.set_ylabel("Radius")
                ax.set_title("Final Merged Segmentation")
                ax.legend()
                plt.pause(1)
        
        if live_plot:
            plt.ioff()
            plt.show()
        
        return self.segments

    def get_breaks(self, indexes=False):
        if not indexes:
            return [seg.data[0][-1] for seg in self.segments]
        else:
            return [seg.nend for seg in self.segments]
     
    def get_breaks_vals(self):
        return [seg.data[1][-1] for seg in self.segments]
     
    def preprocess_duplicates(self, time, radius, tol):
        """
        Collapse runs of consecutive duplicate points unless the run length exceeds allowed_duplicates.
        
        Parameters:
            time (array): 1D array of time values.
            radius (array): 1D array of corresponding radius values.
            tol (float): Tolerance for considering two points equal.
            
        Returns:
            new_time, new_radius (arrays): Processed arrays.
        """
        new_time = []
        new_radius = []
        i = 0
        n = len(time)
        while i < n:
            j = i + 1
            count = 1
            while j < n and np.isclose(radius[j], radius[i], atol=tol):
                count += 1
                j += 1
            if count > self.allowed_duplicates:
                new_time.append(time[i])
                new_radius.append(radius[i])
            else:
                avet = np.mean(time[i:j])
                aver = np.mean(radius[i:j])
                new_time.append(avet)
                new_radius.append(aver)
            i = j
        return np.array(new_time), np.array(new_radius)
    
    def plot_segments(self, ax):
        """
        Plots the original data and the piecewise power-law fits.
        """
        if self.segments is None:
            raise ValueError("No segments found. Please run segment_data() first.")
        colors = ['r', 'b', 'g', 'm', 'c']
        for j, seg in enumerate(self.segments):
            seg_time = self.processed_time[seg.nstart:seg.nend]
            fit_log_radius = seg.slope * np.log10(seg_time) + seg.intercept
            fit_radius = 10 ** fit_log_radius
            label = fr"$\alpha_{j} = {seg.slope:.2f}$",
            ax.plot(seg_time, fit_radius,
                    alpha=0.5, lw=2,
                    linestyle="dashed",
                    label=label,
                    color=colors[j % len(colors)])
    
    def select_segments_to_merge(self, indices):
        """Allow user to select specific segments to merge by providing indices."""
        merge_list = [(indices[i], indices[i+1]) for i in range(len(indices)-1)]
        self.merge_selected_segments(merge_list)
    
    def merge_selected_segments(self, merge_list):
        """Merge specific segments given a list of index pairs."""
        merged_segments = []
        merged_coeffs = []
        skip_indices = set()
        
        for i in range(len(self.segments)):
            if i in skip_indices:
                continue
            for merge_pair in merge_list:
                if i == merge_pair[0]:
                    start = self.segments[merge_pair[0]].nstart
                    end = self.segments[merge_pair[1]].nend
                    merged_slope = (self.coefficients[merge_pair[0]][0] + self.coefficients[merge_pair[1]][0]) / 2
                    merged_intercept = (self.coefficients[merge_pair[0]][1] + self.coefficients[merge_pair[1]][1]) / 2
                    merged_segments.append((start, end))
                    merged_coeffs.append((merged_slope, merged_intercept))
                    skip_indices.update(merge_pair)
                    break
            else:
                merged_segments.append((self.segments[i].nstart, self.segments[i].nend))
                merged_coeffs.append(self.coefficients[i])
        
        # Rebuild segments as objects
        new_segments = []
        for i, seg in enumerate(merged_segments):
            start, end = seg
            new_segments.append(self.segment(start, end, merged_coeffs[i][0], merged_coeffs[i][1],
                                             data=(self.processed_time[start:end], self.processed_radius[start:end]),
                                             index=i))
        self.segments = new_segments
        self.coefficients = merged_coeffs
    
    def merge_segment_with_next(self, i):
        """
        Merge the segment at index i with the segment at index i+1.
        Updates the segments and coefficients accordingly and prints merge information.
        """
        if self.segments is None:
            raise ValueError("No segments available. Run segment_data() first.")
        if i < 0 or i >= len(self.segments) - 1:
            raise IndexError("Invalid segment index. Must be less than len(segments)-1.")
        
        seg1 = self.segments[i]
        seg2 = self.segments[i+1]
        
        new_start = seg1.nstart
        new_end = seg2.nend
        new_data_time = self.processed_time[new_start:new_end]
        new_data_radius = self.processed_radius[new_start:new_end]
        
        # Re-fit the power law on the merged data (log-log linear regression)
        coeffs = np.polyfit(np.log10(new_data_time), np.log10(new_data_radius), 1)
        new_slope = coeffs[0]
        new_intercept = coeffs[1]
        
        # Create a new merged segment object
        new_seg = self.segment(new_start, new_end, new_slope, new_intercept,
                            data=(new_data_time, new_data_radius), index=seg1.index)
        
        # Replace the segment at index i with the merged segment
        self.segments[i] = new_seg
        # Print merging information
        self.logger.info(f"Merged segments at indices {i} and {i+1}:")
        self.logger.trace(f"           |   New boundaries: {new_start} to {new_end}")
        self.logger.trace(f"           |   New slope: {new_slope:.4e}")
        self.logger.trace(f"           |__ New intercept: {new_intercept:.4e}")
        
        # Remove the segment at index i+1 and update coefficients
        del self.segments[i+1]
        self.coefficients[i] = (new_slope, new_intercept)
        del self.coefficients[i+1]
        self.nsegments = len(self.segments)

