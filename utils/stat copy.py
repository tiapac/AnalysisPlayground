import numpy as np
import matplotlib.pyplot as plt
from utils.logger import setup_logger
from utils.utils import printr, printg, printyellow, green, red, reset

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

class PowerLawSegmenter:
    class segment:
        def __init__(self, nstart, nend, slope, intercept, data = None, index = None):
            self.nstart      = nstart
            self.nend        = nend
            self.slope       = slope
            self.intercept   = intercept
            self.data        = data
            self.index       = index
        pass
        @property 
        def size(self): return self.nstart-self.nend
        def set_data(self, data): self.data = data
        def __getitem__(self, index):
            return  (self.nstart,
                     self.nend,
                     self.slope,
                     self.intercept)[index]
    class tmpseg:
        def __init__(self, time, radius):
            self.time   = time
            self.radius = radius
            
            pass
        
        def fit(self):
            slope, intercept = np.polyfit(self.time, self.radius, 1)
 
            return slope, intercept
        
        def add(self, time, radius):
            self.time.append(time)
            self.radius.append(radius)




    def __init__(self,  max_delta          = 0.1,
                        allowed_duplicates = 3,
                        min_segment_length = 3,
                        ntests             = 3):
        """
        Initialize the segmenter with specified parameters.
        
        Parameters:
            max_delta (float): Maximum allowed change in slope (in log–log space) before starting a new segment.
            allowed_duplicates (int): If the number of repeated (duplicate) points is below this, keep only one.
            min_segment_length (int): Minimum number of points per segment.
            ntests (int): Number of subsequent points to test when a candidate outlier is detected.
        """
        self.logger             = setup_logger("PowerLawSegmenter")
        self.logger.setLevel(self.logger.loglevels.TRACE)
        self.max_delta          = max_delta
        self.allowed_duplicates = allowed_duplicates
        self.min_segment_length = min_segment_length
        self.ntests             = ntests
        self.segments           = None
        self.processed_time     = None
        self.processed_radius   = None

    def test_slope(self, icand, time, radius, start, slope_old, ):   
        
        
        condition_exceeded = True  # assume the deviation persists
        # Check the next ntests-1 points (if available)       
        n = len(time)
        passed = 0
        for jtest in range(icand + 1, min(icand + 1 + self.ntests , n)):
            
            seg_log_time_test   = np.concatenate((time   [start:icand], time   [icand+1:jtest+1]))
            seg_log_radius_test = np.concatenate((radius [start:icand], radius [icand+1:jtest+1]))
            
            slope_test, intercept_test = np.polyfit(seg_log_time_test, seg_log_radius_test, 1)
            DeltaSlopeTest = np.abs(slope_test - slope_old)
            
            if DeltaSlopeTest / self.max_delta <= 1.0:
                print(f"      |_____  test {jtest-icand}/{self.ntests} (i = {jtest}): (s_test,s_old) = ({slope_test:.3e},{slope_old:.3e}),  Δs_test/Δs_max  = {DeltaSlopeTest/self.max_delta:.4e}" +
                   " | " + green + " Passed! "+reset)
            
                passed += 1
                if passed > self.ntests//3:
                    printyellow(f"            |_____ Outlier. Including tested points in the segment. Passed: {passed}/{self.ntests}")                
                    condition_exceeded = False                
                    break  # the deviation did not persist; candidate is outlier
            else:
                print(f"      |_____  test {jtest-icand}/{self.ntests} (i = {jtest}): (s_test,s_old) = ({slope_test:.3e},{slope_old:.3e}),  Δs_test/Δs_max  = {DeltaSlopeTest/self.max_delta:.4e}"+
                   " | " +  red + " Failed! "+reset)
            
        return condition_exceeded, jtest, slope_test, intercept_test
    def add_segment(self,*args,**kwargs):
        
        self.segments.append(self.segment(*args, **kwargs))
        self.self.nsegments += 1 

    def segment_data(self, time, radius, tol=1e-1):
        """
        Segments the data into pieces that follow a power law (linear in log-log space).
        
        The algorithm starts with a minimum number of points (min_segment_length) and then adds one point at a time.
        When a new point causes the slope to deviate by more than max_delta from the current segment's slope,
        it does not immediately end the segment. Instead, it checks the next ntests points:
          - If all of them confirm the deviation (i.e. the deviation persists), the segment is broken
            at the first offending point.
          - Otherwise, the candidate point is treated as an outlier and included in the segment.
        
        Parameters:
            time (array): 1D array of time values.
            radius (array): 1D array of corresponding radius values.
            
        Returns:
            segments (list of tuples): Each tuple is (start_index, end_index, slope, intercept) in the processed arrays.
        """
        # Preprocess the data to remove duplicate runs.
        max_segment_length = 10
        new_time, new_radius = self.preprocess_duplicates(time, radius, tol)
        log_time   = np.log10(new_time)
        log_radius = np.log10(new_radius)
        self.segments   = []
        n          = len(new_time)
        seg_start  = 0
        i          = seg_start + self.min_segment_length  # start testing once we have enough points


        trace = self.logger.trace
        info  = self.logger.info
        self.self.nsegments = 0
        newseg = False



        info ("Examining array.")
        info (f"\nLooking for segment #{self.self.nsegments+1}")
        seg_log_time_old         = log_time  [seg_start:i]
        seg_log_radius_old       = log_radius[seg_start:i]
        slope_old, intercept_old = np.polyfit(seg_log_time_old, seg_log_radius_old, 1)
        i += 1

        while i < n:

            if (i-seg_start>max_segment_length):
                self.add_segment(seg_start, i, slope_old, intercept_old, index = self.nsegments)
                i+=1
                print(f"|__   In {i} S is {slope_new:.4e} (old {slope_old:.4e}) | "+
                  f"Δs/Δs_max= {DeltaSlope / self.max_delta:.4e}")
            
                continue

            seg_log_time_new         = log_time  [seg_start:i]
            seg_log_radius_new       = log_radius[seg_start:i]
            



            slope_new, intercept_new = np.polyfit(seg_log_time_new, seg_log_radius_new, 1)
            DeltaSlope = np.abs(slope_new - slope_old)
            
            print(f"|__   In {i} S is {slope_new:.4e} (old {slope_old:.4e}) | "+
                  f"Δs/Δs_max= {DeltaSlope / self.max_delta:.4e}")
            
            if DeltaSlope / self.max_delta > 1.0:
                # A candidate outlier is detected at index i.
                candidate = i
                condition_exceeded, ichecked, slope_test, intercept_test = self.test_slope(candidate, log_time, log_radius, seg_start, slope_old)
                
                if condition_exceeded:

                    printg(f"      |_____ Slope changed. Cutting out a segment.")

                    # The deviation persists: end the segment at the candidate index.
                    SegmentLength = candidate - seg_start

                    if SegmentLength >= self.min_segment_length:

                        

                        self.add_segment(seg_start, candidate, slope_old, intercept_old, index = self.nsegments)
                        
                        seg_log_time_old         = log_time  [seg_start:candidate]
                        seg_log_radius_old       = log_radius[seg_start:candidate]
                        slope_old, intercept_old = float(slope_test), float(intercept_test)
                        
                        info(f"New segment {self.nsegments} found from {seg_start:3d} to {i:3d}. Slope is {slope_old}. ")
                        trace (f"\nLooking for a segment #{self.nsegments+1}")

                    seg_start = candidate
                    i         = seg_start + self.min_segment_length
                    # restart the while  loop with the new segment
                    continue  
                
                else:
                    # The deviation was temporary (an outlier): include this point in the segment.         
                    # remove outlier from stat 

                    printyellow(f"            |_____ Slope did not change. Adding {ichecked-candidate} point to  segment.")

                    i = int(ichecked)

            else:
                # go to the next data
                # print(f"\n checking next...")
                seg_log_time_old    = seg_log_time_new.copy()
                seg_log_radius_old  = seg_log_radius_new.copy()
                slope_old, intercept_old = float(slope_new), float(intercept_new)
                i += 1
                #newseg = False
            
            







        # Add the final segment if it has enough points.
        if n - seg_start >= self.min_segment_length:
            seg_log_time_final           = log_time[seg_start:]
            seg_log_radius_final         = log_radius[seg_start:]
            slope_final, intercept_final = np.polyfit(seg_log_time_final, seg_log_radius_final, 1)
            self.self.nsegments += 1
            self.add_segment(seg_start, n, slope_final, intercept_final, index = self.self.nsegments)
            newseg = True
        
        self.processed_time   = new_time
        self.processed_radius = new_radius
        
        for seg in self.segments:
            self.logger.info(f"  Start index: {seg[0]}, End index: {seg[1]}, Slope: {seg[2]:.2e}, Intercept: {seg[3]:.2e}")
        
        return self.segments
    



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
        new_time   = []
        new_radius = []
        i = 0
        n = len(time)
        while i < n:
            j = i + 1
            count = 1
            while (j < n) and np.isclose(radius[j], radius[i], atol=tol):
                count += 1
                j += 1
            if count > self.allowed_duplicates:
                new_time.append(time[i])
                new_radius.append(radius[i])
            else:
                avet = 0
                aver = 0
                for k in range(i, j):
                    avet += time[k]
                    aver += radius[k]
                if j - i > 0:
                    avet /= (j - i)
                    aver /= (j - i)
                else:
                    raise Exception("Error in duplicate preprocessing")
                new_time.append(avet)
                new_radius.append(aver)
            i = j
        return np.array(new_time), np.array(new_radius)
    
    def plot_segments(self):
        """
        Plots the original data and the piecewise power-law fits.
        """
        if self.segments is None:
            raise ValueError("No segments found. Please run segment_data() first.")
        plt.figure()
        plt.loglog(self.processed_time, self.processed_radius,
                   'kx', ms=2.5, label="Data")
        colors = ['r', 'b', 'g', 'm', 'c']
        for j, (start, end, slope, intercept) in enumerate(self.segments):
            seg_time = self.processed_time[start:end]
            fit_log_radius = slope * np.log10(seg_time) + intercept
            fit_radius = 10 ** fit_log_radius
            plt.plot(seg_time, fit_radius, alpha=0.5,
                     label=fr"$s_{j} = {slope:.4f}$",
                     color=colors[j % len(colors)])
        plt.xlabel("Time")
        plt.ylabel("Radius")
        plt.legend()
        plt.title("Piecewise Power-law Segmentation")
        plt.savefig("YYYYY.png", dpi=256)

# ---------------------
# Example Usage
# ---------------------
# Uncomment the following block to test the segmentation with synthetic data.
# if __name__ == '__main__':
#     # Create synthetic data with three power-law regimes.
#     time = np.linspace(1, 100, 300)
#     radius = np.piecewise(time,
#                           [time < 30, (time >= 30) & (time < 70), time >= 70],
#                           [lambda t: 1 * t**0.5,
#                            lambda t: 2 * t**0.8,
#                            lambda t: 4 * t**1.2])
#     
#     # Create an instance of the segmenter.
#     segmenter = PowerLawSegmenter(max_delta=0.01, allowed_duplicates=3, min_segment_length=5, ntests=3)
#     segments = segmenter.segment_data(time, radius)
#     
#     segmenter.plot_segments()
