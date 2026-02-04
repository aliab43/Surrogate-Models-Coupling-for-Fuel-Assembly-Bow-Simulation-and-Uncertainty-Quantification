import numpy as np
from scipy.interpolate import interp1d
from scipy.optimize import minimize_scalar

class GridClampingUpdater:
    def __init__(self, initial_Grid_Clamping, fluence_increment=0.277e25):
        """
        Initialize the GridClampingUpdater with an initial Grid_Clamping value and fluence increment.
        
        :param initial_Grid_Clamping: The initial grid clamping (Grid_Clamping initiale) value.
        :param fluence_increment: The amount of fluence increment (default is 0.277e25).
        """
        # Data for interpolation
        self.fluence_data = np.array([0, 0.1e25, 0.2e25, 0.3e25, 0.4e25, 0.5e25, 0.7e25, 2e25, 4e25, 1e26])
        self.Grid_Clamping_data = np.array([1, 0.65, 0.49, 0.35, 0.30, 0.24, 0.19, 0.11, 0.04, 0.01]) * 40
        self.interp_func = interp1d(
            self.fluence_data, 
            self.Grid_Clamping_data, 
            kind='linear', 
            bounds_error=False, 
            fill_value=(self.Grid_Clamping_data[0], self.Grid_Clamping_data[-1])
        )
        
        self.initial_Grid_Clamping = initial_Grid_Clamping
        self.fluence_increment = fluence_increment

    def find_fluence_for_Grid_Clamping(self, Grid_Clamping):
        """
        Find the fluence corresponding to a given Grid_Clamping value using interpolation.
        
        :param Grid_Clamping: The clamping force (Grid_Clamping) to find the fluence for.
        :return: The fluence corresponding to the Grid_Clamping.
        """
        cost_function = lambda x: np.abs(self.interp_func(x) - Grid_Clamping)
        result = minimize_scalar(cost_function, bounds=(min(self.fluence_data), max(self.fluence_data)), method='bounded')
        return result.x

    def update_Grid_Clamping(self):
        """
        Update the Grid_Clamping value after applying the fluence increment.
        
        :return: A tuple (updated_fluence, updated_Grid_Clamping).
        """
        # Step 1: Find initial fluence corresponding to the initial Grid_Clamping
        initial_fluence = self.find_fluence_for_Grid_Clamping(self.initial_Grid_Clamping)
        
        # Step 2: Increment the fluence
        updated_fluence = initial_fluence + self.fluence_increment
        
        # Step 3: Find the updated Grid_Clamping corresponding to the updated fluence
        updated_Grid_Clamping = self.interp_func(updated_fluence)
        if updated_Grid_Clamping < 1:
            updated_Grid_Clamping = 1
        return  updated_Grid_Clamping



