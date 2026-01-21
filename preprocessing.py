import numpy as np
import seaborn as sns
import matplotlib.pyplot as plt
# import pandas as pd
from astropy.table import Table


transit_data_file_path: str = "transit_dataset.dat"

try:
  transit_table: Table = Table.read(transit_data_file_path, format="ascii")

except Exception as err_transit:
  print(f' {type(err_transit).__name__} occured!')

# print(transit_table.columns)
# print(transit_table)


def clean_import() -> list(np.ndarray, np.ndarray, np.ndarray): 

    time_container: list[float, ...] = transit_table['Time_(hours)']
    flux_container: list[float, ...] = transit_table['Flux']
    flux_error_container: list[float, ...] = transit_table['Flux_err']

    print(f' The length of the time array is {len(time_container)} ')

    return time_container, flux_container, flux_error_container


def model(time, parameters):
    A_value, b_value, tc_value, w_value = parameters
    # piece_wise_condition = (tc_value - w_value/2) <= self.time &  self.time <= (tc_value + w_value/2)
    transit_function_condition = ((tc_value - w_value/2) <= time) &  (time <= (tc_value + w_value/2))
    mu = np.where(transit_function_condition, A_value - b_value, A_value)

    return mu


time_container, flux_container, flux_error_container = clean_import() # declaring variables


def light_curve_plot() -> None:

  mean_parameters = [1, 0.5, 4 , 2.9]

  transit_mean_fit = model(time_container, mean_parameters)
  
  figure, ax_main = plt.subplots(figsize = (10,8.5))
  ax_main.errorbar(time_container, flux_container, flux_error_container, fmt='ok', label="Mock data", capsize=1.5, elinewidth=1, markersize=4, alpha = 0.5)
  ax_main.plot(time_container, transit_mean_fit, label="Mean model", color='orange', linestyle='-', alpha= 1)
  ax_main.set_xlabel(' time[arbitrary units] ')
  ax_main.set_ylabel(' flux[arbitrary units] ')
  plt.legend(loc="lower left")
  ax_main.grid(True)
  figure.tight_layout()


  # my code here for the model plot. Write out the mechanism for the piecewise function as defined above.


if __name__ == "__main__":
    light_curve_plot()
    plt.show()


# print(time_container.data)
