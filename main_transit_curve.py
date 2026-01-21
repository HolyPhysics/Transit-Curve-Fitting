import numpy as np
import seaborn as sns
import emcee
import corner
import matplotlib.pyplot as plt
# import pandas as pd
from astropy.table import Table
from scipy.stats import norm
from scipy.optimize import minimize
from preprocessing import clean_import


time_container, flux_container, flux_error_container = clean_import() # declaring variables


class transit_object(object):

  def __init__(self, initial_guess: list[float], time: list[float], flux: list[float], flux_error: list[float]) -> None:
    self.initial_guess = initial_guess
    self.time = time
    self.flux = flux
    self.flux_error = flux_error

    self.samples_with_chain = None
    self.samples_without_chain = None

  def model(self, time, parameters):
    A_value, b_value, tc_value, w_value = parameters

    # piece_wise_condition = (tc_value - w_value/2) <= self.time &  self.time <= (tc_value + w_value/2)

    transit_function_condition = ((tc_value - w_value/2) <= time) &  (time <= (tc_value + w_value/2))
    mu = np.where(transit_function_condition, A_value - b_value, A_value)

    return mu

  def log_prior(self, parameters):

    A_value, b_value, tc_value, w_value = parameters
    # since tc_value, b_value, and w_value >=0, then we can choose these parameters to be uniform in the range 0(greater than 0) to x( here we set x to 10 for no special reason at all)
    # so they each have uniform probability equal to 1/10 ( or 1/x )

    if not 0 < A_value < 2:
      return -np.inf
    if not 0 < b_value < 1:
      return -np.inf
    if not min(self.time) < tc_value < max(self.time):
      return -np.inf
    if not 0 < w_value <  max(self.time) - min(self.time) + 1:
      return -np.inf

    return 0 # the weakly informative uniform prior

  def log_likelihood(self, parameters):

    A_value, b_value, tc_value, w_value = parameters

    mu = self.model(self.time, parameters)
    sigma = self.flux_error

    return np.sum(norm.logpdf(self.flux, mu, sigma))

  def transit_maximum_likelihood_estimation(self):

    def mle_log_likelihood(parameters):
      return -self.log_likelihood(parameters)

    maximized_ouputs = minimize(mle_log_likelihood, self.initial_guess, method="Nelder-Mead")

    # if maximized_ouputs.success == True:
    #   print(maximized_ouputs.message)
    # else:
    #   print(' Minimization unsuccesful!')

    return maximized_ouputs.x

  def log_posterior(self,parameters):

    log_prior = self.log_prior(parameters)
    log_likelihood = self.log_likelihood(parameters)
    if np.isinf(log_prior) or np.isinf(log_likelihood):
      return log_prior
    else:
      return log_prior + log_likelihood


  def transit_MCMC_sampler(self, burning_point = 9, number_of_steps=7500, number_of_walkers=50, number_of_dimensons = 4):

    # A_value, b_value, tc_value, w_value = self.transit_maximum_likelihood_estimation()

    def transit_MCMC_log_posterior(parameters):

      return self.log_posterior(parameters)

    # mle_values: list[float, ...] = self.transit_maximum_likelihood_estimation()

    # starting_point = np.array(mle_values) + np.random.randn(number_of_walkers, number_of_dimensons)

    # initialize walkers (specify initial guess range for all walkers)
    # I guess I don't have to always start close to the MLE as estimated from my mle function. Another much convenient method
    starting_guesses1 = np.random.uniform(0,1.9,(number_of_walkers,1)) # A_value
    starting_guesses2 = np.random.uniform(0,0.9,(number_of_walkers,1)) # b_value
    starting_guesses3 = np.random.uniform(min(self.time)+1 ,max(self.time) - 1,(number_of_walkers,1)) # tc_value
    starting_guesses4 = np.random.uniform(0,max(self.time) - min(self.time) ,(number_of_walkers,1)) # w_value

    starting_guesses = np.hstack((starting_guesses1,starting_guesses2,starting_guesses3,starting_guesses4))

    transit_sampler = emcee.EnsembleSampler(number_of_walkers, number_of_dimensons,transit_MCMC_log_posterior)

    transit_sampler.run_mcmc(starting_guesses, number_of_steps, progress=True)

    self.samples_with_chain = transit_sampler.get_chain()
    self.samples_without_chain = transit_sampler.get_chain(discard=burning_point, flat=True)
    # print(starting_guesses)

    # because it doesn't hurt to ensure that the samples_with_chain are big enough for discardgin of the burnt points we have this check
    if len(self.samples_with_chain) <= burning_point:
      print(' \n Careful not to burn all the points, Chidi! \n')
    else:
      print('\n All clear! Proceed! \n')

    return self.samples_with_chain, self.samples_without_chain

  def transit_trace_plot(self, burning_point=9):

    chains = self.samples_with_chain[:]

    figure, ax_main = plt.subplots(4,1, figsize=(10,8.5), sharex=True)

    parameter_names = ['Nontransit Stellar Flux', 'Transit Depth', 'Time at Midpoint', ' Transit Width']

    for values in range(4):

      ax_current = ax_main[values]

      ax_current.plot(chains[:,:,values], '-k', alpha=0.3)
      ax_current.set_ylabel(parameter_names[values])
      ax_current.set_ylim(-5,15)
      ax_current.axvline(burning_point, color='red', linestyle='--', label="Burning point")
      # ax_current.grid(True)

    ax_main[3].set_xlabel(' Iterations')
    plt.legend(loc="best")
    figure.tight_layout()



  def transit_corner_plot(self):

    parameter_names = ['Nontransit Stellar Flux', 'Transit Depth', 'Time at Midpoint', ' Transit Width']

    corner.corner(self.samples_without_chain,labels=parameter_names, show_titles=True)


  def transit_mean_model_fit(self, number_of_points = 300):

    A_value = self.samples_without_chain[:,0]
    b_value = self.samples_without_chain[:,1]
    tc_value = self.samples_without_chain[:,2]
    w_value = self.samples_without_chain[:,3]

    mean_parameters = [np.mean(A_value), np.mean(b_value), np.mean(tc_value), np.mean(w_value)]

    x_range = np.linspace(min(self.time), max(self.time),number_of_points)

    transit_mean_fit = self.model(x_range, mean_parameters)

    figure, ax_main = plt.subplots(figsize = (10,8.5))

    ax_main.errorbar(self.time, self.flux, self.flux_error, fmt='ok', label="Mock data", capsize=1.5, elinewidth=1, markersize=4, alpha = 0.5)

    ax_main.plot(x_range, transit_mean_fit, label="Mean model", color='orange', linestyle='-', alpha= 1)

    ax_main.set_xlabel(' time[arbitrary units] ')
    ax_main.set_ylabel(' flux[arbitrary units] ')
    ax_main.set_title(' \n Bayesian fitting of Transit Model to Data \n')
    plt.legend(loc="lower left")
    ax_main.grid(True)
    plt.show()
    figure.tight_layout()




if __name__ == "__main__":

    # Using emcee, create and initialize a sampler with 4 dimensions, 50 walkers and draw ~8000 samples from the posterior.
    # I found I had to really tweak the starting guesses of the four parameters
    # I don't known if my result is fine, but I didn't havw to tweak it and I made use of the MLE values.
    #running the transit object which contains the MCMC sampler
    transit_test_container = transit_object([1,1.4,3,1], time_container, flux_container, flux_error_container)
    sample_parameters = transit_test_container.transit_MCMC_sampler()


    # Plot the four chains to determine when they stabilize
    transit_test_container.transit_trace_plot()
    plt.show()


    # extract the samples (removing the burn-in)
    # already done in the transit_MCMC_sampler method in transit_object as defined above


    # Use corner.py to visualize the two-dimensional posterior
    transit_test_container.transit_corner_plot()
    plt.show()


    # A plot of the mean line
    transit_test_container.transit_mean_model_fit()
