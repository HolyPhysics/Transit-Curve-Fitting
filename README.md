# Transit Curve Fitting
In this project, I run an MCMC analysis on an exoplanet dataset to determine the transit parameters A, b, t_c, and w.

Most known extrasolar planets (exoplanets) have been discovered using the transit method. A transit occurs when a planet passes between a star and its observer, causing a brief dip in the brightness of the parent star. This dimming can be seen in light curves – graphs showing light received over a period of time. When the exoplanet passes in front of the star, the light curve will show a dip in brightness. Transits can help determine a variety of different exoplanet characteristics. The size of the exoplanet’s orbit can be calculated from how long it takes to orbit once (the period), and the size of the planet itself can be calculated based on how much the star’s brightness lowered.


---

## Methods and Data

### Data Sources
- Data file is included in the file transit_dataset.dat

### Features 
- transit curve
- Parameters burn-in chain plots
- corner plots
- final transit fit

### Parameters
- A -- stellar flux outside transit
- b -- transit depth
- t_c -- time at midpoint
- w -- transit width

### Relationship Between Parameters 

Instead of using a physical model for modeling the exoplanet transit (a model which would have to include parameters such as the period of the planet, the ratio of its radius to the radius of the star, the impact parameter, the eccentricity of the orbit and other physics), we'll use a mock transit model specified by a simple piecwise function defined as follows:

𝜇(𝑡)={𝐴,𝐴−𝑏,(𝑡𝑐−𝑤/2)≤𝑡≤(𝑡𝑐+𝑤/2) otherwise  

where  𝐴  is the stellar flux outside of the transit,  𝑏  is the depth of the transit,  𝑡𝑐  is the time at midpoint and  𝑤  is the width of the transit. The inference procedure is exactly the same using a more realistic transit model except the function we plug into the likelihood is different and a lot more complex.

---

## Installation

### Requirements

- Python 3.x
- numpy
- matplotlib
- astropy
- astroML

### Installation Command

```bash
pip install numpy matplotlib astropy astroML scikit-learn
```

---

## Example Plots
