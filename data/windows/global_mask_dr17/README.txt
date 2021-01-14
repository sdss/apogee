This folder contains the material for making the so called global mask for the possible DR17 synspec-run.
This mask masks out (small) parts of spectra where the linelist tuning did not succeed well for one reason or another. This could for example be that the code couldn't fit the Sun and Arcturus well at the same time, or the quoted uncertainty of the loggf-values is too low, and the code stops adjusting it.

Pixels with a deviation between synthetic and observed spectra greater than 0.03 (in normalized flux) for the sun or a deviation between synthetic and observed spectra greater than 0.05 (in normalized flux) for Arcturus are masked out. These differences are calculated in spectra smoothed to APOGEE-like resolutions. Regions affected by telluric lines in Arcturus are not masked, since there are some telluric residuals in the spectrum used.

In the plots, every page has two panels showing different wavelength regions.
The spectra with continuum around 2 are arcturus observed (black), arcturus synthetic (red) and a telluric observation (green).
The spectra with continuum around 1.2 are the sun observed (black) and synthetic (red).
All spectra have been convolved to APOGEE-like resolution. 
The color-coding in the background is the deviation of the synthetic spectrum from the observation: red means that the synthetic line is too strong and vice versa.
This means that when you see red-blue or blue-red in the arcturus-sun spectra — for example around 15014 Å where the background is blue in arcturus and red in the sun — the code is struggling to even out the fit between the two spectra.
When you see blue-blue as is the case just below 15180, we have a missing line in our list.
When you see red-red as is the case with the line around 16435 the synthetic line is too strong in both the sun and arcturus, and the linestrength is adjusted by Dmitry’s code to the minimum loggf allowed by the quoted uncertainty (in turn meaning that the uncertainty probably is larger than stated). 
At the top are the line-IDs from the Arcturus atlas.
At the bottom are the DR14 windows for the different elements (you can read the elements to the right if you zoom in on the screen. Not if you print the file, though…)
The final, global mask is marked in black around intensities 1.20-1.25 (some problematic regions found by the code, but outside of the APOGEE-region are marked in gray color).
The DR16-mask is marked around intensities 1.25-1.30

Observed Arcturus spectrum:                  arct.20190830.txt
Synthetic Arcturus spectrum:                 arct_nlte_20200921.rspec
Observed Sun spectrum:                       solar.20190819.txt
Synthetic Sun spectrum:                      sun_nlte_20200921.rspec
Plots of mask:                               linelist_20200921_synspec.pdf

