# Analysis of Super-Kamiokande SK-IV dataset

This repository contains ROOT macros designed to read multi-variable tuples and apply cuts on data and Monte Carlo simulation. The sample used in this analysis has been divided according to several selection criteria. For both data and Monte Carlo simulation, the sample been reduced to neutrino interactions where the main vertex is at least 2 m from the ID walls. Events with a single Cherenkov ring originating from the charged lepton
were filtered by 400 MeV/c momentum threshold and separated by the particle identification based on the ring pattern. 

The zenith angle $\theta$ is defined as the angle of incoming neutrino with respect to the vertical direction at the detector. Therefore, $\cos\theta$ < 0 means upward-going neutrinos traveling through the Earth and reaching the detector from the bottom, while $\cos\theta$ > 0 means downward-going neutrinos reaching the detector from the top. 

During the SK-IV period (2008-2018), the Super-Kamiokande detector collected data equivalent to 50% of the pure-water data livetime from SK-I to SK-V (1996-2020). The starting data sample contains only FC events, with a total livetime of 3,244.4 days. The corresponding Monte Carlo (MC) simulation's livetime was 182,625.0 days, which is comparable to 500 years of detector operation. Such extensive statistics are necessary for a reliable MC simulation.

[Zenith angle distributions](https://github.com/mkondzie/analysis/blob/main/zenith.pdf)


A deficit of upward-going $\nu_{\mu}$ events was observed, while $\nu_e$ events are almost the same as predicted without oscillations. The missing muon neutrinos can be explained by $\nu_{\mu} \rightarrow \nu_{\tau}$ oscillations, where $\nu_{\mu}$ change flavor to $\nu_{\tau}$, which are undetectable in SK on an event-by-event basis. This evidence for neutrino oscillations was first published in 1998 by [Super Kamiokande](https://arxiv.org/abs/hep-ex/9807003). 


## Requirements
ROOT-6.30.06 or later
