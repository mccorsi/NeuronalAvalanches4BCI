# NeuronalAvalanches4BCI

---
This repository contains the code and supporting documents associated with the following manuscript:

M.-C. Corsi*, P. Sorrentino*, D. Schwartz, N. George, L. Hugueville, A. E. Kahn, S. Dupont, D. S. Bassett, V. Jirsa, F. De Vico Fallani (2022). Measuring  brain critical dynamics to inform Brain-Computer Interfaces. Biorxiv. https://www.biorxiv.org/content/10.1101/2022.06.14.495887v1


*These authors contributed equally to this work.
 
---
## Authors:
* [Marie-Constance Corsi](https://marieconstance-corsi.netlify.app), Sorbonne Université, Institut du Cerveau
* [Pierpaolo Sorrentino](https://scholar.google.nl/citations?user=T1k8qBsAAAAJ&hl=en), Institut de Neuroscience des Systèmes, Aix-Marseille University
* Denis Schwartz, CERMEP, Lyon
* Nathalie George, Sorbonne Université, Institut du Cerveau
* Laurent Hugueville, Sorbonne Université, Institut du Cerveau
* Ari E. Kahn, University of Pennsylvania, Philadelphia
* Sophie Dupont, Sorbonne Université, Institut du Cerveau
* Danielle S. Bassett, University of Pennsylvania, Philadelphia
* Viktor Jirsa, Institut de Neuroscience des Systèmes, Aix-Marseille University
* Fabrizio De Vico Fallani, Sorbonne Université, Institut du Cerveau


---
## Abstract
Large-scale interactions among multiple brain regions manifest as bursts of activations, called neuronal avalanches, which reconfigure according to the task at hand and, hence, might constitute natural candidates to design brain-computer interfaces (BCI). To test this hypothesis, we compared source-reconstructed magneto/electroencephalography during resting-state and a motor imagery task performed within a BCI protocol. For each condition, we defined an individual avalanche transition matrix, tracking the probability that an avalanche would spread across any two regions. The edges whose transition probabilities significantly differed between conditions hinged selectively on premotor regions in all subjects, defining a topography related to the task. Furthermore, the individual differences in the transition probabilities for edges between pre/motor regions and parietal ones positively related to the individual task performance. Our results show that the patterns of propagation of large-scale perturbations are related to behavior and can be used to inform brain-computer interfaces.


## Code
This repository contains the code used to run the analysis performed and to plot the figures.



---
## Figures

### Figure 1 - Overview of the analysis 
![Fig. 1](./Figures_paper/Fig1.png)

*A. Subject-level analysis. B. Group-level analysis. C. Correlation analysis.*


### Figure 2 - Main results
![Fig. 2](./Figures_paper/Fig2.png)

*A. Edge-wise differences in transition probability. B. Edge-wise correlations. C. Node-wise differences in transition probabilities. D. Node-wise differences in transition probabilities.*


