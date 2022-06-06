# NeuronalAvalanches4BCI

---
This repository contains the code and supporting documents associated with the following manuscript:

M.-C. Corsi*, P. Sorrentino*, D. Schwartz, N. George, L. Hugueville, A. E. Kahn, S. Dupont, D. S. Bassett, V. Jirsa, (2022). Measuring  brain critical dynamics to inform Brain-Computer Interfaces. ArXiv. **TODO-add link**


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
Reconfiguring large-scale interactions among multiple brain regions, which underpin complex behavior, manifest themselves in bursts of activations, called neuronal avalanches, which can be tracked non-invasively as they expand across the brain. Responding to a new task requires brain regions to appropriately reconfigure their interactions, which might in turn affect the path of propagation of neuronal avalanches, providing a readily-accessible read out of the processes sustaining behavior. As such, they constitute  natural candidates to design brain-computer interfaces. To test this hypothesis, we used source-reconstructed magneto/electroencephalography, comparing the resting-state to a motor imagery task. For each experimental condition, we computed an individual avalanche transition matrix, to track the probability that an avalanche would spread across any two regions. Then, we selected those edges whose transition probabilities significantly differed between conditions, at the individual level and in the majority of the participants. We found a robust topography of the edges that were affected by the execution of the task, which mainly hinge upon the premotor regions. Furthermore, we related the individual differences to the task performance, showing that significant correlations are direct and predominantly involve edges connecting pre/motor regions to parietal ones. Our results show that the pattern of propagation of higher-order large-scale perturbations are related to behavior, and that they can be exploited for the design of optimized brain-computer interfaces.


## Code
This repository contains the code used to run the analysis performed and to plot the figures.



---
## Figures

### Figure 1 - Overview of the analysis 
![Fig. 1](./Figures_paper/Fig1.jpg)
*A. Subject-level analysis. B. Group-level analysis. C. Correlation analysis.*
**TODO**

### Figure 2 - Main results
![Fig. 2](./Figures_paper/Fig2.jpg)
*A. Edge-wise differences in transition probability. B. Edge-wise correlations. C. Node-wise differences in transition probabilities. D. Node-wise differences in transition probabilities.*
**TODO**

