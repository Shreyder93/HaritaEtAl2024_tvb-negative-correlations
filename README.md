## The role of inhibition in fMRI resting-state negative correlations 

## Shreyas Harita<sup>1, 2</sup>, Davide Momi<sup>2, 3</sup>, Zheng Wang<sup>2</sup>, Sorenza P. Bastiaens<sup>1, 2</sup>, John D. Griffiths<sup>1, 2, 4, 5, **</sup>

### Affiliations:  

1 = Institute of Medical Science, University of Toronto    
2 = Krembil Centre for Neuroinformatics, Centre for Addiction and Mental Health (CAMH), Toronto    
3 = Unit of Neuroimaging and Neurointervention, Department of Neurological and Neurosensorial Sciences, AOUS, 53100, Siena, Italy  
4 = Department of Physiology, University of Toronto    
5 = Department of Psychiatry, University of Toronto    
** = Corresponding Author  

### Highlights:  

- The reduced Wong-Wang (RWW) neural mass model was used to understand and study the physiological origin of negative correlations in rs-fMRI.
- The simulated data is able to partially capture the dynamics of the empirical rs-fMRI data. 
- Increasing the level of inhibition leads to an increase in the number of negative correlations. 
- The levels of excitation and inhibition have differential effects on model stability and negative correlations.
- Insights from this work can contribute towards fine-tuning model parameters and in future, aim to further validate findings using empirical data from clinical studies.
  

### Keywords:  

Neural-mass model, Negative correlation, Resting-state, Functional Connectivity, Structural Connectivity.

### This paper is currently under review at PLOS Computational Biology. You can access the bioRxiv version here: https://www.biorxiv.org/content/10.1101/2024.03.01.583030v2

### Abstract  

Resting-state brain fluctuations, observed via functional magnetic resonance imaging (fMRI), reveal non-random and structured patterns. These patterns, known as resting-state networks (RSNs), consist of positively correlated spontaneous fluctuations within specific regions. Additionally, RSNs exhibit negative correlations (NCs) between networks, indicating a systematic trade-off in activity. The balance between positive and negative correlations is linked to individual cognitive functions and personality traits, serving as a crucial indicator of healthy brain function and a potential predictor for therapeutic outcomes. Despite their clinical relevance, questions about the physiological origins and emergence of these NCs from the brain's structural connections remain unresolved. To uncover the origins of negative correlations (NCs), we employed the reduced Wong-Wang (RWW) brain network model to simulate fMRI data. By manipulating inhibitory parameters WI and λ, representing the external input scaling weight and feed-forward inhibition variable, respectively, we explored the impact of inhibition levels on NC emergence. Utilizing diffusion-weighted imaging data from the human connectome project for 24 subjects, we incorporated structural connectivity (SC) information in our simulations, determining input based on individualized SC extracted from a brain parcellation covering 200 regions. The model ran for 20 minutes, and we further examined the relationship between SC and NCs by analyzing the correlation between their respective row-wise sums, aiming to elucidate how SC influences NCs. Increasing the overall level of inhibition increased the number and spatial extent of NC in simulated data. In the RWW model, the mean NC across all subjects increased from 50.08 ± 92.77 to 13111.25 ± 5313.50 when WI was increased from 0.7 to 0.9. When λ was increased from 0 to 1, the mean NC increased from 50.08 ± 92.77 to 2161.00 ± 1859.02. However, we observed a reduction in PCs with this increase in NC in our simulated data. Additionally, the RWW simulations showed an inverse relationship between the row-wise sum of the SC weights and the row-wise of NCs. This study introduces a novel hypothesis regarding the effects of increased inhibition in brain network models on the presence of NCs in rs-fMRI data. This novel approach adds to the existing literature on brain network dynamics and provides new insights into the potential mechanisms underlying NCs in brain connectivity. We have identified several areas for future research, including further validation of the findings using empirical data from clinical studies, exploration of other parameters and factors that may influence resting-state brain connectivity, and clarification of the underlying mechanisms and interactions between inhibitory and excitatory dynamics in the brain.

There are 3 main folders in this GitHub repository:  

1. Data: This folder contains all pertinent information regarding the data used for the analysis presented in this study. Details about the HCP dataset can be found here.

2. Scripts: This folder contains all the relevant scripts/code used for the analysis in this study.

3. Figures: This folder contains the figures from the manuscript.  

Please see the README within each folder for more information.  

If you have any further questions, please contact Shreyas Harita at shreyas.harita@mail.utoronto.ca.    
