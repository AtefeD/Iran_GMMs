# Ground-Motion Prediction Models for Iran

This repository provides MATLAB codes to implement a set of consistent GMMs for Iran. The GMMs are developed to predict three components of ground motions, namely, horizontal, vertical, and also vertical-to-horizontal ratio of PGA, PGV, and 5%-damped spectral acceleration (0.01–10 s). Each GMM set (i.e., H, V, or V/H) includes four consistent models based on different distance measures; Repi- and Rhyp-based model are associated with point-source, and RJB- and Rrup-based models are associated with finite-fault sources.
 
## References:
1- Darzi A, Zolfaghari MR, Cauzzi C, Fäh D (2019). An Empirical Ground Motion Model for Horizontal PGV, PGA and 5%-Damped Elastic Response Spectra (0.01-10 s) in Iran. Bulletin of the Seismological Society of America. 10.1785/0120180196. See [this](https://pubs.geoscienceworld.org/ssa/bssa/article/109/3/1041/569755/An-Empirical-Ground-Motion-Model-for-Horizontal)
2- Zolfaghari MR, Darzi A (2019). Ground-Motion Models for Predicting Vertical Components of PGA, PGV and 5%-damped Spectral Acceleration (0.01-10 s) in Iran. Bulletin of Earthquake Engineering, 10.1007/s10518-019-00623-2 . See [this](https://link.springer.com/article/10.1007/s10518-019-00623-2)
3- Zolfaghari MR, Darzi A (2019). A Prediction Model for Vertical-to-Horizontal Ratios of PGA, PGV and 5%-damped Response Spectra (0.01-10 s) for Iran. Journal of Seismology. 10.1007/s10950-019-09836-z. See [this](https://link.springer.com/article/10.1007/s10950-019-09836-z)

If you use the codes, please cite the relevant article appropriately.

In case of any questions or comments, do not hesitate to contact me:  atefe@hi.is
<br/>
<br/>

## Codes Description:  
-   <code>H-GMM/</code> 
    -   <code>H_GMM_Darzi19.m</code> a function that implements the GMM to predict horizontal PGA and spectral accelerations 
    -   <code>text files</code> period-dependent coefficients	
-   <code>V-GMM/</code>
    -   <code>V_GMM_ZD19.m</code>  a function that implements the GMM to predict vertical PGA and spectral accelerations
    -   <code>text files</code> period-dependent coefficients	
-   <code>V2H-GMM/</code>
    -   <code>V2H_GMM_ZD19.m</code> a function that implements the GMM to predict vertical-to-horizontal ratio of PGA and spectral accelerations
    -   <code>text files</code> period-dependent coefficients	
-   <code>visualization/</code>
    -   <code>plot_GMMs.m</code> plot to demonstrate H, V, and V2H GMMs
    -   <code>figures.jpg</code> output spectral accerelation predections derived from H, V, and V2H GMMs 	
    
