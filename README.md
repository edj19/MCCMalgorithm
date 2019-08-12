# MCCMalgorithm
Monte Carlo assisted Classical Methods of electron and gamma radiation

A new research in 2019 about radiation damage in gamma-ray irradiation.And relative formula is:

$$
N_{dpa}=\sum_k(n_k\sum_i N^e_{dpa,k}(E_i)\Phi(E_i,z)\Delta E_i)
$$




## MainCode Contents

This folder includes two Fortran code which to calculate the PKA and dpa cross section for different materials under its different displacement thresholds

$$
\sigma_{dpa}(E)=\sigma_{PKA}(E)\nu(T)
$$

where

$$
\sigma_{PKA}=\frac{\pi Z^2_a r^2_0}{\beta^4\gamma^2}\{(\frac{T_m}{T_d}-1)-\beta^2\ln(\frac{T_m}{T_d})+\pi\alpha\beta[2(\sqrt{\frac{T_m}{T_d}}-1)-\ln(\frac{T_m}{T_d})]\}
$$


 And another code is calculate the 

$N_d(dpa)$ to calculate the $\gamma$-ray radiation damage in the materials.
$$
N^e_{dpa}(E)=\int^E_{E_c}N_a \sigma_{dpa}(E')\frac{1}{(-dE'/dx)}dE'
$$

## Tools Contents

This folder is a post-process or some useful code to learn the theory or some physical phenomena more and more clear,and I will optimize this folder in my research tour.





## post_tools Contents

This folder is my first completely calculate the results and process the materials,I think it can help many researchers a lot.



If you have any questions,please contact me: edj17@lzu.edu.cn
