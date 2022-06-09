# Abstract
This work presents experimental implementation of an improved single-column chromatographic (ISCC) separation process and its optimizing controller. A mixture of guaifenesin enantiomers has been used to evaluate the performance and integrity of (i) the ISCC process and its online monitoring system in open-loop experiments and (ii) the model predictive optimizing controller for closed-loop operation. The open-loop operation has been particularly aimed at assessing the accuracy and precision of the online monitoring system. In the closed-loop operation, the performance of the developed model predictive control (MPC) scheme has been tested for set point tracking and disturbance rejection with an objective function that reflects the process economics. The online optimal operating condition was also compared to the offline optimum condition obtained using a genetic algorithm. Results confirm that the optimizing controller is adequate to operate and maintain the ISCC process at an optimal operating point while fulfilling the product requirements.

N.B.: The codes for this ISCC process were developed by Dr. Bijan Medi. I have amended some parts of this code to implement the offline multi-objective optimization algorithms, which acted as the building block for the online monitoring system for the MPC. 

The details can be found in the following published journals:

https://www.sciencedirect.com/science/article/pii/S0021967312001690

https://pubs.acs.org/doi/full/10.1021/acs.iecr.5b00553

## Conceptual development of fraction collection scheme

![image](https://user-images.githubusercontent.com/62577662/172882857-df806c5c-7e64-438a-8658-1c5dc65e2285.png)

## Process design of ISCC

![image](https://user-images.githubusercontent.com/62577662/172882875-aa4a544a-41b3-4a71-8322-d90d7874f4b6.png)

## Cuttimes for different products collection

![image](https://user-images.githubusercontent.com/62577662/172882955-f60b6713-56dc-49a0-a36e-08180d03d3b6.png)

## MPC strategy

![image](https://user-images.githubusercontent.com/62577662/172883002-b4d2ae91-5de2-40b0-8821-cd53c9cae223.png)
![image](https://user-images.githubusercontent.com/62577662/172883082-9ce515b7-98e0-4634-949f-4583d4e1dcd3.png)


