[//]: # (Image References)

[image1]: ./src/radar_nis.png "nis radar"
[image2]: ./src/laser_nis.png "nis laser"


# Unscented Kalman Filter Project 

In this project an Unscented Kalman Filter was utilized to estimate the state of a moving object of interest with noisy lidar and radar measurements. Passing the project requires obtaining RMSE values that are lower that the tolerance outlined in the project rubric. 

The source code implemented to accomplish the project are src/ukf.cpp, src/ukf.h, src/tools.cpp, and src/tools.h

The program main.cpp handles the communication with the simulator via uWebSocket.

Here is the main protocol that main.cpp uses for uWebSocketIO in communicating with the simulator.


INPUT: values provided by the simulator to the c++ program

["sensor_measurement"] => the measurement that the simulator observed (either lidar or radar)


OUTPUT: values provided by the c++ program to the simulator

["estimate_x"] <= kalman filter estimated position x

["estimate_y"] <= kalman filter estimated position y

["rmse_x"]
["rmse_y"]
["rmse_vx"]
["rmse_vy"]

---

## Summary of Project Results

I fixed on the following set of parameters to achieve my low RMSE results:

* The standard deviation of the longitudinal acceleration std_a_ was set to 0.5 m/s^2.
* The process noise standard deviation of yaw acceleration std_yawdd was set to 1.0 rad/s^2.
* The main diagonal of the initial covariance matrix P_ was initialized to 1.0 and rest to zero.
* When initializing the state vector with a radar measurement I used rho_dot to get a first estimation for vx and vy.

The above parameters combined with an UKF approach for prediction and radar measurement updates, as well as a linear Kalman Filter approach for the lidar updates gets me the following RMSE values when using Dataset1 on the simulator:
["rmse_x"] = 0.0596
["rmse_y"] = 0.0872
["rmse_vx"] = 0.3343
["rmse_vy"] = 0.2210

Which are more than good enough to pass the project rubric. When using "Radar only" or "Laser only" measurements the accuracy of the estimations decreased significantly.

Below I will plot the results of the NIS calculations after each step for Radar and Lidar measurements, which can give an indication if the uncertainty or process noise of the UKF was set up too high or to low (compared to the real noise of the measurements and phenomena observed). The plots are compared to the statistical 95% reference:

![alt text][image1]
![alt text][image2]

As one can appreciate the majority of the of the NIS values stay under the 95% reference line without overestimation happening either. This means the chosen parameters do a good job, at the very least for the encountered datasets, and while further fine tuning for slight improvements might be possible we can be satisfied with the achieved results.


