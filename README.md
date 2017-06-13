#Unscented Kalman Filter Project
Self-Driving Car Engineer Nanodegree Program

Utilization of an Unscented Kalman Filter to estimate the state of a moving object of interest with noisy lidar and radar measurements. 

This project involves the Term 2 Simulator which can be downloaded [here](https://github.com/udacity/self-driving-car-sim/releases)

This repository includes two files that can be used to set up and intall [uWebSocketIO](https://github.com/uWebSockets/uWebSockets) for either Linux or Mac systems. For windows Docker or VMware is recommended. The simulator is not compatible with Ubuntu 16.04. 

Once the install for uWebSocketIO is complete, the main program can be built and ran by doing the following from the project top directory.

1. mkdir build
2. cd build
3. cmake ..
4. make
5. ./UnscentedKF


####Protocol for communication between Simulator and c++ program


INPUT: values provided by the simulator to the c++ program

["sensor_measurement"] => the measurment that the simulator observed (either lidar or radar)


OUTPUT: values provided by the c++ program to the simulator

["estimate_x"] <= kalman filter estimated position x
["estimate_y"] <= kalman filter estimated position y
["rmse_x"]
["rmse_y"]
["rmse_vx"]
["rmse_vy"]

---

## Other Important Dependencies

* cmake >= v3.5
* make >= v4.1
* gcc/g++ >= v5.4

## Basic Build Instructions

1. Clone this repo.
2. Make a build directory: `mkdir build && cd build`
3. Compile: `cmake .. && make`
4. Run it: `./UnscentedKF`.

## Editor Settings

The project compiles with cmake and make. However on `build with Eclipse` folder an Eclipse project is ready with all building instructions.

## Results

The simulator achieves RMSE values for dataset 1: RMSE: [0.0640,0.0879,0.3625,0.2086] for the variables [X,Y,VX,VY] respectively where X: the position in X axis, Y: the position in Y axis, VX: the velocity on X axis direction, VY: the velocity on Y axis direction

![dataset1](http://i.imgur.com/jznBnMF.png)

The RMSE values for the variables [X,Y,VX,VY] using the dataset 2 are: 
RMSE: [0.0615,0.0580,0.4764,0.1791]

![dataset2](http://i.imgur.com/NESOuJZ.png)