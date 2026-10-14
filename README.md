# IOPOS

IMU and Odometer Fusion Positioning for Automatic Parking Scenarios

## Overview

IOPOS is a lightweight positioning library for fusing IMU data with odometer measurements in automatic parking scenarios. It is designed for low-cost IMU systems and high-rate wheel/vehicle odometry, especially in GNSS-denied conditions such as underground garages and parking lots.

The project is intended for compact, practical, and research-friendly localization in constrained driving environments.

## Features

- IMU and odometer fusion for 2D vehicle positioning
- quaternion-based orientation estimation
- compact C/C++ implementation
- odometer-assisted motion updates
- adaptive bias handling and state updates
- trajectory visualization support with OpenCV

## Requirements

- CMake
- OpenCV
- C/C++ compiler
- standard math library support

## Quick Start

```bash
mkdir -p build
cd build
cmake ..
make
```

This builds the sample executable:

```bash
./test_imu_odo
```

## Project Layout

```text
.
├── CMakeLists.txt
├── README.md
├── LICENSE
├── iopos.c
├── iopos.h
├── matrix.c
├── quaternion.c
├── so3.c
├── test_iopos.cpp
├── data/
├── result1.png
├── result2.png
└── ...
```

## Example Workflow

The sample program reads IMU and odometer data, fuses them over time, and visualizes the estimated trajectory. It is meant to demonstrate the algorithm in automatic parking use cases.

## Notes

- The implementation is designed for efficiency and portability.
- It is suitable for embedded or research-oriented navigation experiments.
- This project focuses on low-cost fusion in parking environments rather than full GNSS/INS pipeline integration.

## Contact

- QQ: 1971129844, Su Jinglan
- Email: 1971129844@qq.com

## License

This project is released under the license included in the repository.
