# IOPOS

IMU and Odometer Fusion Positioning for Automatic Parking Scenarios

## Overview

IOPOS is a lightweight positioning library for fusing inertial measurement unit (IMU) data with odometer information in automatic parking scenarios. It is designed for low-cost IMU systems and high-rate wheel-encoder or CAN-based odometry, with the goal of achieving robust and continuous position estimation in constrained environments such as underground parking garages, parking lots, and other structured navigation contexts.

The implementation emphasizes a compact and efficient formulation suitable for embedded or research-oriented applications. It combines IMU orientation estimation, motion-state handling, and odometer-assisted position updates into a unified localization framework.

## Motivation

Automatic parking systems often operate in environments where GNSS is unavailable or unreliable. In such cases, IMU and wheel-odometry data become essential for local positioning. However, raw IMU signals are subject to drift, while odometry alone may be insufficient under complex motion conditions. IOPOS addresses this by combining IMU propagation with odometer constraints and adaptive state updates, providing a practical and robust solution for low-cost navigation.

## Key Features

- IMU and odometer fusion for 2D positioning
- orientation estimation using complementary and quaternion-based filtering
- support for low-cost inertial sensors and high-rate wheel odometry
- adaptive handling of motion states and dynamic updates
- robust estimation in parking and constrained maneuvering scenarios
- compact C implementation suitable for integration into embedded systems
- visualization support for trajectory evaluation

## Core Functional Modules

The project includes the following core implementation files:

- `iopos.c` — main fusion and positioning algorithms
- `iopos.h` — declarations and public interfaces
- `matrix.c` — matrix-related computational utilities
- `quaternion.c` — quaternion operations and transformations
- `so3.c` — SO(3)-related mathematical routines
- `test_iopos.cpp` — example program demonstrating usage and trajectory evaluation

## Repository Layout

```text
.
├── CMakeLists.txt
├── LICENSE
├── README.md
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

## Dependencies

- CMake
- OpenCV
- C/C++ compiler with standard math support

## Build

From the repository root:

```bash
mkdir build
cd build
cmake ..
make
```

This builds the example executable:

- `test_imu_odo`

## Example Usage

The repository includes a sample test program, `test_iopos.cpp`, which demonstrates how to:

1. read IMU and odometer measurements,
2. initialize the fused positioning state,
3. update the system over time, and
4. visualize the estimated trajectory.

The example is intended to show the algorithmic workflow for trajectory estimation in automatic parking scenarios.

## Application Scope

IOPOS is suitable for:

- automatic parking systems
- underground parking localization
- low-cost INS/odometry fusion
- vehicle navigation in GNSS-denied environments
- embedded research prototypes and algorithm validation

## Contact

QQ: 1971129844, Su Jinglan

Email: 1971129844@qq.com

## License

This project is released under the MIT-style license included in the repository.

## Note

IOPOS is a compact and practical implementation for IMU and odometer fusion in parking scenarios, with emphasis on efficiency, portability, and ease of integration into vehicle navigation systems.
