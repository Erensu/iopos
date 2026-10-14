# IOPOS — Developer Guide

IMU and Odometer Fusion Positioning for Automatic Parking — Detailed Developer README

## One-line summary
IOPOS is a compact C implementation that fuses IMU and odometer (wheel/CAN) measurements for robust 2D positioning in parking/low-speed scenarios.

## Goals of this document
Provide developers with precise information about code structure, APIs, data formats, build options, integration points, and how to run/extend the algorithms.

## Repository structure

- `CMakeLists.txt` — build script that compiles the example executable `test_imu_odo`.
- `iopos.h` — public header with data types and API declarations.
- `iopos.c` — main implementation file with the fusion pipeline (complementary filters, orientation filters, position update logic).
- `quaternion.c` — quaternion utilities and conversions.
- `matrix.c` — small matrix utilities used by the implementation.
- `so3.c` — SO(3)-related math routines.
- `test_iopos.cpp` — example program illustrating how to read CSV sensor logs and run the fusion.
- `data/` — placeholder for example IMU/odometry CSV files.
- `result1.png`, `result2.png` — sample outputs showing trajectory visualizations.

## Important types and APIs (inspect `iopos.h`)

Key types (names here are illustrative — check header for exact fields):
- `imu_odo_measurement_t` — container for a single sensor measurement (timestamp, type, value array)
- `imu_measurement_in_2odos_t` — buffer type to hold IMU samples between odometer ticks
- `imu_odo_pos2d_t` — main filter state: position (p), velocity (v), yaw, biases, internal filters (Madgwick/Mahony)

Primary functions:
- `void iopos2dinit(imu_odo_pos2d_t *state);` — initialize the state
- `void iopos2d(imu_odo_pos2d_t *state, const imu_odo_measurement_t *meas);` — feed one measurement (IMU or odometer) into the filter and update state
- `void imu_complementary_filter(...);` — attitude estimation with complementary/quaternion-based filters
- `void madgwick_filter(...);` / `void mayhony_filter(...);` — orientation filters used internally

The design follows a push-based API: push measurements with `iopos2d()` and read the updated state from `imu_odo_pos2d_t`.

## Data input format (as used by `test_iopos.cpp`)

- IMU CSV: timestamp (ns), gyro_x, gyro_y, gyro_z, acc_x, acc_y, acc_z, ...
- Odometer CSV: timestamp (ns), gear, speed_fl, speed_fr, speed_rl, speed_rr, ...

`test_iopos.cpp` parses CSV files where timestamps are in nanoseconds and converts to seconds (`timestamp = atof(val)/1e9`). The exact column order is visible in `test_iopos.cpp` parsing logic.

## Build and run

Requirements:
- CMake >= 3.10
- OpenCV (for visualization used by `test_iopos.cpp`)
- C/C++ compiler

Build:

```bash
mkdir -p build
cd build
cmake ..
make -j$(nproc)
```

Run example (after populating `data/` with IMU and odometer CSVs or adapting paths in `test_iopos.cpp`):

```bash
./test_imu_odo
```

The example will write `meas.txt` and `sols_odo.txt` and open an OpenCV window with trajectory visualization.

## How the algorithm works (brief)

- Measurements are buffered in `imu_measurement_in_2odos_t` until a new odometer tick is found.
- Attitude is estimated from IMU using complementary, Madgwick or Mahony filters.
- When odometer ticks arrive, displacement estimates (`dv`) are computed from IMU deltas and odometer differences, then used to update 2D position and yaw.
- Bias estimation uses stationary windows and simple IIR updates (see `iopos2d_ebias`).

For more details, read `iopos.c` where each step is implemented and commented.

## Extending and integrating

- To adapt to a new odometer format, modify parsing in `test_iopos.cpp` and the places where `imu_odo_measurement_t` is populated.
- To integrate into a larger navigation stack, create a thin wrapper that converts your sensor messages into `imu_odo_measurement_t` and calls `iopos2d()`.
- If you need 3D odometry or multi-axle vehicles, extend the state type and propagation equations accordingly.

## Testing and validation

- Use recorded sequences with ground truth to validate positions written to `sols_odo.txt`.
- Add unit tests for quaternion utilities, integrators and bias estimation routines.

## Performance

- Pure C implementation with few dynamic allocations; suitable for embedded deployment.
- Visualization uses OpenCV and can be disabled for headless embedded builds.

## Contact

- Su Jinglan — sujinglan@whu.edu.cn — QQ: 1971129844

## License

Check the repository `LICENSE` file. If none present, contact the author for licensing terms.
