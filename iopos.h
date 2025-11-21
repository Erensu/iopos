#ifndef IOPOS_H
#define IOPOS_H

#ifdef __cplusplus
extern "C" {
#endif

#include <math.h>
#include <stdio.h>
#include <assert.h>
#include <string.h>

/* constants -----------------------------------------------------------------*/
#define SMALL_EPS         1E-10
#define PI                3.1415926535897932  /* pi */
#define SQRT(x)           ((x)<=0.0||(x)!=(x)?0.0:sqrt(x))
#define SQR(x)            ((x)*(x))
#define quatnorm(q)       quat_normalize_self(q)

#define IOPOS_MEASUREMENT_IMU  0
#define IOPOS_MEASUREMENT_ODO  1
#define IOPOS_MEASUREMENT_GEAR  2

#define IOPOS_NONE       0
#define IOPOS_INIT       1
#define IOPOS_UDODO_OK   2
#define IOPOS_UDODO_FAIL 3
#define IOPOS_REINIT     4
#define IOPOS_UDIMU      5

typedef union {             /* generic 3d vector */
    struct {double x,y,z;};
    double vec[3];
} vec3_t;

typedef union {             /* quaternion type */
    struct {double q0,q1,q2,q3;};
    struct {double w,x,y,z;};
    double vec[4];
} quat_t;

typedef union {             /* euler angle */
    struct {
        double yaw;
        double pitch;
        double roll;
    };
    double vec[3];
} euler_t;

typedef struct {          /* so3 struct data type */
    double unitq[4];      /* unit quaternion for handling so3 */
} so3_t;

typedef struct {
    int type;
    double timestamp;
    double val[16];        /* gears position (D:1,R=-1) */
} imu_odo_measurement_t;

typedef struct imu_complementary_filter {
    double timestamp;
    double gain_acc;      /* gain parameter for the complementary filter, belongs in [0, 1]. */
    double bias_alpha;    /* bias estimation gain parameter, belongs in [0, 1]. */

    int do_bias_estimation;
    int do_adaptive_gain;

    int initialized;
    int steady_state;
    double q0,q1,q2,q3;             /* the orientation as a Hamilton quaternion (q0 is the scalar). Represents
                                     * the orientation of the fixed frame wrt the body frame. */
    double wx_prev,wy_prev,wz_prev; /* bias in angular velocities */
    double wx_bias,wy_bias,wz_bias; /* bias in angular velocities */
    double gravity;
    double acc_iirflt[3];
    double gry_iirflt[3];
} imu_complementary_filter_t;

typedef struct madgwick_filter {
    double timestamp;
    double q0,q1,q2,q3;
    double beta;
} madgwick_filter_t;

typedef struct mayhony_filter {
    double timestamp;
    double q0,q1,q2,q3;
    double twoKp;
    double twoKi;
    double integralFBx;
    double integralFBy;
    double integralFBz;
} mayhony_filter_t;

typedef struct {
    int nimu;
    int uflg;
    double dt;
    imu_odo_measurement_t *odos;
    imu_odo_measurement_t *odoe;
    imu_odo_measurement_t imus[1024];
} imu_measurement_in_2odos_t;

typedef struct imu_odo_pos2d {
    int state;
    int udahrs;
    int direction;
    double timestamp;
    double timestamp_init;
    double p[2];
    double v[2];
    double yaw;
    double yawu;
    double yaw_var;
    double ba[3];
    double bg[3];

    double vehicle_length;
    double vehicle_width;
    double vehicle_wheelbase;
    double vehicle_frontoverhang;
    double vehicle_rearoverhang;

    madgwick_filter_t madgwickflt;
    mayhony_filter_t mayhonyflt;
    imu_measurement_in_2odos_t imu2odos;
} imu_odo_pos2d_t;

extern void iopos2dinit(imu_odo_pos2d_t *iopos2d);
extern void iopos2d(imu_odo_pos2d_t *iopos2d, const imu_odo_measurement_t *meas);

extern double *mat(int n, int m);
extern int *imat(int n, int m);
extern void matcpy(double *A, const double *B, int n, int m);
extern void matmul(const char *tr, int n, int k, int m, double alpha,
                   const double *A, const double *B, double beta, double *C);
extern int matinv(double *A, int n);
extern int filter(double *x, double *P, const double *H, const double *v,
                  const double *R, int n, int m);
extern double *eye(int n);
extern double *zeros(int n, int m);
extern void seteye(double* A,int n);
extern double norm(const double *a, int n);
extern void skewsym3(const double *ang, double *C);
extern void cross3(const double *a, const double *b, double *c);

extern void quat_init(quat_t *q, const vec3_t *acc, const vec3_t *mag);
extern void quat_init_axis(quat_t *q, double x, double y, double z, double a);
extern void quat_init_axis_v(quat_t *q, const vec3_t *v, double a);
extern void quat_to_axis(const quat_t *q, double *x, double *y, double *z, double *a);
extern void quat_to_axis_v(const quat_t *q, vec3_t *v, double *a);
extern void quat_rot_vec_self(vec3_t *v, const quat_t *q);
extern void quat_rot_vec(vec3_t *vo, const vec3_t *vi, const quat_t *q);
extern void quat_copy(quat_t *qo, const quat_t *qi);
extern double quat_len(const quat_t *q);
extern void quat_conj(quat_t *q_out, const quat_t *q_in);
extern void quat_to_euler(euler_t *euler, const quat_t *quat);
extern void quat_mul(quat_t *o, const quat_t *q1, const quat_t *q2);
extern void quat_add(quat_t *o, const quat_t *q1, const quat_t *q2);
extern void quat_scale(quat_t *o, const quat_t *q, double f);
extern void quat_normalize(quat_t *o, const quat_t *q);
extern void quat_normalize_self(quat_t *q);
extern void quat_to_rh_rot_matrix(const quat_t *q, double *m);
extern vec3_t *vec3_rot_axis_self(vec3_t *vo, double x, double y, double z,
                                  double angle);
extern void quat_from_u2v(quat_t *q, const vec3_t *u, const vec3_t *v,
                          const vec3_t *up);
extern double quat_dot(const quat_t *q1, const quat_t *q2);
extern quat_t *quat_nlerp(quat_t *qo, const quat_t *qfrom, const quat_t *qto,
                          double t);
extern quat_t *quat_slerp(quat_t *qo, const quat_t *qfrom, const quat_t *qto,
                          double t);
extern void quat_inv(const quat_t *q,quat_t *qi);
extern quat_t invquat(const quat_t *q);
extern void euler_to_quat(const euler_t *euler,quat_t *q);
extern quat_t *quat_apply_relative_yaw_pitch_roll(quat_t *q, double yaw, double pitch, double roll);
extern void rpy2quat(const double *rpy,quat_t *q);
extern void quat2dcm(const quat_t *q,double *C);
extern void dcm2quat(const double *C,quat_t *q);

extern void so3log(const double *C,double *phi);
extern void so3exp(const double *phi,double *C);
extern void so3jac(const double *phi,double *Jri);
extern void so3_vee(const double *Omg,double *omega);
extern void so3_hat(const double *omega,double *Omg);
extern void so3_log(const double *C,double *omega,double *theta);
extern void so3_exp(const double *omega,double *C);

extern void imu_complementary_filter(imu_complementary_filter_t *filter, const double *accl, const double *gyro, double dt);
extern void imu_complementary_filter2(imu_complementary_filter_t *filter, const double *accl, const double *gyro, double timestamp);

extern void madgwick_filter(madgwick_filter_t *mfliter, const double *accl, const double *gyro, const double *magn, double timestamp);
extern void mayhony_filter(mayhony_filter_t *mfliter, const double *accl, const double *gyro, const double *magn, double timestamp);

extern int iopos_imumeass_in_2odos(const imu_odo_measurement_t *meas, imu_measurement_in_2odos_t *imu2odos);
extern void iopos_imudv_in_2odos(const imu_measurement_in_2odos_t *imu2odos, double *dv);

#ifdef __cplusplus
}
#endif
#endif /* RTKLIB_H */