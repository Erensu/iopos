#include "iopos.h"

static double iopos_vehicle_angular_velocity(imu_odo_pos2d_t *iopos, const imu_odo_measurement_t *meas)
{
    double vfr=meas->val[0];
    double vfl=meas->val[1];
    double vrr=meas->val[2];
    double vrl=meas->val[3];

    double w=(vrr-vrl)/iopos->vehicle_width;
    if (fabs(w)<1E-3) return 0.0;

    if (w<0.0) {
        double R=fabs(vrr/w);
        double L1=sqrt(R*R+iopos->vehicle_wheelbase*iopos->vehicle_wheelbase);
        double L2=sqrt((R+iopos->vehicle_width)*(R+iopos->vehicle_width)+iopos->vehicle_wheelbase*iopos->vehicle_wheelbase);
        double L3=R;
        double L4=R+iopos->vehicle_width;
        double LL=L1*L1+L2*L2+L3*L3+L4*L4;

        double W=L1/LL*vfr+L2/LL*vfl+L3/LL*vrr+L4/LL*vrl;
        return W;
    }
    else {
        double R=fabs(vrl/w);
        double L1=sqrt(R*R+iopos->vehicle_wheelbase*iopos->vehicle_wheelbase);
        double L2=sqrt((R+iopos->vehicle_width)*(R+iopos->vehicle_width)+iopos->vehicle_wheelbase*iopos->vehicle_wheelbase);
        double L3=R;
        double L4=R+iopos->vehicle_width;
        double LL=L1*L1+L2*L2+L3*L3+L4*L4;

        double W=L1/LL*vfl+L2/LL*vfr+L3/LL*vrl+L4/LL*vrr;
        return W;
    }
}
/* q_acc is the quaternion obtained from the acceleration vector
 * representing the orientation of the Global frame wrt the Local frame with
 * arbitrary yaw (intermediary frame). q3_acc is defined as 0.
 * */
static void imu_compleflt_init_q_acc(imu_complementary_filter_t *filter, const double *accl)
{
    double ax=accl[0]/norm(accl,3);
    double ay=accl[1]/norm(accl,3);
    double az=accl[2]/norm(accl,3);

    if (az>=0.0) {
        filter->q0=sqrt((az+1)*0.5);
        filter->q1=-ay/(2.0*filter->q0);
        filter->q2= ax/(2.0*filter->q0);
        filter->q3=0.0;
    }
    else {
        double X=sqrt((1.0-az)*0.5);
        filter->q0=-ay/(2.0*X);
        filter->q1=X;
        filter->q2=0.0;
        filter->q3=ax/(2.0*X);
    }
}

static int imu_compleflt_checkstate(imu_complementary_filter_t *filter, double ax, double ay, double az, double wx, double wy, double wz)
{
    const double kAngularVelocityThreshold=0.2;
    const double kAccelerationThreshold=0.1;
    const double kDeltaAngularVelocityThreshold=0.01;
    double acc_magnitude=sqrt(ax*ax+ay*ay+az*az);

    if (fabs(acc_magnitude-filter->gravity)>kAccelerationThreshold) return 0;

    if (fabs(wx-filter->wx_prev)>kDeltaAngularVelocityThreshold||
        fabs(wy-filter->wy_prev)>kDeltaAngularVelocityThreshold||
        fabs(wz-filter->wz_prev)>kDeltaAngularVelocityThreshold)
        return 0;

    if (fabs(wx-filter->wx_bias)>kAngularVelocityThreshold||
        fabs(wy-filter->wy_bias)>kAngularVelocityThreshold||
        fabs(wz-filter->wz_bias)>kAngularVelocityThreshold)
        return 0;

    return 1;
}

static void imu_compleflt_estiamte_biases(imu_complementary_filter_t *filter, const double *accl, const double *gyro)
{
    filter->steady_state=imu_compleflt_checkstate(filter,accl[0],accl[1],accl[2],gyro[0],gyro[1],gyro[2]);

    if (filter->steady_state) {
        filter->wx_bias+=filter->bias_alpha*(gyro[0]-filter->wx_bias);
        filter->wy_bias+=filter->bias_alpha*(gyro[1]-filter->wy_bias);
        filter->wz_bias+=filter->bias_alpha*(gyro[2]-filter->wz_bias);
    }
    filter->wx_prev=gyro[0];
    filter->wy_prev=gyro[1];
    filter->wz_prev=gyro[2];
}

static void normz_quat(double* q0, double* q1, double* q2, double* q3)
{
    double norm=sqrt(*q0**q0+*q1**q1+*q2**q2+*q3**q3);
    *q0/=norm;
    *q1/=norm;
    *q2/=norm;
    *q3/=norm;
}

static void imu_compleflt_prediction(imu_complementary_filter_t *filter, const double *gyro, double dt,
                                     double* q0_pred, double* q1_pred, double* q2_pred, double* q3_pred)
{
    double wx_unb=gyro[0]-filter->wx_bias;
    double wy_unb=gyro[1]-filter->wy_bias;
    double wz_unb=gyro[2]-filter->wz_bias;

    *q0_pred=filter->q0+0.5*dt*( wx_unb*filter->q1+wy_unb*filter->q2+wz_unb*filter->q3);
    *q1_pred=filter->q1+0.5*dt*(-wx_unb*filter->q0-wy_unb*filter->q3+wz_unb*filter->q2);
    *q2_pred=filter->q2+0.5*dt*( wx_unb*filter->q3-wy_unb*filter->q0-wz_unb*filter->q1);
    *q3_pred=filter->q3+0.5*dt*(-wx_unb*filter->q2+wy_unb*filter->q1-wz_unb*filter->q0);

    double norm=sqrt(*q0_pred**q0_pred+*q1_pred**q1_pred+*q1_pred**q1_pred+*q1_pred**q1_pred);
    *q0_pred/=norm;
    *q1_pred/=norm;
    *q2_pred/=norm;
    *q3_pred/=norm;
}

static void rotvec_quat(double x, double y, double z, double q0, double q1, double q2, double q3,
                        double* vx, double* vy, double* vz)
{
    *vx=(q0*q0+q1*q1-q2*q2-q3*q3)*x+2*(q1*q2-q0*q3)*y+2*(q1*q3+q0*q2)*z;
    *vy=2*(q1*q2+q0*q3)*x+(q0*q0-q1*q1+q2*q2-q3*q3)*y+2*(q2*q3-q0*q1)*z;
    *vz=2*(q1*q3-q0*q2)*x+2*(q2*q3+q0*q1)*y+(q0*q0-q1*q1-q2*q2+q3*q3)*z;
}

static void scale_quat(double gain, double* dq0, double* dq1, double* dq2, double* dq3)
{
    if (*dq0<0.0) {
        /* Slerp (Spherical linear interpolation) */
        double angle=acos(*dq0);
        double A=sin(angle*(1.0-gain))/sin(angle);
        double B=sin(angle*gain)/sin(angle);
        *dq0=A+B** dq0;
        *dq1=B**dq1;
        *dq2=B**dq2;
        *dq3=B**dq3;
    } else {
        /* Lerp (Linear interpolation) */
        *dq0=(1.0-gain)+gain**dq0;
        *dq1=gain**dq1;
        *dq2=gain**dq2;
        *dq3=gain**dq3;
    }
    normz_quat(dq0,dq1,dq2,dq3);
}

static void multi_quat(double p0, double p1, double p2, double p3,
                       double q0, double q1, double q2, double q3,
                       double* r0, double* r1, double* r2, double* r3)
{
    *r0=p0*q0-p1*q1-p2*q2-p3*q3;
    *r1=p0*q1+p1*q0+p2*q3-p3*q2;
    *r2=p0*q2-p1*q3+p2*q0+p3*q1;
    *r3=p0*q3+p1*q2-p2*q1+p3*q0;
}

static void imu_compleflt_acc_correction(const double *accl, double p0, double p1, double p2, double p3,
                                         double* dq0, double* dq1, double* dq2, double* dq3)
{
    /* normalize acceleration vector */
    double norma=sqrt(accl[0]*accl[0]+accl[1]*accl[1]+accl[2]*accl[2]);
    double ax=accl[0]/norma;
    double ay=accl[1]/norma;
    double az=accl[2]/norma;

    /* acceleration reading rotated into the world frame by the inverse
     * predicted quaternion (predicted gravity) */
    double gx,gy,gz;
    rotvec_quat(ax,ay,az,p0,-p1,-p2,-p3,&gx,&gy,&gz);

    /* delta quaternion that rotates the predicted gravity into the real
     * gravity */
    *dq0=sqrt((gz+1)*0.5);
    *dq1=-gy/(2.0**dq0);
    *dq2= gx/(2.0**dq0);
    *dq3=0.0;
}

static double imu_compleflt_adaptive_gain(imu_complementary_filter_t *filter, double alpha, const double *accl)
{
    double a_mag=sqrt(accl[0]*accl[0]+accl[1]*accl[1]+accl[2]*accl[2]);
    double error=fabs(a_mag-filter->gravity)/filter->gravity;
    double factor;
    double error1=0.1;
    double error2=0.2;
    double m=1.0/(error1-error2);
    double b=1.0-m*error1;

    if      (error<error1) factor=1.0;
    else if (error<error2) factor=m*error+b;
    else factor = 0.0;

    return factor*alpha;
}

static void imu_compleflt_autoscale_gravity(imu_complementary_filter_t *filter, const double *accl)
{
    const double alpha=0.1;
    double ga=norm(accl,3);

    if (fabs(ga-filter->gravity)<0.1) {
        filter->gravity=filter->gravity*(1.0-alpha)+alpha*ga;
    }
}

static void imu_compleflt_accflt(imu_complementary_filter_t *filter, const double *accl, double dt)
{
    const double acc_flt_constdt=0.15;
    double iir_alpha=1.0-exp(-dt/acc_flt_constdt);

    if (norm(filter->acc_iirflt,3)==0.0) {
        filter->acc_iirflt[0]=accl[0];
        filter->acc_iirflt[1]=accl[1];
        filter->acc_iirflt[2]=accl[2];
    }
    filter->acc_iirflt[0]=iir_alpha*accl[0]+(1.0-iir_alpha)*filter->acc_iirflt[0];
    filter->acc_iirflt[1]=iir_alpha*accl[1]+(1.0-iir_alpha)*filter->acc_iirflt[1];
    filter->acc_iirflt[2]=iir_alpha*accl[2]+(1.0-iir_alpha)*filter->acc_iirflt[2];
}

extern void imu_complementary_filter(imu_complementary_filter_t *filter, const double *accl, const double *gyro, double dt)
{
    imu_compleflt_accflt(filter,accl,dt);

    if (filter->initialized==0) {
        imu_compleflt_init_q_acc(filter,filter->acc_iirflt);
        filter->initialized=1;
        return;
    }
    /* bias estimation */
    if (filter->do_bias_estimation) {
        imu_compleflt_estiamte_biases(filter,filter->acc_iirflt,gyro);
    }
    if (filter->steady_state) {
        imu_compleflt_autoscale_gravity(filter,filter->acc_iirflt);
    }
    /* update prediction */
    double q0_pred,q1_pred,q2_pred,q3_pred;
    imu_compleflt_prediction(filter,gyro,dt,&q0_pred,&q1_pred,&q2_pred,&q3_pred);

    /* update correction */
    double dq0_acc,dq1_acc,dq2_acc,dq3_acc;
    imu_compleflt_acc_correction(filter->acc_iirflt,q0_pred,q1_pred,q2_pred,q3_pred,&dq0_acc,&dq1_acc,&dq2_acc,&dq3_acc);

    double gain=filter->gain_acc;
    if (filter->do_adaptive_gain) {
        gain=imu_compleflt_adaptive_gain(filter,filter->gain_acc,filter->acc_iirflt);
    }
    /* update attitude */
    scale_quat(gain,&dq0_acc,&dq1_acc,&dq2_acc,&dq3_acc);
    multi_quat(q0_pred,q1_pred,q2_pred,q3_pred,dq0_acc,dq1_acc,dq2_acc,dq3_acc,&filter->q0,&filter->q1,&filter->q2,&filter->q3);
    normz_quat(&filter->q0,&filter->q1,&filter->q2,&filter->q3);
}

extern void imu_complementary_filter2(imu_complementary_filter_t *filter, const double *accl, const double *gyro, double timestamp)
{
    if (filter->timestamp==0) {
        filter->timestamp=timestamp;
    }
    imu_complementary_filter(filter,accl,gyro,timestamp-filter->timestamp);
    filter->timestamp=timestamp;
}

static void madgwick_filter_imu(madgwick_filter_t *mfliter, const double *accl, const double *gyro, double dt)
{
    double gx=gyro[0],gy=gyro[1],gz=gyro[2];
    double ax=accl[0],ay=accl[1],az=accl[2];
    double recipNorm;
    double s0,s1,s2,s3;
    double qDot1,qDot2,qDot3,qDot4;
    double _2q0,_2q1,_2q2,_2q3,_4q0,_4q1,_4q2,_8q1,_8q2,q0q0,q1q1,q2q2,q3q3;
    double q0=mfliter->q0;
    double q1=mfliter->q1;
    double q2=mfliter->q2;
    double q3=mfliter->q3;

    /* Rate of change of quaternion from gyroscope */
    qDot1=0.5*(-q1*gx-q2*gy-q3*gz);
    qDot2=0.5*(q0*gx+q2*gz-q3*gy);
    qDot3=0.5*(q0*gy-q1*gz+q3*gx);
    qDot4=0.5*(q0*gz+q1*gy-q2*gx);

    /* Compute feedback only if accelerometer measurement valid (avoids NaN in accelerometer normalisation) */
    if (!((ax==0.0)&&(ay==0.0)&&(az==0.0))){

        /* Normalise accelerometer measurement */
        recipNorm=1.0/sqrt(ax*ax+ay*ay+az*az);
        ax*=recipNorm;
        ay*=recipNorm;
        az*=recipNorm;

        /* Auxiliary variables to avoid repeated arithmetic */
        _2q0=2.0*q0;
        _2q1=2.0*q1;
        _2q2=2.0*q2;
        _2q3=2.0*q3;
        _4q0=4.0*q0;
        _4q1=4.0*q1;
        _4q2=4.0*q2;
        _8q1=8.0*q1;
        _8q2=8.0*q2;
        q0q0=q0*q0;
        q1q1=q1*q1;
        q2q2=q2*q2;
        q3q3=q3*q3;

        /* Gradient decent algorithm corrective step */
        s0=_4q0*q2q2+_2q2*ax+_4q0*q1q1-_2q1*ay;
        s1=_4q1*q3q3-_2q3*ax+4.0*q0q0*q1-_2q0*ay-_4q1+_8q1*q1q1+_8q1*q2q2+_4q1*az;
        s2=4.0*q0q0*q2+_2q0*ax+_4q2*q3q3-_2q3*ay-_4q2+_8q2*q1q1+_8q2*q2q2+_4q2*az;
        s3=4.0*q1q1*q3-_2q1*ax+4.0*q2q2*q3-_2q2*ay;
        recipNorm=1.0/sqrt(s0*s0+s1*s1+s2*s2+s3*s3);
        s0*=recipNorm;
        s1*=recipNorm;
        s2*=recipNorm;
        s3*=recipNorm;

        /* Apply feedback step */
        qDot1-=mfliter->beta*s0;
        qDot2-=mfliter->beta*s1;
        qDot3-=mfliter->beta*s2;
        qDot4-=mfliter->beta*s3;
    }
    /* Integrate rate of change of quaternion to yield quaternion */
    q0+=qDot1*(dt);
    q1+=qDot2*(dt);
    q2+=qDot3*(dt);
    q3+=qDot4*(dt);

    /* Normalise quaternion */
    recipNorm=1.0/sqrt(q0*q0+q1*q1+q2*q2+q3*q3);
    q0*=recipNorm;
    q1*=recipNorm;
    q2*=recipNorm;
    q3*=recipNorm;
    mfliter->q0=q0;
    mfliter->q1=q1;
    mfliter->q2=q2;
    mfliter->q3=q3;
}

extern void madgwick_filter(madgwick_filter_t *mfliter, const double *accl, const double *gyro, const double *magn, double timestamp)
{
    double gx=gyro[0],gy=gyro[1],gz=gyro[2];
    double ax=accl[0],ay=accl[1],az=accl[2];
    double mx=magn[0],my=magn[1],mz=magn[2];
    double recipNorm;
	double s0,s1,s2,s3;
	double qDot1,qDot2,qDot3,qDot4;
	double hx,hy;
	double _2q0mx,_2q0my,_2q0mz,_2q1mx,_2bx,_2bz,_4bx,_4bz,_2q0,_2q1,_2q2,_2q3,_2q0q2,_2q2q3,q0q0,q0q1,q0q2,q0q3,q1q1,q1q2,q1q3,q2q2,q2q3,q3q3;
    double q0=mfliter->q0;
    double q1=mfliter->q1;
    double q2=mfliter->q2;
    double q3=mfliter->q3;

    if (mfliter->q0==0.0) {
        mfliter->q0=1.0;
        mfliter->q1=0.0;
        mfliter->q2=0.0;
        mfliter->q3=0.0;
    }
    if (mfliter->timestamp==0.0) mfliter->timestamp=timestamp;
    double dt=timestamp-mfliter->timestamp;

	/* Use IMU algorithm if magnetometer measurement invalid (avoids NaN in magnetometer normalisation) */
	if ((mx==0.0)&&(my==0.0)&&(mz==0.0)) {
		madgwick_filter_imu(mfliter,accl,gyro,dt);
	    mfliter->timestamp=timestamp;
		return;
	}
	/* Rate of change of quaternion from gyroscope */
	qDot1=0.5*(-q1*gx-q2*gy-q3*gz);
	qDot2=0.5*( q0*gx+q2*gz-q3*gy);
	qDot3=0.5*( q0*gy-q1*gz+q3*gx);
	qDot4=0.5*( q0*gz+q1*gy-q2*gx);

	/* Compute feedback only if accelerometer measurement valid (avoids NaN in accelerometer normalisation) */
	if(!((ax==0.0)&&(ay==0.0)&&(az==0.0))) {

		/* Normalise accelerometer measurement */
		recipNorm=1.0/sqrt(ax*ax+ay*ay+az*az);
		ax*=recipNorm;
		ay*=recipNorm;
		az*=recipNorm;

		/* Normalise magnetometer measurement */
		recipNorm=1.0/sqrt(mx*mx+my*my+mz*mz);
		mx*=recipNorm;
		my*=recipNorm;
		mz*=recipNorm;

		/* Auxiliary variables to avoid repeated arithmetic */
	    _2q0mx=2.0*q0*mx;
	    _2q0my=2.0*q0*my;
	    _2q0mz=2.0*q0*mz;
	    _2q1mx=2.0*q1*mx;
	    _2q0=2.0*q0;
	    _2q1=2.0*q1;
	    _2q2=2.0*q2;
	    _2q3=2.0*q3;
	    _2q0q2=2.0*q0*q2;
	    _2q2q3=2.0*q2*q3;
	    q0q0=q0*q0;
	    q0q1=q0*q1;
	    q0q2=q0*q2;
	    q0q3=q0*q3;
	    q1q1=q1*q1;
	    q1q2=q1*q2;
	    q1q3=q1*q3;
	    q2q2=q2*q2;
	    q2q3=q2*q3;
	    q3q3=q3*q3;

		/* Reference direction of Earth's magnetic field */
	    hx=mx*q0q0-_2q0my*q3+_2q0mz*q2+mx*q1q1+_2q1*my*q2+_2q1*mz*q3-mx*q2q2-mx*q3q3;
	    hy=_2q0mx*q3+my*q0q0-_2q0mz*q1+_2q1mx*q2-my*q1q1+my*q2q2+_2q2*mz*q3-my*q3q3;
	    _2bx=sqrt(hx*hx+hy*hy);
	    _2bz=-_2q0mx*q2+_2q0my*q1+mz*q0q0+_2q1mx*q3-mz*q1q1+_2q2*my*q3-mz*q2q2+mz*q3q3;
	    _4bx=2.0*_2bx;
	    _4bz=2.0*_2bz;

		/* Gradient decent algorithm corrective step */
	    s0=-_2q2*(2.0*q1q3-_2q0q2-ax)+_2q1*(2.0*q0q1+_2q2q3-ay)-_2bz*q2*(_2bx*(0.5-q2q2-q3q3)+_2bz*(q1q3-q0q2)-mx)+(-_2bx*q3+_2bz*q1)*(_2bx*(q1q2-q0q3)+_2bz*(q0q1+q2q3)-my)+_2bx*q2*(_2bx*(q0q2+q1q3)+_2bz*(0.5-q1q1-q2q2)-mz);
	    s1=_2q3*(2.0*q1q3-_2q0q2-ax)+_2q0*(2.0*q0q1+_2q2q3-ay)-4.0f*q1*(1-2.0*q1q1-2.0*q2q2-az)+_2bz*q3*(_2bx*(0.5-q2q2-q3q3)+_2bz*(q1q3-q0q2)-mx)+(_2bx*q2+_2bz*q0)*(_2bx*(q1q2-q0q3)+_2bz*(q0q1+q2q3)-my)+(_2bx*q3-_4bz*q1)*(_2bx*(q0q2+q1q3)+_2bz*(0.5-q1q1-q2q2)-mz);
	    s2=-_2q0*(2.0*q1q3-_2q0q2-ax)+_2q3*(2.0*q0q1+_2q2q3-ay)-4.0f*q2*(1-2.0*q1q1-2.0*q2q2-az)+(-_4bx*q2-_2bz*q0)*(_2bx*(0.5-q2q2-q3q3)+_2bz*(q1q3-q0q2)-mx)+(_2bx*q1+_2bz*q3)*(_2bx*(q1q2-q0q3)+_2bz*(q0q1+q2q3)-my)+(_2bx*q0-_4bz*q2)*(_2bx*(q0q2+q1q3)+_2bz*(0.5-q1q1-q2q2)-mz);
	    s3=_2q1*(2.0*q1q3-_2q0q2-ax)+_2q2*(2.0*q0q1+_2q2q3-ay)+(-_4bx*q3+_2bz*q1)*(_2bx*(0.5-q2q2-q3q3)+_2bz*(q1q3-q0q2)-mx)+(-_2bx*q0+_2bz*q2)*(_2bx*(q1q2-q0q3)+_2bz*(q0q1+q2q3)-my)+_2bx*q1*(_2bx*(q0q2+q1q3)+_2bz*(0.5-q1q1-q2q2)-mz);
	    recipNorm=1.0/sqrt(s0*s0+s1*s1+s2*s2+s3*s3);
	    s0*=recipNorm;
	    s1*=recipNorm;
	    s2*=recipNorm;
	    s3*=recipNorm;

		/* Apply feedback step */
		qDot1-=mfliter->beta*s0;
		qDot2-=mfliter->beta*s1;
		qDot3-=mfliter->beta*s2;
		qDot4-=mfliter->beta*s3;
	}
	/* Integrate rate of change of quaternion to yield quaternion */
	q0+=qDot1*dt;
	q1+=qDot2*dt;
	q2+=qDot3*dt;
	q3+=qDot4*dt;

	/* Normalise quaternion */
	recipNorm=1.0/sqrt(q0*q0+q1*q1+q2*q2+q3*q3);
	q0*=recipNorm;
	q1*=recipNorm;
	q2*=recipNorm;
	q3*=recipNorm;
    mfliter->q0=q0;
    mfliter->q1=q1;
    mfliter->q2=q2;
    mfliter->q3=q3;
    mfliter->timestamp=timestamp;
}

static void mayhony_filter_imu(mayhony_filter_t *mfliter, const double *accl, const double *gyro, double dt)
{
    double gx=gyro[0],gy=gyro[1],gz=gyro[2];
    double ax=accl[0],ay=accl[1],az=accl[2];
    double recipNorm;
    double halfvx,halfvy,halfvz;
    double halfex,halfey,halfez;
    double qa,qb,qc;
    double q0=mfliter->q0;
    double q1=mfliter->q1;
    double q2=mfliter->q2;
    double q3=mfliter->q3;

    /* Compute feedback only if accelerometer measurement valid (avoids NaN in accelerometer normalisation) */
    if (!((ax==0.0f)&&(ay==0.0f)&&(az==0.0f))){

        /* Normalise accelerometer measurement */
        recipNorm=1.0/sqrt(ax*ax+ay*ay+az*az);
        ax*=recipNorm;
        ay*=recipNorm;
        az*=recipNorm;

        /* Estimated direction of gravity and vector perpendicular to magnetic flux */
        halfvx=q1*q3-q0*q2;
        halfvy=q0*q1+q2*q3;
        halfvz=q0*q0-0.5f+q3*q3;

        /* Error is sum of cross product between estimated and measured direction of gravity */
        halfex=(ay*halfvz-az*halfvy);
        halfey=(az*halfvx-ax*halfvz);
        halfez=(ax*halfvy-ay*halfvx);

        /* Compute and apply integral feedback if enabled */
        if (mfliter->twoKi>0.0){
            mfliter->integralFBx+=mfliter->twoKi*halfex*(dt);
            mfliter->integralFBy+=mfliter->twoKi*halfey*(dt);
            mfliter->integralFBz+=mfliter->twoKi*halfez*(dt);
            gx+=mfliter->integralFBx;
            gy+=mfliter->integralFBy;
            gz+=mfliter->integralFBz;
        }
        else{
            mfliter->integralFBx=0.0f;
            mfliter->integralFBy=0.0f;
            mfliter->integralFBz=0.0f;
        }
        /* Apply proportional feedback */
        gx+=mfliter->twoKp*halfex;
        gy+=mfliter->twoKp*halfey;
        gz+=mfliter->twoKp*halfez;
    }
    /* Integrate rate of change of quaternion */
    gx*=(0.5f*(dt));
    gy*=(0.5f*(dt));
    gz*=(0.5f*(dt));
    qa=q0;
    qb=q1;
    qc=q2;
    q0+=(-qb*gx-qc*gy-q3*gz);
    q1+=(qa*gx+qc*gz-q3*gy);
    q2+=(qa*gy-qb*gz+q3*gx);
    q3+=(qa*gz+qb*gy-qc*gx);

    /* Normalise quaternion */
    recipNorm=1.0/sqrt(q0*q0+q1*q1+q2*q2+q3*q3);
    q0*=recipNorm;
    q1*=recipNorm;
    q2*=recipNorm;
    q3*=recipNorm;
    mfliter->q0=q0;
    mfliter->q1=q1;
    mfliter->q2=q2;
    mfliter->q3=q3;
}
extern void mayhony_filter(mayhony_filter_t *mfliter, const double *accl, const double *gyro, const double *magn, double timestamp)
{
    double gx=gyro[0],gy=gyro[1],gz=gyro[2];
    double ax=accl[0],ay=accl[1],az=accl[2];
    double mx=magn[0],my=magn[1],mz=magn[2];
    double recipNorm;
    double q0q0,q0q1,q0q2,q0q3,q1q1,q1q2,q1q3,q2q2,q2q3,q3q3;
	double hx,hy,bx,bz;
	double halfvx,halfvy,halfvz,halfwx,halfwy,halfwz;
	double halfex,halfey,halfez;
	double qa,qb,qc;
    double q0=mfliter->q0;
    double q1=mfliter->q1;
    double q2=mfliter->q2;
    double q3=mfliter->q3;

    if (mfliter->q0==0.0) {
        mfliter->q0=1.0;
        mfliter->q1=0.0;
        mfliter->q2=0.0;
        mfliter->q3=0.0;
    }
    if (mfliter->timestamp==0.0) mfliter->timestamp=timestamp;
    double dt=timestamp-mfliter->timestamp;

    /* Use IMU algorithm if magnetometer measurement invalid (avoids NaN in magnetometer normalisation) */
	if ((mx==0.0f)&&(my==0.0f)&&(mz==0.0f)){
		mayhony_filter_imu(mfliter,accl,gyro,dt);
	    mfliter->timestamp=timestamp;
		return;
	}
    /* Compute feedback only if accelerometer measurement valid (avoids NaN in accelerometer normalisation) */
	if (!((ax==0.0f)&&(ay==0.0f)&&(az==0.0f))) {

		/* Normalise accelerometer measurement */
		recipNorm=1.0/sqrt(ax*ax+ay*ay+az*az);
		ax*=recipNorm;
		ay*=recipNorm;
		az*=recipNorm;

		/* Normalise magnetometer measurement */
		recipNorm=1.0/sqrt(mx*mx+my*my+mz*mz);
		mx*=recipNorm;
		my*=recipNorm;
		mz*=recipNorm;

	    /* Auxiliary variables to avoid repeated arithmetic */
        q0q0=q0*q0;
        q0q1=q0*q1;
        q0q2=q0*q2;
        q0q3=q0*q3;
        q1q1=q1*q1;
        q1q2=q1*q2;
        q1q3=q1*q3;
        q2q2=q2*q2;
        q2q3=q2*q3;
        q3q3=q3*q3;

	    /* Reference direction of Earth's magnetic field */
        hx=2.0f*(mx*(0.5f-q2q2-q3q3)+my*(q1q2-q0q3)+mz*(q1q3+q0q2));
        hy=2.0f*(mx*(q1q2+q0q3)+my*(0.5f-q1q1-q3q3)+mz*(q2q3-q0q1));
        bx=sqrt(hx*hx+hy*hy);
        bz=2.0f*(mx*(q1q3-q0q2)+my*(q2q3+q0q1)+mz*(0.5f-q1q1-q2q2));

		/* Estimated direction of gravity and magnetic-field */
		halfvx=q1q3-q0q2;
		halfvy=q0q1+q2q3;
		halfvz=q0q0-0.5f+q3q3;
        halfwx=bx*(0.5f-q2q2-q3q3)+bz*(q1q3-q0q2);
        halfwy=bx*(q1q2-q0q3)+bz*(q0q1+q2q3);
        halfwz=bx*(q0q2+q1q3)+bz*(0.5f-q1q1-q2q2);

		/* Error is sum of cross product between estimated direction and measured direction of field vectors */
		halfex=(ay*halfvz-az*halfvy)+(my*halfwz-mz*halfwy);
		halfey=(az*halfvx-ax*halfvz)+(mz*halfwx-mx*halfwz);
		halfez=(ax*halfvy-ay*halfvx)+(mx*halfwy-my*halfwx);

		/* Compute and apply integral feedback if enabled */
		if (mfliter->twoKi>0.0f){
			mfliter->integralFBx+=mfliter->twoKi*halfex*(dt);
			mfliter->integralFBy+=mfliter->twoKi*halfey*(dt);
			mfliter->integralFBz+=mfliter->twoKi*halfez*(dt);
			gx+=mfliter->integralFBx;
			gy+=mfliter->integralFBy;
			gz+=mfliter->integralFBz;
		}
		else{
			mfliter->integralFBx=0.0f;
			mfliter->integralFBy=0.0f;
			mfliter->integralFBz=0.0f;
		}
		/* Apply proportional feedback */
		gx+=mfliter->twoKp*halfex;
		gy+=mfliter->twoKp*halfey;
		gz+=mfliter->twoKp*halfez;
	}
	/* Integrate rate of change of quaternion */
	gx*=(0.5*dt);
	gy*=(0.5*dt);
	gz*=(0.5*dt);
	qa=q0;
	qb=q1;
	qc=q2;
	q0+=(-qb*gx-qc*gy-q3*gz);
	q1+=(qa*gx+qc*gz-q3*gy);
	q2+=(qa*gy-qb*gz+q3*gx);
	q3+=(qa*gz+qb*gy-qc*gx);

	/* Normalise quaternion */
	recipNorm=1.0/sqrt(q0*q0+q1*q1+q2*q2+q3*q3);
	q0*=recipNorm;
	q1*=recipNorm;
	q2*=recipNorm;
	q3*=recipNorm;
    mfliter->q0=q0;
    mfliter->q1=q1;
    mfliter->q2=q2;
    mfliter->q3=q3;
}

extern int iopos_imumeass_in_2odos(const imu_odo_measurement_t *meas, imu_measurement_in_2odos_t *imu2odos)
{
    if (imu2odos->nimu<1024) imu2odos->imus[imu2odos->nimu++]=*meas;
    else {
        mempcpy(imu2odos->imus,imu2odos->imus+1,sizeof(*meas)*(imu2odos->nimu-1));
        imu2odos->imus[imu2odos->nimu-1]=*meas;
    }
    int flag=0;

    imu2odos->uflg=0;

    if (meas->type==IOPOS_MEASUREMENT_ODO) {
        imu2odos->odoe=&imu2odos->imus[imu2odos->nimu-1];

        imu2odos->odos=imu2odos->odoe;
        while (imu2odos->odos!=imu2odos->imus) {
            imu2odos->odos--;
            if (imu2odos->odos->type==IOPOS_MEASUREMENT_ODO) {
                if (imu2odos->odoe->timestamp-imu2odos->odos->timestamp>imu2odos->dt) {flag=1; break;}
            }
        }
    }
    imu2odos->uflg=flag;
    return flag;
}

extern void iopos_imudv_in_2odos(const imu_measurement_in_2odos_t *imu2odos, double *dv)
{
    double dt,dv_odo[2];
    imu_odo_measurement_t *pimu;

    dv[0]=dv[1]=0.0;

    if (!imu2odos->odos) return;
    if (!imu2odos->odoe) return;
    pimu=imu2odos->odos+1;

    dv_odo[0]=imu2odos->odoe->val[2]-imu2odos->odos->val[2];
    dv_odo[1]=imu2odos->odoe->val[3]-imu2odos->odos->val[3];

    for (;;) {
        if (pimu->type==IOPOS_MEASUREMENT_IMU) {
            dt=pimu->timestamp-(pimu-1)->timestamp;
            dv[0]+=pimu->val[3]*dt;
            dv[1]+=pimu->val[4]*dt;
        }
        if (++pimu>imu2odos->odoe) break;
    }
}

extern void iopos2dinit(imu_odo_pos2d_t *iopos2d)
{
    iopos2d->madgwickflt.beta=0.1;
    iopos2d->mayhonyflt.twoKi=(2.0*0.0);
    iopos2d->mayhonyflt.twoKp=(2.0*0.5);
}

static void iopos2d_ebias(imu_odo_pos2d_t *iopos2d, const imu_odo_measurement_t *meas)
{
    imu_measurement_in_2odos_t *imu2odos=&iopos2d->imu2odos;
    if (imu2odos->uflg) {
        if (norm(imu2odos->odos->val,4)<1E-6&&
            norm(imu2odos->odoe->val,4)<1E-6&&
            imu2odos->odoe->timestamp-imu2odos->odos->timestamp<1.0) {

            imu_odo_measurement_t *pimu=imu2odos->odos+1;
            const double bias_alpha=0.2;
            const double g=9.81;

            for (;;) {
                if (pimu->type==IOPOS_MEASUREMENT_IMU) {
                    iopos2d->bg[0]+=bias_alpha*(pimu->val[0]-iopos2d->bg[0]);
                    iopos2d->bg[1]+=bias_alpha*(pimu->val[1]-iopos2d->bg[1]);
                    iopos2d->bg[2]+=bias_alpha*(pimu->val[2]-iopos2d->bg[2]);
                    iopos2d->ba[0]+=bias_alpha*(pimu->val[3]-iopos2d->ba[0]);
                    iopos2d->ba[1]+=bias_alpha*(pimu->val[4]-iopos2d->ba[1]);
                    iopos2d->ba[2]+=bias_alpha*(pimu->val[5]-g-iopos2d->ba[2]);
                }
                if (++pimu>iopos2d->imu2odos.odoe) break;
            }
        }
    }
}

static void iopos2d_inite(imu_odo_pos2d_t *iopos2d, const imu_odo_measurement_t *meas)
{
    if (meas->type!=IOPOS_MEASUREMENT_ODO) return;

    if (iopos2d->timestamp==0.0) {
        iopos2d->timestamp=meas->timestamp;
        iopos2d->timestamp_init=meas->timestamp;
        iopos2d->yaw=0.0;
        iopos2d->yawu=0.0;
        iopos2d->yaw_var=SQR(1.0/180.0*PI);
        iopos2d->p[0]=0.0;
        iopos2d->p[1]=0.0;
        iopos2d->v[0]=(meas->val[2]+meas->val[3])/2.0;
        iopos2d->v[1]=0.0;
        iopos2d->direction=0;
        iopos2d->state=IOPOS_INIT;

        iopos2d->mayhonyflt.timestamp=meas->timestamp;
        iopos2d->madgwickflt.timestamp=meas->timestamp;
        iopos2d->mayhonyflt.q0=1.0;
        iopos2d->mayhonyflt.q1=0.0;
        iopos2d->mayhonyflt.q2=0.0;
        iopos2d->mayhonyflt.q3=0.0;
        iopos2d->madgwickflt.q0=1.0;
        iopos2d->madgwickflt.q1=0.0;
        iopos2d->madgwickflt.q2=0.0;
        iopos2d->madgwickflt.q3=0.0;
    }
}

static double normang(double ang)
{
    ang=fmod(ang,2.0*PI);
    if (ang<0.0) ang+=2.0*PI;
    return ang;
}

static void iopos2d_udodo(imu_odo_pos2d_t *iopos2d, const imu_odo_measurement_t *meas)
{
    if (meas->type!=IOPOS_MEASUREMENT_ODO) return;
    if (iopos2d->state==IOPOS_NONE) return;

    const double thresdt=3.0;
    double dt=meas->timestamp-iopos2d->timestamp;
    double vel=(meas->val[2]+meas->val[3])/2.0;
    double gyro[3]={0};
    int i;

    for (i=iopos2d->imu2odos.nimu-1;i>=0;i--) {
        if (iopos2d->imu2odos.imus[i].type==IOPOS_MEASUREMENT_IMU) {
            if (fabs(meas->timestamp-iopos2d->imu2odos.imus[i].timestamp)<thresdt) {
                gyro[0]=iopos2d->imu2odos.imus[i].val[0]-iopos2d->bg[0];
                gyro[1]=iopos2d->imu2odos.imus[i].val[1]-iopos2d->bg[1];
                gyro[2]=iopos2d->imu2odos.imus[i].val[2]-iopos2d->bg[2];
                break;
            }
        }
    }
    iopos2d->v[0]=vel*cos(iopos2d->yaw);
    iopos2d->v[1]=vel*sin(iopos2d->yaw);
    iopos2d->p[0]+=vel*cos(iopos2d->yaw+0.5*gyro[2]*dt)*dt;
    iopos2d->p[1]+=vel*sin(iopos2d->yaw+0.5*gyro[2]*dt)*dt;
    iopos2d->yaw+=gyro[2]*dt;
    iopos2d->yawu+=gyro[2]*dt;
    iopos2d->state=IOPOS_UDODO_OK;
    iopos2d->timestamp=meas->timestamp;
    iopos2d->yaw=normang(iopos2d->yaw);
}

static void iopos2d_uahrs(imu_odo_pos2d_t *iopos2d, const imu_odo_measurement_t *meas)
{
    if (meas->type!=IOPOS_MEASUREMENT_ODO) return;
    if (iopos2d->state==IOPOS_NONE) return;
    if (iopos2d->udahrs==0) return;

    imu_measurement_in_2odos_t *imu2odos=&iopos2d->imu2odos;
    if (imu2odos->uflg) {
        imu_odo_measurement_t *pimu=imu2odos->odos+1;
        double magn[3]={0};
        double accl[3]={0};
        double gyro[3]={0};
        int i;

        for (;;) {
            if (pimu->type==IOPOS_MEASUREMENT_IMU) {
                for (i=0;i<3;i++) {
                    accl[i]=pimu->val[3+i]-iopos2d->ba[i];
                    gyro[i]=pimu->val[0+i]-iopos2d->bg[i];
                }
                madgwick_filter(&iopos2d->madgwickflt,accl,gyro,magn,pimu->timestamp);
                mayhony_filter (&iopos2d->mayhonyflt ,accl,gyro,magn,pimu->timestamp);
            }
            if (++pimu>iopos2d->imu2odos.odoe) break;
        }
    }
}

static void iopos2d_udyaw(imu_odo_pos2d_t *iopos2d, const imu_odo_measurement_t *meas)
{
    if (iopos2d->udahrs==0) return;
    if (meas->type!=IOPOS_MEASUREMENT_ODO) return;
    if (iopos2d->state==IOPOS_NONE) return;

    quat_t q1;
    q1.q0=iopos2d->madgwickflt.q0;
    q1.q1=iopos2d->madgwickflt.q1;
    q1.q2=iopos2d->madgwickflt.q2;
    q1.q3=iopos2d->madgwickflt.q3;
    quat_t q2;
    q2.q0=iopos2d->mayhonyflt.q0;
    q2.q1=iopos2d->mayhonyflt.q1;
    q2.q2=iopos2d->mayhonyflt.q2;
    q2.q3=iopos2d->mayhonyflt.q3;

    quat_t qo;
    quat_slerp(&qo,&q1,&q2,0.5);

    euler_t euler={0};
    euler.yaw=iopos2d->yaw;
    quat_t qu;
    euler_to_quat(&euler,&qu);

    double yaw_var=SQR(3.0/180.0*PI*exp((meas->timestamp-iopos2d->timestamp_init)/3600.0));
    double yaw_noise=SQR(0.5/180.0*PI);

    iopos2d->yaw_var+=yaw_noise;
    double yaw_var_s=1.0/(1.0/iopos2d->yaw_var+1.0/yaw_var);

    quat_t qq;
    quat_slerp(&qq,&qo,&qu,yaw_var_s/iopos2d->yaw_var);

    euler.yaw=iopos2d->yawu;
    euler_to_quat(&euler,&qu);
    quat_slerp(&qo,&qq,&qu,0.1);

    double x=qo.x,y=qo.y,z=qo.z,w=qo.w;
    double ww=w*w,xx=x*x,yy=y*y,zz=z*z;
    iopos2d->yaw=normang(atan2(2.0*(x*y+z*w),xx-yy-zz+ww));
    iopos2d->yaw_var=yaw_var_s;
}

static void iopos2d_ugear(imu_odo_pos2d_t *iopos2d, const imu_odo_measurement_t *meas)
{
    if (meas->type==IOPOS_MEASUREMENT_GEAR) {
        if (iopos2d->direction==0) {
            iopos2d->direction=meas->val[0];
        }
        else if (iopos2d->direction*meas->val[0]<0) {
            iopos2d->direction*=-1;
            iopos2d->yaw*=1.0;
            iopos2d->yawu*=1.0;
        }
    }
}

extern void iopos2d(imu_odo_pos2d_t *iopos2d, const imu_odo_measurement_t *meas)
{
    iopos_imumeass_in_2odos(meas,&iopos2d->imu2odos);
    iopos2d_ugear(iopos2d,meas);
    iopos2d_ebias(iopos2d,meas);
    iopos2d_inite(iopos2d,meas);
    iopos2d_uahrs(iopos2d,meas);
    iopos2d_udodo(iopos2d,meas);
    iopos2d_udyaw(iopos2d,meas);
}






