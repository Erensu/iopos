
#include "iopos.h"

#include "opencv2/imgcodecs.hpp"
#include "opencv2/highgui.hpp"
#include "opencv2/imgproc.hpp"

using namespace cv;
using namespace std;

static int cmpmeas(const void *p1, const void *p2)
{
    imu_odo_measurement_t *q1=(imu_odo_measurement_t *)p1,*q2=(imu_odo_measurement_t *)p2;
    double tt=q1->timestamp-q2->timestamp;
    if (fabs(tt)>1E-6) return tt<0?-1:1;
    return 0;
}

static int sortmeas(imu_odo_measurement_t *meas, int nmeas)
{
    int i,j;

    if (nmeas<=0) return 0;

    qsort(meas,nmeas,sizeof(imu_odo_measurement_t),cmpmeas);

    for (i=j=0;i<nmeas;i++) {
        if (meas[i].timestamp-meas[j].timestamp!=0.0) {
            meas[++j]=meas[i];
        }
    }
    nmeas=j+1;
    return nmeas;
}

static int readimuodomeas(const char *fileimu, const char *fileodo, imu_odo_measurement_t* meas)
{
    FILE *fpimu=fopen(fileimu,"r");
    FILE *fpodo=fopen(fileodo,"r");
    int n=0;

    if (fpimu) {
        char buff[1024];
        while (fgets(buff,1024,fpimu)) {
            imu_odo_measurement_t t={0};

            char *p,*q,*val[64];
            int m=0;
            for (p=buff;*p&&m<64;p=q+1) {
                if ((q=strchr(p,','))||(q=strchr(p,'\n'))) {
                    val[m++]=p; *q='\0';
                }
                else break;
            }
            if (m<=0) continue;
            t.timestamp=atof(val[0])/1E9;
            t.val[0]=atof(val[1]);
            t.val[1]=atof(val[2]);
            t.val[2]=atof(val[3]);
            t.val[3]=atof(val[4]);
            t.val[4]=atof(val[5]);
            t.val[5]=atof(val[6]);
            t.type=IOPOS_MEASUREMENT_IMU;
            meas[n++]=t;
        }
    }
    if (fpodo) {
        char buff[1024];
        while (fgets(buff,1024,fpodo)) {
            imu_odo_measurement_t t={0};
            char *p,*q,*val[64];
            int m=0;
            for (p=buff;*p&&m<64;p=q+1) {
                if ((q=strchr(p,','))||(q=strchr(p,'\n'))) {
                    val[m++]=p; *q='\0';
                }
                else break;
            }
            if (m<=0) continue;
            t.timestamp=atof(val[0])/1E9;
            t.val[0]=atof(val[1]);
            t.val[1]=atof(val[2]);
            t.val[2]=atof(val[3]);
            t.val[3]=atof(val[4]);
            t.val[4]=atof(val[5]);
            t.type=IOPOS_MEASUREMENT_ODO;
            meas[n++]=t;
        }
    }
    if (fpimu) fclose(fpimu);
    if (fpodo) fclose(fpodo);
    return sortmeas(meas,n);
}

static void drawtrajectory(const double *pts, int npts)
{
    if (npts<=0) return;

    double width=1000;
    double height=1000;
    cv::Mat traj=cv::Mat::zeros(static_cast<int>(width),static_cast<int>(height),CV_8UC3);

    double maxx,maxy;
    double minx,miny;

    int i;
    for (i=0;i<npts;i++) {
        if (i==0||maxx<pts[2*i+0]) maxx=pts[2*i+0];
        if (i==0||maxy<pts[2*i+1]) maxy=pts[2*i+1];
        if (i==0||minx>pts[2*i+0]) minx=pts[2*i+0];
        if (i==0||miny>pts[2*i+1]) miny=pts[2*i+1];
    }
    double sx=fabs((width)/(maxx-minx));
    double sy=fabs((height)/(maxy-miny));

    double s=MIN(MIN(sx,sy),10);

    for (i=0;i<npts;i++) {
        double x=(pts[2*i+0]-minx)*s;
        double y=(pts[2*i+1]-miny)*s;

        cv::circle(traj,cv::Point2d(x,y),i==0?6:2,i==0?Scalar(0,255,0):Scalar(255,0,0),-1);
        cv::imshow("trajectory",traj);
    }
    cv::imshow("trajectory",traj);
    cv::waitKey();
}

int main( int nargv, char** argv )
{
    const char *imufile=R"(..\data\1\IMU\Data\data.csv)";
    const char *odofile=R"(..\data\1\Wheel\Data\data.csv)";
    imu_odo_measurement_t *meas=(imu_odo_measurement_t*)malloc(1024*64*sizeof(imu_odo_measurement_t));
    double *pts_odo=(double*)malloc(2*1024*64*sizeof(double));
    int nmeas=0;
    FILE *fpmeas=fopen("meas.txt","w");
    FILE *fpsols_odo=fopen("sols_odo.txt","w");

    nmeas=readimuodomeas(imufile,odofile,meas);

    for (int i=0;i<nmeas-1;i++) {
        fprintf(fpmeas,"%10.5lf %d %6.3lf  %6.3lf %6.3lf\n",meas[i].timestamp,meas[i].type,meas[i+1].timestamp-meas[i].timestamp,meas[i].type==0?meas[i].val[3]:0.0,meas[i].type==0?meas[i].val[4]:0.0);
        fflush(fpmeas);
    }
    int npts_odo=0;

    imu_measurement_in_2odos_t imu2odos={0};
    imu2odos.dt=2.0;

    imu_odo_pos2d_t iopos2dg={0};

    iopos2dinit(&iopos2dg);
    iopos2dg.udahrs=1;
    iopos2dg.bg[0]=0.00024886440677966102;
    iopos2dg.bg[1]=0.001842457627118643;
    iopos2dg.bg[2]=0.002685135593220338;

    iopos2dg.ba[0]=2.3483144121275168e-05;
    iopos2dg.ba[1]=1.26891666150028e-05;
    iopos2dg.ba[2]=-0.0029452899366244623;

    for (int i=0;i<nmeas;i++) {

        iopos2d(&iopos2dg,&meas[i]);

        if (iopos2dg.state==IOPOS_UDODO_OK) {
            fprintf(fpsols_odo,"%12.8lf %12.4lf %12.4lf %12.4lf\n",iopos2dg.timestamp,iopos2dg.p[0],iopos2dg.p[1],iopos2dg.yaw/PI*180);
            fflush(fpsols_odo);

            pts_odo[2*npts_odo+0]=iopos2dg.p[0];
            pts_odo[2*npts_odo+1]=iopos2dg.p[1];
            npts_odo++;;
        }
    }
    drawtrajectory(pts_odo,npts_odo);
    free(meas);
    fclose(fpmeas);
    fclose(fpsols_odo);
    return 0;

}