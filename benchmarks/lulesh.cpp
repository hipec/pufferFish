#include <cmath>
#include <iostream>
#include <chrono>
#include <bits/stdc++.h>
#include <string>
#include "common.h"

#define CONST_FACTOR 32
using namespace std;

#ifdef USE_REGULAR
#define irregular_recursion1D recursion1D
#endif

class Domain{
    private:
        /* Node-centered */

        /* coordinates */
        double *m_x;
        double *m_y;
        double *m_z;

        /* velocities */
        double *m_xd;
        double *m_yd;
        double *m_zd;

        /* accelerations */
        double *m_xdd;
        double *m_ydd;
        double *m_zdd;

        /* forces */
        double *m_fx;
        double *m_fy;
        double *m_fz;

        /* mass */
        double *m_nodalMass;

        /* symmetry plane nodesets */
        int *m_symmX;
        int *m_symmY;
        int *m_symmZ;

        int *m_nodeElemCount;
        int *m_nodeElemStart;
        int *m_nodeElemCornerList;

        /* Element-centered */

        /* material indexset */
        int *m_matElemlist;
        /* elemToNode connectivity */
        int *m_nodelist; //2D originally

        /* element connectivity across each face */
        int *m_lxim;
        int *m_lxip;
        int *m_letam;
        int *m_letap;
        int *m_lzetam;
        int *m_lzetap;

        /* symmetry/free-surface flags for each elem face */
        int *m_elemBC;

        /* principal strains -- temporary */
        double *m_dxx;
        double *m_dyy;
        double *m_dzz;

        /* velocity gradient -- temporary */
        double *m_delv_xi;
        double *m_delv_eta;
        double *m_delv_zeta;

        /* coordinate gradient -- temporary */
        double *m_delx_xi;
        double *m_delx_eta;
        double *m_delx_zeta;

        double *m_e;   /* energy */

        double *m_p;   /* pressure */
        double *m_q;   /* q */
        double *m_ql;  /* linear term for q */
        double *m_qq;  /* quadratic term for q */

        double *m_v;     /* relative volume */
        double *m_volo;  /* reference volume */
        double *m_vnew;  /* new relative volume -- temporary */
        double *m_delv;  /* m_vnew - m_v */
        double *m_vdov;  /* volume derivative over volume */

        double *m_arealg;  /* characteristic length of an element */

        double *m_ss;      /* "sound speed" */

        double *m_elemMass;  /* mass */

    public:
       ~Domain() {
                hclib::numa_free(m_matElemlist);
                hclib::numa_free(m_nodelist);
                hclib::numa_free(m_lxim);
                hclib::numa_free(m_lxip);
                hclib::numa_free(m_letam);
                hclib::numa_free(m_letap);
                hclib::numa_free(m_lzetam);
                hclib::numa_free(m_lzetap);
                hclib::numa_free(m_elemBC);
                hclib::numa_free(m_e);
                hclib::numa_free(m_p);
                hclib::numa_free(m_q);
                hclib::numa_free(m_ql);
                hclib::numa_free(m_qq);
                hclib::numa_free(m_v);
                hclib::numa_free(m_volo);
                hclib::numa_free(m_delv);
                hclib::numa_free(m_vdov);
                hclib::numa_free(m_arealg);
                hclib::numa_free(m_ss);
                hclib::numa_free(m_elemMass);
                hclib::numa_free(m_dxx);
                hclib::numa_free(m_dyy);
                hclib::numa_free(m_dzz);
                hclib::numa_free(m_delv_xi);
                hclib::numa_free(m_delv_eta);
                hclib::numa_free(m_delv_zeta);
                hclib::numa_free(m_delx_xi);
                hclib::numa_free(m_delx_eta);
                hclib::numa_free(m_delx_zeta);
                hclib::numa_free(m_vnew);
                hclib::numa_free(m_x);
                hclib::numa_free(m_y);
                hclib::numa_free(m_z);
                hclib::numa_free(m_xd);
                hclib::numa_free(m_yd);
                hclib::numa_free(m_zd);
                hclib::numa_free(m_xdd);
                hclib::numa_free(m_ydd);
                hclib::numa_free(m_zdd);
                hclib::numa_free(m_fx);
                hclib::numa_free(m_fy);
                hclib::numa_free(m_fz);
                hclib::numa_free(m_nodalMass);
                hclib::numa_free(m_nodeElemCount);
                hclib::numa_free(m_nodeElemStart);
                hclib::numa_free(m_symmX);
                hclib::numa_free(m_symmY);
                hclib::numa_free(m_symmZ);
        } 

    // ------------------------------ FIELDS ------------------------------

        /* Parameters */
        double m_dtfixed;           /* fixed time increment */
        double m_time;              /* current time */
        double m_deltatime;         /* variable time increment */
        double m_deltatimemultlb;
        double m_deltatimemultub;
        double m_stoptime;          /* end time for simulation */

        double m_u_cut;             /* velocity tolerance */
        double m_hgcoef;            /* hourglass control */
        double m_qstop;             /* excessive q indicator */
        double m_monoq_max_slope;
        double m_monoq_limiter_mult;
        double m_e_cut;             /* energy tolerance */
        double m_p_cut;             /* pressure tolerance */
        double m_ss4o3;
        double m_q_cut;             /* q tolerance */
        double m_v_cut;             /* relative volume tolerance */
        double m_qlc_monoq;         /* linear term coef for q */
        double m_qqc_monoq;         /* quadratic term coef for q */
        double m_qqc;
        double m_eosvmax;
        double m_eosvmin;
        double m_pmin;              /* pressure floor */
        double m_emin;              /* energy floor */
        double m_dvovmax;           /* maximum allowable volume change */
        double m_refdens;           /* reference density */

        double m_dtcourant;         /* courant constraint */
        double m_dthydro;           /* volume change constraint */
        double m_dtmax;             /* maximum allowable time increment */

        int m_cycle;             /* iteration count for simulation */

        int m_sizeX;           /* X,Y,Z extent of this block */
        int m_sizeY;
        int m_sizeZ;

        int m_numElem;         /* Elements/Nodes in this domain */
        int m_numNode;

        void _fill(double *arr, int size, double val);
        void _fill(int *arr, int size, int val);
        void AllocateElemPersistent(int size);
        void AllocateElemTemporary(int size);
        void AllocateNodalPersistent(int size);
        void AllocateNodeElemIndexes();
        int nodeElemCornerList(int idx);
        void nodeElemCornerList(int idx, int val);
        int nodeElemCount(int idx);
        void nodeElemCount(int idx, int val);
        int nodeElemStart(int idx);
        void nodeElemStart(int idx, int val);
        int* nodelist();
        void AllocateNodesets(int size);
        double arealg(int idx);
        void arealg(int idx, double val);
        int cycle();
        void cycle(int val);
        double deltatime();
        void deltatime(double val);
        double deltatimemultlb();
        void deltatimemultlb(double val);
        double deltatimemultub();
        void deltatimemultub(double val);
        double delv(int idx);
        void delv(int idx, double val);
        double delv_eta(int idx);
        void delv_eta(int idx, double val);
        double delv_xi(int idx);
        void delv_xi(int idx, double val);
        double delv_zeta(int idx);
        void delv_zeta(int idx, double val);
        double delx_eta(int idx);
        void delx_eta(int idx, double val);
        double delx_xi(int idx);
        void delx_xi(int idx, double val);
        double delx_zeta(int idx);
        void delx_zeta(int idx, double val);
        double dtcourant();
        void dtcourant(double dtcourant);
        double dtfixed();
        void dtfixed(double val);
        double dthydro();
        void dthydro(double dthydro);
        double dtmax();
        void dtmax(double val);
        double dvovmax();
        void dvovmax(double val);
        double dxx(int idx);
        void dxx(int idx, double val);
        double dyy(int idx);
        void dyy(int idx, double val);
        double dzz(int idx) ;
        void dzz(int idx, double val) ;
        double e(int idx);
        void e(int idx, double val);
        double e_cut();
        void e_cut(double val) ;
        int elemBC(int idx);
        void elemBC(int idx, int val);
        double elemMass(int idx);
        void elemMass(int idx, double val);
        double emin();
        void emin(double val);
        double eosvmax();
        void eosvmax(double val);
        double eosvmin();
        void eosvmin(double val);
        double fx(int idx);
        void fx(int idx, double val);
        void fxFill(double val, int size);
        double fy(int idx);
        void fy(int idx, double val);
        void fyFill(double val, int size);
        double fz(int idx);
        void fz(int idx, double val);
        void fzFill(double val, int size);
        double hgcoef();
        void hgcoef(double val);
        int letam(int idx) ;
        void letam(int idx, int val);
        int letap(int idx);
        void letap(int idx, int val);
        int lxim(int idx);
        void lxim(int idx, int val);
        int lxip(int idx);
        void lxip(int idx, int val);
        int lzetam(int idx);
        void lzetam(int idx, int val);
        int lzetap(int idx);
        void lzetap(int idx, int val);
        int matElemlist(int idx);
        void matElemlist(int idx, int val);
        double monoq_limiter_mult();
        void monoq_limiter_mult(double val);
        double monoq_max_slope();
        void monoq_max_slope(double val);
        double nodalMass(int idx);
        void nodalMass(int idx, double val);
        int numElem();
        void numElem(int val);
        void numNode(int val);
        int numNode();
        double p(int idx);
        void p(int idx, double val);
        double p_cut();
        void p_cut(double val);
        double pmin();
        void pmin(double val);
        double q(int idx);
        void q(int idx, double val);
        double q_cut();
        void q_cut(double val);
        double ql(int idx);
        void ql(int idx, double val);
        double qlc_monoq();
        void qlc_monoq(double val);
        double qq(int idx);
        void qq(int idx, double val);
        double qqc();
        void qqc(double val);
        double qqc_monoq();
        void qqc_monoq(double val);
        double qstop();
        void qstop(double val);
        double refdens();
        void refdens(double val);
        int sizeX();
        void sizeX(int val);
        int sizeY();
        void sizeY(int val);
        int sizeZ();
        void sizeZ(int val);
        double ss(int idx);
        void ss(int idx, double val);
        double ss4o3();
        void ss4o3(double val);
        double stoptime();
        void stoptime(double val);
        int symmX(int idx);
        void symmX(int idx, int val);
        int symmY(int idx);
        void symmY(int idx, int val);
        int symmZ(int idx);
        void symmZ(int idx, int val);
        double time();
        void time(double val);
        double u_cut();
        void u_cut(double val);
        double v(int idx);
        void v(int idx, double val);
        double v_cut();
        void v_cut(double val);
        double vdov(int idx);
        void vdov(int idx, double val);
        double vnew(int idx);
        void vnew(int idx, double val);
        double volo(int idx);
        void volo(int idx, double val);
        double x(int idx);
        void x(int idx, double val);
        double xd(int idx);
        void xd(int idx, double val);
        double xdd(int idx);
        void xdd(int idx, double val);
        double y(int idx);
        void y(int idx, double val);
        double yd(int idx);
        void yd(int idx, double val);
        double ydd(int idx);
        void ydd(int idx, double val);
        double z(int idx) ;
        void z(int idx, double val);
        double zd(int idx);
        void zd(int idx, double val);
        double zdd(int idx);
        void zdd(int idx, double val);
	void* get_m_delv_zeta() {return (void*) m_delv_zeta; }
	void* get_m_delx_eta() {return (void*) m_delx_eta; }
	void* get_m_dxx() {return (void*) m_dxx; }
	void* get_m_dyy() { return (void*) m_dyy; }
	void* get_m_e() {return (void*) m_e; }
	void* get_m_fx() { return (void*) m_fx; }
	void* get_m_p() { return (void*) m_p; }
	void* get_m_v() {return (void*) m_v; }
	void* get_m_xdd() {return (void*) m_xdd; }
	void* get_m_xd() { return (void*)m_xd; }
};


class LuleshConfig{

    protected:
        static void printArgs();

    public:
	    static const int VOLUME_ERROR;
	    static const int Q_STOP_ERROR;

	    /* Stuff needed for boundary conditions, 2 BCs on each of 6 hexahedral faces (12 bits) */
	    static const int XI_M;
	    static const int XI_M_SYMM;
	    static const int XI_M_FREE;

	    static const int XI_P;
	    static const int XI_P_SYMM;
	    static const int XI_P_FREE;

	    static const int ETA_M;
	    static const int ETA_M_SYMM;
	    static const int ETA_M_FREE;

	    static const int ETA_P;
	    static const int ETA_P_SYMM;
	    static const int ETA_P_FREE;

	    static const int ZETA_M;
	    static const int ZETA_M_SYMM;
	    static const int ZETA_M_FREE;

	    static const int ZETA_P;
	    static const int ZETA_P_SYMM;
	    static const int ZETA_P_FREE;

	    static int N;
	    static bool debug;

        static void parseArgs(int, char **);

};


class LuleshParBenchmark{
    public:
        void initialize(int, string[]);
        void printArgInfo();
        void runIteration();
        void cleanupIteration(bool lastIteration, double execTimeMillis);
        static void runLulesh(int edgeElems);
    
    private:
        static double* pHalfStep;
        static double* vnewc;
    	static double *e_old;
    	static double *delvc;
    	static double *p_old;
   	static double *q_old;
    	static double *compression;
    	static double *compHalfStep;
    	static double *qq;
    	static double *ql;
    	static double *work;
    	static double *p_new;
    	static double *e_new;
    	static double *q_new;
    	static double *bvc;
    	static double *pbvc;

	static void initData_evalEOSForElems(int length);
	static void freeData();
        static double calcElemVolume(double x[], double y[], double z[]);
        static void timeIncrement(Domain& domain);
        static void lagrangeLeapFrog(Domain& domain, double sigxx[],
                        double sigyy[], double sigzz[],double determinant[],
                        double dvdx[], double dvdy[], double dvdz[],
                        double x8n[], double y8n[], double z8n[],
                        double fxElem[], double fyElem[], double fzElem[]);
        static void lagrangeNodal(Domain& domain, double sigxx[], double sigyy[],
                        double sigzz[],double determinant[],
                        double dvdx[], double dvdy[], double dvdz[],
                        double x8n[], double y8n[], double z8n[],
                        double fxElem[], double fyElem[], double fzElem[]);
        static void calcForceForNodes(Domain& domain, double sigxx[], double sigyy[],
                        double sigzz[], double determinant[],
                        double dvdx[], double dvdy[], double dvdz[],
                        double x8n[], double y8n[], double z8n[],
                        double fxElem[], double fyElem[], double fzElem[]);
        static void calcVolumeForceForElems(Domain& domain, double sigxx[], double sigyy[], double sigzz[], double determinant[],
                    double dvdx[], double dvdy[], double dvdz[],
                    double x8n[], double y8n[], double z8n[],
                    double fxElem[], double fyElem[], double fzElem[]);
        static void initStressTermsForElems(Domain& domain, int numElem, double sigxx[], double sigyy[], double sigzz[]);
        static void integrateStressForElems(
                Domain& domain, int numElem,
                double sigxx[], double sigyy[], double sigzz[],
                double determinant[],
                double fxElem[], double fyElem[], double fzElem[]);
        static void calcElemNodeNormals(
                double pfx[], double pfy[], double pfz[],
                double x[], double y[], double z[]);
        static void sumElemFaceNormal(
                double nX0[], int nX0index, double nY0[], int nY0index, double nZ0[], int nZ0index,
                double nX1[], int nX1index, double nY1[], int nY1index, double nZ1[], int nZ1index,
                double nX2[], int nX2index, double nY2[], int nY2index, double nZ2[], int nZ2index,
                double nX3[], int nX3index, double nY3[], int nY3index, double nZ3[], int nZ3index,
                double x0, double y0, double z0,
                double x1, double y1, double z1,
                double x2, double y2, double z2,
                double x3, double y3, double z3);
        static void sumElemStressesToNodeForces(
                double B[],
                double xxStress, double yyStress, double zzStress,
                double fx[], int xIndex,
                double fy[], int yIndex,
                double fz[], int zIndex) ;
        static void calcHourglassControlForElems(
                Domain& domain, double determinant[], double hgcoef,
                double dvdx[], double dvdy[], double dvdz[],
                double x8n[], double y8n[], double z8n[],
                double fxElem[], double fyElem[], double fzElem[]);
        static void collectDomainNodesToElemNodes(
                Domain& domain, int elemToNode[],
                double elemX[], double elemY[], double elemZ[]) ;
        static void calcElemVolumeDerivative(
                double dvdx[], double dvdy[], double dvdz[],
                double x[], double y[], double z[]);
        static void voluDer(
                double x0, double x1, double x2,
                double x3, double x4, double x5,
                double y0, double y1, double y2,
                double y3, double y4, double y5,
                double z0, double z1, double z2,
                double z3, double z4, double z5,
                double dvdx[], int dvdxIndex,
                double dvdy[], int dvdyIndex,
                double dvdz[], int dvdzIndex);
        static void calcFBHourglassForceForElems(
                Domain& domain,
                double determinant[],
                double x8n[], double y8n[], double z8n[],
                double dvdx[], double dvdy[], double dvdz[],
                double hourglass, 
                double fxElem[], double fyElem[], double fzElem[]);
        static void calcElemFBHourglassForce(
                double xd[], double yd[], double zd[], double hourglassAm0[],
                double hourglassAm1[], double hourglassAm2[], double hourglassAm3[],
                double hourglassAm4[], double hourglassAm5[], double hourglassAm6[],
                double hourglassAm7[], double coefficient,
                double hgfx[], double hgfy[], double hgfz[]);
        static void calcAccelerationForNodes(Domain& domain) ;
        static void applyAccelerationBoundaryConditionsForNodes(Domain& domain) ;
        static void calcVelocityForNodes(Domain& domain, double dt, double u_cut);
        static void calcPositionForNodes(Domain& domain, double dt) ;
        static void lagrangeElements(Domain& domain);
        static void calcLagrangeElements(Domain& domain, double deltatime);
        static void calcKinematicsForElems(Domain& domain, int numElem, double dt);
        static double calcElemVolume(
                double x0, double x1,
                double x2, double x3,
                double x4, double x5,
                double x6, double x7,
                double y0, double y1,
                double y2, double y3,
                double y4, double y5,
                double y6, double y7,
                double z0, double z1,
                double z2, double z3,
                double z4, double z5,
                double z6, double z7) ;
        static double tripleProduct(
                double x1, double y1, double z1,
                double x2, double y2, double z2,
                double x3, double y3, double z3) ;
        static double calcElemCharacteristicLength(double x[], double y[], double z[], double volume);
        static double areaFace(
                double x0, double x1,
                double x2, double x3,
                double y0, double y1,
                double y2, double y3,
                double z0, double z1,
                double z2, double z3);
        static void calcElemShapeFunctionDerivatives(
                double x[], double y[], double z[], double b[],
                double determinant[], int determinantIndex) ;
        static void calcElemVelocityGrandient(
                double xVel[], double yVel[], double zVel[],
                double b[], double detJ, double d[]) ;
        static void calcQForElems(Domain& domain);
        static void calcMonotonicQGradientsForElems(Domain& domain);
        static void calcMonotonicQForElems(Domain& domain);
        static void calcMonotonicQRegionForElems(
                Domain& domain, double qlcMonoq, double qqcMonoq,
                double monoqLimiterMult, double monoqMaxSlope,
                double ptiny, int elength);

        static void applyMaterialPropertiesForElems(Domain& domain);
        static void evalEOSForElems(Domain& domain, double vnewc[], int length);
        static void calcEnergyForElems(double p_new[], double e_new[], double q_new[],
                                                double bvc[], double pbvc[],
                                                double p_old[], double e_old[], double q_old[],
                                                double compression[], double compHalfStep[],
                                                double vnewc[], double work[], double delvc[], double pmin,
                                                double p_cut, double e_cut, double q_cut, double emin,
                                                double qq[], double ql[],
                                                double rho0,
                                                double eosvmax,
                                                int length) ;
        static void calcPressureForElems(
                double p_new[], double bvc[],
                double pbvc[], double e_old[],
                double compression[], double vnewc[],
                double pmin,
                double p_cut, double eosvmax,
                int length) ;
        static void calcSoundSpeedForElems(
                Domain& domain, double vnewc[], double rho0, double enewc[],
                double pnewc[], double pbvc[], double bvc[], double ss4o3, int nz);
        static void updateVolumesForElems(Domain& domain);
        static void calcTimeConstraintsForElems(Domain& domain);

        static void calcCourantConstraintForElems(Domain& domain);
        static void calcHydroConstraintForElems(Domain& domain);
        static int _compare(double a, double b, double epsilon=0.001);
};

double *LuleshParBenchmark::pHalfStep;
double *LuleshParBenchmark::vnewc;
double *LuleshParBenchmark::e_old;
double *LuleshParBenchmark::delvc;
double *LuleshParBenchmark::p_old;
double *LuleshParBenchmark::q_old;
double *LuleshParBenchmark::compression;
double *LuleshParBenchmark::compHalfStep;
double *LuleshParBenchmark::qq;
double *LuleshParBenchmark::ql;
double *LuleshParBenchmark::work;
double *LuleshParBenchmark::p_new;
double *LuleshParBenchmark::e_new;
double *LuleshParBenchmark::q_new;
double *LuleshParBenchmark::bvc;
double *LuleshParBenchmark::pbvc;

void Domain::_fill(double *arr, int size, double val){
    for (int i=0;i<size; i++){
        arr[i] = val;
    }
}

void Domain::AllocateElemPersistent(int size){
    
    m_matElemlist = hclib::numa_malloc<int>(size);
    m_nodelist = hclib::numa_malloc<int>(size*8);

    m_lxim = hclib::numa_malloc<int>(size);
    m_lxip = hclib::numa_malloc<int>(size);
    m_letam = hclib::numa_malloc<int>(size);
    m_letap = hclib::numa_malloc<int>(size);
    m_lzetam = hclib::numa_malloc<int>(size);
    m_lzetap = hclib::numa_malloc<int>(size);

    m_elemBC = hclib::numa_malloc<int>(size);

    m_e = hclib::numa_malloc<double>(size);
    
    m_p = hclib::numa_malloc<double>(size);
    m_q = hclib::numa_malloc<double>(size);
    m_ql = hclib::numa_malloc<double>(size);
    m_qq = hclib::numa_malloc<double>(size);

    m_v = hclib::numa_malloc<double>(size);
    _fill( m_v, size, 1.0);

    m_volo = hclib::numa_malloc<double>(size);
    m_delv = hclib::numa_malloc<double>(size);
    m_vdov = hclib::numa_malloc<double>(size);

    m_arealg = hclib::numa_malloc<double>(size);

    m_ss = hclib::numa_malloc<double>(size);

    m_elemMass = hclib::numa_malloc<double>(size);
}

/* Temporaries should not be initialized in bulk but */
/* this is a runnable placeholder for now */
void Domain::AllocateElemTemporary(int size) {
    m_dxx = hclib::numa_malloc<double>(size);
    m_dyy = hclib::numa_malloc<double>(size);
    m_dzz = hclib::numa_malloc<double>(size);

    m_delv_xi = hclib::numa_malloc<double>(size);
    m_delv_eta = hclib::numa_malloc<double>(size);
    m_delv_zeta = hclib::numa_malloc<double>(size);

    m_delx_xi = hclib::numa_malloc<double>(size);
    m_delx_eta = hclib::numa_malloc<double>(size);
    m_delx_zeta = hclib::numa_malloc<double>(size);

    m_vnew = hclib::numa_malloc<double>(size);
}

void Domain::AllocateNodalPersistent(int size) {
    m_x = hclib::numa_malloc<double>(size); 
    m_y = hclib::numa_malloc<double>(size);
    m_z = hclib::numa_malloc<double>(size);

    m_xd = hclib::numa_malloc<double>(size);
    m_yd = hclib::numa_malloc<double>(size);
    m_zd = hclib::numa_malloc<double>(size);

    m_xdd = hclib::numa_malloc<double>(size);
    m_ydd = hclib::numa_malloc<double>(size);
    m_zdd = hclib::numa_malloc<double>(size);

    m_fx = hclib::numa_malloc<double>(size);
    m_fy = hclib::numa_malloc<double>(size);
    m_fz = hclib::numa_malloc<double>(size);

    m_nodalMass = hclib::numa_malloc<double>(size);
}


/**
 * <p>AllocateNodeElemIndexes.</p>
 */
void Domain::AllocateNodeElemIndexes() {
    int numElem = Domain::numElem();
    int numNode = Domain::numNode();
    /* set up node-centered indexing of elements */
    m_nodeElemCount = hclib::numa_malloc<int>(numNode);

    for (int i = 0; i < numNode; ++i) {
        nodeElemCount(i, 0);
    }

    for (int i = 0; i < numElem; ++i) {
        int *nl = nodelist();
        for (int j = i*8+0; j < i*8+8; ++j) {
            int index = nl[j];
            nodeElemCount(index, nodeElemCount(index) + 1);
        }
    }

    m_nodeElemStart = hclib::numa_malloc<int>(numNode);
    nodeElemStart(0, 0);

    for (int i = 1; i < numNode; ++i) {
        nodeElemStart(i, nodeElemStart(i - 1) + nodeElemCount(i - 1));
    }

    int clSize = nodeElemStart(numNode - 1) + nodeElemCount(numNode - 1);
    m_nodeElemCornerList = new int[clSize];

    for (int i = 0; i < numNode; ++i) {
        nodeElemCount(i, 0);
    }

    for (int i = 0; i < numElem; ++i) {
        int *nl = nodelist();
        for (int j = i*8+0; j < i*8+8; ++j) {
            int m = nl[j];  
            int offset = nodeElemStart(m) + nodeElemCount(m);

            nodeElemCornerList(offset, j);
            nodeElemCount(m, nodeElemCount(m) + 1);
        }
    }
    for (int i = 0; i < clSize; ++i) {
        int clv = nodeElemCornerList(i);
        if ((clv < 0) || (clv > numElem * 8)) {
            cerr << "AllocateNodeElemIndexes(): nodeElemCornerList entry out of range!";
            exit(1);
        }
    }
}


/**
 * <p>nodeElemCornerList.</p>
 *
 * @param idx a int.
 * @return a int.
 */
int Domain::nodeElemCornerList(int idx) {
    return m_nodeElemCornerList[idx];
}

/**
 * <p>nodeElemCornerList.</p>
 *
 * @param idx a int.
 * @param val a int.
 */
void Domain::nodeElemCornerList(int idx, int val) {
    m_nodeElemCornerList[idx] = val;
}

/**
 * <p>nodeElemCount.</p>
 *
 * @param idx a int.
 * @return a int.
 */
int Domain::nodeElemCount(int idx) {
    return m_nodeElemCount[idx];
}

/**
 * <p>nodeElemCount.</p>
 *
 * @param idx a int.
 * @param val a int.
 */
void Domain::nodeElemCount(int idx, int val) {
    m_nodeElemCount[idx] = val;
}

/**
 * <p>nodeElemStart.</p>
 *
 * @param idx a int.
 * @return a int.
 */
int Domain::nodeElemStart(int idx) {
    return m_nodeElemStart[idx];
}

/**
 * <p>nodeElemStart.</p>
 *
 * @param idx a int.
 * @param val a int.
 */
void Domain::nodeElemStart(int idx, int val) {
    m_nodeElemStart[idx] = val;
}

/**
 * <p>nodelist.</p>
 *
 * @return an array of int.
 */
int* Domain::nodelist() {
    return m_nodelist;
}

/**
 * <p>AllocateNodesets.</p>
 *
 * @param size a int.
 */
void Domain::AllocateNodesets(int size) {
    m_symmX = hclib::numa_malloc<int>(size);
    m_symmY = hclib::numa_malloc<int>(size);
    m_symmZ = hclib::numa_malloc<int>(size);
}

/**
 * <p>arealg.</p>
 *
 * @param idx a int.
 * @return a double.
 */
double Domain::arealg(int idx) {
    return m_arealg[idx];
}

/**
 * <p>arealg.</p>
 *
 * @param idx a int.
 * @param val a double.
 */
void Domain::arealg(int idx, double val) {
    m_arealg[idx] = val;
}

/**
 * <p>cycle.</p>
 *
 * @return a int.
 */
int Domain::cycle() {
    return m_cycle;
}

/**
 * <p>cycle.</p>
 *
 * @param val a int.
 */
void Domain::cycle(int val) {
    m_cycle = val;
}

/**
 * <p>deltatime.</p>
 *
 * @return a double.
 */
double Domain::deltatime() {
    return m_deltatime;
}

/**
 * <p>deltatime.</p>
 *
 * @param val a double.
 */
void Domain::deltatime(double val) {
    m_deltatime = val;
}

/**
 * <p>deltatimemultlb.</p>
 *
 * @return a double.
 */
double Domain::deltatimemultlb() {
    return m_deltatimemultlb;
}

/**
 * <p>deltatimemultlb.</p>
 *
 * @param val a double.
 */
void Domain::deltatimemultlb(double val) {
    m_deltatimemultlb = val;
}

/**
 * <p>deltatimemultub.</p>
 *
 * @return a double.
 */
double Domain::deltatimemultub() {
    return m_deltatimemultub;
}

/**
 * <p>deltatimemultub.</p>
 *
 * @param val a double.
 */
void Domain::deltatimemultub(double val) {
    m_deltatimemultub = val;
}

/**
 * <p>delv.</p>
 *
 * @param idx a int.
 * @return a double.
 */
double Domain::delv(int idx) {
    return m_delv[idx];
}

/**
 * <p>delv.</p>
 *
 * @param idx a int.
 * @param val a double.
 */
void Domain::delv(int idx, double val) {
    m_delv[idx] = val;

}

/**
 * <p>delv_eta.</p>
 *
 * @param idx a int.
 * @return a double.
 */
double Domain::delv_eta(int idx) {
    return m_delv_eta[idx];
}

/**
 * <p>delv_eta.</p>
 *
 * @param idx a int.
 * @param val a double.
 */
void Domain::delv_eta(int idx, double val) {
    m_delv_eta[idx] = val;
}

/**
 * <p>delv_xi.</p>
 *
 * @param idx a int.
 * @return a double.
 */
double Domain::delv_xi(int idx) {
    return m_delv_xi[idx];
}

/**
 * <p>delv_xi.</p>
 *
 * @param idx a int.
 * @param val a double.
 */
void Domain::delv_xi(int idx, double val) {
    m_delv_xi[idx] = val;
}

/**
 * <p>delv_zeta.</p>
 *
 * @param idx a int.
 * @return a double.
 */
double Domain::delv_zeta(int idx) {
    return m_delv_zeta[idx];
}

/**
 * <p>delv_zeta.</p>
 *
 * @param idx a int.
 * @param val a double.
 */
void Domain::delv_zeta(int idx, double val) {
    m_delv_zeta[idx] = val;
}

/**
 * <p>delx_eta.</p>
 *
 * @param idx a int.
 * @return a double.
 */
double Domain::delx_eta(int idx) {
    return m_delx_eta[idx];
}

/**
 * <p>delx_eta.</p>
 *
 * @param idx a int.
 * @param val a double.
 */
void Domain::delx_eta(int idx, double val) {
    m_delx_eta[idx] = val;
}

/**
 * <p>delx_xi.</p>
 *
 * @param idx a int.
 * @return a double.
 */
double Domain::delx_xi(int idx) {
    return m_delx_xi[idx];
}

/**
 * <p>delx_xi.</p>
 *
 * @param idx a int.
 * @param val a double.
 */
void Domain::delx_xi(int idx, double val) {
    m_delx_xi[idx] = val;
}

/**
 * <p>delx_zeta.</p>
 *
 * @param idx a int.
 * @return a double.
 */
double Domain::delx_zeta(int idx) {
    return m_delx_zeta[idx];
}

/**
 * <p>delx_zeta.</p>
 *
 * @param idx a int.
 * @param val a double.
 */
void Domain::delx_zeta(int idx, double val) {
    m_delx_zeta[idx] = val;
}

/**
 * <p>dtcourant.</p>
 *
 * @return a double.
 */
double Domain::dtcourant() {
    return m_dtcourant;
}

/**
 * <p>dtcourant.</p>
 *
 * @param dtcourant a double.
 */
void Domain::dtcourant(double dtcourant) {
    m_dtcourant = dtcourant;
}

/* Params */

double Domain::dtfixed() {
    return m_dtfixed;
}

/**
 * <p>dtfixed.</p>
 *
 * @param val a double.
 */
void Domain::dtfixed(double val) {
    m_dtfixed = val;
}

/**
 * <p>dthydro.</p>
 *
 * @return a double.
 */
double Domain::dthydro() {
    return m_dthydro;
}

/**
 * <p>dthydro.</p>
 *
 * @param dthydro a double.
 */
void Domain::dthydro(double dthydro) {
    m_dthydro = dthydro;
}

/**
 * <p>dtmax.</p>
 *
 * @return a double.
 */
double Domain::dtmax() {
    return m_dtmax;
}

/**
 * <p>dtmax.</p>
 *
 * @param val a double.
 */
void Domain::dtmax(double val) {
    m_dtmax = val;
}

/**
 * <p>dvovmax.</p>
 *
 * @return a double.
 */
double Domain::dvovmax() {
    return m_dvovmax;
}

/**
 * <p>dvovmax.</p>
 *
 * @param val a double.
 */
void Domain::dvovmax(double val) {
    m_dvovmax = val;
}

/**
 * <p>dxx.</p>
 *
 * @param idx a int.
 * @return a double.
 */
double Domain::dxx(int idx) {
    return m_dxx[idx];
}

/**
 * <p>dxx.</p>
 *
 * @param idx a int.
 * @param val a double.
 */
void Domain::dxx(int idx, double val) {
    m_dxx[idx] = val;
}

/**
 * <p>dyy.</p>
 *
 * @param idx a int.
 * @return a double.
 */
double Domain::dyy(int idx) {
    return m_dyy[idx];
}

/**
 * <p>dyy.</p>
 *
 * @param idx a int.
 * @param val a double.
 */
void Domain::dyy(int idx, double val) {
    m_dyy[idx] = val;
}

/**
 * <p>dzz.</p>
 *
 * @param idx a int.
 * @return a double.
 */
double Domain::dzz(int idx) {
    return m_dzz[idx];
}

/**
 * <p>dzz.</p>
 *
 * @param idx a int.
 * @param val a double.
 */
void Domain::dzz(int idx, double val) {
    m_dzz[idx] = val;
}

/**
 * <p>e.</p>
 *
 * @param idx a int.
 * @return a double.
 */
double Domain::e(int idx) {
    return m_e[idx];
}

/**
 * <p>e.</p>
 *
 * @param idx a int.
 * @param val a double.
 */
void Domain::e(int idx, double val) {
    m_e[idx] = val;
}

/**
 * <p>e_cut.</p>
 *
 * @return a double.
 */
double Domain::e_cut() {
    return m_e_cut;
}

/**
 * <p>e_cut.</p>
 *
 * @param val a double.
 */
void Domain::e_cut(double val) {
    m_e_cut = val;
}

/**
 * <p>elemBC.</p>
 *
 * @param idx a int.
 * @return a int.
 */
int Domain::elemBC(int idx) {
    return m_elemBC[idx];
}

/**
 * <p>elemBC.</p>
 *
 * @param idx a int.
 * @param val a int.
 */
void Domain::elemBC(int idx, int val) {
    m_elemBC[idx] = val;
}

/**
 * <p>elemMass.</p>
 *
 * @param idx a int.
 * @return a double.
 */
double Domain::elemMass(int idx) {
    return m_elemMass[idx];
}

/**
 * <p>elemMass.</p>
 *
 * @param idx a int.
 * @param val a double.
 */
void Domain::elemMass(int idx, double val) {
    m_elemMass[idx] = val;
}

/**
 * <p>emin.</p>
 *
 * @return a double.
 */
double Domain::emin() {
    return m_emin;
}

/**
 * <p>emin.</p>
 *
 * @param val a double.
 */
void Domain::emin(double val) {
    m_emin = val;
}

/**
 * <p>eosvmax.</p>
 *
 * @return a double.
 */
double Domain::eosvmax() {
    return m_eosvmax;
}

/**
 * <p>eosvmax.</p>
 *
 * @param val a double.
 */
void Domain::eosvmax(double val) {
    m_eosvmax = val;
}

/**
 * <p>eosvmin.</p>
 *
 * @return a double.
 */
double Domain::eosvmin() {
    return m_eosvmin;
}

/**
 * <p>eosvmin.</p>
 *
 * @param val a double.
 */
void Domain::eosvmin(double val) {
    m_eosvmin = val;
}

/**
 * <p>fx.</p>
 *
 * @param idx a int.
 * @return a double.
 */
double Domain::fx(int idx) {
    return m_fx[idx];
}

/**
 * <p>fx.</p>
 *
 * @param idx a int.
 * @param val a double.
 */
void Domain::fx(int idx, double val) {
    m_fx[idx] = val;
}

/**
 * <p>fxFill.</p>
 *
 * @param val a double.
 * size an int. size of array m_fx.
 */
void Domain::fxFill(double val, int size) {
    _fill(m_fx, size, val);
}

/**
 * <p>fy.</p>
 *
 * @param idx a int.
 * @return a double.
 */
double Domain::fy(int idx) {
    return m_fy[idx];
}

/**
 * <p>fy.</p>
 *
 * @param idx a int.
 * @param val a double.
 */
void Domain::fy(int idx, double val) {
    m_fy[idx] = val;
}

/**
 * <p>fyFill.</p>
 *
 * @param val a double.
 * size a int. size of array m_fy.
 */
void Domain::fyFill(double val, int size) {
    _fill(m_fy, size, val);
}

/**
 * <p>fz.</p>
 *
 * @param idx a int.
 * @return a double.
 */
double Domain::fz(int idx) {
    return m_fz[idx];
}

/**
 * <p>fz.</p>
 *
 * @param idx a int.
 * @param val a double.
 */
void Domain::fz(int idx, double val) {
    m_fz[idx] = val;
}

/**
 * <p>fzFill.</p>
 *
 * @param val a double.
 * size an int. size of array m_fz
 */
void Domain::fzFill(double val, int size) {
    _fill(m_fz, size, val);
}

/**
 * <p>hgcoef.</p>
 *
 * @return a double.
 */
double Domain::hgcoef() {
    return m_hgcoef;
}

/**
 * <p>hgcoef.</p>
 *
 * @param val a double.
 */
void Domain::hgcoef(double val) {
    m_hgcoef = val;
}

/**
 * <p>letam.</p>
 *
 * @param idx a int.
 * @return a int.
 */
int Domain::letam(int idx) {
    return m_letam[idx];
}

/**
 * <p>letam.</p>
 *
 * @param idx a int.
 * @param val a int.
 */
void Domain::letam(int idx, int val) {
    m_letam[idx] = val;
}

/**
 * <p>letap.</p>
 *
 * @param idx a int.
 * @return a int.
 */
int Domain::letap(int idx) {
    return m_letap[idx];
}

/**
 * <p>letap.</p>
 *
 * @param idx a int.
 * @param val a int.
 */
void Domain::letap(int idx, int val) {
    m_letap[idx] = val;
}

/**
 * <p>lxim.</p>
 *
 * @param idx a int.
 * @return a int.
 */
int Domain::lxim(int idx) {
    return m_lxim[idx];
}

/**
 * <p>lxim.</p>
 *
 * @param idx a int.
 * @param val a int.
 */
void Domain::lxim(int idx, int val) {
    m_lxim[idx] = val;
}

/**
 * <p>lxip.</p>
 *
 * @param idx a int.
 * @return a int.
 */
int Domain::lxip(int idx) {
    return m_lxip[idx];
}

/**
 * <p>lxip.</p>
 *
 * @param idx a int.
 * @param val a int.
 */
void Domain::lxip(int idx, int val) {
    m_lxip[idx] = val;
}

/**
 * <p>lzetam.</p>
 *
 * @param idx a int.
 * @return a int.
 */
int Domain::lzetam(int idx) {
    return m_lzetam[idx];
}

/**
 * <p>lzetam.</p>
 *
 * @param idx a int.
 * @param val a int.
 */
void Domain::lzetam(int idx, int val) {
    m_lzetam[idx] = val;
}

/**
 * <p>lzetap.</p>
 *
 * @param idx a int.
 * @return a int.
 */
int Domain::lzetap(int idx) {
    return m_lzetap[idx];
}

/**
 * <p>lzetap.</p>
 *
 * @param idx a int.
 * @param val a int.
 */
void Domain::lzetap(int idx, int val) {
    m_lzetap[idx] = val;
}

/* Element-centered */

int Domain::matElemlist(int idx) {
    return m_matElemlist[idx];
}

/**
 * <p>matElemlist.</p>
 *
 * @param idx a int.
 * @param val a int.
 */
void Domain::matElemlist(int idx, int val) {
    m_matElemlist[idx] = val;
}

/**
 * <p>monoq_limiter_mult.</p>
 *
 * @return a double.
 */
double Domain::monoq_limiter_mult() {
    return m_monoq_limiter_mult;
}

/**
 * <p>monoq_limiter_mult.</p>
 *
 * @param val a double.
 */
void Domain::monoq_limiter_mult(double val) {
    m_monoq_limiter_mult = val;
}

/**
 * <p>monoq_max_slope.</p>
 *
 * @return a double.
 */
double Domain::monoq_max_slope() {
    return m_monoq_max_slope;
}

/**
 * <p>monoq_max_slope.</p>
 *
 * @param val a double.
 */
void Domain::monoq_max_slope(double val) {
    m_monoq_max_slope = val;
}

/**
 * <p>nodalMass.</p>
 *
 * @param idx a int.
 * @return a double.
 */
double Domain::nodalMass(int idx) {
    return m_nodalMass[idx];
}

/**
 * <p>nodalMass.</p>
 *
 * @param idx a int.
 * @param val a double.
 */
void Domain::nodalMass(int idx, double val) {
    m_nodalMass[idx] = val;
}

/**
 * <p>numElem.</p>
 *
 * @return a int.
 */
int Domain::numElem() {
    return m_numElem;
}

/**
 * <p>numElem.</p>
 *
 * @param val a int.
 */
void Domain::numElem(int val) {
    m_numElem = val;
}

/**
 * <p>numNode.</p>
 *
 * @param val a int.
 */
void Domain::numNode(int val) {
    m_numNode = val;
}

/**
 * <p>numNode.</p>
 *
 * @return a int.
 */
int Domain::numNode() {
    return m_numNode;
}

/**
 * <p>p.</p>
 *
 * @param idx a int.
 * @return a double.
 */
double Domain::p(int idx) {
    return m_p[idx];
}

/**
 * <p>p.</p>
 *
 * @param idx a int.
 * @param val a double.
 */
void Domain::p(int idx, double val) {
    m_p[idx] = val;
}

/**
 * <p>p_cut.</p>
 *
 * @return a double.
 */
double Domain::p_cut() {
    return m_p_cut;
}

/**
 * <p>p_cut.</p>
 *
 * @param val a double.
 */
void Domain::p_cut(double val) {
    m_p_cut = val;
}

/**
 * <p>pmin.</p>
 *
 * @return a double.
 */
double Domain::pmin() {
    return m_pmin;
}

/**
 * <p>pmin.</p>
 *
 * @param val a double.
 */
void Domain::pmin(double val) {
    m_pmin = val;
}

/**
 * <p>q.</p>
 *
 * @param idx a int.
 * @return a double.
 */
double Domain::q(int idx) {
    return m_q[idx];
}

/**
 * <p>q.</p>
 *
 * @param idx a int.
 * @param val a double.
 */
void Domain::q(int idx, double val) {
    m_q[idx] = val;
}

/**
 * <p>q_cut.</p>
 *
 * @return a double.
 */
double Domain::q_cut() {
    return m_q_cut;
}

/**
 * <p>q_cut.</p>
 *
 * @param val a double.
 */
void Domain::q_cut(double val) {
    m_q_cut = val;
}

/**
 * <p>ql.</p>
 *
 * @param idx a int.
 * @return a double.
 */
double Domain::ql(int idx) {
    return m_ql[idx];
}

/**
 * <p>ql.</p>
 *
 * @param idx a int.
 * @param val a double.
 */
void Domain::ql(int idx, double val) {
    m_ql[idx] = val;
}

/**
 * <p>qlc_monoq.</p>
 *
 * @return a double.
 */
double Domain::qlc_monoq() {
    return m_qlc_monoq;
}

/**
 * <p>qlc_monoq.</p>
 *
 * @param val a double.
 */
void Domain::qlc_monoq(double val) {
    m_qlc_monoq = val;
}

/**
 * <p>qq.</p>
 *
 * @param idx a int.
 * @return a double.
 */
double Domain::qq(int idx) {
    return m_qq[idx];
}

/**
 * <p>qq.</p>
 *
 * @param idx a int.
 * @param val a double.
 */
void Domain::qq(int idx, double val) {
    m_qq[idx] = val;
}

/**
 * <p>qqc.</p>
 *
 * @return a double.
 */
double Domain::qqc() {
    return m_qqc;
}

/**
 * <p>qqc.</p>
 *
 * @param val a double.
 */
void Domain::qqc(double val) {
    m_qqc = val;
}

/**
 * <p>qqc_monoq.</p>
 *
 * @return a double.
 */
double Domain::qqc_monoq() {
    return m_qqc_monoq;
}

/**
 * <p>qqc_monoq.</p>
 *
 * @param val a double.
 */
void Domain::qqc_monoq(double val) {
    m_qqc_monoq = val;
}

/**
 * <p>qstop.</p>
 *
 * @return a double.
 */
double Domain::qstop() {
    return m_qstop;
}

/**
 * <p>qstop.</p>
 *
 * @param val a double.
 */
void Domain::qstop(double val) {
    m_qstop = val;
}

/**
 * <p>refdens.</p>
 *
 * @return a double.
 */
double Domain::refdens() {
    return m_refdens;
}

/**
 * <p>refdens.</p>
 *
 * @param val a double.
 */
void Domain::refdens(double val) {
    m_refdens = val;
}

/**
 * <p>sizeX.</p>
 *
 * @return a int.
 */
int Domain::sizeX() {
    return m_sizeX;
}

/**
 * <p>sizeX.</p>
 *
 * @param val a int.
 */
void Domain::sizeX(int val) {
    m_sizeX = val;
}

/**
 * <p>sizeY.</p>
 *
 * @return a int.
 */
int Domain::sizeY() {
    return m_sizeY;
}

/**
 * <p>sizeY.</p>
 *
 * @param val a int.
 */
void Domain::sizeY(int val) {
    m_sizeY = val;
}

/**
 * <p>sizeZ.</p>
 *
 * @return a int.
 */
int Domain::sizeZ() {
    return m_sizeZ;
}

/**
 * <p>sizeZ.</p>
 *
 * @param val a int.
 */
void Domain::sizeZ(int val) {
    m_sizeZ = val;
}

/**
 * <p>ss.</p>
 *
 * @param idx a int.
 * @return a double.
 */
double Domain::ss(int idx) {
    return m_ss[idx];
}

/**
 * <p>ss.</p>
 *
 * @param idx a int.
 * @param val a double.
 */
void Domain::ss(int idx, double val) {
    m_ss[idx] = val;
}

/**
 * <p>ss4o3.</p>
 *
 * @return a double.
 */
double Domain::ss4o3() {
    return m_ss4o3;
}

/**
 * <p>ss4o3.</p>
 *
 * @param val a double.
 */
void Domain::ss4o3(double val) {
    m_ss4o3 = val;
}

/**
 * <p>stoptime.</p>
 *
 * @return a double.
 */
double Domain::stoptime() {
    return m_stoptime;
}

/**
 * <p>stoptime.</p>
 *
 * @param val a double.
 */
void Domain::stoptime(double val) {
    m_stoptime = val;
}

/**
 * <p>symmX.</p>
 *
 * @param idx a int.
 * @return a int.
 */
int Domain::symmX(int idx) {
    return m_symmX[idx];
}

/**
 * <p>symmX.</p>
 *
 * @param idx a int.
 * @param val a int.
 */
void Domain::symmX(int idx, int val) {
    m_symmX[idx] = val;
}

/**
 * <p>symmY.</p>
 *
 * @param idx a int.
 * @return a int.
 */
int Domain::symmY(int idx) {
    return m_symmY[idx];
}

/**
 * <p>symmY.</p>
 *
 * @param idx a int.
 * @param val a int.
 */
void Domain::symmY(int idx, int val) {
    m_symmY[idx] = val;
}

/**
 * <p>symmZ.</p>
 *
 * @param idx a int.
 * @return a int.
 */
int Domain::symmZ(int idx) {
    return m_symmZ[idx];
}

/**
 * <p>symmZ.</p>
 *
 * @param idx a int.
 * @param val a int.
 */
void Domain::symmZ(int idx, int val) {
    m_symmZ[idx] = val;
}

/**
 * <p>time.</p>
 *
 * @return a double.
 */
double Domain::time() {
    return m_time;
}

/**
 * <p>time.</p>
 *
 * @param val a double.
 */
void Domain::time(double val) {
    m_time = val;
}

/**
 * <p>u_cut.</p>
 *
 * @return a double.
 */
double Domain::u_cut() {
    return m_u_cut;
}

/**
 * <p>u_cut.</p>
 *
 * @param val a double.
 */
void Domain::u_cut(double val) {
    m_u_cut = val;
}

/**
 * <p>v.</p>
 *
 * @param idx a int.
 * @return a double.
 */
double Domain::v(int idx) {
    return m_v[idx];
}

/**
 * <p>v.</p>
 *
 * @param idx a int.
 * @param val a double.
 */
void Domain::v(int idx, double val) {
    m_v[idx] = val;
}

/**
 * <p>v_cut.</p>
 *
 * @return a double.
 */
double Domain::v_cut() {
    return m_v_cut;
}

/**
 * <p>v_cut.</p>
 *
 * @param val a double.
 */
void Domain::v_cut(double val) {
    m_v_cut = val;
}

/**
 * <p>vdov.</p>
 *
 * @param idx a int.
 * @return a double.
 */
double Domain::vdov(int idx) {
    return m_vdov[idx];
}

/**
 * <p>vdov.</p>
 *
 * @param idx a int.
 * @param val a double.
 */
void Domain::vdov(int idx, double val) {
    m_vdov[idx] = val;
}

/**
 * <p>vnew.</p>
 *
 * @param idx a int.
 * @return a double.
 */
double Domain::vnew(int idx) {
    return m_vnew[idx];
}

/**
 * <p>vnew.</p>
 *
 * @param idx a int.
 * @param val a double.
 */
void Domain::vnew(int idx, double val) {
    m_vnew[idx] = val;
}

/**
 * <p>volo.</p>
 *
 * @param idx a int.
 * @return a double.
 */
double Domain::volo(int idx) {
    return m_volo[idx];
}

/**
 * <p>volo.</p>
 *
 * @param idx a int.
 * @param val a double.
 */
void Domain::volo(int idx, double val) {
    m_volo[idx] = val;
}

/* Node-centered */

double Domain::x(int idx) {
    return m_x[idx];
}

/**
 * <p>x.</p>
 *
 * @param idx a int.
 * @param val a double.
 */
void Domain::x(int idx, double val) {
    m_x[idx] = val;
}

/**
 * <p>xd.</p>
 *
 * @param idx a int.
 * @return a double.
 */
double Domain::xd(int idx) {
    return m_xd[idx];
}

/**
 * <p>xd.</p>
 *
 * @param idx a int.
 * @param val a double.
 */
void Domain::xd(int idx, double val) {
    m_xd[idx] = val;
}

/**
 * <p>xdd.</p>
 *
 * @param idx a int.
 * @return a double.
 */
double Domain::xdd(int idx) {
    return m_xdd[idx];
}

/**
 * <p>xdd.</p>
 *
 * @param idx a int.
 * @param val a double.
 */
void Domain::xdd(int idx, double val) {
    m_xdd[idx] = val;
}

/**
 * <p>y.</p>
 *
 * @param idx a int.
 * @return a double.
 */
double Domain::y(int idx) {
    return m_y[idx];
}

/**
 * <p>y.</p>
 *
 * @param idx a int.
 * @param val a double.
 */
void Domain::y(int idx, double val) {
    m_y[idx] = val;
}

/**
 * <p>yd.</p>
 *
 * @param idx a int.
 * @return a double.
 */
double Domain::yd(int idx) {
    return m_yd[idx];
}

/**
 * <p>yd.</p>
 *
 * @param idx a int.
 * @param val a double.
 */
void Domain::yd(int idx, double val) {
    m_yd[idx] = val;
}

/**
 * <p>ydd.</p>
 *
 * @param idx a int.
 * @return a double.
 */
double Domain::ydd(int idx) {
    return m_ydd[idx];
}

/**
 * <p>ydd.</p>
 *
 * @param idx a int.
 * @param val a double.
 */
void Domain::ydd(int idx, double val) {
    m_ydd[idx] = val;
}

/**
 * <p>z.</p>
 *
 * @param idx a int.
 * @return a double.
 */
double Domain::z(int idx) {
    return m_z[idx];
}

/**
 * <p>z.</p>
 *
 * @param idx a int.
 * @param val a double.
 */
void Domain::z(int idx, double val) {
    m_z[idx] = val;
}

/**
 * <p>zd.</p>
 *
 * @param idx a int.
 * @return a double.
 */
double Domain::zd(int idx) {
    return m_zd[idx];
}

/**
 * <p>zd.</p>
 *
 * @param idx a int.
 * @param val a double.
 */
void Domain::zd(int idx, double val) {
    m_zd[idx] = val;
}

/**
 * <p>zdd.</p>
 *
 * @param idx a int.
 * @return a double.
 */
double Domain::zdd(int idx) {
    return m_zdd[idx];
}

/**
 * <p>zdd.</p>
 *
 * @param idx a int.
 * @param val a double.
 */
void Domain::zdd(int idx, double val) {
    m_zdd[idx] = val;
}

const int LuleshConfig::VOLUME_ERROR = -1;
const int LuleshConfig::Q_STOP_ERROR = -2;

/* Stuff needed for boundary conditions, 2 BCs on each of 6 hexahedral faces (12 bits) */
const int LuleshConfig::XI_M = 0x003;
const int LuleshConfig::XI_M_SYMM = 0x001;
const int LuleshConfig::XI_M_FREE = 0x002;

const int LuleshConfig::XI_P = 0x00c;
const int LuleshConfig::XI_P_SYMM = 0x004;
const int LuleshConfig::XI_P_FREE = 0x008;

const int LuleshConfig::ETA_M = 0x030;
const int LuleshConfig::ETA_M_SYMM = 0x010;
const int LuleshConfig::ETA_M_FREE = 0x020;

const int LuleshConfig::ETA_P = 0x0c0;
const int LuleshConfig::ETA_P_SYMM = 0x040;
const int LuleshConfig::ETA_P_FREE = 0x080;

const int LuleshConfig::ZETA_M = 0x300;
const int LuleshConfig::ZETA_M_SYMM = 0x100;
const int LuleshConfig::ZETA_M_FREE = 0x200;

const int LuleshConfig::ZETA_P = 0xc00;
const int LuleshConfig::ZETA_P_SYMM = 0x400;
const int LuleshConfig::ZETA_P_FREE = 0x800;

int LuleshConfig::N = 32;
bool LuleshConfig::debug = false;

void LuleshConfig::parseArgs(int argc, char** argv) {
    int i = 0;
    while (i < argc) {
        string loopOptionKey = argv[i];
        if (loopOptionKey == "-n"){
            i += 1;
            LuleshConfig::N = stoi(argv[i]);
        }
        else if ((loopOptionKey == "-debug") ||
                    (loopOptionKey == "-verbose")){
            debug = true;
        }
        i += 1;
    }
}

// TODO: change the formatting
void LuleshConfig::printArgs(){
    cout<< N;
    cout<< debug;
}


void runIteration(int argc, char** argv) {
    LuleshConfig::parseArgs(argc, argv);
    LuleshParBenchmark::runLulesh(LuleshConfig::N);
}

int main(int argc, char** argv){
    hclib::launch([&](){
        runIteration(argc, argv);
    });
     
}
void LuleshParBenchmark::runLulesh(int edgeElems){
    int edgeNodes = edgeElems + 1;
    // double ds = 1.125/double(edgeElems) ; /* may accumulate roundoff */
        
    /* get run options to measure various metrics */

    /*   initialize Sedov Mesh  */

    /* construct a uniform box for this processor */
    Domain domain;
    
    domain.sizeX(edgeElems);
    domain.sizeY(edgeElems);
    domain.sizeZ(edgeElems);
    domain.numElem(edgeElems * edgeElems * edgeElems);
    domain.numNode(edgeNodes * edgeNodes * edgeNodes);

    int domElems = domain.numElem();

    /* allocate field memory */
    domain.AllocateElemPersistent(domain.numElem());
    domain.AllocateElemTemporary(domain.numElem());
    domain.AllocateNodalPersistent(domain.numNode());
    domain.AllocateNodesets(edgeNodes * edgeNodes);

    {
        double tx;
        double ty;
        double tz;

        int nIndex = 0;
        tz = 0.0;
        for (int plane = 0; plane < edgeNodes; ++plane) {
            ty = 0.0;
            for (int row = 0; row < edgeNodes; ++row) {
                tx = 0.0;
                for (int col = 0; col < edgeNodes; ++col) {
                    domain.x(nIndex, tx);
                    domain.y(nIndex, ty);
                    domain.z(nIndex, tz);
                    ++nIndex;
                    // tx += ds ; /* may accumulate round-off... */
                    tx = 1.125 * (col + 1) / (edgeElems);
                }
                // ty += ds ;  /* may accumulate round-off... */
                ty = 1.125 * (row + 1) / (edgeElems);
            }
            // tz += ds ;  /* may accumulate round-off... */
            tz = 1.125 * (plane + 1) / (edgeElems);
            }
    }


    {
        int zIndex = 0;
        int nIndex = 0;
        for (int plane = 0; plane < edgeElems; ++plane) {
            for (int row = 0; row < edgeElems; ++row) {
                for (int col = 0; col < edgeElems; ++col) {
                    int *localNode = domain.nodelist();
                    localNode[zIndex*8+0] = nIndex;
                    localNode[zIndex*8+1] = nIndex + 1;
                    localNode[zIndex*8+2] = nIndex + edgeNodes + 1;
                    localNode[zIndex*8+3] = nIndex + edgeNodes;
                    localNode[zIndex*8+4] = nIndex + edgeNodes * edgeNodes;
                    localNode[zIndex*8+5] = nIndex + edgeNodes * edgeNodes + 1;
                    localNode[zIndex*8+6] = nIndex + edgeNodes * edgeNodes + edgeNodes + 1;
                    localNode[zIndex*8+7] = nIndex + edgeNodes * edgeNodes + edgeNodes;
                    ++zIndex;
                    ++nIndex;
                }
                ++nIndex;
            }
            nIndex += edgeNodes;
        }
    }

    domain.AllocateNodeElemIndexes();


    /* Create a material indexSet (entire domain same material for now) */
    for (int i = 0; i < domElems; ++i) {
        domain.matElemlist(i, i);
    }

    /* initialize material parameters */
    domain.dtfixed(-1.0e-7);
    domain.deltatime(1.0e-7);
    domain.deltatimemultlb(1.1);
    domain.deltatimemultub(1.2);
    domain.stoptime(1.0e-2);
    domain.dtcourant(1.0e+20);
    domain.dthydro(1.0e+20);
    domain.dtmax(1.0e-2);
    domain.time(0.0);
    domain.cycle(0);

    domain.e_cut(1.0e-7);
    domain.p_cut(1.0e-7);
    domain.q_cut(1.0e-7);
    domain.u_cut(1.0e-7);
    domain.v_cut(1.0e-10);

    domain.hgcoef(3.0);
    domain.ss4o3(4.0 / 3.0);

    domain.qstop(1.0e+12);
    domain.monoq_max_slope(1.0);
    domain.monoq_limiter_mult(2.0);
    domain.qlc_monoq(0.5);
    domain.qqc_monoq(2.0 / 3.0);
    domain.qqc(2.0);

    domain.pmin(0.0);
    domain.emin(-1.0e+15);

    domain.dvovmax(0.1);

    domain.eosvmax(1.0e+9);
    domain.eosvmin(1.0e-9);

    domain.refdens(1.0);

    /* initialize field data */
    for (int i = 0; i < domElems; ++i) {
        double xLocal[8];
        double yLocal[8];
        double zLocal[8];
        int *elemToNode = domain.nodelist();
        for (int lnode = 0; lnode < 8; ++lnode) {
            int gnode = elemToNode[i*8+lnode];
            xLocal[lnode] = domain.x(gnode);
            yLocal[lnode] = domain.y(gnode);
            zLocal[lnode] = domain.z(gnode);
        }

        // volume calculations
        double volume = calcElemVolume(xLocal, yLocal, zLocal);
        domain.volo(i, volume);
        domain.elemMass(i, volume);
        for (int j = 0; j < 8; ++j) {
            int index = elemToNode[i*8+j];
            domain.nodalMass(index, domain.nodalMass(index) + (volume / 8.0));
        }
    }

    /* deposit energy */
    domain.e(0, 3.948746e+7);

    /* set up symmetry nodesets */
    {
        int nIndex = 0;
        for (int i = 0; i < edgeNodes; ++i) {
            int planeInc = i * edgeNodes * edgeNodes;
            int rowInc = i * edgeNodes;
            for (int j = 0; j < edgeNodes; ++j) {
                domain.symmX(nIndex, planeInc + j * edgeNodes);
                domain.symmY(nIndex, planeInc + j);
                domain.symmZ(nIndex, rowInc + j);
                ++nIndex;
            }
        }
    }

    /* set up element connectivity information */
    domain.lxim(0, 0);
    for (int i = 1; i < domElems; ++i) {
        domain.lxim(i, i - 1);
        domain.lxip(i - 1, i);
    }
    domain.lxip(domElems - 1, domElems - 1);

    for (int i = 0; i < edgeElems; ++i) {
        domain.letam(i, i);
        domain.letap(domElems - edgeElems + i, domElems - edgeElems + i);
    }
    for (int i = edgeElems; i < domElems; ++i) {
        domain.letam(i, i - edgeElems);
        domain.letap(i - edgeElems, i);
    }

    for (int i = 0; i < edgeElems * edgeElems; ++i) {
        domain.lzetam(i, i);
        domain.lzetap(domElems - edgeElems * edgeElems + i, domElems - edgeElems * edgeElems + i);
    }
    for (int i = edgeElems * edgeElems; i < domElems; ++i) {
        domain.lzetam(i, i - edgeElems * edgeElems);
        domain.lzetap(i - edgeElems * edgeElems, i);
    }

    /* set up boundary condition information */
    for (int i = 0; i < domElems; ++i) {
        /* clear BCs by default */
        domain.elemBC(i, 0);
    }

    /* faces on "external" boundaries will be symmetry plane or free surface BCs */
    for (int i = 0; i < edgeElems; ++i) {
        int planeInc = i * edgeElems * edgeElems;
        int rowInc = i * edgeElems;
        for (int j = 0; j < edgeElems; ++j) {
            int j1 = planeInc + j * edgeElems;
            domain.elemBC(j1, domain.elemBC(j1) | LuleshConfig::XI_M_SYMM);
            int j2 = planeInc + j * edgeElems + edgeElems - 1;
            domain.elemBC(j2, domain.elemBC(j2) | LuleshConfig::XI_P_FREE);
            int j3 = planeInc + j;
            domain.elemBC(j3, domain.elemBC(j3) | LuleshConfig::ETA_M_SYMM);
            int j4 = planeInc + j + edgeElems * edgeElems - edgeElems;
            domain.elemBC(j4, domain.elemBC(j4) | LuleshConfig::ETA_P_FREE);
            int j5 = rowInc + j;
            domain.elemBC(j5, domain.elemBC(j5) | LuleshConfig::ZETA_M_SYMM);
            int j6 = rowInc + j + domElems - edgeElems * edgeElems;
            domain.elemBC(j6, domain.elemBC(j6) | LuleshConfig::ZETA_P_FREE);
        }
    }
    /* Initialize arrays*/
    
    double *sigxx = hclib::numa_malloc<double>(domElems);
    double *sigyy = hclib::numa_malloc<double>(domElems);
    double *sigzz =  hclib::numa_malloc<double>(domElems);
    double *determinant = hclib::numa_malloc<double>(domElems);
    double *dvdx = hclib::numa_malloc<double>(domElems*8);
    double *dvdy = hclib::numa_malloc<double>(domElems*8);
    double *dvdz = hclib::numa_malloc<double>(domElems*8);
    double *x8n = hclib::numa_malloc<double>(domElems*8);
    double *y8n = hclib::numa_malloc<double>(domElems*8);
    double *z8n = hclib::numa_malloc<double>(domElems*8);
    double *fxElem = hclib::numa_malloc<double>(domElems*8);
    double *fyElem = hclib::numa_malloc<double>(domElems*8);
    double *fzElem = hclib::numa_malloc<double>(domElems*8);
    
    initData_evalEOSForElems(domain.numElem());

    /* timestep to solution */

    auto startTime = chrono::steady_clock::now();

    hclib::kernel([&]() {
    while (domain.time() < domain.stoptime()) {
        timeIncrement(domain);
        lagrangeLeapFrog(domain, sigxx, sigyy, sigzz, determinant,
                        dvdx, dvdy, dvdz,
                        x8n, y8n, z8n,
                        fxElem, fyElem, fzElem);
    }
    });

    auto endTime = chrono::steady_clock::now();
    std::chrono::duration<double> elapsed_seconds = endTime - startTime;
    cout << "\n\n  Elapsed time = " << elapsed_seconds.count()*1000 <<  " ms \n\n";

    int elemId = 0;
    printf("  Run completed:  \n");
    printf("    Problem size        =  %d \n", edgeElems);
    printf("    Iteration count     =  %d \n", domain.cycle());
    printf("    Origin Energy =  %12.6f \n", domain.e(elemId));

    double maxAbsDiff = 0.0;
    double totalAbsDiff = 0.0;
    double maxRelDiff = 0.0;


    for (int j = 0; j < edgeElems; ++j) {
        for (int k = j + 1; k < edgeElems; ++k) {
            double absDiff = abs(domain.e(j * edgeElems + k) - domain.e(k * edgeElems + j));
            totalAbsDiff += absDiff;

            if (maxAbsDiff < absDiff) {
                maxAbsDiff = absDiff;
            }

            double relDiff = absDiff / domain.e(k * edgeElems + j);

            if (maxRelDiff < relDiff) {
                maxRelDiff = relDiff;
            }
        }
    }

    printf("    Testing Plane 0 of Energy Array:\n");
    printf("      Max Abs Diff   = %12.6e \n", maxAbsDiff);
    printf("      Total Abs Diff = %12.6e \n", totalAbsDiff);
    printf("      Max Rel Diff   = %12.6e \n\n", maxRelDiff);


    /*Delete arrays*/
    {
        hclib::numa_free(sigxx);
        hclib::numa_free(sigyy);
        hclib::numa_free(sigzz);
        hclib::numa_free(determinant);
        hclib::numa_free(dvdx);
        hclib::numa_free(dvdy);
        hclib::numa_free(dvdz);
        hclib::numa_free(x8n);
        hclib::numa_free(y8n);
        hclib::numa_free(z8n);
        hclib::numa_free(fxElem);
        hclib::numa_free(fyElem);
        hclib::numa_free(fzElem);
    }
    freeData();
} 


void LuleshParBenchmark::timeIncrement(Domain& domain){
    double targetDt = domain.stoptime() - domain.time();

    if ((domain.dtfixed() <= 0.0 && (domain.cycle() != 0))) {
        double oldDt = domain.deltatime();
    
        /* This will require a reduction in parallel */
        double newDt = (1.0e+20);
        if (domain.dtcourant() < newDt) {
            newDt = domain.dtcourant() / (2.0);
        }
        if (domain.dthydro() < newDt) {
            newDt = domain.dthydro() * (2.0) / (3.0);
        }
    
        double ratio = newDt / oldDt;
        if (ratio >= (1.0)) {
            if (ratio < domain.deltatimemultlb()) {
                newDt = oldDt;
            } else if (ratio > domain.deltatimemultub()) {
                newDt = oldDt * domain.deltatimemultub();
            }
        }

        if (newDt > domain.dtmax()) {
            newDt = domain.dtmax();
        }
        domain.deltatime(newDt);
    }

    /* TRY TO PREVENT VERY SMALL SCALING ON THE NEXT CYCLE */
    if ((targetDt > domain.deltatime()) && (targetDt < 4.0 * domain.deltatime() / 3.0)) {
        targetDt = 2.0 * domain.deltatime() / 3.0;
    }

    if (targetDt < domain.deltatime()) {
        domain.deltatime(targetDt);
    }

    domain.time(domain.time() + domain.deltatime());

    domain.cycle(domain.cycle() + 1);
}

void LuleshParBenchmark::lagrangeLeapFrog(Domain& domain, 
            double sigxx[], double sigyy[], double sigzz[], double determinant[],
            double dvdx[], double dvdy[], double dvdz[],
            double x8n[], double y8n[], double z8n[],
            double fxElem[], double fyElem[], double fzElem[]){
    /* calculate nodal forces, accelerations, velocities, positions, with
        * applied boundary conditions and slide surface considerations */
        lagrangeNodal(domain, sigxx, sigyy, sigzz, determinant,
                        dvdx, dvdy, dvdz,
                        x8n, y8n, z8n,
                        fxElem, fyElem, fzElem);

	    /* calculate element quantities (i.e. velocity gradient & q), and update
        * material states */
        lagrangeElements(domain);

        calcTimeConstraintsForElems(domain);

        // lagrangeRelease() ;  Creation/destruction of temps may be important to capture
}

void LuleshParBenchmark::lagrangeNodal(Domain& domain,
                            double sigxx[], double sigyy[], double sigzz[], double determinant[],
                            double dvdx[], double dvdy[], double dvdz[],
                            double x8n[], double y8n[], double z8n[],
                            double fxElem[], double fyElem[], double fzElem[]) {
    double delta = domain.deltatime();
    double uCut = domain.u_cut();

    /* time of boundary condition evaluation is beginning of step for force and acceleration boundary conditions. */
    calcForceForNodes(domain, sigxx, sigyy, sigzz, determinant,
                    dvdx, dvdy, dvdz,
                    x8n, y8n, z8n,
                    fxElem, fyElem, fzElem);

    calcAccelerationForNodes(domain);

    applyAccelerationBoundaryConditionsForNodes(domain);

    calcVelocityForNodes(domain, delta, uCut);

    calcPositionForNodes(domain, delta);

    return;
}

void LuleshParBenchmark::calcForceForNodes(Domain& domain,
                    double sigxx[], double sigyy[], double sigzz[], double determinant[], 
                    double dvdx[], double dvdy[], double dvdz[],
                    double x8n[], double y8n[], double z8n[],
                    double fxElem[], double fyElem[], double fzElem[]){
    domain.fxFill(0.0, domain.numNode());
    domain.fyFill(0.0, domain.numNode());
    domain.fzFill(0.0, domain.numNode());

    /* calcforce calls partial, force, hourq */
    calcVolumeForceForElems(domain, sigxx, sigyy, sigzz, determinant,
                            dvdx, dvdy, dvdz,
                            x8n, y8n, z8n,
                            fxElem, fyElem, fzElem);

    /* calculate Nodal Forces at domain boundaries */
    /* problem->commSBN->Transfer(CommSBN::forces); */
}

void LuleshParBenchmark::calcVolumeForceForElems(Domain& domain, 
                                            double sigxx[], double sigyy[], double sigzz[], double determinant[],
                                            double dvdx[], double dvdy[], double dvdz[],
                                            double x8n[], double y8n[], double z8n[],
                                            double fxElem[], double fyElem[], double fzElem[]){
    int numElem = domain.numElem();
        if (numElem != 0) {
            double hgcoef = domain.hgcoef();
            
		    /* sum contributions to total stress tensor */
            initStressTermsForElems(domain, numElem, sigxx, sigyy, sigzz);
            
            // call elemlib stress integration loop to produce nodal forces from
            // material stresses.
            integrateStressForElems(domain, numElem, sigxx, sigyy, sigzz, determinant,
                                fxElem, fyElem, fzElem);

            //check for negative element volume

            /**
             Time taken(in microseconds) for sequential execution:  86.04225

            Problem Size = 32, Number of Threads = 20
            Time taken(in microseconds)for the execution of this loop for different values of 'Factor'

                Factor    Time(microsec)
                1          34.14255
                2          31.4918
                4          39.7053
                8          71.7314
                16         128.646
                32         271.7345
            **/
	    auto start_i = [=](int i) { return i; };
	    auto end_i = [=](int i) { return i-1; };
            HCLIB_FINISH {
                int Factor = 1;
                hclib::loop_domain_t loop = {0, numElem,1, int(numElem/(CONST_FACTOR*Factor))};
                    hclib::irregular_recursion1D(&loop, determinant, start_i, end_i, [=, &domain](int k){
                        if (determinant[k] <= 0.0) {
                            throw "Negative volume (" + to_string(determinant[k]) + ") at position " + to_string(k);
                        }
                    }, FORASYNC_MODE_RECURSIVE);
            }         
            calcHourglassControlForElems(domain, determinant, hgcoef,
                                        dvdx, dvdy, dvdz, x8n, y8n, z8n,
                                        fxElem, fyElem, fzElem);
        }
}

void LuleshParBenchmark::initStressTermsForElems(Domain& domain, int numElem, double sigxx[], double sigyy[], double sigzz[]){
    // pull in the stresses appropriate to the hydro integration
    
    /**
     Time taken(in microseconds) for sequential execution:  412.7545

    Problem Size = 32, Number of Threads = 20
    Time taken(in microseconds)for the execution of this loop for different values of 'Factor'

        Factor    Time(microsec)
        1          140.42495
        2          68.0008
        4          72.50225
        8          89.9324
        16         145.6685
        32         314.5485
    **/
    auto start_i = [=](int i) { return i; };
    auto end_i = [=](int i) { return i-1; };
    HCLIB_FINISH {
        int Factor = 1;
        hclib::loop_domain_t loop = {0, numElem,1,  int(numElem/(CONST_FACTOR*Factor))};
            hclib::irregular_recursion1D(&loop, sigxx, start_i, end_i, [=, &domain](int i){
                    double temp = -domain.p(i) - domain.q(i);
                    sigxx[i] = temp;
                    sigyy[i] = temp;
                    sigzz[i] = temp;
            }, FORASYNC_MODE_RECURSIVE);
    }
}

void LuleshParBenchmark::integrateStressForElems(
            Domain& domain, int numElem,
            double sigxx[], double sigyy[], double sigzz[],
            double determinant[],
            double fxElem[], double fyElem[], double fzElem[]){

    int numElem8 = numElem * 8;
    
    /**
     Time taken(in microseconds) for sequential execution:  21037.7

    Problem Size = 32, Number of Threads = 20
    Time taken(in microseconds)for the execution of this loop for different values of 'Factor'

        Factor    Time(microsec)
        1          1582.745
        2          1648.32
        4          1312.74
        8          1552.165
        16         1281.55
        32         1364.095
    **/
    auto start_i = [=](int i) { return i; };
    auto end_i = [=](int i) { return i-1; };
    HCLIB_FINISH {
        int Factor = 17;
        hclib::loop_domain_t loop = {0, numElem,1,  int(numElem/(CONST_FACTOR*Factor))};
            hclib::irregular_recursion1D(&loop, sigxx, start_i, end_i, [=, &domain](int k){
                double B[3*8];
                double xLocal[8];
                double yLocal[8];
                double zLocal[8];

                int *elemNodes = domain.nodelist();

                // get nodal coordinates from global arrays and copy into local arrays.
                for (int lnode = 0; lnode < 8; ++lnode) {
                    int gnode = elemNodes[k*8+lnode];
                    xLocal[lnode] = domain.x(gnode);
                    yLocal[lnode] = domain.y(gnode);
                    zLocal[lnode] = domain.z(gnode);
                }

                /* Volume calculation involves extra work for numerical consistency. */
                calcElemShapeFunctionDerivatives(xLocal, yLocal, zLocal, B, determinant, k);

                calcElemNodeNormals(&B[8*0], &B[8*1], &B[8*2], xLocal, yLocal, zLocal);

                sumElemStressesToNodeForces(B, sigxx[k], sigyy[k], sigzz[k],
                                    fxElem, k * 8, fyElem, k * 8, fzElem, k * 8);

                
    } , FORASYNC_MODE_RECURSIVE);
    }

    {
        int numNode = domain.numNode();

        /**
         Time taken(in microseconds) for sequential execution:  3774.685

        Problem Size = 32, Number of Threads = 20
        Time taken(in microseconds)for the execution of this loop for different values of 'Factor'

            Factor    Time(microsec)
            1          374.665
            2          371.7625
            4          353.9555
            8          354.8845
            16         378.487
            32         404.9975
        **/
    	auto start_i = [=](int i) { return i; };
    	auto end_i = [=](int i) { return i-1; };
        HCLIB_FINISH {
            int Factor = 2;
            hclib::loop_domain_t loop = {0, numNode,1,  int(numNode/(CONST_FACTOR*Factor))};
            hclib::irregular_recursion1D(& loop, domain.get_m_fx(), start_i, end_i, [=, &domain](int gnode){
                int count = domain.nodeElemCount(gnode);
                int start = domain.nodeElemStart(gnode);
                double fx = 0.0;
                double fy = 0.0;
                double fz = 0.0;
                for (int i = 0; i < count; ++i) {
                    int elem = domain.nodeElemCornerList(start + i);
                    fx += fxElem[elem];
                    fy += fyElem[elem];
                    fz += fzElem[elem];
                }
                domain.fx(gnode, fx);
                domain.fy(gnode, fy);
                domain.fz(gnode, fz);
            }, FORASYNC_MODE_RECURSIVE);
        }
    }
}

void LuleshParBenchmark::calcElemNodeNormals(
        double pfx[], double pfy[], double pfz[],
        double x[], double y[], double z[]){
    
    for (int i = 0; i < 8; ++i) {
            pfx[i] = 0.0;
            pfy[i] = 0.0;
            pfz[i] = 0.0;
        }
    /* evaluate face one: nodes 0, 1, 2, 3 */
    sumElemFaceNormal(
            pfx, 0, pfy, 0, pfz, 0,
            pfx, 1, pfy, 1, pfz, 1,
            pfx, 2, pfy, 2, pfz, 2,
            pfx, 3, pfy, 3, pfz, 3,
            x[0], y[0], z[0], x[1], y[1], z[1],
            x[2], y[2], z[2], x[3], y[3], z[3]);
    /* evaluate face two: nodes 0, 4, 5, 1 */
    sumElemFaceNormal(
            pfx, 0, pfy, 0, pfz, 0,
            pfx, 4, pfy, 4, pfz, 4,
            pfx, 5, pfy, 5, pfz, 5,
            pfx, 1, pfy, 1, pfz, 1,
            x[0], y[0], z[0], x[4], y[4], z[4],
            x[5], y[5], z[5], x[1], y[1], z[1]);
    /* evaluate face three: nodes 1, 5, 6, 2 */
    sumElemFaceNormal(
            pfx, 1, pfy, 1, pfz, 1,
            pfx, 5, pfy, 5, pfz, 5,
            pfx, 6, pfy, 6, pfz, 6,
            pfx, 2, pfy, 2, pfz, 2,
            x[1], y[1], z[1], x[5], y[5], z[5],
            x[6], y[6], z[6], x[2], y[2], z[2]);
    /* evaluate face four: nodes 2, 6, 7, 3 */
    sumElemFaceNormal(
            pfx, 2, pfy, 2, pfz, 2,
            pfx, 6, pfy, 6, pfz, 6,
            pfx, 7, pfy, 7, pfz, 7,
            pfx, 3, pfy, 3, pfz, 3,
            x[2], y[2], z[2], x[6], y[6], z[6],
            x[7], y[7], z[7], x[3], y[3], z[3]);
    /* evaluate face five: nodes 3, 7, 4, 0 */
    sumElemFaceNormal(
            pfx, 3, pfy, 3, pfz, 3,
            pfx, 7, pfy, 7, pfz, 7,
            pfx, 4, pfy, 4, pfz, 4,
            pfx, 0, pfy, 0, pfz, 0,
            x[3], y[3], z[3], x[7], y[7], z[7],
            x[4], y[4], z[4], x[0], y[0], z[0]);
    /* evaluate face six: nodes 4, 7, 6, 5 */
    sumElemFaceNormal(
            pfx, 4, pfy, 4, pfz, 4,
            pfx, 7, pfy, 7, pfz, 7,
            pfx, 6, pfy, 6, pfz, 6,
            pfx, 5, pfy, 5, pfz, 5,
            x[4], y[4], z[4], x[7], y[7], z[7],
            x[6], y[6], z[6], x[5], y[5], z[5]);

}

void LuleshParBenchmark::sumElemFaceNormal(
        double nX0[], int nX0index, double nY0[], int nY0index, double nZ0[], int nZ0index,
        double nX1[], int nX1index, double nY1[], int nY1index, double nZ1[], int nZ1index,
        double nX2[], int nX2index, double nY2[], int nY2index, double nZ2[], int nZ2index,
        double nX3[], int nX3index, double nY3[], int nY3index, double nZ3[], int nZ3index,
        double x0, double y0, double z0,
        double x1, double y1, double z1,
        double x2, double y2, double z2,
        double x3, double y3, double z3){

    double bisectX0 = (0.50) * (x3 + x2 - x1 - x0);
    double bisectY0 = (0.50) * (y3 + y2 - y1 - y0);
    double bisectZ0 = (0.50) * (z3 + z2 - z1 - z0);
    double bisectX1 = (0.50) * (x2 + x1 - x3 - x0);
    double bisectY1 = (0.50) * (y2 + y1 - y3 - y0);
    double bisectZ1 = (0.50) * (z2 + z1 - z3 - z0);
    double areaX = (0.25) * (bisectY0 * bisectZ1 - bisectZ0 * bisectY1);
    double areaY = (0.25) * (bisectZ0 * bisectX1 - bisectX0 * bisectZ1);
    double areaZ = (0.25) * (bisectX0 * bisectY1 - bisectY0 * bisectX1);

    nX0[nX0index] += areaX;
    nX1[nX1index] += areaX;
    nX2[nX2index] += areaX;
    nX3[nX3index] += areaX;

    nY0[nY0index] += areaY;
    nY1[nY1index] += areaY;
    nY2[nY2index] += areaY;
    nY3[nY3index] += areaY;

    nZ0[nZ0index] += areaZ;
    nZ1[nZ1index] += areaZ;
    nZ2[nZ2index] += areaZ;
    nZ3[nZ3index] += areaZ;    

}

void LuleshParBenchmark::sumElemStressesToNodeForces(
        double B[],
        double xxStress, double yyStress, double zzStress,
        double fx[], int xIndex,
        double fy[], int yIndex,
        double fz[], int zIndex) {

    double pfx0 = B[0*8+0];
    double pfx1 = B[0*8+1];
    double pfx2 = B[0*8+2];
    double pfx3 = B[0*8+3];
    double pfx4 = B[0*8+4];
    double pfx5 = B[0*8+5];
    double pfx6 = B[0*8+6];
    double pfx7 = B[0*8+7];

    double pfy0 = B[1*8+0];
    double pfy1 = B[1*8+1];
    double pfy2 = B[1*8+2];
    double pfy3 = B[1*8+3];
    double pfy4 = B[1*8+4];
    double pfy5 = B[1*8+5];
    double pfy6 = B[1*8+6];
    double pfy7 = B[1*8+7];

    double pfz0 = B[2*8+0];
    double pfz1 = B[2*8+1];
    double pfz2 = B[2*8+2];
    double pfz3 = B[2*8+3];
    double pfz4 = B[2*8+4];
    double pfz5 = B[2*8+5];
    double pfz6 = B[2*8+6];
    double pfz7 = B[2*8+7];

    fx[xIndex + 0] = -(xxStress * pfx0);
    fx[xIndex + 1] = -(xxStress * pfx1);
    fx[xIndex + 2] = -(xxStress * pfx2);
    fx[xIndex + 3] = -(xxStress * pfx3);
    fx[xIndex + 4] = -(xxStress * pfx4);
    fx[xIndex + 5] = -(xxStress * pfx5);
    fx[xIndex + 6] = -(xxStress * pfx6);
    fx[xIndex + 7] = -(xxStress * pfx7);

    fy[yIndex + 0] = -(yyStress * pfy0);
    fy[yIndex + 1] = -(yyStress * pfy1);
    fy[yIndex + 2] = -(yyStress * pfy2);
    fy[yIndex + 3] = -(yyStress * pfy3);
    fy[yIndex + 4] = -(yyStress * pfy4);
    fy[yIndex + 5] = -(yyStress * pfy5);
    fy[yIndex + 6] = -(yyStress * pfy6);
    fy[yIndex + 7] = -(yyStress * pfy7);

    fz[zIndex + 0] = -(zzStress * pfz0);
    fz[zIndex + 1] = -(zzStress * pfz1);
    fz[zIndex + 2] = -(zzStress * pfz2);
    fz[zIndex + 3] = -(zzStress * pfz3);
    fz[zIndex + 4] = -(zzStress * pfz4);
    fz[zIndex + 5] = -(zzStress * pfz5);
    fz[zIndex + 6] = -(zzStress * pfz6);
    fz[zIndex + 7] = -(zzStress * pfz7);

}

void LuleshParBenchmark::calcHourglassControlForElems(
        Domain& domain, double determinant[], double hgcoef,
        double dvdx[], double dvdy[], double dvdz[],
        double x8n[], double y8n[], double z8n[],
        double fxElem[], double fyElem[], double fzElem[]){

    int numElem = domain.numElem();
    int numElem8 = numElem * 8;
    
    /* start loop over elements */

    /**
     Time taken(in microseconds) for sequential execution:  19963.95

     Problem Size = 32, Number of Threads = 20
     Time taken(in microseconds)for the execution of this loop for different values of 'Factor'

        Factor    Time(microsec)
        1          1579.13
        2          1545.75
        4          1289.555
        8          1361.12
        16         1289.975
        32         1308.21
    **/
    auto start_i = [=](int i) { return i; };
    auto end_i = [=](int i) { return i-1; };
    HCLIB_FINISH {
        int Factor = 9;
        hclib::loop_domain_t loop = {0, numElem,1,  int(numElem/(CONST_FACTOR*Factor))};
        hclib::irregular_recursion1D(&loop, determinant, start_i, end_i, [&](int i){
            
            double x1[8];
            double y1[8];
            double z1[8];

            double pfx[8];
            double pfy[8];
            double pfz[8];

            int *nodeList = domain.nodelist();
            int *elemToNode = &nodeList[i*8];
            collectDomainNodesToElemNodes(domain, elemToNode, x1, y1, z1);

            calcElemVolumeDerivative(pfx, pfy, pfz, x1, y1, z1);

            /* load into temporary storage for FB Hour Glass control */
            for (int ii = 0; ii < 8; ++ii) {
                int jj = 8 * i + ii;

                dvdx[jj] = pfx[ii];
                dvdy[jj] = pfy[ii];
                dvdz[jj] = pfz[ii];
                x8n[jj] = x1[ii];
                y8n[jj] = y1[ii];
                z8n[jj] = z1[ii];
            }

            determinant[i] = domain.volo(i) * domain.v(i);

            /* Do a check for negative volumes */
            if (domain.v(i) <= 0.0) {
                throw "Negative volume (" + to_string(domain.v(i)) + ") at position " + to_string(i);
            }
        }, FORASYNC_MODE_RECURSIVE);
    }

    if (hgcoef > 0.0) {
        calcFBHourglassForceForElems(domain, determinant, x8n, y8n, z8n, dvdx, dvdy, dvdz, hgcoef,
                                fxElem, fyElem, fzElem);
    }
}
        
void LuleshParBenchmark::collectDomainNodesToElemNodes(
        Domain& domain, int elemToNode[],
        double elemX[], double elemY[], double elemZ[]) {
    
    int nd0i = elemToNode[0];
    int nd1i = elemToNode[1];
    int nd2i = elemToNode[2];
    int nd3i = elemToNode[3];
    int nd4i = elemToNode[4];
    int nd5i = elemToNode[5];
    int nd6i = elemToNode[6];
    int nd7i = elemToNode[7];

    elemX[0] = domain.x(nd0i);
    elemX[1] = domain.x(nd1i);
    elemX[2] = domain.x(nd2i);
    elemX[3] = domain.x(nd3i);
    elemX[4] = domain.x(nd4i);
    elemX[5] = domain.x(nd5i);
    elemX[6] = domain.x(nd6i);
    elemX[7] = domain.x(nd7i);

    elemY[0] = domain.y(nd0i);
    elemY[1] = domain.y(nd1i);
    elemY[2] = domain.y(nd2i);
    elemY[3] = domain.y(nd3i);
    elemY[4] = domain.y(nd4i);
    elemY[5] = domain.y(nd5i);
    elemY[6] = domain.y(nd6i);
    elemY[7] = domain.y(nd7i);

    elemZ[0] = domain.z(nd0i);
    elemZ[1] = domain.z(nd1i);
    elemZ[2] = domain.z(nd2i);
    elemZ[3] = domain.z(nd3i);
    elemZ[4] = domain.z(nd4i);
    elemZ[5] = domain.z(nd5i);
    elemZ[6] = domain.z(nd6i);
    elemZ[7] = domain.z(nd7i);       

}

void LuleshParBenchmark::calcElemVolumeDerivative(
        double dvdx[], double dvdy[], double dvdz[],
        double x[], double y[], double z[]){
    
    voluDer(x[1], x[2], x[3], x[4], x[5], x[7],
            y[1], y[2], y[3], y[4], y[5], y[7],
            z[1], z[2], z[3], z[4], z[5], z[7],
            dvdx, 0, dvdy, 0, dvdz, 0);
    voluDer(x[0], x[1], x[2], x[7], x[4], x[6],
            y[0], y[1], y[2], y[7], y[4], y[6],
            z[0], z[1], z[2], z[7], z[4], z[6],
            dvdx, 3, dvdy, 3, dvdz, 3);
    voluDer(x[3], x[0], x[1], x[6], x[7], x[5],
            y[3], y[0], y[1], y[6], y[7], y[5],
            z[3], z[0], z[1], z[6], z[7], z[5],
            dvdx, 2, dvdy, 2, dvdz, 2);
    voluDer(x[2], x[3], x[0], x[5], x[6], x[4],
            y[2], y[3], y[0], y[5], y[6], y[4],
            z[2], z[3], z[0], z[5], z[6], z[4],
            dvdx, 1, dvdy, 1, dvdz, 1);
    voluDer(x[7], x[6], x[5], x[0], x[3], x[1],
            y[7], y[6], y[5], y[0], y[3], y[1],
            z[7], z[6], z[5], z[0], z[3], z[1],
            dvdx, 4, dvdy, 4, dvdz, 4);
    voluDer(x[4], x[7], x[6], x[1], x[0], x[2],
            y[4], y[7], y[6], y[1], y[0], y[2],
            z[4], z[7], z[6], z[1], z[0], z[2],
            dvdx, 5, dvdy, 5, dvdz, 5);
    voluDer(x[5], x[4], x[7], x[2], x[1], x[3],
            y[5], y[4], y[7], y[2], y[1], y[3],
            z[5], z[4], z[7], z[2], z[1], z[3],
            dvdx, 6, dvdy, 6, dvdz, 6);
    voluDer(x[6], x[5], x[4], x[3], x[2], x[0],
            y[6], y[5], y[4], y[3], y[2], y[0],
            z[6], z[5], z[4], z[3], z[2], z[0],
            dvdx, 7, dvdy, 7, dvdz, 7);       
    
}


void LuleshParBenchmark::voluDer(
        double x0, double x1, double x2,
        double x3, double x4, double x5,
        double y0, double y1, double y2,
        double y3, double y4, double y5,
        double z0, double z1, double z2,
        double z3, double z4, double z5,
        double dvdx[], int dvdxIndex,
        double dvdy[], int dvdyIndex,
        double dvdz[], int dvdzIndex){
    
    double twelfth = (1.0) / (12.0);

    dvdx[dvdxIndex] =
            (y1 + y2) * (z0 + z1) - (y0 + y1) * (z1 + z2) +
                    (y0 + y4) * (z3 + z4) - (y3 + y4) * (z0 + z4) -
                    (y2 + y5) * (z3 + z5) + (y3 + y5) * (z2 + z5);
    dvdy[dvdyIndex] =
            -(x1 + x2) * (z0 + z1) + (x0 + x1) * (z1 + z2) -
                    (x0 + x4) * (z3 + z4) + (x3 + x4) * (z0 + z4) +
                    (x2 + x5) * (z3 + z5) - (x3 + x5) * (z2 + z5);

    dvdz[dvdzIndex] =
            -(y1 + y2) * (x0 + x1) + (y0 + y1) * (x1 + x2) -
                    (y0 + y4) * (x3 + x4) + (y3 + y4) * (x0 + x4) +
                    (y2 + y5) * (x3 + x5) - (y3 + y5) * (x2 + x5);

    dvdx[dvdxIndex] *= twelfth;
    dvdy[dvdyIndex] *= twelfth;
    dvdz[dvdzIndex] *= twelfth;
    
}


void LuleshParBenchmark::calcFBHourglassForceForElems(
        Domain& domain,
        double determinant[],
        double x8n[], double y8n[], double z8n[],
        double dvdx[], double dvdy[], double dvdz[],
        double hourglass,
        double fxElem[], double fyElem[], double fzElem[]){

    /***************************************************************************
     *    FUNCTION: calculates the Flanagan-Belytschko anti-hourglass force.   *
     ***************************************************************************/

    int numElem = domain.numElem();
    int numElem8 = numElem * 8;
    
    double gamma[4*8];

    gamma[0*8+0] = 1.0;
    gamma[0*8+1] = 1.0;
    gamma[0*8+2] = -1.0;
    gamma[0*8+3] = -1.0;
    gamma[0*8+4] = -1.0;
    gamma[0*8+5] = -1.0;
    gamma[0*8+6] = 1.0;
    gamma[0*8+7] = 1.0;
    gamma[1*8+0] = 1.0;
    gamma[1*8+1] = -1.0;
    gamma[1*8+2] = -1.0;
    gamma[1*8+3] = 1.0;
    gamma[1*8+4] = -1.0;
    gamma[1*8+5] = 1.0;
    gamma[1*8+6] = 1.0;
    gamma[1*8+7] = -1.0;
    gamma[2*8+0] = 1.0;
    gamma[2*8+1] = -1.0;
    gamma[2*8+2] = 1.0;
    gamma[2*8+3] = -1.0;
    gamma[2*8+4] = 1.0;
    gamma[2*8+5] = -1.0;
    gamma[2*8+6] = 1.0;
    gamma[2*8+7] = -1.0;
    gamma[3*8+0] = -1.0;
    gamma[3*8+1] = 1.0;
    gamma[3*8+2] = -1.0;
    gamma[3*8+3] = 1.0;
    gamma[3*8+4] = 1.0;
    gamma[3*8+5] = -1.0;
    gamma[3*8+6] = 1.0;
    gamma[3*8+7] = -1.0;

    /*************************************************/
    /*    compute the hourglass modes */

    /**
     Time taken(in microseconds) for sequential execution:  38375.85

    Problem Size = 32, Number of Threads = 20
    Time taken(in microseconds)for the execution of this loop for different values of 'Factor'

        Factor    Time(microsec)
        1          2623.705
        2          2373.395
        4          2137.58
        8          2043.415
        16         2071.075
        32         2098.435
    **/
    auto start_i = [=](int i) { return i; };
    auto end_i = [=](int i) { return i-1; };
    HCLIB_FINISH {
        int Factor = 9;
        hclib::loop_domain_t loop = {0, numElem,1,  int(numElem/(CONST_FACTOR*Factor))};
        hclib::irregular_recursion1D(&loop, determinant, start_i, end_i, [=, &domain](int i2){
            double hgfx[8];
            double hgfy[8];
            double hgfz[8];

            double hourglassAm0[4];
            double hourglassAm1[4];
            double hourglassAm2[4];
            double hourglassAm3[4];
            double hourglassAm4[4];
            double hourglassAm5[4];
            double hourglassAm6[4];
            double hourglassAm7[4];

            double xd1[8];
            double yd1[8];
            double zd1[8];

            int *nodeList = domain.nodelist();
            int *elemToNode = &nodeList[i2*8];
            int i3 = 8 * i2;
            double volinv = (1.0) / determinant[i2];

            for (int i1 = 0; i1 < 4; ++i1) {

                double hourmodx =
                        x8n[i3] * gamma[i1*8+0] + x8n[i3 + 1] * gamma[i1*8+1] +
                                x8n[i3 + 2] * gamma[i1*8+2] + x8n[i3 + 3] * gamma[i1*8+3] +
                                x8n[i3 + 4] * gamma[i1*8+4] + x8n[i3 + 5] * gamma[i1*8+5] +
                                x8n[i3 + 6] * gamma[i1*8+6] + x8n[i3 + 7] * gamma[i1*8+7];

                double hourmody =
                        y8n[i3] * gamma[i1*8+0] + y8n[i3 + 1] * gamma[i1*8+1] +
                                y8n[i3 + 2] * gamma[i1*8+2] + y8n[i3 + 3] * gamma[i1*8+3] +
                                y8n[i3 + 4] * gamma[i1*8+4] + y8n[i3 + 5] * gamma[i1*8+5] +
                                y8n[i3 + 6] * gamma[i1*8+6] + y8n[i3 + 7] * gamma[i1*8+7];

                double hourmodz =
                        z8n[i3] * gamma[i1*8+0] + z8n[i3 + 1] * gamma[i1*8+1] +
                                z8n[i3 + 2] * gamma[i1*8+2] + z8n[i3 + 3] * gamma[i1*8+3] +
                                z8n[i3 + 4] * gamma[i1*8+4] + z8n[i3 + 5] * gamma[i1*8+5] +
                                z8n[i3 + 6] * gamma[i1*8+6] + z8n[i3 + 7] * gamma[i1*8+7];

                hourglassAm0[i1] = gamma[i1*8+0] - volinv * (dvdx[i3] * hourmodx +
                        dvdy[i3] * hourmody +
                        dvdz[i3] * hourmodz);

                hourglassAm1[i1] = gamma[i1*8+1] - volinv * (dvdx[i3 + 1] * hourmodx +
                        dvdy[i3 + 1] * hourmody +
                        dvdz[i3 + 1] * hourmodz);

                hourglassAm2[i1] = gamma[i1*8+2] - volinv * (dvdx[i3 + 2] * hourmodx +
                        dvdy[i3 + 2] * hourmody +
                        dvdz[i3 + 2] * hourmodz);

                hourglassAm3[i1] = gamma[i1*8+3] - volinv * (dvdx[i3 + 3] * hourmodx +
                        dvdy[i3 + 3] * hourmody +
                        dvdz[i3 + 3] * hourmodz);

                hourglassAm4[i1] = gamma[i1*8+4] - volinv * (dvdx[i3 + 4] * hourmodx +
                        dvdy[i3 + 4] * hourmody +
                        dvdz[i3 + 4] * hourmodz);

                hourglassAm5[i1] = gamma[i1*8+5] - volinv * (dvdx[i3 + 5] * hourmodx +
                        dvdy[i3 + 5] * hourmody +
                        dvdz[i3 + 5] * hourmodz);

                hourglassAm6[i1] = gamma[i1*8+6] - volinv * (dvdx[i3 + 6] * hourmodx +
                        dvdy[i3 + 6] * hourmody +
                        dvdz[i3 + 6] * hourmodz);

                hourglassAm7[i1] = gamma[i1*8+7] - volinv * (dvdx[i3 + 7] * hourmodx +
                        dvdy[i3 + 7] * hourmody +
                        dvdz[i3 + 7] * hourmodz);
            }
            /* compute forces */
            /* store forces into h arrays (force arrays) */

            double ss1 = domain.ss(i2);
            double mass1 = domain.elemMass(i2);
            double volume13 = std::cbrt(determinant[i2]);

            int n0si2 = elemToNode[0];
            int n1si2 = elemToNode[1];
            int n2si2 = elemToNode[2];
            int n3si2 = elemToNode[3];
            int n4si2 = elemToNode[4];
            int n5si2 = elemToNode[5];
            int n6si2 = elemToNode[6];
            int n7si2 = elemToNode[7];

            xd1[0] = domain.xd(n0si2);
            xd1[1] = domain.xd(n1si2);
            xd1[2] = domain.xd(n2si2);
            xd1[3] = domain.xd(n3si2);
            xd1[4] = domain.xd(n4si2);
            xd1[5] = domain.xd(n5si2);
            xd1[6] = domain.xd(n6si2);
            xd1[7] = domain.xd(n7si2);

            yd1[0] = domain.yd(n0si2);
            yd1[1] = domain.yd(n1si2);
            yd1[2] = domain.yd(n2si2);
            yd1[3] = domain.yd(n3si2);
            yd1[4] = domain.yd(n4si2);
            yd1[5] = domain.yd(n5si2);
            yd1[6] = domain.yd(n6si2);
            yd1[7] = domain.yd(n7si2);

            zd1[0] = domain.zd(n0si2);
            zd1[1] = domain.zd(n1si2);
            zd1[2] = domain.zd(n2si2);
            zd1[3] = domain.zd(n3si2);
            zd1[4] = domain.zd(n4si2);
            zd1[5] = domain.zd(n5si2);
            zd1[6] = domain.zd(n6si2);
            zd1[7] = domain.zd(n7si2);

            double coefficient = -hourglass * (0.01) * ss1 * mass1 / volume13;

            calcElemFBHourglassForce(xd1, yd1, zd1,
                                        hourglassAm0, hourglassAm1, hourglassAm2, hourglassAm3,
                                        hourglassAm4, hourglassAm5, hourglassAm6, hourglassAm7,
                                        coefficient, hgfx, hgfy, hgfz);

            fxElem[i3 + 0] = hgfx[0];
            fxElem[i3 + 1] = hgfx[1];
            fxElem[i3 + 2] = hgfx[2];
            fxElem[i3 + 3] = hgfx[3];
            fxElem[i3 + 4] = hgfx[4];
            fxElem[i3 + 5] = hgfx[5];
            fxElem[i3 + 6] = hgfx[6];
            fxElem[i3 + 7] = hgfx[7];

            fyElem[i3 + 0] = hgfy[0];
            fyElem[i3 + 1] = hgfy[1];
            fyElem[i3 + 2] = hgfy[2];
            fyElem[i3 + 3] = hgfy[3];
            fyElem[i3 + 4] = hgfy[4];
            fyElem[i3 + 5] = hgfy[5];
            fyElem[i3 + 6] = hgfy[6];
            fyElem[i3 + 7] = hgfy[7];

            fzElem[i3 + 0] = hgfz[0];
            fzElem[i3 + 1] = hgfz[1];
            fzElem[i3 + 2] = hgfz[2];
            fzElem[i3 + 3] = hgfz[3];
            fzElem[i3 + 4] = hgfz[4];
            fzElem[i3 + 5] = hgfz[5];
            fzElem[i3 + 6] = hgfz[6];
            fzElem[i3 + 7] = hgfz[7];
        }, FORASYNC_MODE_RECURSIVE);
    }

    {
        int numNode = domain.numNode();

        /**
         Time taken(in microseconds) for sequential execution:  4019.445

        Problem Size = 32, Number of Threads = 20
        Time taken(in microseconds)for the execution of this loop for different values of 'Factor'

            Factor    Time(microsec)
            1          401.574
            2          388.5215
            4          431.3875
            8          546.3835
            16         414.3725
            32         441.4775
        **/
    	auto start_i = [=](int i) { return i; };
    	auto end_i = [=](int i) { return i-1; };
        HCLIB_FINISH {
            int Factor = 2;
            hclib::loop_domain_t loop = {0, numNode,1,  int(numNode/(CONST_FACTOR*Factor))};
            hclib::irregular_recursion1D(& loop, domain.get_m_fx(), start_i, end_i, [=, &domain](int gnode){
                    int count = domain.nodeElemCount(gnode);
                int start = domain.nodeElemStart(gnode);
                double fx = 0.0;
                double fy = 0.0;
                double fz = 0.0;
                for (int i = 0; i < count; ++i) {
                    int elem = domain.nodeElemCornerList(start + i);
                    fx += fxElem[elem];
                    fy += fyElem[elem];
                    fz += fzElem[elem];
                }
                domain.fx(gnode, domain.fx(gnode) + fx);
                domain.fy(gnode, domain.fy(gnode) + fy);
                domain.fz(gnode, domain.fz(gnode) + fz);
            }, FORASYNC_MODE_RECURSIVE);
        }
    }
}

void LuleshParBenchmark::calcElemFBHourglassForce(
        double xd[], double yd[], double zd[], double hourglassAm0[],
        double hourglassAm1[], double hourglassAm2[], double hourglassAm3[],
        double hourglassAm4[], double hourglassAm5[], double hourglassAm6[],
        double hourglassAm7[], double coefficient,
        double hgfx[], double hgfy[], double hgfz[]){

    int i00 = 0;
    int i01 = 1;
    int i02 = 2;
    int i03 = 3;

    double h00 = hourglassAm0[i00] * xd[0] + hourglassAm1[i00] * xd[1] +
            hourglassAm2[i00] * xd[2] + hourglassAm3[i00] * xd[3] +
            hourglassAm4[i00] * xd[4] + hourglassAm5[i00] * xd[5] +
            hourglassAm6[i00] * xd[6] + hourglassAm7[i00] * xd[7];

    double h01 = hourglassAm0[i01] * xd[0] + hourglassAm1[i01] * xd[1] +
            hourglassAm2[i01] * xd[2] + hourglassAm3[i01] * xd[3] +
            hourglassAm4[i01] * xd[4] + hourglassAm5[i01] * xd[5] +
            hourglassAm6[i01] * xd[6] + hourglassAm7[i01] * xd[7];

    double h02 = hourglassAm0[i02] * xd[0] + hourglassAm1[i02] * xd[1] +
            hourglassAm2[i02] * xd[2] + hourglassAm3[i02] * xd[3] +
            hourglassAm4[i02] * xd[4] + hourglassAm5[i02] * xd[5] +
            hourglassAm6[i02] * xd[6] + hourglassAm7[i02] * xd[7];

    double h03 = hourglassAm0[i03] * xd[0] + hourglassAm1[i03] * xd[1] +
            hourglassAm2[i03] * xd[2] + hourglassAm3[i03] * xd[3] +
            hourglassAm4[i03] * xd[4] + hourglassAm5[i03] * xd[5] +
            hourglassAm6[i03] * xd[6] + hourglassAm7[i03] * xd[7];

    hgfx[0] = coefficient *
            (hourglassAm0[i00] * h00 + hourglassAm0[i01] * h01 +
                    hourglassAm0[i02] * h02 + hourglassAm0[i03] * h03);

    hgfx[1] = coefficient *
            (hourglassAm1[i00] * h00 + hourglassAm1[i01] * h01 +
                    hourglassAm1[i02] * h02 + hourglassAm1[i03] * h03);

    hgfx[2] = coefficient *
            (hourglassAm2[i00] * h00 + hourglassAm2[i01] * h01 +
                    hourglassAm2[i02] * h02 + hourglassAm2[i03] * h03);

    hgfx[3] = coefficient *
            (hourglassAm3[i00] * h00 + hourglassAm3[i01] * h01 +
                    hourglassAm3[i02] * h02 + hourglassAm3[i03] * h03);

    hgfx[4] = coefficient *
            (hourglassAm4[i00] * h00 + hourglassAm4[i01] * h01 +
                    hourglassAm4[i02] * h02 + hourglassAm4[i03] * h03);

    hgfx[5] = coefficient *
            (hourglassAm5[i00] * h00 + hourglassAm5[i01] * h01 +
                    hourglassAm5[i02] * h02 + hourglassAm5[i03] * h03);

    hgfx[6] = coefficient *
            (hourglassAm6[i00] * h00 + hourglassAm6[i01] * h01 +
                    hourglassAm6[i02] * h02 + hourglassAm6[i03] * h03);

    hgfx[7] = coefficient *
            (hourglassAm7[i00] * h00 + hourglassAm7[i01] * h01 +
                    hourglassAm7[i02] * h02 + hourglassAm7[i03] * h03);

    h00 = hourglassAm0[i00] * yd[0] + hourglassAm1[i00] * yd[1] +
            hourglassAm2[i00] * yd[2] + hourglassAm3[i00] * yd[3] +
            hourglassAm4[i00] * yd[4] + hourglassAm5[i00] * yd[5] +
            hourglassAm6[i00] * yd[6] + hourglassAm7[i00] * yd[7];

    h01 = hourglassAm0[i01] * yd[0] + hourglassAm1[i01] * yd[1] +
            hourglassAm2[i01] * yd[2] + hourglassAm3[i01] * yd[3] +
            hourglassAm4[i01] * yd[4] + hourglassAm5[i01] * yd[5] +
            hourglassAm6[i01] * yd[6] + hourglassAm7[i01] * yd[7];

    h02 = hourglassAm0[i02] * yd[0] + hourglassAm1[i02] * yd[1] +
            hourglassAm2[i02] * yd[2] + hourglassAm3[i02] * yd[3] +
            hourglassAm4[i02] * yd[4] + hourglassAm5[i02] * yd[5] +
            hourglassAm6[i02] * yd[6] + hourglassAm7[i02] * yd[7];

    h03 = hourglassAm0[i03] * yd[0] + hourglassAm1[i03] * yd[1] +
            hourglassAm2[i03] * yd[2] + hourglassAm3[i03] * yd[3] +
            hourglassAm4[i03] * yd[4] + hourglassAm5[i03] * yd[5] +
            hourglassAm6[i03] * yd[6] + hourglassAm7[i03] * yd[7];

    hgfy[0] = coefficient *
            (hourglassAm0[i00] * h00 + hourglassAm0[i01] * h01 +
                    hourglassAm0[i02] * h02 + hourglassAm0[i03] * h03);

    hgfy[1] = coefficient *
            (hourglassAm1[i00] * h00 + hourglassAm1[i01] * h01 +
                    hourglassAm1[i02] * h02 + hourglassAm1[i03] * h03);

    hgfy[2] = coefficient *
            (hourglassAm2[i00] * h00 + hourglassAm2[i01] * h01 +
                    hourglassAm2[i02] * h02 + hourglassAm2[i03] * h03);

    hgfy[3] = coefficient *
            (hourglassAm3[i00] * h00 + hourglassAm3[i01] * h01 +
                    hourglassAm3[i02] * h02 + hourglassAm3[i03] * h03);

    hgfy[4] = coefficient *
            (hourglassAm4[i00] * h00 + hourglassAm4[i01] * h01 +
                    hourglassAm4[i02] * h02 + hourglassAm4[i03] * h03);

    hgfy[5] = coefficient *
            (hourglassAm5[i00] * h00 + hourglassAm5[i01] * h01 +
                    hourglassAm5[i02] * h02 + hourglassAm5[i03] * h03);

    hgfy[6] = coefficient *
            (hourglassAm6[i00] * h00 + hourglassAm6[i01] * h01 +
                    hourglassAm6[i02] * h02 + hourglassAm6[i03] * h03);

    hgfy[7] = coefficient *
            (hourglassAm7[i00] * h00 + hourglassAm7[i01] * h01 +
                    hourglassAm7[i02] * h02 + hourglassAm7[i03] * h03);

    h00 = hourglassAm0[i00] * zd[0] + hourglassAm1[i00] * zd[1] +
            hourglassAm2[i00] * zd[2] + hourglassAm3[i00] * zd[3] +
            hourglassAm4[i00] * zd[4] + hourglassAm5[i00] * zd[5] +
            hourglassAm6[i00] * zd[6] + hourglassAm7[i00] * zd[7];

    h01 = hourglassAm0[i01] * zd[0] + hourglassAm1[i01] * zd[1] +
            hourglassAm2[i01] * zd[2] + hourglassAm3[i01] * zd[3] +
            hourglassAm4[i01] * zd[4] + hourglassAm5[i01] * zd[5] +
            hourglassAm6[i01] * zd[6] + hourglassAm7[i01] * zd[7];

    h02 = hourglassAm0[i02] * zd[0] + hourglassAm1[i02] * zd[1] +
            hourglassAm2[i02] * zd[2] + hourglassAm3[i02] * zd[3] +
            hourglassAm4[i02] * zd[4] + hourglassAm5[i02] * zd[5] +
            hourglassAm6[i02] * zd[6] + hourglassAm7[i02] * zd[7];

    h03 = hourglassAm0[i03] * zd[0] + hourglassAm1[i03] * zd[1] +
            hourglassAm2[i03] * zd[2] + hourglassAm3[i03] * zd[3] +
            hourglassAm4[i03] * zd[4] + hourglassAm5[i03] * zd[5] +
            hourglassAm6[i03] * zd[6] + hourglassAm7[i03] * zd[7];

    hgfz[0] = coefficient *
            (hourglassAm0[i00] * h00 + hourglassAm0[i01] * h01 +
                    hourglassAm0[i02] * h02 + hourglassAm0[i03] * h03);

    hgfz[1] = coefficient *
            (hourglassAm1[i00] * h00 + hourglassAm1[i01] * h01 +
                    hourglassAm1[i02] * h02 + hourglassAm1[i03] * h03);

    hgfz[2] = coefficient *
            (hourglassAm2[i00] * h00 + hourglassAm2[i01] * h01 +
                    hourglassAm2[i02] * h02 + hourglassAm2[i03] * h03);

    hgfz[3] = coefficient *
            (hourglassAm3[i00] * h00 + hourglassAm3[i01] * h01 +
                    hourglassAm3[i02] * h02 + hourglassAm3[i03] * h03);

    hgfz[4] = coefficient *
            (hourglassAm4[i00] * h00 + hourglassAm4[i01] * h01 +
                    hourglassAm4[i02] * h02 + hourglassAm4[i03] * h03);

    hgfz[5] = coefficient *
            (hourglassAm5[i00] * h00 + hourglassAm5[i01] * h01 +
                    hourglassAm5[i02] * h02 + hourglassAm5[i03] * h03);

    hgfz[6] = coefficient *
            (hourglassAm6[i00] * h00 + hourglassAm6[i01] * h01 +
                    hourglassAm6[i02] * h02 + hourglassAm6[i03] * h03);

    hgfz[7] = coefficient *
            (hourglassAm7[i00] * h00 + hourglassAm7[i01] * h01 +
                        hourglassAm7[i02] * h02 + hourglassAm7[i03] * h03);
    
}


void LuleshParBenchmark::calcAccelerationForNodes(Domain& domain) {
    int numNode = domain.numNode();

    /**
     Time taken(in microseconds) for sequential execution:  1407.375

    Problem Size = 32, Number of Threads = 20
    Time taken(in microseconds)for the execution of this loop for different values of 'Factor'

        Factor    Time(microsec)
        1          279.2695
        2          146.2335
        4          134.3815
        8          144.8535
        16         180.262
        32         232.2245
    **/
    auto start_i = [=](int i) { return i; };
    auto end_i = [=](int i) { return i-1; };
    HCLIB_FINISH {
        int Factor = 1;
        hclib::loop_domain_t loop = {0, numNode,1,  int(numNode/(CONST_FACTOR*Factor))};
        hclib::irregular_recursion1D(& loop, domain.get_m_xdd(), start_i, end_i, [=, &domain](int i){
   
            domain.xdd(i, domain.fx(i) / domain.nodalMass(i));
            domain.ydd(i, domain.fy(i) / domain.nodalMass(i));
            domain.zdd(i, domain.fz(i) / domain.nodalMass(i));
        }, FORASYNC_MODE_RECURSIVE);
    }
}

void LuleshParBenchmark::applyAccelerationBoundaryConditionsForNodes(Domain& domain) {
    
    int numNodeBC = (domain.sizeX() + 1) * (domain.sizeX() + 1);

    {
        for (int i = 0; i < numNodeBC; ++i) {
            domain.xdd(domain.symmX(i), 0.0);
        }

        for (int i = 0; i < numNodeBC; ++i) {
            domain.ydd(domain.symmY(i), 0.0);
        }

        for (int i = 0; i < numNodeBC; ++i) {
            domain.zdd(domain.symmZ(i), 0.0);
        }
    }  
}

void LuleshParBenchmark::calcVelocityForNodes(Domain& domain, double dt, double u_cut){
    int numNode = domain.numNode();

    /**
     Time taken(in microseconds) for sequential execution:  2818.18

    Problem Size = 32, Number of Threads = 20
    Time taken(in microseconds)for the execution of this loop for different values of 'Factor'

        Factor    Time(microsec)
        1          235.883
        2          233.079
        4          212.609
        8          215.24
        16         236.285
        32         288.873
    **/
    auto start_i = [=](int i) { return i; };
    auto end_i = [=](int i) { return i-1; };
    HCLIB_FINISH {
        int Factor = 1;
        hclib::loop_domain_t loop = {0, numNode,1,  int(numNode/(CONST_FACTOR*Factor))};
        hclib::irregular_recursion1D(& loop, domain.get_m_xd(), start_i, end_i, [=, &domain](int i){

            double xdtmp, ydtmp, zdtmp;

            xdtmp = domain.xd(i) + domain.xdd(i) * dt;
            if (_compare(std::abs(xdtmp), u_cut) < 0) {
                xdtmp = 0.0;
            }
            domain.xd(i, xdtmp);

            ydtmp = domain.yd(i) + domain.ydd(i) * dt;
            if (_compare(abs(ydtmp), u_cut) < 0) {
                ydtmp = 0.0;
            }
            domain.yd(i, ydtmp);

            zdtmp = domain.zd(i) + domain.zdd(i) * dt;
            if (_compare(abs(zdtmp), u_cut) < 0) {
                zdtmp = 0.0;
            }
            domain.zd(i, zdtmp);
        }, FORASYNC_MODE_RECURSIVE);
    }
}

void LuleshParBenchmark::calcPositionForNodes(Domain& domain, double dt) {
    int numNode = domain.numNode();

    /**
     Time taken(in microseconds) for sequential execution:  1260.095

    Problem Size = 32, Number of Threads = 20
    Time taken(in microseconds)for the execution of this loop for different values of 'Factor'

        Factor    Time(microsec)
        1          153.06
        2          131.0875
        4          128.327
        8          141.7085
        16         166.6355
        32         326.173
    **/
    auto start_i = [=](int i) { return i; };
    auto end_i = [=](int i) { return i-1; };
    HCLIB_FINISH {
        int Factor = 1;
        hclib::loop_domain_t loop = {0, numNode,1,  int(numNode/(CONST_FACTOR*Factor))};
        hclib::irregular_recursion1D(& loop, domain.get_m_xd(), start_i, end_i, [=, &domain](int i){
            
            domain.x(i, domain.x(i) + (domain.xd(i) * dt));
            domain.y(i, domain.y(i) + (domain.yd(i) * dt));
            domain.z(i, domain.z(i) + (domain.zd(i) * dt));
        }, FORASYNC_MODE_RECURSIVE);
        }
}

void LuleshParBenchmark::lagrangeElements(Domain& domain){
    double deltatime = domain.deltatime();

    calcLagrangeElements(domain, deltatime);

    /* calculate Q.  (Monotonic q option requires communication) */
    calcQForElems(domain);

    applyMaterialPropertiesForElems(domain);

    updateVolumesForElems(domain);
}

void LuleshParBenchmark::calcLagrangeElements(Domain& domain, double deltatime){
   int numElem = domain.numElem();
    if (numElem > 0) {
        calcKinematicsForElems(domain, numElem, deltatime);

        // element loop to do some stuff not included in the elemlib function.

        /**
         Time taken(in microseconds) for sequential execution:  1778.835

        Problem Size = 32, Number of Threads = 20
        Time taken(in microseconds)for the execution of this loop for different values of 'Factor'

            Factor    Time(microsec)
            1          153.868
            2          153.7705
            4          147.6245
            8          156.4835
            16         175.4845
            32         311.9375
        **/
    	auto start_i = [=](int i) { return i; };
    	auto end_i = [=](int i) { return i-1; };
        HCLIB_FINISH {
            int Factor = 1;
            hclib::loop_domain_t loop = {0, numElem,1,  int(numElem/(CONST_FACTOR*Factor))};
            hclib::irregular_recursion1D(& loop, domain.get_m_dxx(), start_i, end_i, [=, &domain](int k){
                // calc strain rate and apply as constraint (only done in FB element)
                double vdov = domain.dxx(k) + domain.dyy(k) + domain.dzz(k);
                double vdovthird = vdov / 3.0;

                // make the rate of deformation tensor deviatoric
                domain.vdov(k, vdov);
                domain.dxx(k, domain.dxx(k) - vdovthird);
                domain.dyy(k, domain.dyy(k) - vdovthird);
                domain.dzz(k, domain.dzz(k) - vdovthird);

                // See if any volumes are negative, and take appropriate action.
                if (_compare(domain.vnew(k), 0.0) <= 0) {
                    throw "Volume at position " + to_string(k) + " is negative with value " + to_string(domain.vnew(k));
                }
            }, FORASYNC_MODE_RECURSIVE);
        }
    } 
}

void LuleshParBenchmark::calcKinematicsForElems(Domain& domain, int numElem, double dt){
    // loop over all elements

    /**
     Time taken(in microseconds) for sequential execution:  30588.2

    Problem Size = 32, Number of Threads = 20
    Time taken(in microseconds)for the execution of this loop for different values of 'Factor'

        Factor    Time(microsec)
        1          2214.66
        2          2038.37
        4          1966.965
        8          1902.115
        16         1817.98
        32         1926.785
    **/
    auto start_i = [=](int i) { return i; };
    auto end_i = [=](int i) { return i-1; };
    HCLIB_FINISH {
        int Factor = 9;
        hclib::loop_domain_t loop = {0, numElem,1,  int(numElem/(CONST_FACTOR*Factor))};
        hclib::irregular_recursion1D(& loop, domain.get_m_dyy(), start_i, end_i, [=, &domain](int k){
            double B[3*8];
            double D[8];
            double xLocal[8];
            double yLocal[8];
            double zLocal[8];
            double xdLocal[8];
            double ydLocal[8];
            double zdLocal[8];
            double detJ[] = {0.0};

            int *nodeList = domain.nodelist();
            int *elemToNode = &nodeList[k*8];

            // get nodal coordinates from global arrays and copy into local arrays.
            for (int lnode = 0; lnode < 8; ++lnode) {
                int gnode = elemToNode[lnode];
                xLocal[lnode] = domain.x(gnode);
                yLocal[lnode] = domain.y(gnode);
                zLocal[lnode] = domain.z(gnode);
            }
            ;

            // volume calculations
            double volume = calcElemVolume(xLocal, yLocal, zLocal);
            double relativeVolume = volume / domain.volo(k);
            domain.vnew(k, relativeVolume);
            domain.delv(k, relativeVolume - domain.v(k));


            // set characteristic length
            domain.arealg(k, calcElemCharacteristicLength(xLocal, yLocal, zLocal, volume));

            // get nodal velocities from global array and copy into local arrays.
            for (int lnode = 0; lnode < 8; ++lnode) {
                int gnode = elemToNode[lnode];
                xdLocal[lnode] = domain.xd(gnode);
                ydLocal[lnode] = domain.yd(gnode);
                zdLocal[lnode] = domain.zd(gnode);
            }

            double dt2 = 0.50 * dt;
            for (int j = 0; j < 8; ++j) {
                xLocal[j] -= dt2 * xdLocal[j];
                yLocal[j] -= dt2 * ydLocal[j];
                zLocal[j] -= dt2 * zdLocal[j];
            }

            calcElemShapeFunctionDerivatives(xLocal, yLocal, zLocal, B, detJ, 0);

            calcElemVelocityGrandient(xdLocal, ydLocal, zdLocal, B, detJ[0], D);

            // put velocity gradient quantities into their global arrays.
            domain.dxx(k, D[0]);
            domain.dyy(k, D[1]);
            domain.dzz(k, D[2]);

        }, FORASYNC_MODE_RECURSIVE);
    }
}

double LuleshParBenchmark::calcElemVolume(double x[], double y[], double z[]) {
    return calcElemVolume(x[0], x[1], x[2], x[3], x[4], x[5], x[6], x[7],
                            y[0], y[1], y[2], y[3], y[4], y[5], y[6], y[7],
                            z[0], z[1], z[2], z[3], z[4], z[5], z[6], z[7]);
}

double LuleshParBenchmark::calcElemVolume(
        double x0, double x1,
        double x2, double x3,
        double x4, double x5,
        double x6, double x7,
        double y0, double y1,
        double y2, double y3,
        double y4, double y5,
        double y6, double y7,
        double z0, double z1,
        double z2, double z3,
        double z4, double z5,
        double z6, double z7) {
    
    double twelveth = 1.0 / 12.0;

        double dx61 = x6 - x1;
        double dy61 = y6 - y1;
        double dz61 = z6 - z1;

        double dx70 = x7 - x0;
        double dy70 = y7 - y0;
        double dz70 = z7 - z0;

        double dx63 = x6 - x3;
        double dy63 = y6 - y3;
        double dz63 = z6 - z3;

        double dx20 = x2 - x0;
        double dy20 = y2 - y0;
        double dz20 = z2 - z0;

        double dx50 = x5 - x0;
        double dy50 = y5 - y0;
        double dz50 = z5 - z0;

        double dx64 = x6 - x4;
        double dy64 = y6 - y4;
        double dz64 = z6 - z4;

        double dx31 = x3 - x1;
        double dy31 = y3 - y1;
        double dz31 = z3 - z1;

        double dx72 = x7 - x2;
        double dy72 = y7 - y2;
        double dz72 = z7 - z2;

        double dx43 = x4 - x3;
        double dy43 = y4 - y3;
        double dz43 = z4 - z3;

        double dx57 = x5 - x7;
        double dy57 = y5 - y7;
        double dz57 = z5 - z7;

        double dx14 = x1 - x4;
        double dy14 = y1 - y4;
        double dz14 = z1 - z4;

        double dx25 = x2 - x5;
        double dy25 = y2 - y5;
        double dz25 = z2 - z5;

        double volume = twelveth *
                (tripleProduct(dx31 + dx72, dx63, dx20, dy31 + dy72, dy63, dy20, dz31 + dz72, dz63, dz20) +
                        tripleProduct(dx43 + dx57, dx64, dx70, dy43 + dy57, dy64, dy70, dz43 + dz57, dz64, dz70) +
                        tripleProduct(dx14 + dx25, dx61, dx50, dy14 + dy25, dy61, dy50, dz14 + dz25, dz61, dz50));

        return volume;
}

double LuleshParBenchmark::tripleProduct(
        double x1, double y1, double z1,
        double x2, double y2, double z2,
        double x3, double y3, double z3) {
    return ((x1) * ((y2) * (z3) - (z2) * (y3)) +
                (x2) * ((z1) * (y3) - (y1) * (z3)) +
                (x3) * ((y1) * (z2) - (z1) * (y2)));
}

double LuleshParBenchmark::calcElemCharacteristicLength(double x[], double y[], double z[], double volume){
    double a;
    double charLength = 0.0;

    a = areaFace(x[0], x[1], x[2], x[3],
                    y[0], y[1], y[2], y[3],
                    z[0], z[1], z[2], z[3]);
    charLength = max(a, charLength);

    a = areaFace(x[4], x[5], x[6], x[7],
                    y[4], y[5], y[6], y[7],
                    z[4], z[5], z[6], z[7]);
    charLength = max(a, charLength);

    a = areaFace(x[0], x[1], x[5], x[4],
                    y[0], y[1], y[5], y[4],
                    z[0], z[1], z[5], z[4]);
    charLength = max(a, charLength);

    a = areaFace(x[1], x[2], x[6], x[5],
                    y[1], y[2], y[6], y[5],
                    z[1], z[2], z[6], z[5]);
    charLength = max(a, charLength);

    a = areaFace(x[2], x[3], x[7], x[6],
                    y[2], y[3], y[7], y[6],
                    z[2], z[3], z[7], z[6]);
    charLength = max(a, charLength);

    a = areaFace(x[3], x[0], x[4], x[7],
                    y[3], y[0], y[4], y[7],
                    z[3], z[0], z[4], z[7]);
    charLength = max(a, charLength);

    charLength = (4.0) * volume / std::sqrt(charLength);

    return charLength;
}

double LuleshParBenchmark::areaFace(
        double x0, double x1,
        double x2, double x3,
        double y0, double y1,
        double y2, double y3,
        double z0, double z1,
        double z2, double z3){
    double fx = (x2 - x0) - (x3 - x1);
    double fy = (y2 - y0) - (y3 - y1);
    double fz = (z2 - z0) - (z3 - z1);
    double gx = (x2 - x0) + (x3 - x1);
    double gy = (y2 - y0) + (y3 - y1);
    double gz = (z2 - z0) + (z3 - z1);
    double area =
            (fx * fx + fy * fy + fz * fz) *
                    (gx * gx + gy * gy + gz * gz) -
                    (fx * gx + fy * gy + fz * gz) *
                            (fx * gx + fy * gy + fz * gz);
    return area;
}

void LuleshParBenchmark::calcElemShapeFunctionDerivatives(
        double x[], double y[], double z[], double b[],
        double determinant[], int determinantIndex) {
    
    double x0 = x[0];
    double x1 = x[1];
    double x2 = x[2];
    double x3 = x[3];
    double x4 = x[4];
    double x5 = x[5];
    double x6 = x[6];
    double x7 = x[7];

    double y0 = y[0];
    double y1 = y[1];
    double y2 = y[2];
    double y3 = y[3];
    double y4 = y[4];
    double y5 = y[5];
    double y6 = y[6];
    double y7 = y[7];

    double z0 = z[0];
    double z1 = z[1];
    double z2 = z[2];
    double z3 = z[3];
    double z4 = z[4];
    double z5 = z[5];
    double z6 = z[6];
    double z7 = z[7];

    double fjxxi;
    double fjxet;
    double fjxze;
    double fjyxi;
    double fjyet;
    double fjyze;
    double fjzxi;
    double fjzet;
    double fjzze;
    double cjxxi;
    double cjxet;
    double cjxze;
    double cjyxi;
    double cjyet;
    double cjyze;
    double cjzxi;
    double cjzet;
    double cjzze;

    fjxxi = 0.125 * ((x6 - x0) + (x5 - x3) - (x7 - x1) - (x4 - x2));
    fjxet = 0.125 * ((x6 - x0) - (x5 - x3) + (x7 - x1) - (x4 - x2));
    fjxze = 0.125 * ((x6 - x0) + (x5 - x3) + (x7 - x1) + (x4 - x2));

    fjyxi = 0.125 * ((y6 - y0) + (y5 - y3) - (y7 - y1) - (y4 - y2));
    fjyet = 0.125 * ((y6 - y0) - (y5 - y3) + (y7 - y1) - (y4 - y2));
    fjyze = 0.125 * ((y6 - y0) + (y5 - y3) + (y7 - y1) + (y4 - y2));

    fjzxi = 0.125 * ((z6 - z0) + (z5 - z3) - (z7 - z1) - (z4 - z2));
    fjzet = 0.125 * ((z6 - z0) - (z5 - z3) + (z7 - z1) - (z4 - z2));
    fjzze = 0.125 * ((z6 - z0) + (z5 - z3) + (z7 - z1) + (z4 - z2));

    /* compute cofactors */
    cjxxi = (fjyet * fjzze) - (fjzet * fjyze);
    cjxet = -(fjyxi * fjzze) + (fjzxi * fjyze);
    cjxze = (fjyxi * fjzet) - (fjzxi * fjyet);

    cjyxi = -(fjxet * fjzze) + (fjzet * fjxze);
    cjyet = (fjxxi * fjzze) - (fjzxi * fjxze);
    cjyze = -(fjxxi * fjzet) + (fjzxi * fjxet);

    cjzxi = (fjxet * fjyze) - (fjyet * fjxze);
    cjzet = -(fjxxi * fjyze) + (fjyxi * fjxze);
    cjzze = (fjxxi * fjyet) - (fjyxi * fjxet);

    /* calculate partials: this need only be done for l = 0,1,2,3. Since, by symmetry, (6,7,4,5) = - (0,1,2,3).
    */
    b[0*8+0] = -cjxxi - cjxet - cjxze;
    b[0*8+1] = cjxxi - cjxet - cjxze;
    b[0*8+2] = cjxxi + cjxet - cjxze;
    b[0*8+3] = -cjxxi + cjxet - cjxze;
    b[0*8+4] = -b[0*8+2];
    b[0*8+5] = -b[0*8+3];
    b[0*8+6] = -b[0*8+0];
    b[0*8+7] = -b[0*8+1];

    b[1*8+0] = -cjyxi - cjyet - cjyze;
    b[1*8+1] = cjyxi - cjyet - cjyze;
    b[1*8+2] = cjyxi + cjyet - cjyze;
    b[1*8+3] = -cjyxi + cjyet - cjyze;
    b[1*8+4] = -b[1*8+2];
    b[1*8+5] = -b[1*8+3];
    b[1*8+6] = -b[1*8+0];
    b[1*8+7] = -b[1*8+1];

    b[2*8+0] = -cjzxi - cjzet - cjzze;
    b[2*8+1] = cjzxi - cjzet - cjzze;
    b[2*8+2] = cjzxi + cjzet - cjzze;
    b[2*8+3] = -cjzxi + cjzet - cjzze;
    b[2*8+4] = -b[2*8+2];
    b[2*8+5] = -b[2*8+3];
    b[2*8+6] = -b[2*8+0];
    b[2*8+7] = -b[2*8+1];

    /* calculate jacobian determinant (volume) */
    determinant[determinantIndex] = 8.0 * (fjxet * cjxet + fjyet * cjyet + fjzet * cjzet);
}


void LuleshParBenchmark::calcElemVelocityGrandient(
        double xVel[], double yVel[], double zVel[],
        double b[], double detJ, double d[]) {
    
    double inv_detJ = 1.0 / detJ;
    double dyddx;
    double dxddy;
    double dzddx;
    double dxddz;
    double dzddy;
    double dyddz;
    double *pfx = &b[0*8];
    double *pfy = &b[1*8];
    double *pfz = &b[2*8];

    d[0] = inv_detJ * (pfx[0] * (xVel[0] - xVel[6])
            + pfx[1] * (xVel[1] - xVel[7])
            + pfx[2] * (xVel[2] - xVel[4])
            + pfx[3] * (xVel[3] - xVel[5]));

    d[1] = inv_detJ * (pfy[0] * (yVel[0] - yVel[6])
            + pfy[1] * (yVel[1] - yVel[7])
            + pfy[2] * (yVel[2] - yVel[4])
            + pfy[3] * (yVel[3] - yVel[5]));

    d[2] = inv_detJ * (pfz[0] * (zVel[0] - zVel[6])
            + pfz[1] * (zVel[1] - zVel[7])
            + pfz[2] * (zVel[2] - zVel[4])
            + pfz[3] * (zVel[3] - zVel[5]));

    dyddx = inv_detJ * (pfx[0] * (yVel[0] - yVel[6])
            + pfx[1] * (yVel[1] - yVel[7])
            + pfx[2] * (yVel[2] - yVel[4])
            + pfx[3] * (yVel[3] - yVel[5]));

    dxddy = inv_detJ * (pfy[0] * (xVel[0] - xVel[6])
            + pfy[1] * (xVel[1] - xVel[7])
            + pfy[2] * (xVel[2] - xVel[4])
            + pfy[3] * (xVel[3] - xVel[5]));

    dzddx = inv_detJ * (pfx[0] * (zVel[0] - zVel[6])
            + pfx[1] * (zVel[1] - zVel[7])
            + pfx[2] * (zVel[2] - zVel[4])
            + pfx[3] * (zVel[3] - zVel[5]));

    dxddz = inv_detJ * (pfz[0] * (xVel[0] - xVel[6])
            + pfz[1] * (xVel[1] - xVel[7])
            + pfz[2] * (xVel[2] - xVel[4])
            + pfz[3] * (xVel[3] - xVel[5]));

    dzddy = inv_detJ * (pfy[0] * (zVel[0] - zVel[6])
            + pfy[1] * (zVel[1] - zVel[7])
            + pfy[2] * (zVel[2] - zVel[4])
            + pfy[3] * (zVel[3] - zVel[5]));

    dyddz = inv_detJ * (pfz[0] * (yVel[0] - yVel[6])
            + pfz[1] * (yVel[1] - yVel[7])
            + pfz[2] * (yVel[2] - yVel[4])
            + pfz[3] * (yVel[3] - yVel[5]));

    d[5] = 0.50 * (dxddy + dyddx);
    d[4] = 0.50 * (dxddz + dzddx);
    d[3] = 0.50 * (dzddy + dyddz);
}

void LuleshParBenchmark::calcQForElems(Domain& domain){
    double qstop = domain.qstop();
    int numElem = domain.numElem();

    // MONOTONIC Q option

    /* calculate velocity gradients */
    calcMonotonicQGradientsForElems(domain);

    /* Transfer veloctiy gradients in the first order elements */
    /* problem->commElements->Transfer(CommElements::monoQ) ; */
    calcMonotonicQForElems(domain);

    /* Don't allow excessive artificial viscosity */
    if (numElem != 0) {
        int index = -1;
        for (int i = 0; i < numElem; ++i) {
            if (domain.q(i) > qstop) {
                index = i;
                break;
            }
        }

        if (index >= 0) {
            throw "Q_STOP_ERROR: index = " + to_string(index);
        }
    }
}

void LuleshParBenchmark::calcMonotonicQGradientsForElems(Domain& domain){
    int numElem = domain.numElem();

    /**
     Time taken(in microseconds) for sequential execution:  10084.35

    Problem Size = 32, Number of Threads = 20
    Time taken(in microseconds)for the execution of this loop for different values of 'Factor'

        Factor    Time(microsec)
        1          1031.7425
        2          726.9615
        4          999.734
        8          638.0745
        16         648.4075
        32         648.8595
    **/
    auto start_i = [=](int i) { return i; };
    auto end_i = [=](int i) { return i-1; };
    HCLIB_FINISH {
        int Factor = 8;
        hclib::loop_domain_t loop = {0, numElem,1,  int(numElem/(CONST_FACTOR*Factor))};
        hclib::irregular_recursion1D(&loop, domain.get_m_delx_eta(), start_i, end_i, [=, &domain](int i){
            double ptiny = 1.e-36;
            double ax, ay, az;
            double dxv, dyv, dzv;

            int *nodeList = domain.nodelist();
            int *elemToNode = &nodeList[i*8];

            int n0 = elemToNode[0];
            int n1 = elemToNode[1];
            int n2 = elemToNode[2];
            int n3 = elemToNode[3];
            int n4 = elemToNode[4];
            int n5 = elemToNode[5];
            int n6 = elemToNode[6];
            int n7 = elemToNode[7];

            double x0 = domain.x(n0);
            double x1 = domain.x(n1);
            double x2 = domain.x(n2);
            double x3 = domain.x(n3);
            double x4 = domain.x(n4);
            double x5 = domain.x(n5);
            double x6 = domain.x(n6);
            double x7 = domain.x(n7);

            double y0 = domain.y(n0);
            double y1 = domain.y(n1);
            double y2 = domain.y(n2);
            double y3 = domain.y(n3);
            double y4 = domain.y(n4);
            double y5 = domain.y(n5);
            double y6 = domain.y(n6);
            double y7 = domain.y(n7);

            double z0 = domain.z(n0);
            double z1 = domain.z(n1);
            double z2 = domain.z(n2);
            double z3 = domain.z(n3);
            double z4 = domain.z(n4);
            double z5 = domain.z(n5);
            double z6 = domain.z(n6);
            double z7 = domain.z(n7);

            double xv0 = domain.xd(n0);
            double xv1 = domain.xd(n1);
            double xv2 = domain.xd(n2);
            double xv3 = domain.xd(n3);
            double xv4 = domain.xd(n4);
            double xv5 = domain.xd(n5);
            double xv6 = domain.xd(n6);
            double xv7 = domain.xd(n7);

            double yv0 = domain.yd(n0);
            double yv1 = domain.yd(n1);
            double yv2 = domain.yd(n2);
            double yv3 = domain.yd(n3);
            double yv4 = domain.yd(n4);
            double yv5 = domain.yd(n5);
            double yv6 = domain.yd(n6);
            double yv7 = domain.yd(n7);

            double zv0 = domain.zd(n0);
            double zv1 = domain.zd(n1);
            double zv2 = domain.zd(n2);
            double zv3 = domain.zd(n3);
            double zv4 = domain.zd(n4);
            double zv5 = domain.zd(n5);
            double zv6 = domain.zd(n6);
            double zv7 = domain.zd(n7);

            double vol = domain.volo(i) * domain.vnew(i);
            double norm = 1.0 / (vol + ptiny);

            double dxj = (-0.25) * ((x0 + x1 + x5 + x4) - (x3 + x2 + x6 + x7));
            double dyj = (-0.25) * ((y0 + y1 + y5 + y4) - (y3 + y2 + y6 + y7));
            double dzj = (-0.25) * ((z0 + z1 + z5 + z4) - (z3 + z2 + z6 + z7));

            double dxi = (0.25) * ((x1 + x2 + x6 + x5) - (x0 + x3 + x7 + x4));
            double dyi = (0.25) * ((y1 + y2 + y6 + y5) - (y0 + y3 + y7 + y4));
            double dzi = (0.25) * ((z1 + z2 + z6 + z5) - (z0 + z3 + z7 + z4));

            double dxk = (0.25) * ((x4 + x5 + x6 + x7) - (x0 + x1 + x2 + x3));
            double dyk = (0.25) * ((y4 + y5 + y6 + y7) - (y0 + y1 + y2 + y3));
            double dzk = (0.25) * ((z4 + z5 + z6 + z7) - (z0 + z1 + z2 + z3));

            /* find delvk and delxk ( i cross j ) */

            ax = dyi * dzj - dzi * dyj;
            ay = dzi * dxj - dxi * dzj;
            az = dxi * dyj - dyi * dxj;

            domain.delx_zeta(i, vol / std::sqrt(ax * ax + ay * ay + az * az + ptiny));

            ax *= norm;
            ay *= norm;
            az *= norm;

            dxv = (0.25) * ((xv4 + xv5 + xv6 + xv7) - (xv0 + xv1 + xv2 + xv3));
            dyv = (0.25) * ((yv4 + yv5 + yv6 + yv7) - (yv0 + yv1 + yv2 + yv3));
            dzv = (0.25) * ((zv4 + zv5 + zv6 + zv7) - (zv0 + zv1 + zv2 + zv3));

            domain.delv_zeta(i, ax * dxv + ay * dyv + az * dzv);

            /* find delxi and delvi ( j cross k ) */

            ax = dyj * dzk - dzj * dyk;
            ay = dzj * dxk - dxj * dzk;
            az = dxj * dyk - dyj * dxk;

            domain.delx_xi(i, vol / std::sqrt(ax * ax + ay * ay + az * az + ptiny));

            ax *= norm;
            ay *= norm;
            az *= norm;

            dxv = (0.25) * ((xv1 + xv2 + xv6 + xv5) - (xv0 + xv3 + xv7 + xv4));
            dyv = (0.25) * ((yv1 + yv2 + yv6 + yv5) - (yv0 + yv3 + yv7 + yv4));
            dzv = (0.25) * ((zv1 + zv2 + zv6 + zv5) - (zv0 + zv3 + zv7 + zv4));

            domain.delv_xi(i, ax * dxv + ay * dyv + az * dzv);

            /* find delxj and delvj ( k cross i ) */

            ax = dyk * dzi - dzk * dyi;
            ay = dzk * dxi - dxk * dzi;
            az = dxk * dyi - dyk * dxi;

            domain.delx_eta(i, vol / std::sqrt(ax * ax + ay * ay + az * az + ptiny));

            ax *= norm;
            ay *= norm;
            az *= norm;

            dxv = (-0.25) * ((xv0 + xv1 + xv5 + xv4) - (xv3 + xv2 + xv6 + xv7));
            dyv = (-0.25) * ((yv0 + yv1 + yv5 + yv4) - (yv3 + yv2 + yv6 + yv7));
            dzv = (-0.25) * ((zv0 + zv1 + zv5 + zv4) - (zv3 + zv2 + zv6 + zv7));

            domain.delv_eta(i, ax * dxv + ay * dyv + az * dzv);
        }, FORASYNC_MODE_RECURSIVE);
    }
}

void LuleshParBenchmark::calcMonotonicQForElems(Domain& domain){

    // initialize parameters
    // 
    double ptiny = 1.e-36;
    double monoqMaxSlope = domain.monoq_max_slope();
    double monoqLimiterMult = domain.monoq_limiter_mult();

    //
    // calculate the monotonic q for pure regions
    //
    int elength = domain.numElem();
    if (elength > 0) {
        double qlcMonoq = domain.qlc_monoq();
        double qqcMonoq = domain.qqc_monoq();
        calcMonotonicQRegionForElems(
                domain,
                qlcMonoq,
                qqcMonoq,
                monoqLimiterMult,
                monoqMaxSlope,
                ptiny,
                elength);
    }
}

void LuleshParBenchmark::calcMonotonicQRegionForElems(
        Domain& domain, double qlcMonoq, double qqcMonoq,
        double monoqLimiterMult, double monoqMaxSlope,
        double ptiny, int elength){
    
    // TODO #pragma omp parallel for firstprivate(elength, qlcMonoq, qqcMonoq, monoqLimiterMult, monoqMaxSlope, ptiny)

    /**
     Time taken(in microseconds) for sequential execution:  5385.445

    Problem Size = 32, Number of Threads = 20
    Time taken(in microseconds)for the execution of this loop for different values of 'Factor'

        Factor    Time(microsec)
        1          463.8875
        2          409.6305
        4          378.821
        8          396.7955
        16         392.118
        32         411.0545
    **/
    auto start_i = [=](int i) { return i; };
    auto end_i = [=](int i) { return i-1; };
    HCLIB_FINISH {
        int Factor = 8;
        hclib::loop_domain_t loop = {0, elength,1,  int(elength/(CONST_FACTOR*Factor))};
        hclib::irregular_recursion1D(& loop, domain.get_m_delv_zeta(), start_i, end_i, [=, &domain](int ielem){
            double qlin;
            double qquad;
            double phixi, phieta, phizeta;
            int i = domain.matElemlist(ielem);
            int bcMask = domain.elemBC(i);
            double delvm, delvp;

            /*  phixi     */
            double norm = 1.0 / (domain.delv_xi(i) + ptiny);

            switch (bcMask & LuleshConfig::XI_M) {
                case 0:
                    delvm = domain.delv_xi(domain.lxim(i));
                    break;
                case 0x001: //**XI_M_SYMM
                    delvm = domain.delv_xi(i);
                    break;
                case 0x002: //**XI_M_FREE:
                    delvm = 0.0;
                    break;
                default:        /* ERROR */
                    throw "Error in computing delmv";
            }
            switch (bcMask & LuleshConfig::XI_P) {
                case 0:
                    delvp = domain.delv_xi(domain.lxip(i));
                    break;
                case 0x004: //** XI_P_SYMM:
                    delvp = domain.delv_xi(i);
                    break;
                case 0x008: //** XI_P_FREE:
                    delvp = 0.0;
                    break;
                default:        /* ERROR */
                    throw "Error in computing delvp";
            }

            delvm = delvm * norm;
            delvp = delvp * norm;

            phixi = 0.50 * (delvm + delvp);

            delvm *= monoqLimiterMult;
            delvp *= monoqLimiterMult;

            if (delvm < phixi) {
                phixi = delvm;
            }
            if (delvp < phixi) {
                phixi = delvp;
            }
            if (phixi < 0.0) {
                phixi = 0.0;
            }
            if (phixi > monoqMaxSlope) {
                phixi = monoqMaxSlope;
            }

        
            /*  phieta     */
            norm = 1.0 / (domain.delv_eta(i) + ptiny);

            switch (bcMask & LuleshConfig::ETA_M) {
                case 0:
                    delvm = domain.delv_eta(domain.letam(i));
                    break;
                case 0x010: //** ETA_M_SYMM:
                    delvm = domain.delv_eta(i);
                    break;
                case 0x020: //** ETA_M_FREE:
                    delvm = 0.0;
                    break;
                default:         /* ERROR */
                    break;
            }
            switch (bcMask & LuleshConfig::ETA_P) {
                case 0:
                    delvp = domain.delv_eta(domain.letap(i));
                    break;
                case 0x040: //** ETA_P_SYMM:
                    delvp = domain.delv_eta(i);
                    break;
                case 0x080: //** ETA_P_FREE:
                    delvp = 0.0;
                    break;
                default:         /* ERROR */
                    break;
            }

            delvm = delvm * norm;
            delvp = delvp * norm;

            phieta = 0.50 * (delvm + delvp);

            delvm *= monoqLimiterMult;
            delvp *= monoqLimiterMult;

            if (delvm < phieta) {
                phieta = delvm;
            }
            if (delvp < phieta) {
                phieta = delvp;
            }
            if (phieta < 0.0) {
                phieta = 0.0;
            }
            if (phieta > monoqMaxSlope) {
                phieta = monoqMaxSlope;
            }
        
            /*  phizeta     */
            norm = 1.0 / (domain.delv_zeta(i) + ptiny);

            switch (bcMask & LuleshConfig::ZETA_M) {
                case 0:
                    delvm = domain.delv_zeta(domain.lzetam(i));
                    break;
                case 0x100: //**ZETA_M_SYMM:
                    delvm = domain.delv_zeta(i);
                    break;
                case 0x200: //** ZETA_M_FREE:
                    delvm = 0.0;
                    break;
                default:          /* ERROR */
                    break;
            }
            switch (bcMask & LuleshConfig::ZETA_P) {
                case 0:
                    delvp = domain.delv_zeta(domain.lzetap(i));
                    break;
                case 0x400: //** ZETA_P_SYMM:
                    delvp = domain.delv_zeta(i);
                    break;
                case 0x800: //** ZETA_P_FREE:
                    delvp = 0.0;
                    break;
                default:          /* ERROR */
                    break;
            }

            delvm = delvm * norm;
            delvp = delvp * norm;

            phizeta = 0.50 * (delvm + delvp);

            delvm *= monoqLimiterMult;
            delvp *= monoqLimiterMult;

            if (delvm < phizeta) {
                phizeta = delvm;
            }
            if (delvp < phizeta) {
                phizeta = delvp;
            }
            if (phizeta < 0.0) {
                phizeta = 0.0;
            }
            if (phizeta > monoqMaxSlope) {
                phizeta = monoqMaxSlope;
            }

            /* Remove length scale */

            if (domain.vdov(i) > 0.0) {
                qlin = 0.0;
                qquad = 0.0;
            } else {
                double delvxxi = domain.delv_xi(i) * domain.delx_xi(i);
                double delvxeta = domain.delv_eta(i) * domain.delx_eta(i);
                double delvxzeta = domain.delv_zeta(i) * domain.delx_zeta(i);

                if (delvxxi > 0.0) {
                    delvxxi = 0.0;
                }
                if (delvxeta > 0.0) {
                    delvxeta = 0.0;
                }
                if (delvxzeta > 0.0) {
                    delvxzeta = 0.0;
                }

                double rho = domain.elemMass(i) / (domain.volo(i) * domain.vnew(i));

                qlin = -qlcMonoq * rho * (delvxxi * (1.0 - phixi) + delvxeta * (1.0 - phieta) + delvxzeta * (1.0 - phizeta));

                qquad = qqcMonoq * rho *
                        (delvxxi * delvxxi * (1.0 - phixi * phixi) +
                                delvxeta * delvxeta * (1.0 - phieta * phieta) +
                                delvxzeta * delvxzeta * (1.0 - phizeta * phizeta));
            }

            domain.qq(i, qquad);
            domain.ql(i, qlin);
        }, FORASYNC_MODE_RECURSIVE);
    }
}

void LuleshParBenchmark::applyMaterialPropertiesForElems(Domain& domain){
    int length = domain.numElem();

    if (length != 0) {
        /* Expose all of the variables needed for material evaluation */
        double eosvmin = domain.eosvmin();
        double eosvmax = domain.eosvmax();
        
        /**
         Time taken(in microseconds) for sequential execution:  1515.22

        Problem Size = 32, Number of Threads = 20
        Time taken(in microseconds)for the execution of this loop for different values of 'Factor'

            Factor    Time(microsec)
            1          272.074
            2          140.143
            4          153.4745
            8          162.0615
            16         230.8245
            32         371.9965
        **/
    	auto start_i = [=](int i) { return i; };
    	auto end_i = [=](int i) { return i-1; };
        HCLIB_FINISH {
            int Factor = 1;
            hclib::loop_domain_t loop = {0, length,1,  int(length/(CONST_FACTOR*Factor))};
            hclib::irregular_recursion1D(& loop, domain.get_m_v(), start_i, end_i, [&](int i){
                int zn1 = domain.matElemlist(i);
                vnewc[i] = domain.vnew(zn1);

                if (_compare(eosvmin, 0) != 0) {
                    if (vnewc[i] < eosvmin) {
                        vnewc[i] = eosvmin;
                    }
                }

                if (_compare(eosvmax, 0) != 0) {
                    if (vnewc[i] > eosvmax) {
                        vnewc[i] = eosvmax;
                    }
                }

                int zn2 = domain.matElemlist(i);
                double vc = domain.v(zn2);

                if (_compare(eosvmin, 0) != 0) {
                    if (vc < eosvmin) {
                        vc = eosvmin;
                    }
                }
                if (_compare(eosvmax, 0) != 0) {
                    if (vc > eosvmax) {
                        vc = eosvmax;
                    }
                }
                if (vc <= 0.) {
                    throw "Negative volume: " + to_string(vc);
                }
            }, FORASYNC_MODE_RECURSIVE);
        }
        evalEOSForElems(domain, vnewc, length);
    }
}

void LuleshParBenchmark::freeData() {
    hclib::numa_free(e_old);
    hclib::numa_free(delvc);
    hclib::numa_free(p_old);
    hclib::numa_free(q_old);
    hclib::numa_free(compression);
    hclib::numa_free(compHalfStep);
    hclib::numa_free(qq);
    hclib::numa_free(ql);
    hclib::numa_free(work);
    hclib::numa_free(p_new);
    hclib::numa_free(e_new);
    hclib::numa_free(q_new);
    hclib::numa_free(bvc);
    hclib::numa_free(pbvc);
    hclib::numa_free(pHalfStep);
    hclib::numa_free(vnewc);
}
 
void LuleshParBenchmark::initData_evalEOSForElems(int length) {
    pHalfStep = hclib::numa_malloc<double>(length);
    vnewc = hclib::numa_malloc<double>(length);
    e_old = hclib::numa_malloc<double>(length);
    delvc = hclib::numa_malloc<double>(length);
    p_old = hclib::numa_malloc<double>(length);
    q_old = hclib::numa_malloc<double>(length);
    compression = hclib::numa_malloc<double>(length);
    compHalfStep = hclib::numa_malloc<double>(length);
    qq = hclib::numa_malloc<double>(length);
    ql = hclib::numa_malloc<double>(length);
    work = hclib::numa_malloc<double>(length);
    p_new = hclib::numa_malloc<double>(length);
    e_new = hclib::numa_malloc<double>(length);
    q_new = hclib::numa_malloc<double>(length);
    bvc = hclib::numa_malloc<double>(length);
    pbvc = hclib::numa_malloc<double>(length);
}

void LuleshParBenchmark::evalEOSForElems(Domain& domain, double vnewc[], int length){
    double e_cut = domain.e_cut();
    double p_cut = domain.p_cut();
    double ss4o3 = domain.ss4o3();
    double q_cut = domain.q_cut();

    double eosvmax = domain.eosvmax();
    double eosvmin = domain.eosvmin();
    double pmin = domain.pmin();
    double emin = domain.emin();
    double rho0 = domain.refdens();

    /* compress data, minimal set */
    bool eosvminIsZero = _compare(eosvmin, 0) != 0;
    bool eosvmaxIsZero = _compare(eosvmax, 0) != 0;

    // Note: Loop Fusion

    /**
     Time taken(in microseconds) for sequential execution:  1363.985

    Problem Size = 32, Number of Threads = 20
    Time taken(in microseconds)for the execution of this loop for different values of 'Factor'

        Factor    Time(microsec)
        1          243.4195
        2          241.3515
        4          208.618
        8          240.667
        16         227.018
        32         360.1005
    **/
    auto start_i = [=](int i) { return i; };
    auto end_i = [=](int i) { return i-1; };
    HCLIB_FINISH {
        int Factor = 9;
        hclib::loop_domain_t loop = {0, length,1,  int(length/(CONST_FACTOR*Factor))};
        hclib::irregular_recursion1D(& loop, domain.get_m_e(), start_i, end_i, [=,&domain](int i){
            int zIndex = domain.matElemlist(i);
            e_old[i] = domain.e(zIndex);
            delvc[i] = domain.delv(zIndex);
            p_old[i] = domain.p(zIndex);
            q_old[i] = domain.q(zIndex);

            compression[i] = 1.0 / vnewc[i] - 1.0;
            double vchalf = vnewc[i] - delvc[i] * 0.50;
            compHalfStep[i] = 1.0 / vchalf - 1.0;

            /* Check for v > eosvmax or v < eosvmin */
            if (eosvminIsZero) {
                if (vnewc[i] <= eosvmin) { /* impossible due to calling func? */
                    compHalfStep[i] = compression[i];
                }
            }
            if (eosvmaxIsZero) {
                if (vnewc[i] >= eosvmax) { /* impossible due to calling func? */
                    p_old[i] = 0.0;
                    compression[i] = 0.0;
                    compHalfStep[i] = 0.0;
                }
            }

            qq[i] = domain.qq(zIndex);
            ql[i] = domain.ql(zIndex);
            work[i] = 0.0;
        }, FORASYNC_MODE_RECURSIVE);
    }

    calcEnergyForElems(p_new, e_new, q_new, bvc, pbvc,
                        p_old, e_old, q_old, compression, compHalfStep,
                        vnewc, work, delvc, pmin,
                        p_cut, e_cut, q_cut, emin,
                        qq, ql, rho0, eosvmax, length);

    // Note: Loop Fusion

    /**
     Time taken(in microseconds) for sequential execution:  724.2495

    Problem Size = 32, Number of Threads = 20
    Time taken(in microseconds)for the execution of this loop for different values of 'Factor'

        Factor    Time(microsec)
        1          89.7827
        2          93.7881
        4          119.7275
        8          109.205
        16         279.153
        32         318.5055
    **/
    HCLIB_FINISH {
        int Factor = 1;
        hclib::loop_domain_t loop = {0, length,1,  int(length/(CONST_FACTOR*Factor))};
        hclib::irregular_recursion1D(& loop, domain.get_m_p(), start_i, end_i, [=,&domain](int i){
            int zIndex = domain.matElemlist(i);
            domain.p(zIndex, p_new[i]);
            domain.e(zIndex, e_new[i]);
            domain.q(zIndex, q_new[i]);
        }, FORASYNC_MODE_RECURSIVE);
    }

    calcSoundSpeedForElems(domain, vnewc, rho0, e_new, p_new, pbvc, bvc, ss4o3, length);

}

void LuleshParBenchmark::calcEnergyForElems(double p_new[], double e_new[], double q_new[],
                                        double bvc[], double pbvc[],
                                        double p_old[], double e_old[], double q_old[],
                                        double compression[], double compHalfStep[],
                                        double vnewc[], double work[], double delvc[], double pmin,
                                        double p_cut, double e_cut, double q_cut, double emin,
                                        double qq[], double ql[],
                                        double rho0,
                                        double eosvmax,
                                        int length) {

    /**
     Time taken(in microseconds) for sequential execution:  276.0055

    Problem Size = 32, Number of Threads = 20
    Time taken(in microseconds)for the execution of this loop for different values of 'Factor'

        Factor    Time(microsec)
        1          66.1703
        2          61.78685
        4          76.5225
        8          429.86915
        16         146.321
        32         317.879
    **/
    auto start_i = [=](int i) { return i; };
    auto end_i = [=](int i) { return i-1; };
    HCLIB_FINISH {
        int Factor = 1;
        hclib::loop_domain_t loop = {0, length,1,  int(length/(CONST_FACTOR*Factor))};
        hclib::irregular_recursion1D(&loop, e_new, start_i, end_i, [=](int i){
            e_new[i] = e_old[i] - 0.50 * delvc[i] * (p_old[i] + q_old[i]) + 0.50 * work[i];

            if (e_new[i] < emin) {
                e_new[i] = emin;
            }
        
        }, FORASYNC_MODE_RECURSIVE);
    }
    

    calcPressureForElems(pHalfStep, bvc, pbvc, e_new, compHalfStep, vnewc,
                            pmin, p_cut, eosvmax, length);

    /**
     Time taken(in microseconds) for sequential execution:  1426.155

    Problem Size = 32, Number of Threads = 20
    Time taken(in microseconds)for the execution of this loop for different values of 'Factor'

        Factor    Time(microsec)
        1          170.613
        2          169.007
        4          151.235
        8          163.856
        16         208.9725
        32         314.6415
    **/
    HCLIB_FINISH {
        int Factor = 1;
        hclib::loop_domain_t loop = {0, length,1,  int(length/(CONST_FACTOR*Factor))};
        hclib::irregular_recursion1D(&loop, e_new, start_i, end_i, [=](int i){
            double vhalf = 1.0 / (1.0 + compHalfStep[i]);

            if (delvc[i] > 0.0) {
                q_new[i] /* = qq[i] = ql[i] */ = 0.0;
            } else {
                double ssc = (pbvc[i] * e_new[i]
                        + vhalf * vhalf * bvc[i] * pHalfStep[i]) / rho0;

                if (_compare(ssc, 0.111111e-36) <= 0) {
                    ssc = 0.333333e-18;
                } else {
                    ssc = std::sqrt(ssc);
                }

                q_new[i] = (ssc * ql[i] + qq[i]);
            }

            e_new[i] = e_new[i] + 0.50 * delvc[i] * (3.0 * (p_old[i] + q_old[i]) - 4.0 * (pHalfStep[i] + q_new[i]));

            e_new[i] += 0.50 * work[i];

            if (_compare(abs(e_new[i]), e_cut) < 0) {
                e_new[i] = 0.0;
            }
            if (e_new[i] < emin) {
                e_new[i] = emin;
            }
        }, FORASYNC_MODE_RECURSIVE);
    }

    calcPressureForElems(p_new, bvc, pbvc, e_new, compression, vnewc,
                            pmin, p_cut, eosvmax, length);

    /**
     Time taken(in microseconds) for sequential execution:  1227.15

    Problem Size = 32, Number of Threads = 20
    Time taken(in microseconds)for the execution of this loop for different values of 'Factor'

        Factor    Time(microsec)
        1          150.2665
        2          120.6625
        4          119.9535
        8          132.8875
        16         158.3035
        32         318.6925
    **/
    HCLIB_FINISH {
        int Factor = 1;
        hclib::loop_domain_t loop = {0, length,1,  int(length/(CONST_FACTOR*Factor))};
        hclib::irregular_recursion1D(&loop, e_new, start_i, end_i, [=](int i){
            double sixth = 1.0 / 6.0;
            double q_tilde;

            if (delvc[i] > 0.0) {
                q_tilde = 0.0;
            } else {
                double ssc = ((pbvc[i] * e_new[i]) + (vnewc[i] * vnewc[i] * bvc[i] * p_new[i])) / rho0;

                if (_compare(ssc, 0.111111e-36) <= 0) {
                    ssc = 0.333333e-18;
                } else {
                    ssc = std::sqrt(ssc);
                }

                q_tilde = (ssc * ql[i] + qq[i]);
            }

            e_new[i] = e_new[i] - ((((7.0 * (p_old[i] + q_old[i])) -
                    (8.0 * (pHalfStep[i] + q_new[i]))) + (p_new[i] + q_tilde)) * delvc[i] * sixth);

            if (_compare(abs(e_new[i]), e_cut) < 0) {
                e_new[i] = 0.0;
            }
            if (e_new[i] < emin) {
                e_new[i] = emin;
            }
        }, FORASYNC_MODE_RECURSIVE);
    }

    calcPressureForElems(p_new, bvc, pbvc, e_new, compression, vnewc, pmin, p_cut, eosvmax, length);

    /**
     Time taken(in microseconds) for sequential execution:  662.134

    Problem Size = 32, Number of Threads = 20
    Time taken(in microseconds)for the execution of this loop for different values of 'Factor'

        Factor    Time(microsec)
        1          92.33645
        2          85.32685
        4          84.55035
        8          98.84255
        16         161.6005
        32         320.008
    **/
    HCLIB_FINISH {
        int Factor = 1;
        hclib::loop_domain_t loop = {0, length,1,  int(length/(CONST_FACTOR*Factor))};
        hclib::irregular_recursion1D(&loop, e_new, start_i, end_i, [=](int i){
            if (delvc[i] <= 0.0) {
                double ssc = ((pbvc[i] * e_new[i]) + (vnewc[i] * vnewc[i] * bvc[i] * p_new[i])) / rho0;

                if (_compare(ssc, 0.111111e-36) <= 0) {
                    ssc = 0.333333e-18;
                } else {
                    ssc = std::sqrt(ssc);
                }

                q_new[i] = (ssc * ql[i] + qq[i]);

                if (_compare(abs(q_new[i]), q_cut) < 0) {
                    q_new[i] = 0.0;
                }
            }
        }, FORASYNC_MODE_RECURSIVE);
    }

    return;
                                      
    
}

void LuleshParBenchmark::calcPressureForElems(
        double p_new[], double bvc[],
        double pbvc[], double e_old[],
        double compression[], double vnewc[],
        double pmin,
        double p_cut, double eosvmax,
        int length) {
    

    
    double c1s = 2.0 / 3.0;
    // Note: loop fusion

    /**
     Time taken(in microseconds) for sequential execution:  766.9885

    Problem Size = 32, Number of Threads = 20
    Time taken(in microseconds)for the execution of this loop for different values of 'Factor'

        Factor    Time(microsec)
        1          91.65615
        2          88.83725
        4          92.86305
        8          103.403
        16         147.6055
        32         317.426
    **/
    auto start_i = [=](int i) { return i; };
    auto end_i = [=](int i) { return i-1; };
    HCLIB_FINISH {
        int Factor = 8;
        hclib::loop_domain_t loop = {0, length,1,  int(length/(CONST_FACTOR*Factor))};
        hclib::irregular_recursion1D(&loop, p_new, start_i, end_i, [=](int i){
            bvc[i] = c1s * (compression[i] + 1.0);
            pbvc[i] = c1s;

            p_new[i] = bvc[i] * e_old[i];

            if (_compare(abs(p_new[i]), p_cut) < 0) {
                p_new[i] = 0.0;
            }

            if (vnewc[i] >= eosvmax) /* impossible condition here? */ {
                p_new[i] = 0.0;
            }

            if (p_new[i] < pmin) {
                p_new[i] = pmin;
            }
        }, FORASYNC_MODE_RECURSIVE);
    }
}

void LuleshParBenchmark::calcSoundSpeedForElems(
        Domain& domain, double vnewc[], double rho0, double enewc[],
        double pnewc[], double pbvc[], double bvc[], double ss4o3, int nz){
    
    /**
     Time taken(in microseconds) for sequential execution:  978.2935

    Problem Size = 32, Number of Threads = 20
    Time taken(in microseconds)for the execution of this loop for different values of 'Factor'

        Factor    Time(microsec)
        1          96.2756
        2          110.20245
        4          95.58215
        8          105.7615
        16         154.55
        32         316.106
    **/
    auto start_i = [=](int i) { return i; };
    auto end_i = [=](int i) { return i-1; };
    HCLIB_FINISH {
        int Factor = 1;
        hclib::loop_domain_t loop = {0, nz,1,  int(nz/(CONST_FACTOR*Factor))};
        hclib::irregular_recursion1D(&loop, pbvc, start_i, end_i, [=, &domain](int i){        
  
            int iz = domain.matElemlist(i);
            double ssTmp = ((pbvc[i] * enewc[i]) + (vnewc[i] * vnewc[i] * bvc[i] * pnewc[i])) / rho0;
            if (_compare(ssTmp, 0.111111e-36) <= 0) {
                ssTmp = 0.333333e-18;
            } else {
                ssTmp = std::sqrt(ssTmp);
            }
            domain.ss(iz, ssTmp);
        }, FORASYNC_MODE_RECURSIVE);
    }
}


void LuleshParBenchmark::updateVolumesForElems(Domain& domain){
    int numElem = domain.numElem();
    if (numElem != 0) {
        double v_cut = domain.v_cut();

        /**
         Time taken(in microseconds) for sequential execution:  505.83

        Problem Size = 32, Number of Threads = 20
        Time taken(in microseconds)for the execution of this loop for different values of 'Factor'

            Factor    Time(microsec)
            1          59.5184
            2          61.833
            4          64.22145
            8          171.54755
            16         159.543
            32         311.731
        **/
	auto start_i = [=](int i) { return i; };
	auto end_i = [=](int i) { return i-1; };
        HCLIB_FINISH {
            int Factor = 1;
            hclib::loop_domain_t loop = {0, numElem,1,  int(numElem/(CONST_FACTOR*Factor))};
            hclib::irregular_recursion1D(& loop, domain.get_m_v(), start_i, end_i, [=, &domain](int i){
                double tmpV = domain.vnew(i);

                if (abs(tmpV - 1.0) < v_cut) {
                    tmpV = 1.0;
                }
                domain.v(i, tmpV);
            }, FORASYNC_MODE_RECURSIVE);
        }
    }

    return;

}

void LuleshParBenchmark::calcTimeConstraintsForElems(Domain& domain){
    /* evaluate time constraint */
    calcCourantConstraintForElems(domain);

    /* check hydro constraint */
    calcHydroConstraintForElems(domain);
}

void LuleshParBenchmark::calcCourantConstraintForElems(Domain& domain){
    double dtCourant = 1.0e+20;
    int courantElem = -1;
    double qqc = domain.qqc();
    int length = domain.numElem();

    double qqc2 = 64.0 * qqc * qqc;

    int threads = 1;

    int courantElemPerThread[threads];
    double dtCourantPerThread[threads];

    for (int i = 0; i < threads; i++) {
        courantElemPerThread[i] = -1;
        dtCourantPerThread[i] = 1.0e+20;
    }

    for (int i = 0; i < length; ++i) {
        int indx = domain.matElemlist(i);

        double dtf = domain.ss(indx) * domain.ss(indx);

        if (_compare(domain.vdov(indx), 0.0) < 0) {

            dtf = dtf + (qqc2 * domain.arealg(indx) * domain.arealg(indx) * domain.vdov(indx) * domain.vdov(indx));
        }

        dtf = std::sqrt(dtf);

        dtf = domain.arealg(indx) / dtf;

        /* determine minimum timestep with its corresponding elem */
        if (_compare(domain.vdov(indx), 0) != 0) {
            int thread_num = 0;
            if (dtf < dtCourantPerThread[thread_num]) {
                dtCourantPerThread[thread_num] = dtf;
                courantElemPerThread[thread_num] = indx;
            }
        }
    }

    for (int i = 0; i < threads; i++) {
        if (dtCourantPerThread[i] < dtCourant) {
            dtCourant = dtCourantPerThread[i];
            courantElem = courantElemPerThread[i];
        }
    }

    /* Don't try to register a time constraint if none of the elements were active */
    if (courantElem != -1) {
        domain.dtcourant(dtCourant);
    }

    return;
}

void LuleshParBenchmark::calcHydroConstraintForElems(Domain& domain){
    double dtHydro = 1.0e+20;
    int hydroElem = -1;
    double dvovmax = domain.dvovmax();
    int length = domain.numElem();

    int threads = 1;

    double dtHydroPerThread[threads];
    int hydroElemPerThread[threads];

    for (int i = 0; i < threads; i++) {
        hydroElemPerThread[i] = hydroElem;
        dtHydroPerThread[i] = dtHydro;
    }

    for (int i = 0; i < length; ++i) {
        int indx = domain.matElemlist(i);

        if (_compare(domain.vdov(indx), 0) != 0) {
            double dtdvov = dvovmax / (abs(domain.vdov(indx)) + 1.e-20);

            int thread_num = 0;
            if (dtHydroPerThread[thread_num] > dtdvov) {
                dtHydroPerThread[thread_num] = dtdvov;
                hydroElemPerThread[thread_num] = indx;
            }
        }
    }

    for (int i = 0; i < threads; i++) {
        if (dtHydroPerThread[i] < dtHydro) {
            dtHydro = dtHydroPerThread[i];
            hydroElem = hydroElemPerThread[i];
        }
    }

    if (hydroElem != -1) {
        domain.dthydro(dtHydro);
    }

    return;
}

int LuleshParBenchmark::_compare(double a, double b, double epsilon)
{   
    if (abs(a - b) < epsilon){
        return 0;
    }
    if (a < b){
        return -1;
    }
    else{
        return 1;
    }
}

