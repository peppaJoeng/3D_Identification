/*
 Created by kiki on 2020-02-15.
*/

#include "../thirdparty/clopdemo-09/src/mixture.hpp"


using namespace cp;

#include "iostream"
#include "mex.h"
#include <cstdlib>
#include "../thirdparty/clopdemo-09/src/nanooff.hpp"

using namespace std;
#define IN_X prhs[0]
#define IN_LEVEL prhs[1]
#define IN_DIR prhs[2]
#define IN_ALPHA prhs[3]
#define OUT_X plhs[0]

/*transform matrix matlab to c++ pointSet*/
template<typename T1, typename T2>
void m2c(T1 *x, T2 *y, int N, int D) {
    for (int row = 0; row < N; ++row) {
        for (int col = 0; col < D; ++col) {
            *y++ = (T2) x[row + col * N];
        }
    }
}

/*transform c++ pointSet to matrix matlab*/
template<typename T1, typename T2>
void c2m(T1 *x, T2 *y, int N, int D) {
    for (int col = 0; col < D; ++col) {
        for (int row = 0; row < N; ++row) {
            *y++ = (T2) x[row * D + col];
        }
    }
}

PointSet *
HEM(const double *x, const int N, const int D, const int &nLevels, const string &dir, float alpha) {
    string filename_O = dir + "_pointcloud_O.off";
    string filename_H = dir + "_pointcloud_H.off";
    cout << "Now we begin to cluster !!! " << endl;

    /*define PointSet and change double array to PointSet*/
    auto *points = new PointSet((uint) N);
    m2c(x, (float *) points->data(), N, D);
    cout<<"points size : "<<points->size()<<endl;
    
    nanooff::savePointCloud(filename_O, (float *) points->data(), points->size());

    // compute bounding box (with epsilon space border)
    vec3 mBBmin = vec3(FLT_MAX, FLT_MAX, FLT_MAX);
    vec3 mBBmax = vec3(-FLT_MAX, -FLT_MAX, -FLT_MAX);
    for (vector<vec3>::const_iterator p = points->begin(); p != points->end(); ++p) {
        const vec3 pt = *p;
        mBBmin = min(mBBmin, pt);
        mBBmax = max(mBBmax, pt);
    }
    float radius = sqrt(length(mBBmax - mBBmin) / points->size());
    
    /* hierarchical gaussian mixture model parameters*/
    Mixture::Params mParams;
    
    mParams.globalInitRadius = radius * alpha;
    cout << "mParams.globalInitRadius = " << mParams.globalInitRadius << endl;
    mParams.useGlobalInitRadius = true;  /* use global initialization radius instead of NNDist sampling (more stable for non-uniform point sampling)*/
    mParams.useWeightedPotentials = true;               /* if true, performs WLOP-like balancing of the initial Gaussian potentials*/
    mParams.alpha0 = 2.5f;                            /* multiple of nearest neighbor distance to use for initial query radius (<= 2.5f recommended)*/
    mParams.alpha = 2.3f;                               /* multiple of cluster maximum std deviation to use for query radius (<= 2.5f recommended)*/
    mParams.nLevels = (uint) nLevels;


    /*reduce the number of points*/
    auto Mix = new Mixture(points, mParams);
    auto *HGMMPoints = new PointSet(Mix->size());
    auto *tmp = (float *) HGMMPoints->data();
    for (auto g : *Mix) {
        *tmp++ = g.u.x;
        *tmp++ = g.u.y;
        *tmp++ = g.u.z;
    }
    cout << "After clustering, Mixture's number: " << Mix->size() << endl;
    nanooff::savePointCloud(filename_H, (float *) HGMMPoints->data(), HGMMPoints->size());
    return HGMMPoints;
}


void mexFunction(int nlhs, mxArray *plhs[], int nrhs, const mxArray *prhs[]) {

    double *X, *Y;
    int N, D;
    int nlevel;
    float globalInitRadius = 0.9f;
    float alpha = 4.0f;
    /* input*/
    N = (int) mxGetM(IN_X);
    D = (int) mxGetN(IN_X);
    X = mxGetPr(IN_X);
    nlevel = (int) mxGetScalar(IN_LEVEL);
    string dir(mxArrayToString(IN_DIR));
    alpha = (float) mxGetScalar(IN_ALPHA);
   
    cout << "piont cloud store in  " << dir.substr(0, dir.size() - 2) << endl;

    auto *points = HEM(X, N, D, nlevel, dir, alpha);

    OUT_X = mxCreateDoubleMatrix((mwSize) points->size(), (mwSize) D, mxREAL);
    Y = mxGetPr(OUT_X);
    c2m((float *) points->data(), Y, points->size(), D);
}

