#ifndef MAINLOOP_H
#define MAINLOOP_H

#include "General.h"
#include "Param.h"
#include "Arrays.h"
#include "Forcing.h"
#include "Mesh.h"
#include "Write_netcdf.h"
#include "InitialConditions.h"
#include "MemManagement.h"
#include "Boundary.h"
#include "FlowGPU.h"
#include "FlowCPU.h"
#include "Meanmax.h"
#include "Updateforcing.h"
#include "FlowMLGPU.h"

template <class T> void MainLoop(Param& XParam, Forcing<float> XForcing, Model<T>& XModel, Model<T>& XModel_g);

template <class T> void DebugLoop(Param& XParam, Forcing<float> XForcing, Model<T>& XModel, Model<T>& XModel_g);

template <class T> __host__ double initdt(Param XParam, Loop<T> XLoop, Model<T> XModel);

template <class T> Loop<T> InitLoop(Param& XParam, Model<T>& XModel);

template <class T> void printstatus(T totaltime, T dt);

// Volume balance of the run so far (Param massbalance). Volumes in m3.
struct VolBalance
{
	bool valid = false;
	bool closed = false; // all boundary segments are walls: the residual is the mass error
	double t = 0.0;
	double storage = 0.0; // water on the grid now
	double start = 0.0; // water on the grid at the first balance
	double rain = 0.0; // rain actually added by the model
	double rainseries = -1.0; // rain in the uniform rain series over the rain area (-1: rain maps or no rain)
	double rainarea = 0.0; // area receiving rain (m2)
	double river = 0.0; // river inflow
	double loss = 0.0; // initial + continuing losses (hgw) where losses act
	double residual = 0.0; // start + rain + river - loss - storage: out through open boundaries + mass error
};
template <class T> VolBalance VolumeBalance(Param& XParam, Loop<T>& XLoop, Forcing<float>& XForcing, Model<T>& XModel, Model<T>& XModel_g);
extern VolBalance LastVolBalance; // the latest balance written (read by test 21)


template <class T> __global__ void storeTSout(Param XParam, int noutnodes, int outnode, int istep, int blknode, int inode, int jnode, int* blkTS, EvolvingP<T> XEv, T* store);


// End of global definition
#endif
