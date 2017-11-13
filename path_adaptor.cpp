/* CHOMP class implementation
 *
 * Copyright (C) 2016 Rafael Valencia. All rights reserved.
 * License (3-Cluase BSD): https://github.com/rafaelvalencia
 * 
 * This code uses and is based on code from:
 *   Project: trychomp https://github.com/poftwaresatent/trychomp
 *   Copyright (C) 2014 Roland Philippsen. All rights reserved.
 *   License (3-Clause BSD) : https://github.com/poftwaresatent/trychomp
 * **
 * \file chomp.cpp
 *
 * CHOMP for point vehicles (x,y) moving holonomously in the plane. It will
 * plan a trajectory (xi) connecting start point (qs) to end point (qe) while
 * avoiding obstacles (obs)
 */
#include "path_adaptor.hpp"
#include <iostream>
#include <vector>
#include <Eigen/Dense>
#include <stdlib.h>
#include <sys/time.h>
#include <err.h>
#include <algorithm>

typedef Eigen::VectorXd Vector;
typedef Eigen::MatrixXd Matrix;
typedef Eigen::Isometry3d Transform;
static size_t const obs_dim(3); 	// obstacle dimensions (x,y,Radius)

using namespace std;

#define INF 1E20
#define epsilon 2.0

static double *dt(double *f, int n) {
  double *d = new double[n];
  int *v = new int[n];
  double *z = new double[n+1];
  int k = 0;
  v[0] = 0;
  z[0] = -INF;
  z[1] = +INF;
  for (int q = 1; q <= n-1; q++) {        // Compute lower envelope
    double s  = ((f[q]+pow(q,2.0))-(f[v[k]]+pow(v[k],2.0)))/(2*q-2*v[k]);
    while (s <= z[k]) {
      k--;
      s  = ((f[q]+pow(q,2.0))-(f[v[k]]+pow(v[k],2.0)))/(2*q-2*v[k]);
    }
    k++;
    v[k] = q;
    z[k] = s;
    z[k+1] = +INF;
  }

  k = 0;
  for (int q = 0; q <= n-1; q++) {
    while (z[k+1] < q)
      k++;
    d[q] = pow(q-v[k],2.0) + f[v[k]];
  }

  delete [] v;
  delete [] z;
  return d;
}



Matrix ComputeDistanceTranform(Matrix cost_init){

	int width = cost_init.rows();
	int height = cost_init.cols();
	double *f = new double[max(width, height)];

	// transfrom over rows
	for (int x=0; x < width; x++){
		for (int y=0; y<height; y++){
			f[y] = cost_init(x,y);
		}
		double *d = dt(f,height);
		for (int y=0; y< height; y++){
			cost_init(x,y) = d[y];
		}
		delete [] d;
	}

	//cout << cost_init << endl;

	// transfrom over columns
	for (int y=0; y < height; y++){
		for (int x=0; x < width; x++ ){
			f[x] = cost_init(x,y);
		}
		double *d = dt(f, width);
		for (int x=0; x < width; x++){
			cost_init(x,y) = d[x];
		}
		delete [] d;
	}

	delete f;

	//cout << cost_init << endl;
	return cost_init;

}


CHOMP::CHOMP(double dt_input, double eta_input, double lambda_input, size_t nq_input, size_t cdim_input, size_t numIt_input, double gain, double gamma_input) 
{
	//Sets basic parameters.	
	nq_ = nq_input;			// number of poses q in xi
	cdim_ = cdim_input;		// dimension of config space
	xidim_ = nq_ * cdim_; 	// dimension of trajectory, xidim = nq * cdim
	dt_ =  dt_input;	    // time step
	eta_ = eta_input; 		// >= 1, regularization factor for gradient descent
	lambda_ = lambda_input; // weight of smoothness objective	
	numIt_  = numIt_input; 	// Number of iterations
	costGain_ = gain;		// Gain inside cost function (usually 10) 
	gamma_ = gamma_input;
	
	res_ = 0.0001;			// Residual from optimization
	cter_ = 0;				// Zero iterations so far
	
	PATH_INIT_ = false;     // path is not initialized yet
	
	OBS_ = Matrix::Zero (obs_dim, 1); //initialize obstacle matrix
	
	cout << "-----------------------------------------------------:"<< endl;  		  	
	cout << "CHOMP has the following parameters:"<< endl;  
	cout << "-----------------------------------------------------:"<< endl;  	  	
	cout << "Time step: " << dt_ << endl;   
	cout << "Eta: " << eta_ << endl;  
	cout << "Lambda: " << lambda_ << endl;  	
	cout << "Number of poses in xi: " << nq_ << endl;
	cout << "Dimensions of config. space: " << cdim_ << endl;	
	cout << "Number of iterations: " << numIt_ << endl;
	cout << "Cost function gain: " << costGain_ << endl;
	cout << "-----------------------------------------------------:"<< endl;  
		  	
}

void CHOMP::makeGrid(double xmin, double ymin, double xmax, double ymax)
{
	int xdim = (int) xmax - xmin;
	int ydim = (int) ymax - ymin;
	grid_ = Matrix::Zero (xdim, ydim);
}

double CHOMP::chompIteration(Vector  &xi, Vector &ti)
{  
	
	// Before performing the iteration check if a path has been given
	
	if (PATH_INIT_==false)
	{
		cout << "A path was not initialized. Leaving CHOMP iteration! " << endl;
		return NAN;	
	}
	
	//////////////////////////////////////////////////
	// beginning of "the" CHOMP iteration
	
	Vector nabla_smooth (AA_ * xi_ + bb_);
	Vector nabla_smooth_time(gamma_* (BB_ * ti_ + tt_));

	Vector const & xidd (- nabla_smooth); // indeed, it is the same in this formulation...
	Vector const & tidd ( - nabla_smooth_time/gamma_);  // t''

	Vector nabla_obs (Vector::Zero (xidim_));  // xidim_ : dimension of trajectory
	Vector nabla_obs_time (Vector::Zero (nq_));

	for (size_t iq (0); iq < nq_; ++iq) 
	{
		Vector const qq (xi_.block (iq * cdim_, 0, cdim_, 1));
		double currentTime = ti_(iq);

		Vector qd;
		double tprime;

		if (0 == iq) {
			//cout << "block:  " << xi_.block ((iq+1) * cdim_, 0, cdim_, 1) << endl;
		  qd = 0.5 * (xi_.block ((iq+1) * cdim_, 0, cdim_, 1) - qs_);
		  tprime =  1/(2*dt_) * (ti_(iq+1) - ts_); 
		}
		else if (iq == nq_ - 1) {
		  qd = 0.5 * (qe_ - xi_.block ((iq-1) * cdim_, 0, cdim_, 1));
		  tprime = 1/(2*dt_) * (te_ - ti_(iq-1));
		}
		else {
		  qd = 0.5 * (xi_.block ((iq+1) * cdim_, 0, cdim_, 1) - xi_.block ((iq-1) * cdim_, 0, cdim_, 1));
		  tprime = 1/(2*dt_) * (ti_(iq+1) - ti_(iq-1));
		}

		//cout << "qs_  " << qs_ << endl;
		//cout << "qe_  " << qe_ << endl;
		//cout << "xi_  "  << xi_ << endl;
		//cout << "qd_  " << qd << endl;
		//cout << "xidd_  " << xidd << endl;
		// In this case, C and W are the same, Jacobian is identity.  We
		// still write more or less the full-fledged CHOMP expressions
		// (but we only use one body point) to make subsequent extension
		// easier.
		//
		Vector const & xx (qq);          // currentTime equivalent to this
		Vector const & xd (qd);
		Matrix const JJ (Matrix::Identity (2, 2)); // a little silly here, as noted above.
		double const vel (xd.norm());         // tprime equivalent to this
		if (vel < 1.0e-3 || tprime < 1.0e-3) 
		{
			// avoid div by zero further down
			continue;
		}
		Vector const xdn (xd / vel);
		Vector const xdd (JJ * xidd.block (iq * cdim_, 0, cdim_ , 1));
		Matrix const prj (Matrix::Identity (2, 2) - xdn * xdn.transpose()); // hardcoded planar case
		Vector const kappa (prj * xdd / pow (vel, 2.0));
		Matrix delta = Matrix::Zero(2,1);
		double cost;

		//Add obstacles		 
		for (int ii = 0; ii < OBS_.cols(); ii++) 
		{
			
			/*
			Vector delta(xx - OBS_.block(0, ii, 2, 1));
			double const dist(delta.norm());
			if ((dist >= OBS_(2, ii)) || (dist < 1e-9))   // Maxdist = radius*2
				continue;
			//double const cost(costGain_ * OBS_(2, ii) * pow(1.0 - dist / OBS_(2, ii), 3.0) / 3.0); 
			double const cost(OBS_(2, ii) * pow(1.0 - dist / OBS_(2, ii), 3.0) / 3.0);   
			delta *= - costGain_ *pow(1.0 - dist / OBS_(2, ii), 2.0) / dist;
			*/
			
			//cout << "xx: " << xx(0) << " " << xx(1) << endl;

			int a = (int)xx(0);
			int b = (int)xx(1);

			// cout << "a " << a <<  "  "  << "b " << b << endl;

			double distanceField = Dx_(a, b);

			delta(0,0) =  Dx_(a+1,b) - Dx_(a-1,b);
			delta(1,0) =  Dx_(a,b+1) - Dx_(a,b-1) ;  

			//cout << distanceField << endl;
			
			if (distanceField > epsilon)
				continue;

			
			if (distanceField < 0)
			{
				 cost = -costGain_*distanceField + 0.5*epsilon;
				 delta = delta*-1;
			}

			else if (distanceField <= epsilon)
			{
				cost = 0.5*epsilon*pow(distanceField - epsilon,2);
				delta = delta*(Dx_(a,b) - epsilon)/epsilon;
			}  
			   

			

			//cout << "JJ transpose  " << JJ.transpose() << endl;
			//cout << "velocity " << vel << endl;
			//cout << "projection matrix " << prj << endl;
			//cout << "cost " << cost << endl;
			//cout << "kappa " << kappa << endl;
			//cout << "delta  " << delta << endl;


			nabla_obs.block(iq * cdim_, 0, cdim_, 1) += JJ.transpose() * vel * (prj * delta - cost * kappa);
		} 
		
	}

	Vector dxi (Ainv_ * (nabla_obs + lambda_ * nabla_smooth));

 	//cout << "change " << dxi / eta_ << endl;

	xi_ -= dxi / eta_;
	
	xi = xi_; //updated path
	
	res_ = dxi.norm() / eta_;
	
	return res_;
	// end of "the" CHOMP iteration
	//////////////////////////////////////////////////
}


double CHOMP::chompUpdate(Vector  &xi, double &U_cost, double &curv)
{  
	
	// Before performing the iteration check if a path has been given
	
	if (PATH_INIT_==false)
	{
		cout << "A path was not initialized. Leaving CHOMP iteration! " << endl;
		return NAN;	
	}
	
	//////////////////////////////////////////////////
	// beginning of "the" CHOMP iteration
	
	double U(NAN); //cost functional value
	double F_obs(0);//obstacle functional value
	double F_smooth(0);//obstacle functional value
	double curvature(0); //curvature	

	Vector nabla_smooth (AA_ * xi_ + bb_);
	Vector const & xidd (nabla_smooth); // indeed, it is the same in this formulation...

	Vector nabla_obs (Vector::Zero (xidim_));
	for (size_t iq (0); iq < nq_; ++iq) 
	{
		Vector const qq (xi_.block (iq * cdim_, 0, cdim_, 1));
		Vector qd;
		if (0 == iq) {
		  qd = 0.5 * (xi_.block ((iq+1) * cdim_, 0, cdim_, 1) - qs_);
		}
		else if (iq == nq_ - 1) {
		  qd = 0.5 * (qe_ - xi_.block ((iq-1) * cdim_, 0, cdim_, 1));
		}
		else {
		  qd = 0.5 * (xi_.block ((iq+1) * cdim_, 0, cdim_, 1) - xi_.block ((iq-1) * cdim_, 0, cdim_, 1));;
		}

		// In this case, C and W are the same, Jacobian is identity.  We
		// still write more or less the full-fledged CHOMP expressions
		// (but we only use one body point) to make subsequent extension
		// easier.
		//
		Vector const & xx (qq);
		Vector const & xd (qd);
		Matrix const JJ (Matrix::Identity (2, 2)); // a little silly here, as noted above.
		double const vel (xd.norm());
		if (vel < 1.0e-3) 
		{
			// avoid div by zero further down
			continue;
		}
		Vector const xdn (xd / vel);
		Vector const xdd (JJ * xidd.block (iq * cdim_, 0, cdim_ , 1));
		Matrix const prj (Matrix::Identity (2, 2) - xdn * xdn.transpose()); // hardcoded planar case
		Vector const kappa (prj * xdd / pow (vel, 2.0));
		
		//curvature
		curvature += kappa.norm();
		
		double acc_cost(0); //sum of cost function values from all obstacles at a given robot pose
				
		//Add obstacles		 
		for (int ii = 0; ii < OBS_.cols(); ii++) 
		{
			Vector delta(xx - OBS_.block(0, ii, 2, 1));
			double const dist(delta.norm());
			if ((dist >= OBS_(2, ii)) || (dist < 1e-9))
				continue;
			double const cost(costGain_ * OBS_(2, ii) * pow(1.0 - dist / OBS_(2, ii), 3.0) / 3.0);  
			delta *= - costGain_ *pow(1.0 - dist / OBS_(2, ii), 2.0) / dist;                        
			nabla_obs.block(iq * cdim_, 0, cdim_, 1) += JJ.transpose() * vel * (prj * delta - cost * kappa);
			acc_cost += cost;
		} 
		
		//smoothness and obstacle costs
		F_smooth +=    pow(vel / dt_, 2.0); 
		F_obs +=    acc_cost * (vel / dt_); 		
	}
	
	//compute cost functional (from smoothness and obstacle costs)
	U =  (F_obs + 0.5 * lambda_ * F_smooth ) / (nq_ + 1);
    U_cost = U;

    //normalize curvature to size
    curv =   curvature / (nq_ + 1);
    
	Vector dxi (Ainv_ * (nabla_obs + lambda_ * nabla_smooth));
	xi_ -= dxi / eta_;
	
	xi = xi_; //updated path
	
	res_ = dxi.norm() / eta_;
	
	return res_;
	// end of "the" CHOMP iteration
	//////////////////////////////////////////////////
}


double CHOMP::chompUpdateWithSearchRegion(Vector  &xi, double theta, double &U_cost, double &curv)
{  
	
	// Before performing the iteration check if a path has been given
	
	if (PATH_INIT_==false)
	{
		cout << "A path was not initialized. Leaving CHOMP iteration! " << endl;
		return NAN;	
	}
	
	//////////////////////////////////////////////////
	//Compute transformation between global and starting vehicle's reference frames (qs_). 
	Eigen::Matrix4d T_global2local, T_local2global;
	T_local2global   <<   cos(theta), -sin(theta), 	0, 	qs_(0),
						  sin(theta),  cos(theta), 	0, 	qs_(1),
						  0, 					0, 	1, 	0,
						  0,					0,	0,	1;
	T_global2local = 	T_local2global.inverse();			  
	
	//Slope for lines defining search region	
	double m  = 1.; // lines at qs_ with 45 deg of slope 
		
	//////////////////////////////////////////////////
	// beginning of "the" CHOMP iteration
	
	double U(NAN); //cost functional value
	double F_obs(0);//obstacle functional value
	double F_smooth(0);//obstacle functional value
	double curvature(0); //curvature

	Vector nabla_smooth (AA_ * xi_ + bb_);
	Vector const & xidd (nabla_smooth); // indeed, it is the same in this formulation...

	Vector nabla_obs (Vector::Zero (xidim_));
	for (size_t iq (0); iq < nq_; ++iq) 
	{
		Vector const qq (xi_.block (iq * cdim_, 0, cdim_, 1));
		Vector qd;
		if (0 == iq) {
		  qd = 0.5 * (xi_.block ((iq+1) * cdim_, 0, cdim_, 1) - qs_);
		}
		else if (iq == nq_ - 1) {
		  qd = 0.5 * (qe_ - xi_.block ((iq-1) * cdim_, 0, cdim_, 1));
		}
		else {
		  qd = 0.5 * (xi_.block ((iq+1) * cdim_, 0, cdim_, 1) - xi_.block ((iq-1) * cdim_, 0, cdim_, 1));;
		}

		// In this case, C and W are the same, Jacobian is identity.  We
		// still write more or less the full-fledged CHOMP expressions
		// (but we only use one body point) to make subsequent extension
		// easier.
		//
		Vector const & xx (qq);
		Vector const & xd (qd);
		Matrix const JJ (Matrix::Identity (2, 2)); // a little silly here, as noted above.
		double const vel (xd.norm());
		if (vel < 1.0e-3) 
		{
			// avoid div by zero further down
			continue;
		}
		Vector const xdn (xd / vel);
		Vector const xdd (JJ * xidd.block (iq * cdim_, 0, cdim_ , 1));
		Matrix const prj (Matrix::Identity (2, 2) - xdn * xdn.transpose()); // hardcoded planar case
		Vector const kappa (prj * xdd / pow (vel, 2.0));
		double gain;
	
		//curvature
		curvature += kappa.norm();
		
		double acc_cost(0); //sum of cost function values from all obstacles at a given robot pose
		//Consider obstacles		 
		for (int ii = 0; ii < OBS_.cols(); ii++) 
		{
			//Check search region, if it is outside the search region set a high obs. function gain
			Eigen::Vector4d xx_l, xx_g;
			xx_g << xx(0), xx(1), 0, 1;
			
			xx_l =  T_global2local *  xx_g;
			
			//allow only poses infront of initial position (avoids reverse motion paths)  
			//if( ( xx_l(0) > 0 ) ) 
			//no reverse motion paths & inside a region defined by 2 lines centered at qs_
			if( ( xx_l(0) > 0 ) && ( xx_l(1) < m * xx_l(0) ) && (  xx_l(1) > - m * xx_l(0) ) )
			{
				//Inside desired search region
				gain = costGain_;
			}
			else
			{
				//Outside desired search region
				gain = 1000 * costGain_;	
				cout << " Outside desired search region!!! " << endl;			
			}
			
			Vector delta(xx - OBS_.block(0, ii, cdim_, 1));
			double const dist(delta.norm());
			if ((dist >= OBS_(2, ii)) || (dist < 1e-9))
				continue;
			double const cost(gain * OBS_(2, ii) * pow(1.0 - dist / OBS_(2, ii), 3.0) / 3.0);  
			delta *= - gain *pow(1.0 - dist / OBS_(2, ii), 2.0) / dist;                        
			nabla_obs.block(iq * cdim_, 0, cdim_, 1) += JJ.transpose() * vel * (prj * delta - cost * kappa);
			acc_cost += cost;
		} 
		
		//smoothness and obstacle costs
		F_smooth +=    pow(vel / dt_, 2.0); 
		F_obs +=    acc_cost * (vel / dt_); 
		

	}

	//compute cost functional (from smoothness and obstacle costs)
	U =  (F_obs + 0.5 * lambda_ * F_smooth ) / (nq_ + 1);
    U_cost = U;
    
    //normalize curvature to size
    curv =   curvature / (nq_ + 1);
	
	Vector dxi (Ainv_ * (nabla_obs + lambda_ * nabla_smooth));
	xi_ -= dxi / eta_;
	
	xi = xi_; //updated path
	
	++cter_; // increase iteration counter	
	
	res_ = dxi.norm() / eta_;
	cout << "It No. " << cter_ << " Res.: " << dxi.norm() / eta_<< " Cost: " << U << " Curvature:" << curv << endl;
	
	//Detect local minimum	
	if(cter_ > 1000 && res_ > 0.1)
	{
		cout << "LOCAL MINIMUM DETECTED! " << endl;
		for (size_t ii (0); ii < xidim_; ++ii) 
		{
			 double noise = 0.1 * (- 50 + rand() % 100); //random noise between -5 to 5
			 xi_[ii] = xi_[ii] + noise;
			 cter_ = 0;
			 
		}	
	}
	else if (cter_ > 1000 && res_ < 0.01)
	{
		cout << "CONVERGENCE, RESETING COUNTER! " << endl;
		cter_ = 0;
	}
		
	
	return res_;
	
	// end of "the" CHOMP iteration
	//////////////////////////////////////////////////
}

void CHOMP::generatePath(Vector  &xi, Vector &ti)
{
	double err;
	for (size_t ii(0); ii < numIt_; ++ii)  
	{
		err = CHOMP::chompIteration(xi, ti);
		cout << "err " << err  << "ii " << ii << endl;
		if (err < 0.01)
		{
			//it converged
			cter_ = 0;
			break;
		}
	}	
	
}

//optimize path using with numIt_ iterations with serch region.	
void CHOMP::generatePathWithSearchReg(Vector  &xi, double orientation, double &U_cost, double &curv)	
{
	double err;
	for (size_t ii(0); ii < numIt_; ++ii)  
	{
		err = CHOMP::chompUpdateWithSearchRegion(xi, orientation,U_cost,curv);
		if (err < 0.01)
		{
			//it converged
			cter_ = 0;
			break;
		}
	}	
	
}

void CHOMP::addObstacle(double px, double py, double radius)
{
	OBS_.conservativeResize(obs_dim, OBS_.cols() + 1);
	OBS_.block(0, OBS_.cols() - 1, obs_dim, 1) << px, py, radius;
}

void CHOMP::setObstacles(Matrix obs)
{
	OBS_.resize (obs_dim, obs.cols()); 
	OBS_ = obs; 

	for (size_t ii(0); ii < OBS_.cols(); ++ii)
	{
		int ox = (int) OBS_(0,ii);
		int oy = (int) OBS_(1,ii);
		int ro = (int) OBS_(2,ii);

		//cout << "ox: " << ox << " " << "oy: " << oy << " " << "ro: " << ro << endl; 

		for (int row = (ox - ro); row <= ox + ro; row++)
		{
			for (int col = (oy - ro); col <= oy + ro; col++)
			{
				grid_(row, col) = 1;
			}
		}
	}

	//cout << "grid_  " << grid_ << endl;

	Matrix cost_init = Matrix::Zero(grid_.rows(), grid_.cols());
	for (int i=0; i< grid_.rows(); i++)
	{
		for (int j=0; j < grid_.cols(); j++)
		{
			if (!grid_(i,j))
				cost_init(i,j) = +INF;
		}
	}

	//cout << "cost_init  " << cost_init << endl;
	cost_init = ComputeDistanceTranform(cost_init);

			// Taking square root
	for (int x=0; x < grid_.rows(); x++)
	{
		for (int y=0; y < grid_.cols(); y++)
		{
				cost_init(x,y) = pow(cost_init(x,y), 0.5);
		}
	}

	Matrix posd_x = cost_init;
	//cout << posd_x << endl;

	// EDT of obstacle field Complement 
	for (int i=0; i< grid_.rows(); i++)
	{
		for (int j=0; j < grid_.cols(); j++)
		{
			if (grid_(i,j))
				cost_init(i,j) = +INF;
			else
				cost_init(i,j) = 0;
		}
	}

	cost_init = ComputeDistanceTranform(cost_init);

		for (int x=0; x < grid_.rows(); x++)
	{
		for (int y=0; y < grid_.cols(); y++)
		{
				cost_init(x,y) = pow(cost_init(x,y), 0.5);
		}
	}

	Matrix negd_x = cost_init;
	//cout << negd_x << endl;

	Dx_ = posd_x - negd_x; // Discretized Signed Distance Field of obstacle
	//cout << Dx_ << endl;  

	//cout << "Final Euclidean Distance Cost " << cost_init << endl;
	//cout << "Current Obstacles: \n" << OBS_ << endl;
}

// Sets an aribitrary value for xi_, qs_, qe_		
void CHOMP::setPath(Vector  &qs, Vector &qe, double ts, double te, Vector &xi)
{
	xi_ = xi;
	qs_ = qs;
	qe_ = qe;

	ts_ = ts;
	te_ = te;
	
	//Sets gradient descent vectors and matrices with the initialized path
	CHOMP::initCHOMP();		
}

//Sets qs_ and qe_ and initializes all points xi_ to qs_ (stacked to qs_)
void CHOMP::initStackedPath(Vector  &qs, Vector &qe, double ts, double te)
{
	qs_ = qs;
	qe_ = qe;

	ts_ = ts;
	te_ = te;

	xi_ = Vector::Zero (xidim_);
	ti_ = Vector::Zero (nq_);
	for (size_t ii (0); ii < nq_; ++ii) 
	{
		xi_.block (cdim_ * ii, 0, cdim_, 1) = qs_ + (ii+1)*(qe_ - qs_)/(nq_ + 1);
		ti_(ii) = ts_ + (ii+1)*(te_ - ts_)/(nq_ + 1);
	}

	//Sets gradient descent vectors and matrices with the initialized path
	CHOMP::initCHOMP();		

}

//Initializes path xi as a straight line conecting the starting point qs_ and the ending point qe_
void CHOMP::initStraightLinePath(Vector  &qs, Vector &qe)
{
	//initalize a new trajectory based on a direct line connecting qs to qe
	qs_ = qs;
	qe_ = qe;
	
	xi_ = Vector::Zero(xidim_);
	Vector dxi(cdim_);
	dxi << (qe_(0) - qs_(0)) / (nq_ - 1), (qe_(0) - qs_(0)) / (nq_ - 1);
	for (size_t ii(0); ii < nq_; ++ii)
	{
		xi_.block(cdim_ * ii, 0, cdim_, 1) = qs_ + ii * dxi;
	}
	
	//Sets gradient descent vectors and matrices with the initialized path
	CHOMP::initCHOMP();			
}

// Gets xi_ 
void CHOMP::getPath(Vector &xi, Vector &ti)
{
	xi = xi_;
	ti = ti_;
}


void CHOMP::initCHOMP(void)
{

	// Initializes gradient descent vectors and matrices
	AA_ = Matrix::Zero (xidim_, xidim_);

	// For timing of trajectory
	BB_ = Matrix::Zero(nq_, nq_);

	for (size_t ii(0); ii < nq_; ++ii) 
	{
		AA_.block (cdim_ * ii, cdim_ * ii, cdim_ , cdim_) = 2.0 * Matrix::Identity (cdim_, cdim_);
		BB_(ii,ii) = 2.0;
		if (ii > 0) 
		{
			AA_.block (cdim_ * (ii-1), cdim_ * ii, cdim_ , cdim_) = -1.0 * Matrix::Identity (cdim_, cdim_);
			AA_.block (cdim_ * ii, cdim_ * (ii-1), cdim_ , cdim_) = -1.0 * Matrix::Identity (cdim_, cdim_);
			BB_(ii-1, ii) = -1.0;
			BB_(ii, ii-1) = -1.0;
		}
	}

	AA_ /= dt_ * dt_ * (nq_ + 1);
	BB_ /= dt_ * dt_ * (nq_ + 1);

	bb_ = Vector::Zero (xidim_);
	bb_.block (0,            0, cdim_, 1) = qs_;
	bb_.block (xidim_ - cdim_, 0, cdim_, 1) = qe_;
	bb_ /= - dt_ * dt_ * (nq_ + 1);

	Ainv_ = AA_.inverse();


	tt_ = Vector::Zero (nq_);
	tt_(0) = ts_;  
	tt_(nq_-1) = te_;
	tt_ /= -dt_ * dt_ * (nq_ + 1);

	Binv_ = BB_.inverse();
	
	PATH_INIT_= true;
}
