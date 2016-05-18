#ifndef _NUMEROV_
#define _NUMEROV_
#include <iostream>
#include <fstream>
#include <stdlib.h>
#include <math.h>
#include <vector>
#include <string>
#include <functional>
#include "mpreal.h"

#define DEBUG true
#define DEBUG_STDERR(x) (std::cerr << "[ERROR][numerov] " << (x) << endl)
#define DEBUG_STDOUT(x) (std::cout << "[DEBUG][numerov] " << (x) << endl)

using namespace std;
using mpfr::mpreal;

struct Point {
    mpreal x = 0;
    mpreal y = 0;
};

struct Point2 {
    double x = 0;
    double y = 0;
};

vector<Point2> wrap2(vector<Point> v) {
    vector<Point2> vWrap;
    vWrap.resize(v.size());
    for(int i = 0; i < v.size(); i++){
        vWrap[i].x = v[i].x.toDouble();
        vWrap[i].y = v[i].y.toDouble();
    }
    
    return vWrap;
}

vector<double> wrap(vector<mpreal> v) {
    vector<double> vWrap;
    vWrap.resize(v.size());
    for(int i = 0; i < v.size(); i++){
        vWrap[i] = v[i].toDouble();
    }
    
    return vWrap;
}

struct Range {
	mpreal eMin = 0;
	mpreal eMax = 1;
};

struct Range2 {
    double eMin = 0;
	double eMax = 1;
};

struct Mode {
	mpreal energy;
	vector<Point> wavefunction;
	Mode(mpreal e, vector<Point> f){
		energy = e;
		wavefunction = f;
	}
};

struct Spectrum {
	vector<Mode>  modes;
	vector<Point> potential;

	void addMode(Mode m) {
	   	modes.push_back(m); 
	}
	
	void clear() {
	   	modes.clear();
		potential.clear();
   	}
       
    vector<mpreal> getEnergies() {
        vector<mpreal> energies;
        for(int i = 0; i < modes.size(); i++) 
            energies.push_back(modes[i].energy);
        
        return energies;
    }
    
    vector<double> getEnergies2() {
        return wrap(getEnergies());
    }
    
    vector<vector<Point>> getWavefunctions() {
        vector<vector<Point>> wfs;
        for(int i = 0; i < modes.size(); i++) 
            wfs.push_back(modes[i].wavefunction);
        
        return wfs;
    }
    
    vector<vector<Point2>> getWavefunctions2() {
        vector<vector<Point2>> wfs;
        for(int i = 0; i < modes.size(); i++) 
            wfs.push_back(wrap2(modes[i].wavefunction));
        
        return wfs;
    }
    
};


class Numerov {
	private:
		vector<Point> potential;
		vector<Range> zeros;
		vector<Point> solLR;
		vector<Point> solRL;
		vector<Point> sol;
		Point minPot;
		Spectrum spectrum;
		std::function<mpreal(mpreal)> potFunc;
		function<mpreal(mpreal)> diffFunc = [this](mpreal E){
			return diff(E);
		};

		mpreal bisection(function<mpreal(mpreal)> diffFunc, mpreal min, mpreal max);
		void buildSol() ;
		mpreal diff(mpreal E);
		void displaySol(vector<Point> sol);
		void findMinPot();
		//void plotXY(std::vector<mpreal> x, std::vector<mpreal> y);
		void savePotential() ;
		void scanForZeroRanges(int nZeros) ;
		void wait_for_key();
	 	mpreal zbrent(function<mpreal(mpreal)>& func, mpreal x1, mpreal x2, mpreal tol);
	
	public:
		// these values can be overrided later
		mpreal h     = "0.001";
		mpreal xMin  = "-15";
		mpreal xMax  = "15";
		mpreal tol   = "1e-9";
		mpreal dEmin = "0.1";
		int nPoints  = 2000;
        void set_h(double dh);
        void set_tol(double dtol);
        void set_dEmin(double ddEmin);

		void setPotential(vector<Point>);
		void setPotential2(vector<Point2>);
		vector<Point> getPotential();
		void dumpEnergies();
		void dumpEnergies2();
		vector<Point2> getPotential2();
		void findSpectrum(int nEigen);
		Spectrum getSpectrum();

};

void Numerov::set_h(double dh) {
    h = mpreal(h);
}

void Numerov::set_tol(double dtol) {
    tol = mpreal(tol);
}

void Numerov::set_dEmin(double ddEmin) {
    dEmin = mpreal(dEmin);
}

void Numerov::dumpEnergies (){
	cout << "ENERGIES:" << endl;
	vector<mpreal> energies = spectrum.getEnergies();
	for (auto e: energies) 
		cout << e << endl;
	cout << "--------" << endl;
}

void Numerov::dumpEnergies2 (){
	cout << "ENERGIES(using double):" << endl;
	vector<mpreal> energies = spectrum.getEnergies();
	for (auto e: wrap(energies)) 
		cout << e << endl;
	cout << "--------" << endl;
}

vector<Point> Numerov::getPotential() {
    return potential;
}

vector<Point2> Numerov::getPotential2() {
    return wrap2(getPotential());
}

Spectrum Numerov::getSpectrum() {
	return spectrum;
}

// Van Wijngaarden–Dekker–Brent Method for finding root, from NR
mpreal EPS = "1e-12";
#define ITMAX 800
#define SIGN(a,b) ((b) >= 0.0 ? fabs(a) : -fabs(a))
mpreal Numerov::zbrent(function<mpreal(mpreal)>& func, mpreal x1, mpreal x2, mpreal tol)
{
	int iter;
	mpreal a=x1,b=x2,c=x2,d,e,min1,min2;
	mpreal fa=func(a),fb=func(b),fc,p,q,r,s,tol1,xm;
	
	if ((fa > "0.0" && fb > "0.0") || (fa < "0.0" && fb < "0.0"))
		cout << "[WARN] Root must be bracketed in zbrent"<<endl;
	fc=fb;
	for (iter=1;iter<=ITMAX;iter++) {
		if ((fb > "0.0" && fc > "0.0") || (fb < "0.0" && fc < "0.0")) {
			c  = a;			//Rename a, b, c and adjust bounding interval d.
			fc = fa;
			e  = d = b-a;
		}
		if (fabs(fc) < fabs(fb)) {
			a  = b;
			b  = c;
			c  = a;
			fa = fb;
			fb = fc;
			fc = fa;
		}

		tol1 = "2.0"*EPS*fabs(b)+"0.5"*tol; //Convergence check.
		xm = "0.5"*(c-b);

		if (fabs(xm) <= tol1 || fb == "0.0") return b;
		if (fabs(e) >= tol1 && fabs(fa) > fabs(fb)) {
			s=fb/fa;		//Attempt inverse quadratic interpolation.
			if (a == c) {
				p = "2.0"*xm*s;
				q = "1.0"-s;
			} else {
				q = fa/fc;
				r = fb/fc;
				p = s*("2.0"*xm*q*(q-r)-(b-a)*(r-"1.0"));
				q = (q-"1.0")*(r-"1.0")*(s-"1.0");
			}
			if (p > "0.0") q = -q;		//Check whether in bounds.
			p=fabs(p);
			min1="3.0"*xm*q-fabs(tol1*q);
			min2=fabs(e*q);
			if ("2.0"*p < (min1 < min2 ? min1 : min2)) {
				e = d;
				//Accept interpolation.
				d = p/q;
			} else {
				d = xm;
				//Interpolation failed, use bisection.
				e = d;
			}
		} else {		//Bounds decreasing too slowly, use bisection.
			d = xm;
			e = d;
		}
		a = b;			//Move last best guess to a.
		fa = fb;
		if (fabs(d) > tol1)	//Evaluate new trial root.
			b += d;
		else
			b += SIGN(tol1,xm);
		fb = func(b);
	}
	cout<<"[WARN] Maximum number of iterations exceeded in zbrent"<<endl;
	return "0.0";
}

// just returns the overall minima of the potential
void Numerov::findMinPot(){
    #ifdef DEBUG
        DEBUG_STDOUT("Finding the minimum of the potential");
    #endif
	minPot = potential[0];

	for (int i = 0; i < potential.size(); i++) 
		if (minPot.y > potential[i].y) 
			minPot = potential[i];

}

void Numerov::wait_for_key()
{
#if defined(WIN32) || defined(_WIN32) || defined(__WIN32__) || defined(__TOS_WIN__)  // every keypress registered, also arrow keys
    cout << endl << "Press any key to continue..." << endl;

    FlushConsoleInputBuffer(GetStdHandle(STD_INPUT_HANDLE));
    _getch();
#elif defined(unix) || defined(__unix) || defined(__unix__) || defined(__APPLE__)
    cout << endl << "Press ENTER to continue..." << endl;

    std::cin.clear();
    std::cin.ignore(std::cin.rdbuf()->in_avail());
    std::cin.get();
#endif
    return;
}

void Numerov::setPotential(vector<Point> v = vector<Point>(1)){
    #ifdef DEBUG
        DEBUG_STDOUT("Setting potential");
    #endif
	// the potential can be read at this point from a file or similar, or can be hard coded.
    if(v.size() > 1) {
        potential = v;
    } else { // show an harmonic oscillator data by default
    	potential.resize(nPoints);
    	for (int i = 0; i < nPoints; i++) {
    		potential[i].x = xMin + i * (xMax - xMin) / nPoints;
    		potential[i].y = potential[i].x * potential[i].x;  // harmonic oscillator for testing
    	}
        cout << "[DEBUG] Using harmonic oscillator potential..." << endl;
    }
    
    nPoints = potential.size();
    xMin    = potential.front().x;
    xMax    = potential.back().x;
	h = (xMax - xMin) / nPoints;

	// here we use the nice feature of c++ 11, the lambda functions
	// oh dear, this is soo cool!
	potFunc = [this] (mpreal x){ 
		// here we have to interpolate the potential since all the testing points are not in the grid!
		// a linear interpolation is used
		// first get the closest left point
		int index = 0;
		for (int i = 0; i < potential.size(); i++){ 
			if (potential[i].x > x) {
				index = i - 1;
				break;
			}
			// if we arrived to the end then it because x is the end point.
			index = potential.size() - 2;
		}

		mpreal x0 = potential[index].x;
		mpreal x1 = potential[index + 1].x;
		mpreal y0 = potential[index].y;
		mpreal y1 = potential[index + 1].y;

		mpreal m = (y1 - y0) / (x1 - x0);
		return y0 + m * (x - x0);
	};
    
}

// just a wrapper for using double as type
void Numerov::setPotential2(vector<Point2> v){
    vector<Point> vWrap;
    vWrap.resize(v.size());
    for(int i = 0; i < v.size(); i++){
        vWrap[i].x = mpreal(v[i].x);
        vWrap[i].y = mpreal(v[i].y);
    }
    
    setPotential(vWrap);
}

void Numerov::displaySol(vector<Point> sol){
	for (int i = 0; i < sol.size(); ++i) {
		cout << i <<" "<< sol[i].x<< ", " << sol[i].y << endl;
	}
}

mpreal Numerov::diff(mpreal E){

	// first we find the matching point
	function<mpreal(mpreal)> shiftedPot = [this, &E](mpreal x) {
		//cout<<"E=" << E << " shiftedPot(" << x <<")=" << E - potFunc(x) << endl;
		return E - potFunc(x);
	};
	
	// since we want to find preferently the right turning point we
	// look from the minimum position to the right
	mpreal matchPoint = zbrent(shiftedPot, minPot.x, xMax, tol);
	int matchPointIndex = int((matchPoint - xMin) / h);
	// now we have to propagate left and right solutions until the matching point
	solLR.clear();
    solLR.resize(matchPointIndex + 1);
    solLR[0].x = xMin;
    solLR[0].y = "0";
    solLR[1].x = xMin + h;
    solLR[1].y = "0.0001";  // this number should not be important, but it has to be small to avoid overflow

	//cout << "matchPointIndex = " << matchPointIndex << " solLR.size() = " << solLR.size() << endl;
	mpreal h2 = h * h / "12";
    for(int i = 2; i < solLR.size(); i++)
    {
        mpreal x2  = xMin + (i - 2) * h; 
        mpreal x1  = xMin + (i - 1) * h;
        mpreal x   = xMin + i * h;     
        mpreal p1  = "2" - "10"  * h2 * shiftedPot(x1);
        mpreal p2  = "1" +  h2 * shiftedPot(x2);
        mpreal p3  = "1" +  h2 * shiftedPot(x);
        solLR[i].y = (p1 * solLR[i-1].y - p2 * solLR[i-2].y) / p3;
        solLR[i].x = x;
    }
 
	solRL.clear(); 
	solRL.resize(potential.size() - solLR.size());
    //fill the array, propagating the solRLution from the Right to the Left
    solRL[solRL.size() - 1].x = xMax;
    solRL[solRL.size() - 1].y = "0";
    solRL[solRL.size() - 2].x = xMax - h;
    solRL[solRL.size() - 2].y = "0.0001";

    for(int i = solRL.size() - 3; i > -1; i--)
    {
        mpreal x2  = matchPoint + (i + 2) * h;   
        mpreal x1  = matchPoint + (i + 1) * h;  
        mpreal x   = matchPoint + i * h;       
        mpreal p1  = "2" - "10" * h2 * shiftedPot(x1);
        mpreal p2  = "1" + h2 * shiftedPot(x2);
        mpreal p3  = "1" + h2 * shiftedPot(x);
        solRL[i].y = (p1 * solRL[i+1].y - p2*solRL[i+2].y) / p3;
        solRL[i].x = x;
    }

	// now we have to find the log derivative at the matching point
	mpreal logDer = ((25/12)*solLR[matchPointIndex].y-4*solLR[matchPointIndex-1].y+3*solLR[matchPointIndex-2].y-(4/3)* solLR[matchPointIndex-3].y+(1/4)* solLR[matchPointIndex-4].y)*solRL[0].y+solLR[matchPointIndex].y*((25/12)*solRL[0].y-4*solRL[1].y+3*solRL[2].y-(4/3)*solRL[3].y+(1/4) * solRL[4].y);
	
	return logDer;
}



/**
 * Starting from the potential minima with a step of dEmin
 * computes diff(E), if a change of sign appears it means
 * that there is a zero thus store the range in the zeros vector.
 * These ranges are later used to find the root in a more refined phase.
 */
void Numerov::scanForZeroRanges(int nZeros) {
    #ifdef DEBUG
        DEBUG_STDOUT("Scanning for zeros");
    #endif
    
	findMinPot();
	mpreal E = minPot.y;
	zeros.clear();
	mpreal lastDiff = diff(E), newDiff;
	while(zeros.size() < nZeros){
		newDiff = diff(E);
		if(newDiff * lastDiff < 0){
			Range range;
			range.eMin = E - dEmin;
			range.eMax = E + dEmin;
			zeros.push_back(range);
			//cout << "zero in [" << range.eMin << ", " << range.eMax << "]" << endl;
		}
		lastDiff = newDiff;
		E += dEmin;
	}
}

void Numerov::buildSol() {
    #ifdef DEBUG
        DEBUG_STDOUT("Building solution");
    #endif
	sol.clear();
	mpreal scale = solRL[0].y / solLR[solLR.size() - 1].y;
	for (int i = 0; i < solLR.size(); i++) {
		Point p;
		p.x = solLR[i].x;
		p.y = scale * solLR[i].y;
		sol.push_back(p);
	}
	for (int i = 0; i < solRL.size(); i++) {
		sol.push_back(solRL[i]);
	}

	// finally we have to normalize it
	mpreal c = 0;
	for (int i = 0; i < sol.size(); i++) 
		c += sol[i].y * sol[i].y;
	
	for (int i = 0; i < sol.size(); i++) 
		sol[i].y /= sqrt(c);
}

void Numerov::findSpectrum(int nEigen){
    #ifdef DEBUG
        DEBUG_STDOUT("Finding spectrum");
    #endif
    
    if(potential.size() < 2) {
        cout << "Please use the setPotential() function before using this one." << endl;
        return;
    }
    
	spectrum.clear();
	spectrum.potential = potential;
	scanForZeroRanges(nEigen);
	for (int i = 0; i < nEigen; i++) {
		Range range = zeros[i];
		mpreal E = zbrent(diffFunc, range.eMin, range.eMax, tol);
		buildSol();
		Mode mode(E, sol);
		spectrum.addMode(mode);
		//cout << "E[" << i << "] = " << E << endl;		
	}
}

void Numerov::savePotential() {
    ofstream f("potential.dat");
	for (mpreal x = xMin; x < xMax; x += 0.01) {
		f << x << " " << potFunc(x) << endl;
	}
	f.close();
}

mpreal Numerov::bisection(function<mpreal(mpreal)> diffFunc, mpreal min, mpreal max){
	mpreal E;
	int i=0;			     /* counter for iterations */
	do{
		i++;  
		E=(max+min)/2.0;		     /* divide energy range */
		if (diffFunc(max)*diffFunc(E)>0) 
			max=E;   /* the bisection algorithm */
		else 
			min=E;
	}while(fabs(diffFunc(E))>EPS);
	return E;
}

#endif
