// Author: Keitaro Yamashita
//
// This software is released under GNU General Public License v3.

/*

g++ -fPIC -O3 -c peakfinders.cpp
cd boost_python
g++ -O3 -c cheetah_ext.cpp -shared -fPIC  -I../ -I/oys/python/python-2.7.6/include/python2.7/ -I/oys/python/python-2.7.6/lib/python2.7/site-packages/PyUblas-2013.1-py2.7-linux-x86_64.egg/pyublas/include -I/oys/python/python-2.7.6/lib/python2.7/site-packages/numpy/core/include
g++ -shared cheetah_ext.o  -o cheetah_ext.so -L/oys/xtal/cctbx/build/lib -lboost_python ../peakfinders.o 

*/
#include <boost/python.hpp>
#include <boost/python/args.hpp>
#include <vector>
#include <cmath>
#include <pyublas/numpy.hpp>

#include "peakfinders.h"
#include "cheetahmodules.h"

namespace bp = boost::python;

struct CheetahSpot {
	float peak_maxintensity;             // Maximum intensity in peak
	float peak_totalintensity;   // Integrated intensity in peak
	float peak_sigma;                    // Signal-to-noise ratio of peak
	float peak_snr;                              // Signal-to-noise ratio of peak
	float peak_npix;                             // Number of pixels in peak
	float peak_com_x;                    // peak center of mass x (in assembled layout)
	float peak_com_y;                    // peak center of mass y (in assembled layout)
	float peak_com_r;  // peak center of mass r (in assembled layout)
	float peak_com_res;                  // REsolution of this peak
};

class CheetahSinglePanel {
public:
	void set_params(float ADCthresh, float MinSNR, long MinPixCount,
					long MaxPixCount, long LocalBGRadius, float MinPeakSeparation,
					float Dmin, float Dmax) {
		hitfinderADCthresh = ADCthresh;
		hitfinderMinSNR = MinSNR;
		hitfinderMinPixCount = MinPixCount;
		hitfinderMaxPixCount = MaxPixCount;
		hitfinderLocalBGRadius = LocalBGRadius;
		hitfinderMinPeakSeparation = MinPeakSeparation;
		hitfinderDmin = Dmin;
		hitfinderDmax = Dmax;
	}

	template<class T>
	boost::python::list run(int asic_nx, int asic_ny, 
							pyublas::numpy_vector<T> data,
							float beam_x, float beam_y, float wavelength,
							float distance, float pixel_size, int algorithm) {
        const long pix_nn = asic_nx * asic_ny;
        const long nasics_x = 1, nasics_y = 1;
		assert(algorithm == 6 || algorithm == 8);
		tPeakList peaklist;
		allocatePeakList(&peaklist, 500);

		// Setup resolution mask
		const double r_min = distance * std::tan(2.*std::asin(wavelength/2./hitfinderDmax)) / pixel_size;
		const double r_max = distance * std::tan(2.*std::asin(wavelength/2./hitfinderDmin)) / pixel_size;

		float *pix_r = (float *) calloc(pix_nn, sizeof(float));
        char *mask = (char*) calloc(pix_nn, sizeof(char));
		for (long iy = 0; iy < asic_ny; ++iy)
			for (long ix = 0; ix < asic_nx; ++ix) {
				const long i = iy * asic_nx + ix;
				pix_r[i] = std::sqrt((ix-beam_x)*(ix-beam_x) + (iy-beam_y)*(iy-beam_y));
				mask[i] = r_min <= pix_r[i] && pix_r[i] <= r_max ? 1 : 0;
			}

		// Reference: source/libcheetah/src/peakfinders.cpp
		
		long nPeaks = algorithm == 8 ?  peakfinder8<T>(&peaklist, &data[0], mask, pix_r,
													   asic_nx, asic_ny, nasics_x, nasics_y,
													   hitfinderADCthresh, hitfinderMinSNR, hitfinderMinPixCount,
													   hitfinderMaxPixCount, hitfinderLocalBGRadius)
			: peakfinder6<T>(&peaklist, &data[0], mask, asic_nx, asic_ny, nasics_x, nasics_y,
							 hitfinderADCthresh, hitfinderMinSNR, hitfinderMinPixCount,
							 hitfinderMaxPixCount, hitfinderLocalBGRadius, hitfinderMinPeakSeparation);

        /*
         *      Too many peaks for the peaklist counter?
         */
        if (nPeaks > peaklist.nPeaks_max) {
			nPeaks = peaklist.nPeaks_max;
			peaklist.nPeaks = peaklist.nPeaks_max;
        }

        /*
         *      Find physical position of peaks on assembled detector
         */
        const long np = nPeaks > peaklist.nPeaks_max ? peaklist.nPeaks_max : nPeaks;

        for(long k=0; k<nPeaks && k<peaklist.nPeaks_max; ++k) {
			const long e = peaklist.peak_com_index[k];
			peaklist.peak_com_x_assembled[k] = e % asic_nx;
			peaklist.peak_com_y_assembled[k] = static_cast<float>(e / asic_nx);
			peaklist.peak_com_r_assembled[k] = pix_r[e];
        }

		/*
         *      eliminate closely spaced peaks
         */
        if(hitfinderMinPeakSeparation > 0 )
			nPeaks = killNearbyPeaks(&peaklist, hitfinderMinPeakSeparation);

		/*
		 * Prepare return value and calculate resolution in Angstrom
		 */
		boost::python::list ret;
		//std::vector<CheetahSpot> ret;
		//ret.reserver(np);
		for(long k=0; k < np; ++k) {
			//ret.push_back(CheetahSpot());
			CheetahSpot s;// = ret.back();
			s.peak_maxintensity = peaklist.peak_maxintensity[k];
			s.peak_totalintensity = peaklist.peak_totalintensity[k];
			s.peak_sigma = peaklist.peak_sigma[k];
			s.peak_snr = peaklist.peak_snr[k];
			s.peak_npix = peaklist.peak_npix[k];
			s.peak_com_x = peaklist.peak_com_x_assembled[k];
			s.peak_com_y = peaklist.peak_com_y_assembled[k];
			s.peak_com_r = peaklist.peak_com_r_assembled[k];

			const double r = peaklist.peak_com_r_assembled[k];
			const float r2 = std::sqrt(distance*distance+(pixel_size*pixel_size*r*r));
			const double twotheta = std::asin(pixel_size*r/r2);
			const double sintheta = std::sin(twotheta/2.0);
			const float A = wavelength/(2*sintheta);
			s.peak_com_res = A;
			ret.append(&s);
		}

        // Release memory
        free(mask);
		free(pix_r);
		freePeakList(peaklist);
		return ret;
	}

private:
	float hitfinderADCthresh;
	float hitfinderMinSNR;
	long hitfinderMinPixCount;
	long hitfinderMaxPixCount;
	long hitfinderLocalBGRadius;
	float hitfinderMinPeakSeparation;
	float hitfinderDmin, hitfinderDmax;
};

BOOST_PYTHON_MODULE(cheetah_ext)
{
	bp::class_<CheetahSpot>("CheetahSpot")
		.add_property("peak_maxintensity", &CheetahSpot::peak_maxintensity)
		.add_property("peak_totalintensity", &CheetahSpot::peak_totalintensity)
		.add_property("peak_sigma", &CheetahSpot::peak_sigma)
		.add_property("peak_snr", &CheetahSpot::peak_snr)
		.add_property("peak_npix", &CheetahSpot::peak_npix)
		.add_property("peak_com_x", &CheetahSpot::peak_com_x)
		.add_property("peak_com_y", &CheetahSpot::peak_com_y)
		.add_property("peak_com_r", &CheetahSpot::peak_com_r)
		.add_property("peak_com_res", &CheetahSpot::peak_com_res)
		;

	bp::class_<CheetahSinglePanel>("CheetahSinglePanel")
		.def("set_params", &CheetahSinglePanel::set_params,
			 (bp::arg("ADCthresh")=5, bp::arg("MinSNR")=8, bp::arg("MinPixCount")=2,
			  bp::arg("MaxPixCount")=40, bp::arg("LocalBGRadius")=2, bp::arg("MinPeakSeparation")=0,
			  bp::arg("Dmin")=5, bp::arg("Dmax")=30))
		.def("run_float", &CheetahSinglePanel::run<float>,
			 (bp::arg("asic_nx"), bp::arg("asic_ny"), bp::arg("data"),
			  bp::arg("beam_x"), bp::arg("beam_y"), bp::arg("wavelength"),
			  bp::arg("distance"), bp::arg("pixel_size"),
			  bp::arg("algorithm")))
		.def("run_uint16", &CheetahSinglePanel::run<uint16_t>,
			 (bp::arg("asic_nx"), bp::arg("asic_ny"), bp::arg("data"),
			  bp::arg("beam_x"), bp::arg("beam_y"), bp::arg("wavelength"),
			  bp::arg("distance"), bp::arg("pixel_size"),
			  bp::arg("algorithm")))
		.def("run_uint32", &CheetahSinglePanel::run<uint32_t>,
			 (bp::arg("asic_nx"), bp::arg("asic_ny"), bp::arg("data"),
			  bp::arg("beam_x"), bp::arg("beam_y"), bp::arg("wavelength"),
			  bp::arg("distance"), bp::arg("pixel_size"),
			  bp::arg("algorithm")))


		;
}
