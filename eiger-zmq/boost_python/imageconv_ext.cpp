// Original code is yamtbx/dataproc/myspotfinder/core_toolbox/image_conv.hpp (New BSD License)
// The original of original is iotbx/detectors/display.h. (New BSD License)

/*
g++ -O3 -c imageconv_ext.cpp -shared -fPIC  -I../ -I/oys/python/python-2.7.6/include/python2.7/ -I/oys/python/python-2.7.6/lib/python2.7/site-packages/PyUblas-2013.1-py2.7-linux-x86_64.egg/pyublas/include -I/oys/python/python-2.7.6/lib/python2.7/site-packages/numpy/core/include
g++ -shared imageconv_ext.o -o imageconv_ext.so -L/oys/xtal/cctbx/build/lib -lboost_python -ljpeg
 */
#include <stdio.h>
#include <vector>
#include <functional> 
#include <algorithm>
#include <jpeglib.h>
#include <boost/python.hpp>
#include <boost/python/args.hpp>
#include <pyublas/numpy.hpp>

namespace bp = boost::python;

template<class T>
class MyImage {

public:
	typedef float data_t;
	int raw_size1, raw_size2;

	double global_bright_contrast() const {
		return local_bright_contrast(0, 0, size2(), size1());
	}

	double local_bright_contrast(int x0, int y0, int w, int h) const {
		const double adjlevel = 0.4;

		std::vector<data_t> roi(w*h, 0);

		for (int y_ = 0; y_ < h; ++y_)
			for (int x_ = 0; x_ < w; ++x_) {
				const int x = x_ + x0, y = y_ + y0;
				const int i = y * size2() + x;
				const int i_ = y_*w + x_;

				roi[i_] = rawdata[i];
			}
    
		const bool is_pilatus = (vendortype=="Pilatus-6M"||vendortype=="Pilatus-2M"||vendortype=="Pilatus-300K");

		// count dead area. should be count_if (<0)
		const int dead_area = is_pilatus ? std::count(roi.begin(), roi.end(), -1) : 0;
		const int saturated_num = std::count_if(roi.begin(), roi.end(), 
												std::bind2nd(std::greater_equal<int>(), saturation));
		const int nth_offset = 0.9 * (w*h - saturated_num) + 0.1 * dead_area;

		std::nth_element(roi.begin(), roi.begin() + nth_offset, roi.end());
		const double percentile = roi[nth_offset];

		return (percentile>0.) ? brightness * adjlevel/percentile : 1.0;
	}

	//data_t *d_rawdata;
	pyublas::numpy_vector<T> rawdata;
	std::string vendortype;
	double brightness;
	data_t saturation;
public:

	MyImage(pyublas::numpy_vector<T> rawdata, 
			int width, int height,
			const std::string& vendortype, const double& brightness = 1.0,
			const data_t& saturation = 65535):
		brightness(brightness),
		saturation(saturation),
		rawdata(rawdata),
		raw_size2(width),
		raw_size1(height),
		vendortype(vendortype)
		{ }

	inline int size1() const {return raw_size1;}
	inline int size2() const {return raw_size2;}

  void prep_string(){
    const int s_size = size1()*size2();
    export_s.resize(s_size);

    const double correction = global_bright_contrast();

    if (vendortype=="Pilatus-6M" || vendortype=="Pilatus-2M" || vendortype=="Pilatus-300K") {
      for (int i = 0; i < size1()*size2(); ++i) {
	prep_str_pilatus(i, export_s, correction, saturation, size1()*size2(), size2());
      }
    }
    else if (vendortype=="EIGER") {
      for (int i = 0; i < size1()*size2(); ++i) {
	prep_str_eiger(i, export_s, correction, saturation, size1()*size2(), size2());
      }
    }
    else {
      for (int i = 0; i < size1()*size2(); ++i) {
	prep_str_inner(export_s, i, i, correction, saturation);
      }
    }
  }

  void prep_string_cropped(int x0, int y0, int w, int h) {
    const int s_size = w * h * 3;
    export_s.resize(s_size);
    assert (x0 >= 0);
    assert (y0 >= 0);
    assert (x0+w <= size2());
    assert (y0+h <= size1());

    const double correction = local_bright_contrast(x0, y0, w, h);

    if (vendortype=="Pilatus-6M" || vendortype=="Pilatus-2M" || vendortype=="Pilatus-300K") {
      // TODO check if the size2() is valid??
      for (int j = 0; j < h; ++j)
	for (int i = 0; i < w; ++i)
		prep_str_cropped_pilatus(i, j, export_s, correction, saturation, x0, y0, w, h, size2());
    }
    else if (vendortype=="EIGER") {
      for (int j = 0; j < h; ++j)
	for (int i = 0; i < w; ++i)
	  prep_str_cropped_eiger(i, j, export_s, correction, saturation, x0, y0, w, h, size2());
    }
    else {
      for (int y_ = 0; y_ < h; ++y_) 
	for (int x_ = 0; x_ < w; ++x_) {
	  const int x = x_ + x0, y = y_ + y0;
	  const int i = y * size2() + x;
	  const int i_ = y_*w + x_;
	  prep_str_inner(export_s, i, i_, correction, saturation);
	}
    }
  }

  void prep_string_cropped_and_scaled(int x0, int y0, int w, int h, int im_w, int im_h) {
    if (w==im_w && h==im_h) {
      prep_string_cropped(x0, y0, w, h);
      return;
    }

    const int s_size = im_w * im_h * 3;
    export_s.resize(s_size);
    assert (x0 >= 0);
    assert (y0 >= 0);
    assert (x0+w <= size2());
    assert (y0+h <= size1());
    assert (im_w >= 0);
    assert (im_h >= 0);

    // Origin of this idea: wxPython-src-2.8.12.1/src/common/image.cpp
    const int delta_x = (w<<16)/im_w, delta_y = (h<<16)/im_h;

    const double correction = local_bright_contrast(x0, y0, w, h);

    if (vendortype=="Pilatus-6M" || vendortype=="Pilatus-2M" || vendortype=="Pilatus-300K") {
      // TODO check if the size2() is valid??
      for (int j = 0; j < im_h; ++j)
	for (int i = 0; i < im_w; ++i)
	  prep_str_cropped_and_scaled_pilatus(i, j, export_s, correction,
					      saturation, x0, y0, w, h, size2(),
					      im_w, im_h, delta_x, delta_y);
    }
    else if (vendortype=="EIGER") {
      for (int j = 0; j < im_h; ++j)
	for (int i = 0; i < im_w; ++i)
	  prep_str_cropped_and_scaled_eiger(i, j, export_s, correction,
					    saturation, x0, y0, w, h, size2(),
					    im_w, im_h, delta_x, delta_y);
    }
    else {
      for (int y_ = 0; y_ < im_h; ++y_)
	for (int x_ = 0; x_ < im_w; ++x_) {
	  const int x = ((x_*delta_x)>>16)+x0, y = ((y_*delta_y)>>16)+y0;
	  const int i = y * size2() + x;
	  const int i_ = y_*im_w + x_;
	  prep_str_inner(export_s, i, i_, correction, saturation);
	}
    }
  }

  std::string export_s; //export data in string form; make public for readonly

  void save_as_jpeg(const std::string &filename, int width, int height, int quality) const {
    assert (export_s.size() == width*height*3);

    struct jpeg_compress_struct cinfo;
    struct jpeg_error_mgr jerr;

    /* this is a pointer to one row of image data */
    FILE *outfile = fopen(filename.c_str(), "wb");

    cinfo.err = jpeg_std_error(&jerr);
    jpeg_create_compress(&cinfo);
    jpeg_stdio_dest(&cinfo, outfile);

    /* Setting the parameters of the output file here */
    cinfo.image_width = width;
    cinfo.image_height = height;
    cinfo.input_components = 3;
    cinfo.in_color_space = JCS_RGB;

    jpeg_set_defaults(&cinfo);
    jpeg_set_quality(&cinfo, quality, TRUE);

    /* Now do the compression .. */
    jpeg_start_compress( &cinfo, TRUE );

    for(int i = 0; i < height; ++i) {
      JSAMPROW row_ptr =  (JSAMPROW)(export_s.c_str() + i*width*3);
      jpeg_write_scanlines( &cinfo, &row_ptr, 1 );
    }

    jpeg_finish_compress( &cinfo );
    jpeg_destroy_compress( &cinfo );
    fclose( outfile );
  }

private:
  inline void prep_str_inner(std::string &export_s, int i_raw, int i_s, float correction, int saturation)
  {
    const float outscale = 256;
    const float val = rawdata[i_raw];
    float corrected = val * correction;
    float outvalue  = outscale * ( 1.0 - corrected );

    if (val == -2) { // inactive pixels (-2) on Pilatus
      export_s[i_s*3+0] = (char)245;
      export_s[i_s*3+1] = (char)1;
      export_s[i_s*3+2] = (char)1;
    } else if (val >= saturation) {
      export_s[i_s*3+0] = (char)254;
      export_s[i_s*3+1] = (char)254;
      export_s[i_s*3+2] = (char)1;
    } else if (outvalue < 0.0) {
      export_s[i_s*3+0] = 0;
      export_s[i_s*3+1] = 0;
      export_s[i_s*3+2] = 0;
    } else if (outvalue >= outscale) {
      export_s[i_s*3+0] = char(int(outscale)-1);
      export_s[i_s*3+1] = char(int(outscale)-1);
      export_s[i_s*3+2] = char(int(outscale)-1);
    } else {
      export_s[i_s*3+0] = int(outvalue);
      export_s[i_s*3+1] = int(outvalue);
      export_s[i_s*3+2] = int(outvalue);
    }
  }

  inline void prep_str_pilatus(int i, std::string &export_s, float correction, int saturation, int size, int nfast)
  {
    if (i >= size)
      return;

    const int y = i/nfast, x = i%nfast;

    if ( (y%212<195) && (x%494<487) )
      prep_str_inner(export_s, i, i, correction, saturation);
    else {
      export_s[i*3+0] = (char)255;
      export_s[i*3+1] = (char)228;
      export_s[i*3+2] = (char)228;
    }
  }

  inline void prep_str_cropped_pilatus(int x_, int y_, std::string &export_s, float correction, int saturation,
				int x0, int y0, int w, int h, int width)
  {
    if (x_>=w) return;
    if (y_>=h) return;

    const int x = x_ + x0, y = y_ + y0;
    const int i = y * width + x;
    const int i_ = y_*w + x_;

    if ( (y%212<195) && (x%494<487) )
      prep_str_inner(export_s, i, i_, correction, saturation);
    else {
      export_s[i_*3+0] = (char)255;
      export_s[i_*3+1] = (char)228;
      export_s[i_*3+2] = (char)228;
    }
  }

  inline void prep_str_cropped_and_scaled_pilatus(int x_, int y_, std::string &export_s, float correction, int saturation,
					   int x0, int y0, int w, int h, int width, int im_w, int im_h,
					   int delta_x, int delta_y)
  {
    if (x_>=im_w) return;
    if (y_>=im_h) return;

    const int x = ((x_*delta_x)>>16)+x0, y = ((y_*delta_y)>>16)+y0;
    const int i = y * width + x;
    const int i_ = y_*im_w + x_;

    if ( (y%212<195) && (x%494<487) )
      prep_str_inner(export_s, i, i_, correction, saturation);
    else {
      export_s[i_*3+0] = (char)255;
      export_s[i_*3+1] = (char)228;
      export_s[i_*3+2] = (char)228;
    }
  }

  inline void prep_str_eiger(int i, std::string &export_s, float correction, int saturation, int size, int nfast)
  {
    if (i >= size)
      return;

    const int y = i/nfast, x = i%nfast;

    if ( (y%551<514) && (x%1040<1030) )
      prep_str_inner(export_s, i, i, correction, saturation);
    else {
      export_s[i*3+0] = (char)255;
      export_s[i*3+1] = (char)228;
      export_s[i*3+2] = (char)228;
    }
  }

  inline void prep_str_cropped_eiger(int x_, int y_, std::string &export_s, float correction, int saturation,
				int x0, int y0, int w, int h, int width)
  {
    if (x_>=w) return;
    if (y_>=h) return;

    const int x = x_ + x0, y = y_ + y0;
    const int i = y * width + x;
    const int i_ = y_*w + x_;

    if ( (y%551<514) && (x%1040<1030) )
      prep_str_inner(export_s, i, i_, correction, saturation);
    else {
      export_s[i_*3+0] = (char)255;
      export_s[i_*3+1] = (char)228;
      export_s[i_*3+2] = (char)228;
    }
  }

  inline void prep_str_cropped_and_scaled_eiger(int x_, int y_, std::string &export_s, float correction, int saturation,
					   int x0, int y0, int w, int h, int width, int im_w, int im_h,
					   int delta_x, int delta_y)
  {
    if (x_>=im_w) return;
    if (y_>=im_h) return;

    const int x = ((x_*delta_x)>>16)+x0, y = ((y_*delta_y)>>16)+y0;
    const int i = y * width + x;
    const int i_ = y_*im_w + x_;

    if ( (y%551<514) && (x%1040<1030) )
      prep_str_inner(export_s, i, i_, correction, saturation);
    else {
      export_s[i_*3+0] = (char)255;
      export_s[i_*3+1] = (char)228;
      export_s[i_*3+2] = (char)228;
    }
  }
};

template<class T>
struct my_image_wrapper {

	typedef MyImage<T> w_t;
	typedef pyublas::numpy_vector<T> array_t;

	static void wrap (const char* python_name){
		bp::class_<w_t >(python_name, bp::no_init)
			.def(bp::init<array_t, int, int, const std::string&,
				 double const&, int const& >(
					 (
						 bp::arg("rawdata"),
						 bp::arg("width"), bp::arg("height"),
						 bp::arg("vendortype"),
						 bp::arg("brightness"),
						 bp::arg("saturation")
						 )
					 ))
			.def("size1", &w_t::size1)
			.def("size2", &w_t::size2)
			.def("prep_string",&w_t::prep_string)
			.def("prep_string_cropped",&w_t::prep_string_cropped, (bp::arg_("x0"),bp::arg_("y0"),bp::arg_("w"),bp::arg_("h")))
			.def("prep_string_cropped_and_scaled",&w_t::prep_string_cropped_and_scaled, (bp::arg_("x0"),bp::arg_("y0"),bp::arg_("w"),bp::arg_("h"),bp::arg_("im_w"),bp::arg_("im_h")))
			.def("save_as_jpeg", &w_t::save_as_jpeg, (bp::arg_("width"),bp::arg_("height"),bp::arg_("quality")))

			.def_readonly("export_string",&w_t::export_s)
			;
	}
};

BOOST_PYTHON_MODULE(imageconv_ext)
{
	my_image_wrapper<uint16_t>::wrap("MyImage_uint16");
	my_image_wrapper<uint32_t>::wrap("MyImage_uint32");
}
