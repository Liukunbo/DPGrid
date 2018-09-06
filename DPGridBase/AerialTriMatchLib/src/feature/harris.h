#ifndef FEA_HARRIS_HEADER
#define FEA_HARRIS_HEADER

#include <string>
#include <vector>

#include <random>

#include "math/vector.h"
#include "image/image.h"
#include "feature/defines.h"

FEA_NAMESPACE_BEGIN

// 1.
//    Ix = (p[c-1])-(p[c+1]) // dX
//    Iy = (p[r-1])-(p[r+1]) // dY
// 2.
//       Ix*Ix   Ix*Iy
// M = [               ]
//       Ix*Iy   Iy*Iy
// 3.
//                       (x*x+y*y)
// Gauss = exp( (-1) * --------------- , sigm = 0.7
//                      2* sigm*sigm
// 4.
//  cim = det(M) - r(trace(M)        
//   = ((Ix*Ix)*(Iy*Iy) - (Ix*Iy)*(Ix*Iy)) - 0.04*(Ix*Ix+Iy*Iy)*(Ix*Ix+Iy*Iy);
// 

class Harris
{
public:
	struct Options
	{
		Options(void);

		/* should do wallis filte befor harris feature extact */
		bool do_wallis_filte;

		/* the window size used to extract harris feature */
		int window_size;

		/* the overlap between adjacent window, smaller than window size */
		int overlap_size;

		/* parameters used to do wallis, default is the probably the best */
		double meanValue;
		double sigmaValue;
		double CValue;
		double BValue;
		double gridWZ;
		double filterWZ;

		/**
		* Produce even more messages on the console.
		*/
		bool debug_output;
	};

	struct Keypoint
	{
		/** Keypoint x-coordinate. */
		short x;
		/** Keypoint y-coordinate. */
		short y;
	};

	struct Descriptor
	{
		/** The sub-pixel x-coordinate of the image keypoint. */
		short x;
		/** The sub-pixel y-coordinate of the image keypoint. */
		short y;
		/** The scale (or sigma value) of the keypoint. */
		float scale;
		/** The descriptor data, elements are unsigned char in [0, 255].  Add only the first 121 was used */
		math::Vector<uint8_t, 128> data;
	};

public:
	typedef std::vector<Keypoint> Keypoints;
	typedef std::vector<Descriptor> Descriptors;

public:
	explicit Harris(Options const& options);

	/** Sets the input image. */
	void set_image(dpgrid::ByteImage::ConstPtr img);

	/** Starts the Harris keypoint detection and descriptor extraction. */
	void process(void);

	/** Returns the list of keypoints. */
	Keypoints const& get_keypoints(void) const;

	/** Returns the list of descriptors. */
	Descriptors const& get_descriptors(void) const;

private:
	bool Harris_ExtrFeatPt(unsigned char *pImg, float *pIx2, float *pIy2, float *pIxy, float *pM, int grdSz, float *px, float *py, float * xvalue, float mxv = 2);
	bool Harris_ExtrFeat(unsigned char *pImg, int cols, int rows, int grdSz0, int ovlp, int sC = 0, int eC = -1, int sR = 0, int eR = -1, float mxv = 2);
private:
	Options opt_;
	dpgrid::ByteImage::Ptr orig_; // Original input image
	Keypoints keypoints_; // Detected keypoints
	Descriptors descriptors_; // Final Harris descriptors
};

/* ---------------------------------------------------------------- */
inline 
Harris::Options::Options(void)
	:do_wallis_filte(true)
	, window_size(11)
	, overlap_size(0)
	, meanValue(137)
	, sigmaValue(190)
	, CValue(0.8)
	, BValue(0.9)
	, gridWZ(16)
	, filterWZ(32)
	, debug_output(false)
{
}

inline Harris::Keypoints const&
Harris::get_keypoints(void) const
{
	return this->keypoints_;
}

inline Harris::Descriptors const&
Harris::get_descriptors(void) const
{
	return this->descriptors_;
}

FEA_NAMESPACE_END

#endif
