#pragma once

#include <string>

#if _WIN32
# ifndef MAKE_DLL
#   define FeatureExtractorAPI __declspec(dllimport)
# else
#	define FeatureExtractorAPI __declspec(dllexport)
# endif
#else
#	define FeatureExtractorAPI
#endif

//The Default Feature Point's Coordinate System is Showned Below:
//	o------> (x pixels)
//  \
//  \
//  \
// \/ (y pixels)

class FeatureExtractorAPI FeatureExtractor
{
public:
	enum FeatureTypes
	{
		FEATURE_SIFT = 0x1,
		FEATURE_SURF = 0x2,
		FEATURE_HARRIS = 0x3, /** Harris will always be extracted at the origin lvl.*/
		FEATURE_SIFTSURF = 0x4
	};

	enum FeatureLevel
	{
		LEVEL_1 = 0x0,
		LEVEL_2 = 0x1,
		LEVEL_4 = 0x2,
		LEVEL_8 = 0x3,
		LEVEL_AUTO = 0x4
	};

	struct Options
	{
		/* feature type. default is sift */
		FeatureTypes feature_type = FEATURE_SIFT;

		/* the image scale that will be used to extract feature. Default is LEVEL_1.
			if set to LEVEL_AUTO, then the maximagesize will be used to calculate the scale.
		*/
		FeatureLevel feature_lvl = LEVEL_1;

		/* Save debug information. Default is false */
		bool debug_output = false;

		/* Do we need to reverse the y coordinates (from up-down to  down-up).
		    default is false. Any change to this will not work!
		*/
		bool reverse_y_coord = false; 

		/* Do we need to normalize the coordinates.
		   default is true. Any change to this will not work!
		*/
		bool normalize_coord = true;
	};

public:
	explicit FeatureExtractor(FeatureExtractor::Options const & options);
	~FeatureExtractor();

	void SetOptions(FeatureExtractor::Options const & options);

	void SetMaxImageSize(std::size_t maxSize);

	int FeatureDetecHarris(std::string const szImagePath, std::string const szFeaturePath);

	int FeatureDetec(std::string const szImagePath, std::string const szFeaturePath);

private:
	Options opt_;
	std::size_t max_image_size_;
};

