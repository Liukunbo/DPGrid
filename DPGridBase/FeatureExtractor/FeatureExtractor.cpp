#include "FeatureExtractor.h"

#include <iostream>
#include <iomanip>
#include <fstream>

#include "AerialTriMatchLib.h"

#define DEFAULT_MAX_IMAGE_SIZE 6000000
#define MIN_IMAGE_SIZE 64 * 64

FeatureExtractor::FeatureExtractor(FeatureExtractor::Options const & options)
	: opt_(options), max_image_size_(DEFAULT_MAX_IMAGE_SIZE)
{
}

FeatureExtractor::~FeatureExtractor()
{
}

void FeatureExtractor::SetOptions(FeatureExtractor::Options const & options)
{
	opt_ = options;

	/* If we extract Harris, then we should extract on the origin level. */
	if (opt_.feature_type == FEATURE_HARRIS)
	{
		opt_.feature_lvl = LEVEL_1;
	}
}

void FeatureExtractor::SetMaxImageSize(std::size_t maxSize)
{
	if (maxSize < MIN_IMAGE_SIZE)
	{
		std::cerr << "*** Error:The image sizes cannot be smaller than 64 x 64" << std::endl;
	}
	max_image_size_ = maxSize;
}

int FeatureExtractor::FeatureDetecHarris(std::string const szImagePath, std::string const szFeaturePath)
{
	if (this->opt_.feature_type != FEATURE_HARRIS)
	{
		return EXIT_FAILURE;
	}

	/* load image */
	dpgrid::ByteImage::Ptr image;
	try { image = dpgrid::image::load_file(szImagePath); }
	catch (const std::exception&) { std::cerr << "*** Error: Can't open file: " << szImagePath << std::endl;	return EXIT_FAILURE; }
	std::cout << "\nOpen file: " << szImagePath << " Success" << std::endl;

	int nOriWid = image->width(), nOriHei = image->height();
	std::cout << "Wid: " << nOriWid << "; Hei: " << nOriHei << std::endl;

	/* extract harris on origin lvl. */
	fea::FeatureSet::Options feature_options;
	fea::Harris::Options harris_option; harris_option.debug_output = this->opt_.debug_output;
	harris_option.window_size = 50;
	feature_options.harris_opts = harris_option;
	feature_options.feature_types = fea::FeatureSet::FEATURE_HARRIS;

	fea::FeatureSet harris_features;
	harris_features.set_options(feature_options);
	harris_features.compute_features(image);
	harris_features.scale = 1;

	// extract sift
	//int nOriWid = image->width();
	while (image->width() * image->height() > max_image_size_)
	{
		image = dpgrid::image::rescale_half_size<uint8_t>(image);
	}
	int nScale = nOriWid / image->width();
	std::cout << "Max Image Size: " << max_image_size_ << std::endl;

	feature_options.feature_types = fea::FeatureSet::FEATURE_SIFT;

	fea::FeatureSet features;
	features.set_options(feature_options);
	features.compute_features(image);
	features.scale = nScale;

	/* feature num is too small */
	if (features.positions.size() <= 10)
	{
		std::cout << "Cannot detec features!" << std::endl;
		return EXIT_FAILURE;
	}

	// combine feature
	fea::CombineFeature(features, harris_features);

	std::cout << "Number of Features: " << features.positions.size() << std::endl;
	std::cout << "Number of Sift: " << features.sift_descriptors.size() << std::endl;
	std::cout << "Number of Harris: " << features.harris_descriptors.size() << std::endl;

	/* Normalize feature coordinates. */
	int sift_fea_size = features.sift_descriptors.size();
	int harris_fea_size = features.harris_descriptors.size();
	float fwidth = static_cast<float>(features.width);
	float fheight = static_cast<float>(features.height);
	float fnorm = std::max(fwidth, fheight);
	if (this->opt_.reverse_y_coord)
	{
		for (std::size_t j = 0; j < sift_fea_size; ++j)
		{
			math::Vec2f& pos = features.positions[j];
			pos[1] = fheight - 1 - pos[1]; /* reverse the y coordinate */
		}
		for (std::size_t j = 0; j < harris_fea_size; ++j)
		{
			math::Vec2f& pos = features.positions[j + sift_fea_size];
			pos[1] = fheight * features.scale - 1 - pos[1]; /* reverse the y coordinate */
		}

	}
	if (this->opt_.normalize_coord)
	{
		for (std::size_t j = 0; j < sift_fea_size; ++j)
		{
			math::Vec2f& pos = features.positions[j];
			pos[0] = (pos[0] + 0.5f - fwidth / 2.0f) / fnorm;
			pos[1] = (pos[1] + 0.5f - fheight / 2.0f) / fnorm;
		}

		fwidth = static_cast<float>(features.width * features.scale);
		fheight = static_cast<float>(features.height * features.scale);
		fnorm = std::max(fwidth, fheight);
		for (std::size_t j = 0; j < harris_fea_size; ++j)
		{
			math::Vec2f& pos = features.positions[j + sift_fea_size];
			pos[0] = (pos[0] + 0.5f - fwidth / 2.0f) / fnorm;
			pos[1] = (pos[1] + 0.5f - fheight / 2.0f) / fnorm;
		}
	}

	/* save feature file */
	int ret = fea::fea_mch::saveFeatureFile(szFeaturePath, nScale, features.width, features.height, features);
	if (this->opt_.debug_output)
	{
		fea::fea_mch::saveHarrisFeatureTxt(szFeaturePath + ".txt", nScale, features.width, features.height, features);
	}

	features.positions.clear();
	features.colors.clear();
	features.clear_descriptors();

	return ret;
}

int FeatureExtractor::FeatureDetec(std::string const szImagePath, std::string const szFeaturePath)
{
	if (this->opt_.feature_type == FEATURE_HARRIS)
	{
		return FeatureDetecHarris(szImagePath, szFeaturePath);
	}

	/* load image */
	dpgrid::ByteImage::Ptr image;
	try { image = dpgrid::image::load_file(szImagePath); }
	catch (const std::exception&) { std::cerr << "*** Error: Can't open file: " << szImagePath << std::endl;	return EXIT_FAILURE; }
	std::cout << "\nOpen file: " << szImagePath << " Success" << std::endl;

	int nOriWid = image->width(), nOriHei = image->height();
	std::cout << "Wid: " << nOriWid << "; Hei: " << nOriHei << std::endl;

	/* scale image */
	int nScale = 1;
	if (opt_.feature_lvl == LEVEL_AUTO)
	{
		while (image->width() * image->height() > max_image_size_)
		{
			image = dpgrid::image::rescale_half_size<uint8_t>(image);
		}
		nScale = nOriWid / image->width();

		std::cout << "Max Image Size: " << max_image_size_ << std::endl;	
	}
	else
	{
		for (std::size_t i = 0; i < (std::size_t)opt_.feature_lvl; i++)
		{
			image = dpgrid::image::rescale_half_size<uint8_t>(image);
		}
		nScale = nOriWid / image->width();
	}

	std::cout << "Scaled : " << nScale << std::endl;

	/* feature detec */
	fea::FeatureSet::Options feature_options;
	fea::Sift::Options sift_option;
	//sift_option.max_octave = 3;
	//sift_option.num_samples_per_octave = 2;
	feature_options.feature_types = (fea::FeatureSet::FeatureTypes)opt_.feature_type;
	feature_options.sift_opts = sift_option;

	fea::FeatureSet features;
	features.set_options(feature_options);
	features.compute_features(image);
	image.reset();

	features.scale = nScale;

	/* feature num is too small */
	if (features.positions.size() <= 10)
	{
		std::cout << "Cannot detec features!" << std::endl;
		return EXIT_FAILURE;
	}

	std::cout << "Number of Features: " << features.positions.size() << std::endl;
	std::cout << "Number of Sift: " << features.sift_descriptors.size() << std::endl;
	std::cout << "Number of Surf: " << features.surf_descriptors.size() << std::endl;

	/* Normalize feature coordinates. */		
	float const fwidth = static_cast<float>(features.width);
	float const fheight = static_cast<float>(features.height);
	float const fnorm = std::max(fwidth, fheight);
	if (this->opt_.reverse_y_coord)
	{
		for (std::size_t j = 0; j < features.positions.size(); ++j)
		{
			math::Vec2f& pos = features.positions[j];
			pos[1] = fheight - 1 - pos[1]; /* reverse the y coordinate */
		}
	}
	if (this->opt_.normalize_coord)
	{
		/* normalize */
		for (std::size_t j = 0; j < features.positions.size(); ++j)
		{
			math::Vec2f& pos = features.positions[j];
			pos[0] = (pos[0] + 0.5f - fwidth / 2.0f) / fnorm;
			pos[1] = (pos[1] + 0.5f - fheight / 2.0f) / fnorm; 
		}
	}	
	
	/* save feature file */
	int ret = fea::fea_mch::saveFeatureFile(szFeaturePath, nScale, features.width, features.height, features);
	if (this->opt_.debug_output)
	{
		fea::fea_mch::saveFeatureFileTxt(szFeaturePath + ".txt", nScale, features.width, features.height, features);
	}	

	features.positions.clear();
	features.colors.clear();
	features.clear_descriptors();

	return ret;
}


