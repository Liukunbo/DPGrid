#include <iostream>

#include "AerialTriMatchLib.h"
#if _DEBUG
#pragma comment(lib, "AerialTriMatchLibD.lib")
#else
#pragma comment(lib, "AerialTriMatchLib.lib")
#endif

int load_image(std::string szImagePath, dpgrid::ByteImage::Ptr & image)
{
	//dpgrid::ByteImage::Ptr image;
	try { image = dpgrid::image::load_file(szImagePath); }
	catch (const std::exception&) { std::cerr << "*** Error:Can't open file: " << szImagePath << std::endl;	return EXIT_FAILURE; }
	std::cout << "\nOpen file: " << szImagePath << " Success" << std::endl;

	int nOriWid = image->width(), nOriHei = image->height();
	std::cout << "Wid: " << nOriWid << "; Hei: " << nOriHei << std::endl;

	return EXIT_SUCCESS;
}

int extract_feature(std::string szImagePath, fea::FeatureSet & feature)
{
	dpgrid::ByteImage::Ptr img;
	if (load_image(szImagePath, img) != EXIT_SUCCESS)
	{
		return EXIT_FAILURE;
	}
	// extract harris
	fea::FeatureSet::Options feature_options;
	fea::Harris::Options harris_option; harris_option.debug_output = true;
	harris_option.window_size = 50;
	feature_options.harris_opts = harris_option;
	feature_options.feature_types = fea::FeatureSet::FEATURE_HARRIS;
	
	fea::FeatureSet harris_features;
	harris_features.set_options(feature_options);
	harris_features.compute_features(img);
	harris_features.scale = 1;

	// extract sift
	int max_image_size = 6000000;
	int nOriWid = img->width();
	while (img->width() * img->height() > max_image_size)
	{
		img = dpgrid::image::rescale_half_size<uint8_t>(img);
	}
	int nScale = nOriWid / img->width();
	std::cout << "Max Image Size: " << max_image_size << std::endl;

	feature_options.feature_types = fea::FeatureSet::FEATURE_SIFT;
	feature.set_options(feature_options);
	feature.compute_features(img);
	feature.scale = nScale;

	// combine feature
	fea::CombineFeature(feature, harris_features);

	// normalize
	float fwidth = static_cast<float>(feature.width);
	float fheight = static_cast<float>(feature.height);
	float fnorm = std::max(fwidth, fheight);
	for (std::size_t j = 0; j < feature.sift_descriptors.size(); j++)
	{
		math::Vec2f& pos = feature.positions[j];
		pos[0] = (pos[0] + 0.5f - fwidth / 2.0f) / fnorm;
		pos[1] = (pos[1] + 0.5f - fheight / 2.0f) / fnorm;
	}

	fwidth = static_cast<float>(feature.width * feature.scale);
	fheight = static_cast<float>(feature.height * feature.scale);
	fnorm = std::max(fwidth, fheight);
	for (std::size_t j = 0; j < feature.harris_descriptors.size(); j++)
	{
		math::Vec2f& pos = feature.positions[j + feature.sift_descriptors.size()];
		pos[0] = (pos[0] + 0.5f - fwidth / 2.0f) / fnorm;
		pos[1] = (pos[1] + 0.5f - fheight / 2.0f) / fnorm;
	}

	return EXIT_SUCCESS;
}

int match_harris(std::string strImage1, std::string strImage2, std::string strMchFile)
{
	/* extract feat */
	fea::FeatureSet fst1, fst2;
	extract_feature(strImage1, fst1);
	extract_feature(strImage2, fst2);

	/* match */
	/* Set Options */
	fea::fea_mch::dpg_mch::Matching::Options matching_opts;
	matching_opts.ransac_opts.verbose_output = false;
	matching_opts.match_num_previous_frames = 0; // disabled 
	matching_opts.min_feature_matches = 24;
	matching_opts.min_matching_inliers = (int)(24 / 5.0);
	matching_opts.use_lowres_matching = false;
	matching_opts.num_lowres_features = 1000;
	matching_opts.matcher_type = fea::fea_mch::dpg_mch::Matching::MatcherType::MATCHER_HARRIS;

	/* Compute match */
	fea::FeatureSets features;
	features.emplace_back(fst1);
	features.emplace_back(fst2);
	fea::fea_mch::dpg_mch::Matching dpg_matching(matching_opts);
	dpg_matching.init(&features);
	fea::PairwiseMatching pairwise_matching;
	dpg_matching.compute2(&pairwise_matching);

	std::vector<std::string> szFeaList; szFeaList.emplace_back(strImage1); szFeaList.emplace_back(strImage2);
	fea::fea_mch::saveMatchFileTxt(strMchFile + ".txt", pairwise_matching, szFeaList, &features);

	return EXIT_SUCCESS;
}

int main(int argc, char* argv[])
{
	match_harris("F:\\Data_In\\20170328\\image\\DJI_0001.JPG", "F:\\Data_In\\20170328\\image\\DJI_0002.JPG", "F:\\Data_Work\\match.mch");
	return 0;
}