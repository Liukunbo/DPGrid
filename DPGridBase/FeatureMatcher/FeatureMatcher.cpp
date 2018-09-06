#include "FeatureMatcher.h"

#include <iostream>
#include <iomanip>
#include <fstream>

#include "AerialTriMatchLib.h"

#include "vector\tri_define.hpp"
#include "vector\vector_function.hpp"
#include "io_file.hpp"

typedef fea::CorrespondenceIndices remap_feature_id;


FeatureMatcher::FeatureMatcher(FeatureMatcher::Options const & options)
	:opt_(options)
{
}

FeatureMatcher::~FeatureMatcher()
{
}

void FeatureMatcher::SetOptions(FeatureMatcher::Options const & options)
{
	opt_ = options;
}

int FeatureMatcher::Match(std::string szFea1, std::string szFea2, std::string szMchFile)
{
	std::cout << std::endl;

	/* Read Features */
	fea::FeatureSet features1;
	int Scale1 = 0, FeaW1 = 0, FeaH1 = 0;
	int ret1 = fea::fea_mch::readFeatureFile(szFea1, Scale1, FeaW1, FeaH1, features1);
	if (ret1 == EXIT_FAILURE) { std::cerr << "*** Error: Can't open file: " << szFea1 << std::endl; return EXIT_FAILURE; }
	std::cout << "Open file: " << szFea1 << " Success" << std::endl;
	fea::FeatureSet features2;
	int Scale2 = 0, FeaW2 = 0, FeaH2 = 0;
	int ret2 = fea::fea_mch::readFeatureFile(szFea2, Scale2, FeaW2, FeaH2, features2);
	if (ret2 == EXIT_FAILURE) { std::cerr << "*** Error: Can't open file: " << szFea2 << std::endl; return EXIT_FAILURE; }
	std::cout << "Open file: " << szFea2 << " Success" << std::endl;

	/* Organize features */
	fea::FeatureSet::Options feature_options;
	feature_options.feature_types = fea::FeatureSet::FEATURE_SIFT;
	fea::FeatureSets features;
	features.clear();	features.resize(2);
	features1.set_options(feature_options);
	features.at(0) = features1;
	std::cout << "Number of Features One: " << features.at(0).get_feature_num() << std::endl;
	features2.set_options(feature_options);
	features.at(1) = features2;
	std::cout << "Number of Features Two: " << features.at(1).get_feature_num() << std::endl;
	std::size_t minFea = std::min(features.at(0).get_feature_num(), features.at(1).get_feature_num());
	if (minFea <= 0)
	{
		std::cerr << "*** Error: Sift Features do not Exist! " << std::endl; return EXIT_FAILURE;
	}

	/* Set Options */
	fea::fea_mch::dpg_mch::Matching::Options matching_opts;
	matching_opts.ransac_opts.verbose_output = false;
	matching_opts.match_num_previous_frames = 0; // disabled 
	matching_opts.min_feature_matches = this->opt_.min_feature_matches;
	matching_opts.min_matching_inliers = (int)(this->opt_.min_feature_matches / 5.0);
	matching_opts.use_lowres_matching = this->opt_.use_lowres_matching;
	matching_opts.num_lowres_features = (int)(minFea / 10.0);;
	matching_opts.min_lowres_matches = (int)(matching_opts.num_lowres_features / 5.0);
	matching_opts.matcher_type = (fea::fea_mch::dpg_mch::Matching::MatcherType)this->opt_.matcher_type;

	/* Compute match */
	fea::fea_mch::dpg_mch::Matching dpg_matching(matching_opts);
	dpg_matching.init(&features);
	fea::PairwiseMatching pairwise_matching;
	dpg_matching.compute(&pairwise_matching);
	if (pairwise_matching.size() < 1) { std::cerr << "*** Error: Wrong matching image pairs. Exiting." << std::endl; return EXIT_FAILURE; }

	/* Save result */
	std::vector<std::string> feat_list; feat_list.emplace_back(szFea1); feat_list.emplace_back(szFea2);
	fea::fea_mch::saveMatchFile(szMchFile, pairwise_matching, feat_list);
	if (this->opt_.debug_output)
	{
		fea::fea_mch::saveMatchFileTxt(szMchFile + ".txt", pairwise_matching, feat_list, &features);
	}

	/* Clear */
	features.clear();
	pairwise_matching.clear();

	return EXIT_SUCCESS;
}

void filter_features(fea::FeatureSet const & feature, int const nScale, tri::contour_list const & poly,
	fea::FeatureSet & feature_filter, remap_feature_id & id_map)
{
	/* init the filter feat*/
	feature_filter.width = feature.width;
	feature_filter.height = feature.height;
	feature_filter.scale = feature.scale;

	/* get some para will be used */
	std::size_t sift_fea_num = feature.sift_descriptors.size();
	std::size_t surf_fea_num = feature.surf_descriptors.size();
	std::size_t harris_fea_num = feature.harris_descriptors.size();
	float const fwidth = static_cast<float>(feature.width);
	float const fheight = static_cast<float>(feature.height);
	float const fnorm = std::max(fwidth, fheight);

	/* get poly boundry */
	vertexd pmin, pmax;
	tri::get_boundary(poly, pmin, pmax);

	/* filter the sift features */
	vertexd pt;
	int cur_filtered_id = 0;
	for (std::size_t i = 0; i < sift_fea_num; i++)
	{
		math::Vec2f const & pos = feature.positions[i];
		pt.x = pos[0] * nScale; // *fnorm + (fwidth - 1.0) / 2.0;  pt.x *= nScale;
		pt.y = pos[1] * nScale; // *fnorm + (fheight - 1.0) / 2.0; pt.y *= nScale;

		if (pt.x < pmin.x || pt.x > pmax.x || pt.y < pmin.y || pt.y > pmax.y) continue;
		if (tri::is_pt_in_poly(poly, pt))
		{
			/* emplace back the postion, colors, descriptors */
			feature_filter.positions.emplace_back(pos);
			feature_filter.colors.emplace_back(feature.colors[i]);
			feature_filter.sift_descriptors.emplace_back(feature.sift_descriptors[i]);

			/* save the id map */
			id_map.emplace_back(cur_filtered_id, i);
			cur_filtered_id++;
		}
	}

	/* filter the surf features */
	int nOffset = sift_fea_num;
	for (std::size_t i = 0; i < surf_fea_num; i++)
	{
		math::Vec2f const & pos = feature.positions[i + nOffset];
		pt.x = pos[0] * nScale; // * fnorm + (fwidth - 1.0) / 2.0;  pt.x *= nScale;
		pt.y = pos[1] * nScale; // * fnorm + (fheight - 1.0) / 2.0; pt.y *= nScale;

		if (pt.x < pmin.x || pt.x > pmax.x || pt.y < pmin.y || pt.y > pmax.y) continue;
		if (tri::is_pt_in_poly(poly, pt))
		{
			/* emplace back the postion, colors, descriptors */
			feature_filter.positions.emplace_back(pos);
			feature_filter.colors.emplace_back(feature.colors[i + nOffset]);
			feature_filter.surf_descriptors.emplace_back(feature.surf_descriptors[i]);

			/* save the id map */
			id_map.emplace_back(cur_filtered_id, i + nOffset);
			cur_filtered_id++;
		}
	}

	/* filter the harris features. */
	nOffset = sift_fea_num + surf_fea_num;
	for (std::size_t i = 0; i < harris_fea_num; i++)
	{
		math::Vec2f const & pos = feature.positions[i + nOffset];
		pt.x = pos[0]; //* nScale; // * fnorm + (fwidth - 1.0) / 2.0;  pt.x *= nScale;
		pt.y = pos[1]; //* nScale; // * fnorm + (fheight - 1.0) / 2.0; pt.y *= nScale;

		if (pt.x < pmin.x || pt.x > pmax.x || pt.y < pmin.y || pt.y > pmax.y) continue;
		if (tri::is_pt_in_poly(poly, pt))
		{
			/* emplace back the postion, colors, descriptors */
			feature_filter.positions.emplace_back(pos);
			feature_filter.colors.emplace_back(feature.colors[i + nOffset]);
			feature_filter.harris_descriptors.emplace_back(feature.harris_descriptors[i]);

			/* save the id map */
			id_map.emplace_back(cur_filtered_id, i + nOffset);
			cur_filtered_id++;
		}
	}
}

int FeatureMatcher::Match(std::vector<std::string> szFeaList, std::string szRegionFile, std::string szMchFile, std::size_t expandSize)
{
	std::cout << std::endl;
	/* Read region file */
	tri::triinfo region_poly;
	bool bHasValidRegion = true;
	if (!tri::read_tri_idx(szRegionFile.c_str(), region_poly) || (int)expandSize == -1)
	{
		bHasValidRegion = false;
		std::cout << "** NOTICE: The region file is not valid and will not be used !" << std::endl;
	}
	if (bHasValidRegion)
	{
		if (region_poly.size() != szFeaList.size())
		{
			std::cerr << "*** Error: The region file is not correspond to the feature list !" << std::endl; return EXIT_FAILURE;
		}

		/* expand the region poly */
		tri::triinfo tmp;
		tmp.assign(region_poly.begin(), region_poly.end());
		for (std::size_t i = 0; i < region_poly.size(); i++)
		{
			tri::poly_expand(tmp[i].contour, region_poly[i].contour, expandSize);
		}
	}

	/* Read feature file */
	fea::FeatureSets features;
	std::vector<int> Widths, Heights;
	std::vector<int> Scales;
	features.resize(szFeaList.size());
	Widths.resize(szFeaList.size()); Heights.resize(szFeaList.size());
	Scales.resize(szFeaList.size());

#pragma omp parallel for schedule(dynamic)
#ifdef _MSC_VER
	for (int64_t i = 0; i < szFeaList.size(); i++)
#else
	for (std::size_t i = 0; i < szFeaList.size(); i++)
#endif
	{
		if (EXIT_FAILURE == fea::fea_mch::readFeatureFile(szFeaList[i], Scales[i], Widths[i], Heights[i], features[i]))
		{
			std::cerr << "*** Error: Read feature file :" << szFeaList[i].c_str() << " error !" << std::endl;
			return EXIT_FAILURE;
		}
	}

	/* Construct region if needed */
	if (!bHasValidRegion)
	{
		region_poly.clear();
		region_poly.resize(szFeaList.size());

#pragma omp parallel for schedule(dynamic)
#ifdef _MSC_VER
		for (int64_t i = 0; i < szFeaList.size(); i++)
#else
		for (std::size_t i = 0; i < szFeaList.size(); i++)
#endif
		{
			tri::tri_contour contour;

			contour.id = i;

			std::string name = szFeaList[i];
			int len = name.length();
			int base = name.rfind('\\');
			name = name.substr(base + 1, len - base);		
			strcpy(contour.name, name.c_str());

			contour.vertex_num = 4;
			contour.contour.emplace_back(0.0, 0.0);
			contour.contour.emplace_back((double)Widths[i] * Scales[i], 0.0);
			contour.contour.emplace_back( (double)Widths[i] * Scales[i], (double)Heights[i] * Scales[i]);
			contour.contour.emplace_back(0.0, (double)Heights[i] * Scales[i] );

			region_poly[i] = contour;
		}
	}
	
	/* do filter */
	fea::FeatureSets features_filter;
	std::vector<remap_feature_id> fea_id_maps;
	features_filter.resize(features.size());
	fea_id_maps.resize(features.size());

	std::size_t minFea = 99999999;
#pragma omp parallel for schedule(dynamic)
#ifdef _MSC_VER
	for (int64_t i = 0; i < features.size(); i++)
#else
	for (std::size_t i = 0; i < features.size(); i++)
#endif
	{
		util::WallTimer timer;
		
		/* do anti coordinate transe. */
		{
			float fwidth = static_cast<float>(features[i].width);
			float fheight = static_cast<float>(features[i].height);
			float fnorm = std::max(fwidth, fheight);
			int scale = Scales[i];
			int sift_fea_num = features[i].sift_descriptors.size();
			int harris_fea_num = features[i].harris_descriptors.size();

			if (this->opt_.anti_normalize_coord)
			{			
				for (std::size_t j = 0; j < sift_fea_num; j++)
				{
					math::Vec2f& pos = features[i].positions[j];
					pos[0] = pos[0] * fnorm + (fwidth - 1.0) / 2.0;  // pos[0] *= Scales[i];
					pos[1] = pos[1] * fnorm + (fheight - 1.0) / 2.0; // pos[1] *= Scales[i];
				}

				for (std::size_t j = 0; j < harris_fea_num; j++)
				{
					math::Vec2f& pos = features[i].positions[j + sift_fea_num];
					pos[0] = pos[0] * fnorm * scale + (fwidth * scale - 1.0) / 2.0;  // pos[0] *= Scales[i];
					pos[1] = pos[1] * fnorm * scale + (fheight * scale - 1.0) / 2.0; // pos[1] *= Scales[i];
				}
			}

			if (this->opt_.reverse_y_coord)
			{
				for (std::size_t j = 0; j < sift_fea_num; j++)
				{
					math::Vec2f& pos = features[i].positions[j];
					pos[1] = fheight - 1 - pos[1];
				}

				for (std::size_t j = 0; j < harris_fea_num; j++)
				{
					math::Vec2f& pos = features[i].positions[j + sift_fea_num];
					pos[1] = fheight * scale - 1 - pos[1];
				}
			}
		}	

		filter_features(features[i], Scales[i], region_poly[i].contour,
			features_filter[i], fea_id_maps[i]);

		/* Must do Normalizztion, or it will fail in outlier detection process */	
		float fwidth = static_cast<float>(features_filter[i].width);
		float fheight = static_cast<float>(features_filter[i].height);
		float fnorm = std::max(fwidth, fheight);
		int scale = Scales[i];
		int sift_fea_num = features_filter[i].sift_descriptors.size();
		int harris_fea_num = features_filter[i].harris_descriptors.size();
		for (std::size_t j = 0; j < sift_fea_num; j++)
		{
			math::Vec2f& pos = features_filter[i].positions[j];
			pos[0] = (pos[0] + 0.5f - fwidth / 2.0f) / fnorm;
			pos[1] = (pos[1] + 0.5f - fheight / 2.0f) / fnorm;
		}
		for (std::size_t j = 0; j < harris_fea_num; j++)
		{
			math::Vec2f& pos = features_filter[i].positions[j + sift_fea_num];
			pos[0] = (pos[0] + 0.5f - fwidth * scale / 2.0f) / (fnorm * scale);
			pos[1] = (pos[1] + 0.5f - fheight * scale / 2.0f) / (fnorm * scale);
		}
#pragma omp critical
		{
			std::cout << "\rGet"
				<< util::string::get_filled(features_filter[i].get_feature_num(), 5, ' ') << " ( / "
				<< util::string::get_filled(features[i].get_feature_num(), 5, ' ') << " ) features"
				<< ", took " << timer.get_elapsed() << " ms." << std::endl;

			minFea = std::min(minFea, features_filter[i].get_feature_num());
		}
	}

	/* empty some vector to save memory */
	if (!this->opt_.debug_output) fea::fea_mch::ClearVector<fea::FeatureSet>(features);
	fea::fea_mch::ClearVector<int>(Widths); fea::fea_mch::ClearVector<int>(Heights); fea::fea_mch::ClearVector<int>(Scales);

	/* do match */
	{
		/* Set Options */
		fea::fea_mch::dpg_mch::Matching::Options matching_opts;
		matching_opts.ransac_opts.verbose_output = false;
		matching_opts.match_num_previous_frames = 0; // disabled 
		matching_opts.min_feature_matches = this->opt_.min_feature_matches;
		matching_opts.min_matching_inliers = (int)(this->opt_.min_feature_matches / 5.0);
		matching_opts.use_lowres_matching = !bHasValidRegion;
		matching_opts.num_lowres_features = 1000;
		matching_opts.matcher_type = (fea::fea_mch::dpg_mch::Matching::MatcherType)this->opt_.matcher_type;

		/* Compute match */
		fea::fea_mch::dpg_mch::Matching dpg_matching(matching_opts);
		dpg_matching.init(&features_filter);
		fea::PairwiseMatching pairwise_matching;
		dpg_matching.compute2(&pairwise_matching);
		if (pairwise_matching.size() < 1) { std::cerr << "*** Error: Wrong matching image pairs. Exiting." << std::endl; return EXIT_FAILURE; }

		/* Remap the Feature ID */
		{
#pragma omp parallel for schedule(dynamic)
#ifdef _MSC_VER
			for (int64_t i = 0; i < pairwise_matching.size(); i++)
#else
			for (std::size_t i = 0; i < pairwise_matching.size(); i++)
#endif
			{
				fea::TwoViewMatching & two_mch = pairwise_matching[i];
				if (two_mch.view_1_id != 0)
				{
					continue;
				}

				fea::CorrespondenceIndices & cor_indi = two_mch.matches;
				remap_feature_id & id_map1 = fea_id_maps[two_mch.view_1_id];
				remap_feature_id & id_map2 = fea_id_maps[two_mch.view_2_id];

				for (fea::CorrespondenceIndices::iterator it = cor_indi.begin(); it != cor_indi.end(); )
				{
					const int* val1 = math::algo::binary_search(id_map1, it->first); /* Remap the first id */
					const int* val2 = math::algo::binary_search(id_map2, it->second); /* Remap the second id */
					if (val1 != nullptr && val2 != nullptr)
					{
						it->first = *(int*)val1;
						it->second = *(int*)val2;
						it++;
					}
					else
					{
						it = cor_indi.erase(it); /* if not found, then erase this*/
					}
				}

				/* remap the view id according to the region file */
				/*if (!this->opt_.debug_output)
				{
					two_mch.view_1_id = region_poly[two_mch.view_1_id].id;
					two_mch.view_2_id = region_poly[two_mch.view_2_id].id;
				}*/		
			}
		}

		/* Save result */		
		if (this->opt_.debug_output)
		{
			if (this->opt_.matcher_type == MATCHER_HARRIS)
			{
				fea::fea_mch::saveHarrisMatchTxt(szMchFile + ".txt", pairwise_matching, szFeaList, &features);
			}
			else
			{
				fea::fea_mch::saveMatchFileTxt(szMchFile + ".txt", pairwise_matching, szFeaList, &features);
			}
			
			/* clear */
			fea::fea_mch::ClearVector<fea::FeatureSet>(features); // features.clear();
		}	

		//if (this->opt_.debug_output)
		//{
		//	/* We need remap the id here */
		//	for (int64_t i = 0; i < pairwise_matching.size(); i++)
		//	{
		//		fea::TwoViewMatching & two_mch = pairwise_matching[i];
		//		if (two_mch.view_1_id != 0)
		//		{
		//			continue;
		//		}

		//		two_mch.view_1_id = region_poly[two_mch.view_1_id].id;
		//		two_mch.view_2_id = region_poly[two_mch.view_2_id].id;
		//	}
		//}
		fea::fea_mch::saveMatchFile(szMchFile, pairwise_matching, szFeaList);

		/* Clear */
		pairwise_matching.clear();
	}
	
	/* Clear */
	features_filter.clear();
	fea_id_maps.clear();

	return EXIT_SUCCESS;
}

