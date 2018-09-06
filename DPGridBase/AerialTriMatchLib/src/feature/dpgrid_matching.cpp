#include "dpgrid_matching.h"

#include "feature\exhaustive_matching.h"
#include "feature\cascade_hashing.h"
#include "feature\harris_matching.h"

FEA_NAMESPACE_BEGIN
FEA_MATCH_NAMESPACE_BEGIN
DPG_MATCH_NAMESPACE_BEGIN

Matching::Matching(Options const& options, Progress* progress)
	: opts_(options)
	, progress_(progress)
{
	switch (this->opts_.matcher_type)
	{
	case MATCHER_EXHAUSTIVE:
		this->matcher_.reset(new ExhaustiveMatching());
		break;
	case MATCHER_CASCADE_HASHING:
		this->matcher_.reset(new CascadeHashing());
		break;
	case MATCHER_HARRIS:
		this->matcher_.reset(new HarrisMatching());
		break;
	default:
		throw std::runtime_error("Unhandled matcher type");
	}
}

void
Matching::init(FeatureSets* featuresets)
{
	if (featuresets == nullptr)
		throw std::invalid_argument("Viewports must not be null");

	this->featuresets_ = featuresets;
	this->matcher_->init(featuresets);

	/* Free descriptors. */
	for (std::size_t i = 0; i < featuresets->size(); i++)
		featuresets->at(i).clear_descriptors();
}

void
Matching::compute(PairwiseMatching* pairwise_matching)
{
	if (this->featuresets_ == nullptr)
		throw std::runtime_error("Viewports must not be null");

	std::size_t num_viewports = this->featuresets_->size();
	std::size_t num_pairs = num_viewports * (num_viewports - 1) / 2;
	std::size_t num_done = 0;

	if (this->progress_ != nullptr)
	{
		this->progress_->num_total = num_pairs;
		this->progress_->num_done = 0;
	}

#pragma omp parallel for schedule(dynamic)
#ifdef _MSC_VER
	for (int64_t i = 0; i < num_pairs; ++i)
#else
	for (std::size_t i = 0; i < num_pairs; ++i)
#endif
	{
#pragma omp critical
		{
			num_done += 1;
			if (this->progress_ != nullptr)
				this->progress_->num_done += 1;

			float percent = (num_done * 1000 / num_pairs) / 10.0f;
			std::cout << "\rMatching pair " << num_done << " of "
				<< num_pairs << " (" << percent << "%)..." << std::flush;
		}

		int const view_1_id = (int)(0.5 + std::sqrt(0.25 + 2.0 * i));
		int const view_2_id = (int)i - view_1_id * (view_1_id - 1) / 2;
		if (this->opts_.match_num_previous_frames != 0
			&& view_2_id + this->opts_.match_num_previous_frames < view_1_id)
			continue;

		FeatureSet const& view_1 = this->featuresets_->at(view_1_id);
		FeatureSet const& view_2 = this->featuresets_->at(view_2_id);
		if (view_1.positions.empty() || view_2.positions.empty())
			continue;

		/* Match the views. */
		util::WallTimer timer;
		std::stringstream message;
		CorrespondenceIndices matches;
		this->two_view_matching(view_1_id, view_2_id, &matches, message);
		std::size_t matching_time = timer.get_elapsed();

		if (matches.empty())
		{
#pragma omp critical
			std::cout << "\rPair (" << view_1_id << ","
				<< view_2_id << ") rejected, "
				<< message.str() << std::endl;
			continue;
		}

		/* Successful two view matching. Add the pair. */
		TwoViewMatching matching;
		matching.view_1_id = view_1_id;
		matching.view_2_id = view_2_id;
		std::swap(matching.matches, matches);

#pragma omp critical
		{
			pairwise_matching->push_back(matching);
			std::cout << "\rPair (" << view_1_id << ","
				<< view_2_id << ") matched, " << matching.matches.size()
				<< " inliers, took " << matching_time << " ms." << std::endl;
		}
	}

	std::cout << "\rFound a total of " << pairwise_matching->size()
		<< " matching image pairs." << std::endl;
}

void Matching::compute2(PairwiseMatching * pairwise_matching)
{
	if (this->featuresets_ == nullptr)
		throw std::runtime_error("Viewports must not be null");

	std::size_t num_featuress = this->featuresets_->size();
	std::size_t num_pairs = num_featuress - 1;
	std::size_t num_done = 0;

	if (this->progress_ != nullptr)
	{
		this->progress_->num_total = num_pairs;
		this->progress_->num_done = 0;
	}


#pragma omp parallel for schedule(dynamic)
#ifdef _MSC_VER
	for (int64_t i = 0; i < num_pairs; ++i)
#else
	for (std::size_t i = 0; i < num_pairs; ++i)
#endif
	{
#pragma omp critical
		{
			num_done += 1;
			if (this->progress_ != nullptr)
				this->progress_->num_done += 1;

			float percent = (num_done * 1000 / num_pairs) / 10.0f;
			std::cout << "\rMatching pair " << num_done << " of "
				<< num_pairs << " (" << percent << "%)..." << std::flush;
		}

		/* Matching the first featureset with the rest one */
		int const view_1_id = 0;
		int const view_2_id = (int)(i + 1);

		FeatureSet const& view_1 = this->featuresets_->at(view_1_id);
		FeatureSet const& view_2 = this->featuresets_->at(view_2_id);
		if (view_1.positions.empty() || view_2.positions.empty())
			continue;

		/* Match the views. */
		util::WallTimer timer;
		std::stringstream message;
		CorrespondenceIndices matches;
		this->two_view_matching(view_1_id, view_2_id, &matches, message);
		std::size_t matching_time = timer.get_elapsed();

		if (matches.empty())
		{
#pragma omp critical
			std::cout << "\rPair (" << view_1_id << ","
				<< view_2_id << ") rejected, "
				<< message.str() << std::endl;
			continue;
		}

		/* Successful two view matching. Add the pair. */
		TwoViewMatching matching;
		matching.view_1_id = view_1_id;
		matching.view_2_id = view_2_id;
		std::swap(matching.matches, matches);

#pragma omp critical
		{
			pairwise_matching->push_back(matching);
			std::cout << "\rPair (" << view_1_id << ","
				<< view_2_id << ") matched, " << matching.matches.size()
				<< " inliers, took " << matching_time << " ms." << std::endl;
		}
	}

	std::cout << "\rFound a total of " << pairwise_matching->size()
		<< " matching image pairs." << std::endl;
}

void
Matching::two_view_matching(int view_1_id, int view_2_id,
	CorrespondenceIndices* matches, std::stringstream& message)
{
	FeatureSet const& view_1 = this->featuresets_->at(view_1_id);
	FeatureSet const& view_2 = this->featuresets_->at(view_2_id);

	/* Low-res matching if number of features is large. */
	if (this->opts_.use_lowres_matching
		&& view_1.positions.size() * view_2.positions.size() > 1000000)
	{
		int const num_matches = this->matcher_->pairwise_match_lowres(view_1_id,
			view_2_id, this->opts_.num_lowres_features);
		if (num_matches < this->opts_.min_lowres_matches)
		{
			message << "only " << num_matches
				<< " of " << this->opts_.min_lowres_matches
				<< " low-res matches.";
			return;
		}
	}

	/* Perform two-view descriptor matching. */
	fea::fea_mch::Matching::Result matching_result;
	this->matcher_->pairwise_match(view_1_id, view_2_id, &matching_result);
	int num_matches = fea::fea_mch::Matching::count_consistent_matches(matching_result);

	//std::cout << "Match Num: " << num_matches << std::endl;
	/* Require at least 8 matches. Check threshold. */
	int const min_matches_thres = std::max(8, this->opts_.min_feature_matches);
	if (num_matches < min_matches_thres)
	{
		message << "matches below threshold of "
			<< min_matches_thres << ".";
		return;
	}

	/* Build correspondences from feature matching result. */
	fea::Correspondences2D2D unfiltered_matches;
	fea::CorrespondenceIndices unfiltered_indices;
	{
		std::vector<int> const& m12 = matching_result.matches_1_2;
		for (std::size_t i = 0; i < m12.size(); ++i)
		{
			if (m12[i] < 0)
				continue;

			fea::Correspondence2D2D match;
			match.p1[0] = view_1.positions[i][0];
			match.p1[1] = view_1.positions[i][1];
			match.p2[0] = view_2.positions[m12[i]][0];
			match.p2[1] = view_2.positions[m12[i]][1];
			unfiltered_matches.push_back(match);
			unfiltered_indices.push_back(std::make_pair(i, m12[i]));
		}
	}

	/* gravity-lize the position */
	{
		/*double gx1 = 0.0, gy1 = 0.0, gx2 = 0.0, gy2 = 0.0;
		int unfiltered_matches_num = unfiltered_matches.size();

		for (int i = 0; i < unfiltered_matches_num; i++)
		{
			fea::Correspondence2D2D const match = unfiltered_matches[i];
			gx1 += (match.p1[0] / (double)unfiltered_matches_num);
			gy1 += (match.p1[1] / (double)unfiltered_matches_num);
			gx2 += (match.p2[0] / (double)unfiltered_matches_num);
			gy2 += (match.p2[1] / (double)unfiltered_matches_num);
		}
		
		for (int i = 0; i < unfiltered_matches_num; i++)
		{
			fea::Correspondence2D2D & match = unfiltered_matches[i];
			match.p1[0] -= gx1;
			match.p1[1] -= gy1;
			match.p2[0] -= gx2;
			match.p2[1] -= gy2;
		}*/
	}

	/* Compute fundamental matrix using RANSAC. */
	fea::RansacFundamental::Result ransac_result;
	int num_inliers = 0;
	{
		fea::RansacFundamental ransac(this->opts_.ransac_opts);
		ransac.estimate(unfiltered_matches, &ransac_result);
		num_inliers = ransac_result.inliers.size();
	}

	/* Require at least 8 inlier matches. */
	int const min_inlier_thres = std::max(8, this->opts_.min_matching_inliers);
	if (num_inliers < min_inlier_thres)
	{
		message << "inliers below threshold of "
			<< min_inlier_thres << ".";
		return;
	}

	/* Create Two-View matching result. */
	matches->clear();
	matches->reserve(num_inliers);
	for (int i = 0; i < num_inliers; ++i)
	{
		int const inlier_id = ransac_result.inliers[i];
		matches->push_back(unfiltered_indices[inlier_id]);
	}
}

DPG_MATCH_NAMESPACE_END
FEA_MATCH_NAMESPACE_END
FEA_NAMESPACE_END
