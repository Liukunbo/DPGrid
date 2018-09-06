#ifndef DPGRID_MATCHING_HEADER
#define DPGRID_MATCHING_HEADER

#include "feature\defines.h"
#include "feature\matching_base.h"
#include "feature\dpgrid_common.h"
#include "feature\ransac_fundamental.h"

FEA_NAMESPACE_BEGIN
FEA_MATCH_NAMESPACE_BEGIN
DPG_MATCH_NAMESPACE_BEGIN

class Matching
{
public:
	enum MatcherType
	{
		MATCHER_EXHAUSTIVE = 0x0,
		MATCHER_CASCADE_HASHING = 0x1,
		MATCHER_HARRIS = 0x2
	};

	/** Options for feature matching. */
	struct Options
	{
		/** Options for RANSAC computation of the fundamental matrix. */
		fea::RansacFundamental::Options ransac_opts;
		/** Minimum number of matching features before RANSAC. */
		int min_feature_matches = 24;
		/** Minimum number of matching features after RANSAC. */
		int min_matching_inliers = 12;
		/** Perform low-resolution matching to reject unlikely pairs. */
		bool use_lowres_matching = false;
		/** Number of features used for low-res matching. */
		int num_lowres_features = 500;
		/** Minimum number of matches from low-res matching. */
		int min_lowres_matches = 5;
		/** Only match to a few previous frames. Disabled by default. */
		int match_num_previous_frames = 0;
		/** Matcher type. Exhaustive by default. */
		MatcherType matcher_type = MATCHER_EXHAUSTIVE;
	};

	struct Progress
	{
		std::size_t num_total;
		std::size_t num_done;
	};

public:
	explicit Matching(Options const& options, Progress* progress = nullptr);

	/**
	 * Initialize matching by passing features to the matcher for
	 * preprocessing.
	 */
	void init(FeatureSets* featuresets);

	/**
	 * Computes the pairwise matching between all pairs of views.
	 * Computation requires both descriptor data and 2D feature positions
	 * in the viewports.
	 */
	void compute(PairwiseMatching* pairwise_matching);

	/**
	 * Computes the pairwise matching between the first view and the remained views.
	 * Computation requires both descriptor data and 2D feature positions
	 * in the viewports.
	 */
	void compute2(PairwiseMatching* pairwise_matching);


private:
	void two_view_matching(int view_1_id, int view_2_id,
		CorrespondenceIndices* matches, std::stringstream& message);

private:
	Options opts_;
	Progress* progress_;
	std::unique_ptr<MatchingBase> matcher_;
	FeatureSets const* featuresets_;
};

DPG_MATCH_NAMESPACE_END
FEA_MATCH_NAMESPACE_END
FEA_NAMESPACE_END


#endif // DPGRID_MATCHING_HEADER

