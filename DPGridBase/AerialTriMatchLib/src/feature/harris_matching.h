#pragma once

#ifndef FEA_HARRIS_MATCHING_HEADER
#define FEA_HARRIS_MATCHING_HEADER

#include "feature/matching_base.h"
#include "feature/sift.h"
#include "feature/harris.h"

/* Whether to use floating point or 8-bit descriptors for matching. */
#define DISCRETIZE_DESCRIPTORS 1

FEA_NAMESPACE_BEGIN
FEA_MATCH_NAMESPACE_BEGIN

class HarrisMatching :public MatchingBase
{
public:
	enum SiftMatcherType
	{
		SIFT_MATCH_EXHAUSTIVE = 0x0,
		SIFT_MATCH_CASCADE_HASHING = 0x1
	};

	struct Options
	{
		SiftMatcherType sift_matcher_type = SIFT_MATCH_CASCADE_HASHING;

		/** For hashing. */
		/** Number of bucket groups. */
		uint8_t num_bucket_groups = 6;

		/** Number of buckets (2^num_bucket_bits) per bucket group. */
		uint8_t num_bucket_bits = 8;

		/** Minimum number of top ranked candidates to collect. */
		uint16_t min_num_candidates = 6;

		/** Maximum number of top ranked candidates to collect. */
		uint16_t max_num_candidates = 10;

		/** For Harris Matching. */
		/** The num of neighbor features used to do harris matching. */
		int near_feat_num = 3;

		/*% The max distance between neighbor and estimated location. */
		float max_near_distance = 25.0f; 

		bool debug_output = true;
	};

public:
	//HarrisMatching(Options option);
	~HarrisMatching(void) override = default;

	void set_option(Options option);

	/** Initialize matcher by preprocessing given SIFT/SURF features. */
	void init(fea::FeatureSets* featuresets) override;

	/** Matches all feature types yielding a single matching result. */
	void pairwise_match(int view_1_id, int view_2_id,
		Matching::Result* result) const override;

	/**
	* Matches the N lowest resolution features and returns the number of
	* matches. Can be used as a guess for full matchability. Useful values
	* are at most 3 matches for 500 features, or 2 matches with 300 features.
	*/
	int pairwise_match_lowres(int view_1_id, int view_2_id,
		std::size_t num_features) const override;

protected:
	typedef util::AlignedMemory<math::Vec128us, 16> HarrisDescriptors;
	typedef util::AlignedMemory<math::Vec2f, 16> SiftFeaturePositions;
	typedef util::AlignedMemory<math::Vec2f, 16> HarrisFeaturePositions;
#if DISCRETIZE_DESCRIPTORS
	typedef util::AlignedMemory<math::Vec128us, 16> SiftDescriptors;
#else
	typedef util::AlignedMemory<math::Vec128f, 16> SiftDescriptors;
#endif

	/** Internal initialization methods for SIFT/SURF features. */
	void init_sift(SiftDescriptors* dst, SiftFeaturePositions* posi, Sift::Descriptors const& src, math::Vec2f const * pos);
	void init_harris(HarrisDescriptors* dst, HarrisFeaturePositions* posi, Harris::Descriptors const& src, math::Vec2f const * pos);

	struct ProcessedFeatureSet
	{
		SiftDescriptors sift_descr;
		SiftFeaturePositions sift_position;
		HarrisDescriptors harris_descr;
		HarrisFeaturePositions harris_position;
	};
	typedef std::vector<ProcessedFeatureSet> ProcessedFeatureSets;
	ProcessedFeatureSets processed_feature_sets;

	int sift_matching(int view_1_id, int view_2_id, Matching::Result * result) const;

	int get_affine_trans(SiftFeaturePositions const& pfs1, SiftFeaturePositions const& pfs2,
		std::vector<int> const& match_result,
		math::Vec6d * affine_trans_param) const;

	int affine_transition(math::Vec6d const & affine_trans_param,
		math::Vec2f const * src, math::Vec2f * dst, int size = 1) const;

	float ncc_harris(math::Vec128us const & des1, math::Vec128us const & des2, int window_len = 7) const;

	void harris_match_one_way(ProcessedFeatureSet const& pfs1, ProcessedFeatureSet const& pfs2,		
		std::vector<int> const& sift_result,
		std::vector<int>* harris_matches) const;

	void harris_match_two_way(ProcessedFeatureSet const& pfs1, ProcessedFeatureSet const& pfs2,
		Matching::Result const& sift_result,
		Matching::Result* matches) const;
private:
	struct LocalData;

	/**
	* Matches all elements in set 1 to all elements in set 2 and vice versa.
	* Works similarly to sfm::Matching::twoway_match().
	*/
	template <typename D>
	void twoway_match(Matching::Options const& matching_opts,
		LocalData const& set_1, LocalData const& set_2,
		D const& set_1_descs, D const& set_2_descs,
		Matching::Result* matches, Options const& cashash_opts) const;

	/**
	* Matches all elements in set 1 to all elements in set 2.
	* Works similarly to sfm::Matching::oneway_match().
	*/
	template <typename D>
	void oneway_match(Matching::Options const& matching_opts,
		LocalData const& set_1, LocalData const& set_2,
		D const& set_1_descs, D const& set_2_descs,
		std::vector<int>* result, Options const& cashash_opts) const;

	/**
	* Cascade hashing projection matrices. T shall be a vector with the same
	* number of elements as the feature descriptor, i.e. Vec128f for SIFT and
	* Vec64f for SURF.
	*/
	template <typename T>
	struct ProjMats
	{
		/**
		* Primary projection matrix. The number of vectors/Ts determines the
		* size of the bit vector that is computed and used for computing the
		* hamming distance later on. The size of prim_proj_mat is set to the
		* size of the feature descriptors
		* (see GlobalData::generate_proj_matrices()), i.e. 128 for SIFT and 64
		* for SURF.
		*/
		std::vector<T> prim_proj_mat;

		/**
		* Secondary projection matrices. The number of Ts is the number of
		* buckets per bucket group (see Options::num_bucket_bits).
		* The size of sec_proj_mats is the number of bucket groups.
		*/
		std::vector<std::vector<T>> sec_proj_mats;
	};

	typedef ProjMats<math::Vec128f> ProjMatsSift;

	/**
	* GlobalData contains the primary and secondary projection matrices
	* for each descriptor type. These are used to compute the hashes of each
	* feature and to which buckets a feature is assigned. It is global in the
	* sense that it is the same and used for all images.
	*/
	class GlobalData
	{
	public:
		void generate_proj_matrices(Options const& cashash_opts);

		ProjMatsSift sift;

	private:
		template <typename T>
		void generate_proj_matrices(std::vector<T>* prim_hash,
			std::vector<std::vector<T>>* sec_hash, Options const& cashash_opts);
	};

	typedef std::vector<std::size_t> BucketFeatureIDs;
	typedef std::vector<BucketFeatureIDs> BucketGroupFeatures;
	typedef std::vector<BucketGroupFeatures> BucketGroupsFeatures;

	typedef std::vector<uint16_t> BucketIDs;
	typedef std::vector<BucketIDs> BucketGroupsBuckets;

	/** Per-image local data. */
	struct LocalData
	{
		/**
		* Compressed hash data. 2x uint64_t for SIFT and 1x uint64_t for SURF.
		* Each uint64_t is a bit vector which represents the inverted signs of
		* the dot products of the feature vector with the vectors in
		* ProjMats::prim_proj_mat. It is compressed in the sense that the bits
		* are tightly packed into uint64_ts.
		*/
		std::vector<uint64_t> comp_hash_data;

		/**
		* Vector of bucket IDs to which a feature is assigned for each bucket
		* group, i.e. bucket_grps_bucket_ids[2][64] = 8 means that in bucket
		* group 2 feature 64 is assigned to bucket 8.
		*/
		BucketGroupsBuckets bucket_grps_bucket_ids;

		/**
		* Vector of feature IDs for each bucket of a bucket group, i.e.
		* bucket_grps_feature_ids[2][8][4] = 6 means that the 4th feature in
		* bucket 8 of bucket group 2 is the feature with the ID 6.
		*/
		BucketGroupsFeatures bucket_grps_feature_ids;
	};

	/** Compute cascade hashing data from SIFT and SURF descriptors. */
	void compute(LocalData* ld_sift,
		std::vector<math::Vec128f> const& sift_zero_mean_descs,
		GlobalData const& cashash_global_data, Options const& cashash_opts);

	/** Compute average descriptors of given feature sets. */
	void compute_avg_descriptors(ProcessedFeatureSets const& pfs,
		math::Vec128f* sift_avg);

	/** Compute zero mean descriptors using given average vectors. */
	void compute_zero_mean_descs(
		std::vector<math::Vec128f>* sift_zero_mean_descs,
		SiftDescriptors const& sift_descs,
		math::Vec128f const& sift_avg);

	/** Compute cascade hashes of zero mean SIFT/SURF descriptors. */
	template <typename T>
	void compute_cascade_hashes(std::vector<T> const& zero_mean_descs,
		std::vector<uint64_t>* comp_hash_data,
		BucketGroupsBuckets* bucket_grps_bucket_ids,
		std::vector<T> const& prim_proj_mat,
		std::vector<std::vector<T>> const& sec_proj_mats, Options const& cashash_opts);

	/** Assign features to a bucket for each bucket group. */
	void build_buckets(BucketGroupsFeatures* bucket_grps_feature_ids,
		BucketGroupsBuckets const& bucket_grps_bucket_ids, size_t num_descs,
		Options const& cashash_opts);

	/**
	* Collect features from the bucket from each bucket group to which
	* feature_id is assigned and group them according to their hamming
	* distance to the feature.
	*/
	template <int DIMHASH>
	void collect_features_from_buckets(
		std::vector<std::vector<uint32_t>>* grouped_features,
		size_t feature_id, std::vector<bool>* data_index_used,
		BucketGroupsBuckets const& bucket_grps_bucket_ids,
		BucketGroupsFeatures const& bucket_grps_feature_ids,
		uint64_t const* comp_hash_data1,
		std::vector<uint64_t> const& comp_hash_data2) const;

	/**
	* Collect the features with the smallest hamming distance from
	* grouped_features into top_candidates.
	*/
	void collect_top_ranked_candidates(std::vector<uint32_t>* top_candidates,
		std::vector<std::vector<uint32_t>> const& grouped_features,
		uint8_t dim_hash_data, uint16_t min_num_candidates,
		uint16_t max_num_candidates) const;

private:
	Options opt_;

	GlobalData global_data;
	std::vector<LocalData> local_data_sift;	

};

template <typename T>
void HarrisMatching::GlobalData::generate_proj_matrices(
	std::vector<T>* prim_hash, std::vector<std::vector<T>>* sec_hash,
	Options const& cashash_opts)
{
	int const dim_desc = T::dim;
	int const dim_hash_data = dim_desc;

	prim_hash->resize(dim_hash_data);

	/* Setup PRNG. */
	std::mt19937 prng_mt(0);
	std::normal_distribution<> dis(0, 1);

	/* Generate values for primary hashing function. */
	for (uint8_t i = 0; i < dim_hash_data; i++)
		for (uint8_t j = 0; j < dim_desc; j++)
			(*prim_hash)[i][j] = dis(prng_mt);

	uint8_t const num_bucket_groups = cashash_opts.num_bucket_groups;
	uint8_t const num_bucket_bits = cashash_opts.num_bucket_bits;

	sec_hash->resize(num_bucket_groups, std::vector<T>(num_bucket_bits));

	/* Generate values for secondary hashing function. */
	for (uint8_t group_idx = 0; group_idx < num_bucket_groups; group_idx++)
		for (uint8_t i = 0; i < num_bucket_bits; i++)
			for (uint8_t j = 0; j < dim_desc; j++)
				(*sec_hash)[group_idx][i][j] = dis(prng_mt);
}

/* ---------------------------------------------------------------- */
template<typename T>
inline void
HarrisMatching::compute_cascade_hashes(std::vector<T> const & zero_mean_descs, std::vector<uint64_t>* comp_hash_data, BucketGroupsBuckets * bucket_grps_bucket_ids, std::vector<T> const & prim_proj_mat, std::vector<std::vector<T>> const & sec_proj_mats, Options const & cashash_opts)
{
	int const dim_desc = T::dim;
	int const dim_hash_data = dim_desc;
	uint8_t const dim_comp_hash_data = dim_hash_data / 64;
	uint8_t const num_bucket_bits = cashash_opts.num_bucket_bits;
	uint8_t const num_bucket_grps = cashash_opts.num_bucket_groups;

	/* Allocate memory. */
	size_t num_descs = zero_mean_descs.size();
	comp_hash_data->resize(num_descs * dim_comp_hash_data);
	bucket_grps_bucket_ids->resize(num_bucket_grps, std::vector<uint16_t>(num_descs));

	for (size_t i = 0; i < num_descs; i++)
	{
		T const& desc = zero_mean_descs[i];

		/* Compute hash code. */
		for (uint8_t j = 0; j < dim_comp_hash_data; j++)
		{
			uint64_t comp_hash = 0;
			uint8_t data_start = j * 64;
			uint8_t data_end = (j + 1) * 64;
			for (uint8_t k = data_start; k < data_end; k++)
			{
				T const& proj_vec = prim_proj_mat[k];
				float sum = desc.dot(proj_vec);
				comp_hash = (comp_hash << 1) | (sum > 0.0f);
			}
			(*comp_hash_data)[i * dim_comp_hash_data + j] = comp_hash;
		}

		/* Determine the descriptor's bucket index for each bucket group. */
		for (uint8_t grp_idx = 0; grp_idx < num_bucket_grps; grp_idx++)
		{
			uint16_t bucket_id = 0;
			for (uint8_t bit_idx = 0; bit_idx < num_bucket_bits; bit_idx++)
			{
				T const& proj_vec = sec_proj_mats[grp_idx][bit_idx];
				float sum = desc.dot(proj_vec);
				bucket_id = (bucket_id << 1) | (sum > 0.0f);
			}

			(*bucket_grps_bucket_ids)[grp_idx][i] = bucket_id;
		}
	}
}

template <typename D>
inline void
HarrisMatching::twoway_match(Matching::Options const& matching_opts,
	LocalData const& set_1, LocalData const& set_2,
	D const& set_1_descs, D const& set_2_descs,
	Matching::Result* matches, Options const& cashash_opts) const
{
	oneway_match(matching_opts, set_1, set_2, set_1_descs, set_2_descs,
		&matches->matches_1_2,
		cashash_opts);
	oneway_match(matching_opts, set_2, set_1, set_2_descs, set_1_descs,
		&matches->matches_2_1,
		cashash_opts);
}

template <typename D>
inline void
HarrisMatching::oneway_match(Matching::Options const& matching_opts,
	LocalData const& set_1, LocalData const& set_2,
	D const& set_1_descs, D const& set_2_descs,
	std::vector<int>* result, Options const& cashash_opts) const
{
	typedef typename D::value_type V;
	typedef typename V::ValueType T;

	size_t set_1_size = set_1_descs.size();
	size_t set_2_size = set_2_descs.size();

	if (set_1_size == 0 || set_2_size == 0)
		return;

	float const square_lowe_thres = MATH_POW2(matching_opts.lowe_ratio_threshold);
	float const square_dist_thres = MATH_POW2(matching_opts.distance_threshold);

	uint16_t const min_num_candidates = cashash_opts.min_num_candidates;
	uint16_t const max_num_candidates = cashash_opts.max_num_candidates;
	int const descriptor_length = V::dim;
	uint8_t const dim_hash_data = descriptor_length;
	uint32_t const dim_comp_hash_data = dim_hash_data / 64;

	result->resize(set_1_size, -1);
	std::vector<bool> data_index_used(set_2_size);
	std::vector<std::vector<uint32_t> > grouped_features(dim_hash_data + 1);
	std::vector<uint32_t> top_candidates;

	top_candidates.reserve(max_num_candidates);

	std::unique_ptr<T> tmp(new T[max_num_candidates * descriptor_length]);
	NearestNeighbor<T> nn;
	nn.set_elements(tmp.get());
	nn.set_element_dimensions(descriptor_length);

	for (size_t i = 0; i < set_1_size; i++)
	{
		std::fill(data_index_used.begin(), data_index_used.end(), false);

		for (size_t j = 0; j < grouped_features.size(); j++)
			grouped_features[j].clear();

		top_candidates.clear();

		/* Fetch candidate features from the buckets in each group. */
		collect_features_from_buckets<dim_hash_data>(
			&grouped_features,
			i,
			&data_index_used,
			set_1.bucket_grps_bucket_ids,
			set_2.bucket_grps_feature_ids,
			&set_1.comp_hash_data[i * dim_comp_hash_data],
			set_2.comp_hash_data);

		/* Add closest candidates by Hamming distance to top_candidates vector. */
		collect_top_ranked_candidates(
			&top_candidates,
			grouped_features,
			dim_hash_data,
			min_num_candidates,
			max_num_candidates);

		/* Copy top candidates' descriptors into a contiguous array. */
		for (size_t j = 0; j < top_candidates.size(); j++)
		{
			uint32_t candidate_id = top_candidates[j];
			std::memcpy(tmp.get() + j * descriptor_length,
				set_2_descs[candidate_id].begin(),
				sizeof(V));
		}

		typename NearestNeighbor<T>::Result nn_result;
		nn.set_num_elements(top_candidates.size());
		nn.find(set_1_descs[i].begin(), &nn_result);

		if (nn_result.dist_1st_best > square_dist_thres)
			continue;

		if (static_cast<float>(nn_result.dist_1st_best)
			/ static_cast<float>(nn_result.dist_2nd_best)
		> square_lowe_thres)
			continue;

		result->at(i) = top_candidates[nn_result.index_1st_best];
	}
}

template <int DIMHASH>
inline void
HarrisMatching::collect_features_from_buckets(
	std::vector<std::vector<uint32_t>>* grouped_features,
	size_t feature_id, std::vector<bool>* data_index_used,
	BucketGroupsBuckets const& bucket_grps_bucket_ids,
	BucketGroupsFeatures const& bucket_grps_feature_ids,
	uint64_t const* comp_hash_data1,
	std::vector<uint64_t> const& comp_hash_data2) const
{
	size_t num_bucket_grps = bucket_grps_bucket_ids.size();
	for (size_t grp_idx = 0; grp_idx < num_bucket_grps; grp_idx++)
	{
		uint8_t bucket_id = bucket_grps_bucket_ids[grp_idx][feature_id];
		BucketFeatureIDs const& bucket_feature_ids = bucket_grps_feature_ids[grp_idx][bucket_id];

		for (size_t j = 0; j < bucket_feature_ids.size(); j++)
		{
			size_t candidate_id = bucket_feature_ids[j];
			if ((*data_index_used)[candidate_id])
				continue;

			uint64_t const *ptr2 = &comp_hash_data2[candidate_id * DIMHASH / 64];
			int hamming_dist = 0;
			for (int k = 0; k < (DIMHASH / 64); k++)
				hamming_dist += _mm_popcnt_u64(comp_hash_data1[k] ^ ptr2[k]);

			(*grouped_features)[hamming_dist].emplace_back(candidate_id);
			(*data_index_used)[candidate_id] = true;
		}
	}
}

inline void
HarrisMatching::collect_top_ranked_candidates(
	std::vector<uint32_t>* top_candidates,
	std::vector<std::vector<uint32_t>> const& grouped_features,
	uint8_t dim_hash_data, uint16_t min_num_candidates,
	uint16_t max_num_candidates) const
{
	for (uint16_t hash_dist = 0; hash_dist <= dim_hash_data; hash_dist++)
	{
		for (size_t j = 0; j < grouped_features[hash_dist].size(); j++)
		{
			uint32_t idx2 = grouped_features[hash_dist][j];
			top_candidates->emplace_back(idx2);

			if (top_candidates->size() >= max_num_candidates)
				break;
		}

		if (top_candidates->size() >= min_num_candidates)
			break;
	}
}

FEA_MATCH_NAMESPACE_END
FEA_NAMESPACE_END

#endif

