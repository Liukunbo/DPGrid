#include "harris_matching.h"

#include <iostream>
#include "util/timer.h"
#include "feature/nearest_neighbor.h"
#include "flann/flann.hpp"

#define ROBUST_SIFT
#ifdef ROBUST_SIFT
#include "feature/ransac_fundamental.h"
#endif

FEA_NAMESPACE_BEGIN
FEA_MATCH_NAMESPACE_BEGIN

namespace
{
	void
		convert_descriptor(Harris::Descriptor const& descr, unsigned short* data)
	{
		for (int i = 0; i < 128; ++i)
		{
			uint8_t value = descr.data[i];
			data[i] = static_cast<unsigned char>(value);
		}
	}

#if DISCRETIZE_DESCRIPTORS
	void
		convert_descriptor(Sift::Descriptor const& descr, unsigned short* data)
	{
		for (int i = 0; i < 128; ++i)
		{
			float value = descr.data[i];
			value = math::clamp(value, 0.0f, 1.0f);
			value = math::round(value * 255.0f);
			data[i] = static_cast<unsigned char>(value);
		}
	}

#else // DISCRETIZE_DESCRIPTORS
	void
		convert_descriptor(Sift::Descriptor const& descr, float* data)
	{
		std::copy(descr.data.begin(), descr.data.end(), data);
	}

#endif // DISCRETIZE_DESCRIPTORS
}


void HarrisMatching::set_option(Options option)
{
	this->opt_ = option;
}

void HarrisMatching::init(fea::FeatureSets * featuresets)
{
	this->processed_feature_sets.clear();
	this->processed_feature_sets.resize(featuresets->size());

#pragma omp parallel for schedule(dynamic)
#ifdef _MSC_VER
	for (int64_t i = 0; i < featuresets->size(); i++)
#else
	for (size_t i = 0; i < featuresets->size(); i++)
#endif
	{
		FeatureSet const& fs = (*featuresets)[i];
		ProcessedFeatureSet& pfs = this->processed_feature_sets[i];

		this->init_sift(&pfs.sift_descr, &pfs.sift_position, fs.sift_descriptors, &fs.positions[0]);
		this->init_harris(&pfs.harris_descr, &pfs.harris_position, fs.harris_descriptors, &fs.positions[fs.sift_descriptors.size()]);
	}

	if (this->opt_.sift_matcher_type == SIFT_MATCH_CASCADE_HASHING)
	{
		util::WallTimer timer;
		this->local_data_sift.clear();
		this->local_data_sift.resize(featuresets->size());

		this->global_data.generate_proj_matrices(this->opt_);

		/* Convert feature descriptors to zero mean. */
		math::Vec128f sift_avg;
		compute_avg_descriptors(this->processed_feature_sets, &sift_avg);

#pragma omp parallel for schedule(dynamic)
#ifdef _MSC_VER
		for (int64_t i = 0; i < featuresets->size(); i++)
#else
		for (std::size_t i = 0; i < featuresets->size(); i++)
#endif
		{
			LocalData* ld_sift = &this->local_data_sift[i];

			ProcessedFeatureSet const& pfs = this->processed_feature_sets[i];
			std::vector<math::Vec128f> sift_zero_mean_descs;
			this->compute_zero_mean_descs(&sift_zero_mean_descs,
				pfs.sift_descr, sift_avg);

			this->compute(ld_sift,sift_zero_mean_descs,
				this->global_data, this->opt_);
		}

		if (this->opt_.debug_output)
		{
			std::cout << "Computing cascade hashes took " << timer.get_elapsed()
				<< " ms" << std::endl;
		}	
	}
}

int HarrisMatching::sift_matching(int view_1_id, int view_2_id, Matching::Result * result) const
{
	ProcessedFeatureSet const& pfs_1 = this->processed_feature_sets[view_1_id];
	ProcessedFeatureSet const& pfs_2 = this->processed_feature_sets[view_2_id];
	Matching::Result sift_result;

	/* SIFT matching. */
	if (this->opt_.sift_matcher_type == SIFT_MATCH_EXHAUSTIVE)
	{
		if (pfs_1.sift_descr.size() > 0)
		{
			Matching::twoway_match(this->opts.sift_matching_opts,
				pfs_1.sift_descr.data()->begin(), pfs_1.sift_descr.size(),
				pfs_2.sift_descr.data()->begin(), pfs_2.sift_descr.size(),
				&sift_result);
			Matching::remove_inconsistent_matches(&sift_result);
		}
	}
	else
	{
		LocalData const& ld_sift_1 = this->local_data_sift[view_1_id];
		LocalData const& ld_sift_2 = this->local_data_sift[view_2_id];
		if (pfs_1.sift_descr.size() > 0)
		{
			this->twoway_match(this->opts.sift_matching_opts,
				ld_sift_1, ld_sift_2,
				pfs_1.sift_descr, pfs_2.sift_descr,
				&sift_result, this->opt_);
			Matching::remove_inconsistent_matches(&sift_result);
		}
	}
	int num_matches_old = fea::fea_mch::Matching::count_consistent_matches(sift_result);
	if (num_matches_old < 8)
	{
		if (this->opt_.debug_output)
		{
			std::cout << "  Inital sift matching inliers below threshold of 8 .";
		}

		return EXIT_FAILURE;
	}

#ifdef ROBUST_SIFT
	/* */
	/* Build correspondences from feature matching result. */
	fea::Correspondences2D2D unfiltered_matches;
	fea::CorrespondenceIndices unfiltered_indices;
	{
		std::vector<int> const& m12 = sift_result.matches_1_2;
		for (std::size_t i = 0; i < m12.size(); ++i)
		{
			if (m12[i] < 0)
				continue;

			fea::Correspondence2D2D match;
			match.p1[0] = pfs_1.sift_position[i][0];
			match.p1[1] = pfs_1.sift_position[i][1];
			match.p2[0] = pfs_2.sift_position[m12[i]][0];
			match.p2[1] = pfs_2.sift_position[m12[i]][1];
			unfiltered_matches.push_back(match);
			unfiltered_indices.push_back(std::make_pair(i, m12[i]));
		}
	}

	/* Compute fundamental matrix using RANSAC. */
	fea::RansacFundamental::Result ransac_result;
	int num_inliers = 0;
	{
		fea::RansacFundamental::Options option;
		fea::RansacFundamental ransac(option);
		ransac.estimate(unfiltered_matches, &ransac_result);
		num_inliers = ransac_result.inliers.size();
	}

	/* Require at least 8 inlier matches. */
	int const min_inlier_thres = 4;
	if (num_inliers < min_inlier_thres)
	{
		if (this->opt_.debug_output)
		{
			std::cout << " Roubst sift matching inliers below threshold of 4 .";
		}
		
		return EXIT_FAILURE;
	}

	/* Create Two-View matching result. */
	int sift_result_num = sift_result.matches_1_2.size();
	sift_result.matches_1_2.resize(0);
	sift_result.matches_1_2.resize(sift_result_num, -1);
	sift_result_num = sift_result.matches_2_1.size();
	sift_result.matches_2_1.resize(0);
	sift_result.matches_2_1.resize(sift_result_num, -1);
	for (int i = 0; i < num_inliers; ++i)
	{
		int const inlier_id = ransac_result.inliers[i];
		fea::CorrespondenceIndex in = unfiltered_indices[inlier_id];
		sift_result.matches_1_2[in.first] = in.second;
		sift_result.matches_2_1[in.second] = in.first;
	}

#endif

	result->matches_1_2.swap(sift_result.matches_1_2);
	result->matches_2_1.swap(sift_result.matches_2_1);

	int num_matches = fea::fea_mch::Matching::count_consistent_matches(*result);
	if (this->opt_.debug_output)
	{
		std::cout << "\n  Get sift seed matches of " << num_matches 
			<< " from ( " << num_matches_old  << " ). ";
	}
	return EXIT_SUCCESS;
}

void combine_results(Matching::Result const& sift_result,
	Matching::Result const& harris_result, Matching::Result* result)
{
	/* Determine size of combined matching result. */
	std::size_t num_matches_1 = sift_result.matches_1_2.size()
		+ harris_result.matches_1_2.size();
	std::size_t num_matches_2 = sift_result.matches_2_1.size()
		+ harris_result.matches_2_1.size();

	/* Combine results. */
	result->matches_1_2.clear();
	result->matches_1_2.reserve(num_matches_1);
	result->matches_1_2.insert(result->matches_1_2.end(),
		sift_result.matches_1_2.begin(), sift_result.matches_1_2.end());
	result->matches_1_2.insert(result->matches_1_2.end(),
		harris_result.matches_1_2.begin(), harris_result.matches_1_2.end());

	result->matches_2_1.clear();
	result->matches_2_1.reserve(num_matches_2);
	result->matches_2_1.insert(result->matches_2_1.end(),
		sift_result.matches_2_1.begin(), sift_result.matches_2_1.end());
	result->matches_2_1.insert(result->matches_2_1.end(),
		harris_result.matches_2_1.begin(), harris_result.matches_2_1.end());

	/* Fix offsets. */
	std::size_t harris_offset_1 = sift_result.matches_1_2.size();
	std::size_t harris_offset_2 = sift_result.matches_2_1.size();

	if (harris_offset_2 > 0)
		for (std::size_t i = harris_offset_1; i < result->matches_1_2.size(); ++i)
			if (result->matches_1_2[i] >= 0)
				result->matches_1_2[i] += harris_offset_2;

	if (harris_offset_1 > 0)
		for (std::size_t i = harris_offset_2; i < result->matches_2_1.size(); ++i)
			if (result->matches_2_1[i] >= 0)
				result->matches_2_1[i] += harris_offset_1;

}

void HarrisMatching::pairwise_match(int view_1_id, int view_2_id, Matching::Result * result) const
{
	ProcessedFeatureSet const& pfs_1 = this->processed_feature_sets[view_1_id];
	ProcessedFeatureSet const& pfs_2 = this->processed_feature_sets[view_2_id];
	Matching::Result sift_result;

	if (EXIT_SUCCESS != this->sift_matching(view_1_id, view_2_id, &sift_result))
	{
		return;
	}

	/*Matching::Result harris_result;
	int harris_result_num = harris_result.matches_1_2.size();
	harris_result.matches_1_2.resize(0);
	harris_result.matches_1_2.resize(harris_result_num, -1);
	harris_result_num = harris_result.matches_2_1.size();
	harris_result.matches_2_1.resize(0);
	harris_result.matches_2_1.resize(harris_result_num, -1);*/

	/* Do harris mathching. */
	Matching::Result harris_result;
	this->harris_match_two_way(pfs_1, pfs_2, sift_result, &harris_result);
	//Matching::remove_inconsistent_matches(&harris_result);
	int num_matches = fea::fea_mch::Matching::count_consistent_matches(harris_result);
	if (this->opt_.debug_output)
	{
		std::cout << "\n  Get inital harris matches paris of " << num_matches
			<< " from ( " << harris_result.matches_1_2.size() << " ). ";
	}

	/* combine matching result. (sift matching are set to -1 )*/
	int sift_result_num = sift_result.matches_1_2.size();
	sift_result.matches_1_2.resize(0);
	sift_result.matches_1_2.resize(sift_result_num, -1);
	sift_result_num = sift_result.matches_2_1.size();
	sift_result.matches_2_1.resize(0);
	sift_result.matches_2_1.resize(sift_result_num, -1);

	combine_results(sift_result, harris_result, result);
}

int HarrisMatching::pairwise_match_lowres(int view_1_id, int view_2_id, std::size_t num_features) const
{
	/* If loweres matching with sift cannot be done then, we cannot do harris matching */
	ProcessedFeatureSet const& pfs_1 = this->processed_feature_sets[view_1_id];
	ProcessedFeatureSet const& pfs_2 = this->processed_feature_sets[view_2_id];

	/* SIFT lowres matching. */
	if (pfs_1.sift_descr.size() > 0)
	{
		Matching::Result sift_result;
		Matching::twoway_match(this->opts.sift_matching_opts,
			pfs_1.sift_descr.data()->begin(),
			std::min(num_features, pfs_1.sift_descr.size()),
			pfs_2.sift_descr.data()->begin(),
			std::min(num_features, pfs_2.sift_descr.size()),
			&sift_result);
		return Matching::count_consistent_matches(sift_result);
	}

	return 0;
}

/* ---------------------------------------------------------------- */

void HarrisMatching::init_sift(SiftDescriptors * dst, SiftFeaturePositions* posi, Sift::Descriptors const & src, math::Vec2f const * pos)
{
	/* Prepare and copy to data structures. */
	dst->resize(src.size());
	posi->resize(src.size());

#if DISCRETIZE_DESCRIPTORS
	uint16_t* ptr = dst->data()->begin();
#else
	float* ptr = dst->data()->begin();
#endif

	float* ppos = posi->data()->begin();

	for (std::size_t i = 0; i < src.size(); ++i, ptr += 128, ppos += 2)
	{
		Sift::Descriptor const& d = src[i];
		convert_descriptor(d, ptr);
		ppos[0] = pos[i][0]; ppos[1] = pos[i][1];
	}
}

void HarrisMatching::init_harris(HarrisDescriptors * dst, HarrisFeaturePositions* posi, Harris::Descriptors const & src, math::Vec2f const * pos)
{
	/* Prepare and copy to data structures. */
	dst->resize(src.size());
	posi->resize(src.size());

	uint16_t* ptr = dst->data()->begin();
	float* ppos = posi->data()->begin();

	for (std::size_t i = 0; i < src.size(); ++i, ptr += 128, ppos += 2)
	{
		Harris::Descriptor const& d = src[i];
		convert_descriptor(d, ptr);
		//ppos[0] = d.x; ppos[1] = d.y;
		ppos[0] = pos[i][0]; ppos[1] = pos[i][1];
	}
}
/* ---------------------------------------------------------------- */

namespace
{
	///////////////   法化   法化   法化   法化   法化
	static 
	void nrml(double *aa, int n, double bb, double *a, double *b)
	{
		int  i, j;

		for (i = 0; i < n; i++)
		{
			for (j = 0; j < n - i; j++)
			{
				*a += *aa * *(aa + j);
				a++;
			}
			*b += *aa * bb;
			b++;
			aa++;
		}
	}

	///////////////   solve   solve   solve   solve   solve
	static
	void zero(double *x, int n)
	{
		memset(x, 0, sizeof(double)*n);
	}

	static
	void ldltban1(double *a, double *d, double *l, int n, int wide)
	{
		int i, j, k, kk, km, m;
		double *ao = nullptr, *aa = nullptr, *co = nullptr, *c = nullptr;

		m = wide*(2 * n + 1 - wide) / 2;
		c = new double[m - wide]; // (double *)calloc((m - wide), sizeof(double));

		ao = a; co = c; a += wide;
		for (i = 0; i < m - wide; i++) *c++ = *a++;
		c = co; a = ao;

		for (k = 1; k < n; k++) {
			if (k < n - wide + 2) kk = wide - 1;
			else kk--;

			*d = *a++; aa = a;  a += kk;

			if (k < n - wide + 1) km = wide;
			else km = n - k + 1;

			for (i = 1; i < kk + 1; i++) {
				*l = *aa++ / *d;
				for (j = 0; j<kk - i + 1; j++) *(a + j) -= *l * *(aa + j - 1);
				l++;

				if (k + i>n - wide + 1) km--;
				a += km;
			}

			a = aa; d++;
			if (k == n - 1)  *d = *a;
		}

		a = ao;  a += wide;
		for (i = 0; i < m - wide; i++) *a++ = *c++;
		c = co;
		delete[]c; c = nullptr;
	}

	static
	void ldltban2(double*l, double*d, double*b, double*x, int n, int wide)
	{
		int i, j, kk, m;
		double *bo = nullptr, *lo = nullptr, *xx = nullptr;
		double *bb = nullptr, *bbo = nullptr;

		bb = new double[n]; // (double*)calloc(n, sizeof(double));
		bbo = bb;

		bo = b; lo = l;

		for (i = 0; i < n; i++)*bb++ = *b++;
		b = bo;  bb = bbo;
		m = wide*(2 * n + 1 - wide) / 2;

		for (i = 1; i < n; i++) {
			if (i < n - wide + 2) kk = wide;
			else kk--;

			b = bo + i;
			for (j = 1; j < kk; j++) *b++ -= *(b - j) * *l++;
		}

		kk = 0;
		b = bo + n - 1;  l = lo + m - n - 1;
		x += n - 1;  xx = x;  d += n - 1;
		*x-- = *b-- / *d--;

		for (i = 1; i < n; i++) {
			if (i < wide) kk++;
			else { kk = wide - 1;  xx--; }

			*x = *b-- / *d--;
			for (j = 1; j < kk + 1; j++) *x -= *l-- * *(xx - j + 1);
			x--;
		}

		b = bo;
		for (i = 0; i < n; i++) *b++ = *bb++;
		bb = bbo;
		delete[]bb; bb = nullptr;
	}

	static
	void solve(double *a, double *b, double *x, int n, int wide)
	{
		int      m;
		double   *d = nullptr, *l = nullptr;

		m = n*(n + 1) / 2;
		d = new double[n]; // (double *)malloc(n * sizeof(double));
		l = new double[m - n]; // (double *)malloc((m - n) * sizeof(double));
		memset(d, 0, sizeof(double)*n);
		memset(l, 0, sizeof(double)*(m - n));

		ldltban1(a, d, l, n, wide);
		ldltban2(l, d, b, x, n, wide);

		delete[]d; d = nullptr;
		delete[]l; l = nullptr;
	}

}

int HarrisMatching::get_affine_trans(SiftFeaturePositions const & sift_fea_1, SiftFeaturePositions const & sift_fea_2, std::vector<int> const & match_result, math::Vec6d * affine_trans_param) const
{
	affine_trans_param->fill(0.0);
	(*affine_trans_param)[0] = 1.0f; (*affine_trans_param)[3] = 1.0f;
	std::vector<int> const& m12 = match_result;
	if (m12.size() < 3)
	{
		if (this->opt_.debug_output)
		{
			std::cout << "At least 3 matches required" << std::endl;
		}

		return -1;
	}

	int size = m12.size();

	/* Create N * 3 matrix A and N * 2 Matrix B. Each match result create one row in A and one row in B*/
	double a[6] = { 0.0 }; // 2 * 6 
	double aa[36] = { 0.0 }; // 6 * 6
	double l[2] = { 0.0 }; // 2 * 1
	double al[6] = { 0.0 }; // 6 * 1
	float srcx, srcy, dstx, dsty;
	for (std::size_t i = 0; i < size; i++)
	{
		if (m12[i] < 0)
			continue;
		
		srcx = sift_fea_1[i][0];
		srcy = sift_fea_1[i][1];

		dstx = sift_fea_2[m12[i]][0];
		dsty = sift_fea_2[m12[i]][1];

		memset(a, 0, sizeof(double) * 6);
		a[0] = 1; a[1] = srcx; a[2] = srcy; a[3] = 0; a[4] = 0; a[5] = 0;
		l[0] = dstx;
		nrml(a, 6, l[0], aa, al);

		memset(a, 0, sizeof(double) * 6);
		a[3] = 1; a[4] = srcx; a[5] = srcy; a[0] = 0; a[1] = 0; a[2] = 0;
		l[1] = dsty;
		nrml(a, 6, l[1], aa, al);
	}

	memset(affine_trans_param->begin(), 0, sizeof(double) * 6);
	solve(aa, al, affine_trans_param->begin(), 6, 6);

	int cnt = 0;
	for (int i = 0; i < 6; i++)
	{
		if (_finite(*(affine_trans_param->begin() + i)) == 0)
		{
			cnt++;
		}
	}
	if (cnt == 6)
	{
		if (this->opt_.debug_output)
		{
			std::cout << "Wrong Calculated Affine Parameters! " << std::endl;
		}

		memset(affine_trans_param->begin(), 0, sizeof(double) * 6);
		(*affine_trans_param)[0] = 1.0f; (*affine_trans_param)[3] = 1.0f;

		return -1;
	}

	return 0;
}

int HarrisMatching::affine_transition(math::Vec6d const & affine_trans_param, math::Vec2f const * src, math::Vec2f * dst, int size) const
{
//#pragma omp parallel for schedule(dynamic)
#ifdef _MSC_VER
	for (int64_t i = 0; i < size; i++)
#else
	for (std::size_t i = 0; i < size; i++)
#endif
	{
		math::Vec2f const & p1 = src[i];
		math::Vec2f & p2 = dst[i];

		p2[0] = (float)(affine_trans_param[0] + affine_trans_param[1] * p1[0] + affine_trans_param[2] * p1[1]);
		p2[1] = (float)(affine_trans_param[3] + affine_trans_param[4] * p1[0] + affine_trans_param[5] * p1[1]);
	}

	return 0;
}

int descriptor_affine_trans(math::Vec6d const & affine_trans_param, math::Vec128us const & src, math::Vec128us * dst, int window_len = 7)
{
	int des_window_len = 11;
	int sub_des_window = des_window_len / 2;
	int sub_window = window_len / 2;
	if (sub_des_window < sub_window) return -1;

	dst->fill(0);
	//memset(dst->begin(), 0, sizeof(unsigned short) * 128);
	double x, y;
	int ix, iy; double dx, dy;

	//int dstx, dsty;
	int srcx, srcy;
	const unsigned short* src_data = nullptr;
	unsigned short* dst_data = nullptr;
	for (int j = -sub_window; j <= sub_window; j++)
	{	
		dst_data = dst->begin() + (sub_des_window + j) * des_window_len + (sub_des_window - sub_window);
		for (int i = -sub_window; i <= sub_window; i++)
		{
			/* Just do rotate */
			x = i * affine_trans_param[1] + j * affine_trans_param[2];
			y = i * affine_trans_param[4] + j * affine_trans_param[5];

			ix = (int)std::floor(x); iy = (int)std::floor(y);
			dx = x - ix; dy = y - iy;

			/* Do resampleing */
			//dstx = i + sub_des_window; dsty = j + sub_des_window;
			srcx = ix + sub_des_window; srcy = iy + sub_des_window;
			if (srcx < 0 || srcx + 1 >= des_window_len || srcy < 0 || srcy + 1 >= des_window_len)
			{
				return -1;
			}

			unsigned short resample_v = 0;
			src_data = src.begin() + srcy * des_window_len + srcx;
			resample_v += unsigned short(*src_data * (1.0 - dx) * (1.0 - dy));
			resample_v += unsigned short(*(src_data + 1) * dx * (1.0 - dy));
			resample_v += unsigned short(*(src_data + des_window_len) * (1.0 - dx) * dy);
			resample_v += unsigned short(*(src_data + des_window_len + 1) * dx * dy);
			*dst_data = resample_v;
			dst_data++;
		}
	}
	src_data = nullptr;
	dst_data = nullptr;

	return 0;
}

float HarrisMatching::ncc_harris(math::Vec128us const & des1, math::Vec128us const & des2, int window_len) const
{
	int des_window_len = 11;
	int sub_des_window = des_window_len / 2;
	int sub_window = window_len / 2;
	if (sub_des_window < sub_window) return 0.0f;

	double w1 = 0.0f, w2 = 0.0f, w3 = 0.0f;
	const unsigned short* ds1 = nullptr;
	const unsigned short* ds2 = nullptr;
	for (int j = -sub_window; j <= sub_window; j++)
	{
		ds1 = des1.begin() + (sub_des_window + j) * des_window_len + (sub_des_window - sub_window);
		ds2 = des2.begin() + (sub_des_window + j) * des_window_len + (sub_des_window - sub_window);
		for (int i = -sub_window; i <= sub_window; i++)
		{
			w1 += (*ds1 * *ds1);
			w2 += (*ds2 * *ds2);
			w3 += (*ds1 * *ds2);
			ds1++;
			ds2++;
		}
	}

	double ncc_val = w3 / (std::sqrt(w1) * std::sqrt(w2));
	return float(ncc_val);
}

void HarrisMatching::harris_match_one_way(ProcessedFeatureSet const & pfs1, ProcessedFeatureSet const & pfs2, std::vector<int> const & sift_result, std::vector<int> * harris_matches) const
{
	harris_matches->clear();
	harris_matches->resize(pfs1.harris_descr.size(), -1);

	/* Compute the Affine trans param using sift matching results. */
	math::Vec6d affine_param(0.0);
	if (this->get_affine_trans(pfs1.sift_position, pfs2.sift_position, sift_result, &affine_param) != 0)
	{		
		return;
	}

	/* Do Harris match */
	int harris_fea_num = pfs1.harris_descr.size();

	/*NearestNeighbor<float> nn;
	nn.set_elements(pfs2.harris_position.data()->begin());
	nn.set_num_elements(pfs2.harris_position.size());
	nn.set_element_dimensions(2);*/
	flann::Matrix<float> data((float*)(pfs2.harris_position.data()->begin()), pfs2.harris_position.size(), 2);
	flann::Index<flann::L2<float>> index(data, flann::KDTreeIndexParams(1));
	index.buildIndex();

#pragma omp parallel for schedule(dynamic)
#ifdef _MSC_VER
	for (int64_t i = 0; i < harris_fea_num; i++)
#else
	for (std::size_t i = 0; i < harris_fea_num; i++)
#endif
	{
		/* estimate the init location */
		math::Vec2f pest;
		this->affine_transition(affine_param, &pfs1.harris_position[i], &pest);
	
		/* find neighbor harris feature */
		int neighbor_num = this->opt_.near_feat_num;
		std::vector<int> neighbor_index;
		std::vector<float> neighbor_dis;
		neighbor_index.resize(neighbor_num, -1);
		neighbor_dis.resize(neighbor_num, -1);
		{
			///* using this we will find 2 neighbor */
			//typename NearestNeighbor<float>::Result nn_result;
			//float const * query_pointer = pest.begin();
			//nn.find(query_pointer, &nn_result);
			//if (nn_result.dist_1st_best > this->opt_.max_near_distance)
			//{
			//	continue;
			//}	
			//neighbor_index.emplace_back(nn_result.index_1st_best);
			//neighbor_index.emplace_back(nn_result.index_2nd_best);
		}
		{
			// using flann
			flann::Matrix<int> NeiIndex(&neighbor_index[0], 1, neighbor_num);
			flann::Matrix<float> NeiDis(&neighbor_dis[0], 1, neighbor_num);
			flann::Matrix<float> query(pest.begin(), 1, 2);
			index.knnSearch(query, NeiIndex, NeiDis, neighbor_num, flann::SearchParams(32));
		}
		
		/* calculate ncc and find the max location. */
		int maxloc = -1;
		float maxNcc = -1;	
		for (std::size_t j = 0; j < neighbor_index.size(); j++)
		{		
			/* Get affine transited descriptor */
			math::Vec128us des_affine; des_affine.fill(0);
			if (0 != descriptor_affine_trans(affine_param, pfs2.harris_descr[neighbor_index[j]], &des_affine, 7))
			{
				continue;
			}

			/* Calculate NCC */
			float ncc_val = ncc_harris(pfs1.harris_descr[i], des_affine, 7);
			if (ncc_val > maxNcc)
			{
				maxNcc = ncc_val;
				maxloc = neighbor_index[j];
			}
		}

		/* add */
		if (maxNcc > 0.75)
		{
			harris_matches->at(i) = maxloc;
		}
	}
}

void HarrisMatching::harris_match_two_way(ProcessedFeatureSet const & pfs1, ProcessedFeatureSet const & pfs2, Matching::Result const & sift_result, Matching::Result * matches) const
{
	this->harris_match_one_way(pfs1, pfs2, sift_result.matches_1_2, &matches->matches_1_2);
	this->harris_match_one_way(pfs2, pfs1, sift_result.matches_2_1, &matches->matches_2_1);
}

/* ---------------------------------------------------------------- */

void
HarrisMatching::GlobalData::generate_proj_matrices(Options const& opts)
{
	generate_proj_matrices(
		&this->sift.prim_proj_mat,
		&this->sift.sec_proj_mats,
		opts);
}

void 
HarrisMatching::compute(LocalData * ld_sift, std::vector<math::Vec128f> const & sift_zero_mean_descs, GlobalData const & cashash_global_data, Options const & cashash_opts)
{
	/* Compute cascade hashes of zero mean SIFT/SURF descriptors. */
	compute_cascade_hashes(
		sift_zero_mean_descs,
		&ld_sift->comp_hash_data,
		&ld_sift->bucket_grps_bucket_ids,
		cashash_global_data.sift.prim_proj_mat,
		cashash_global_data.sift.sec_proj_mats,
		cashash_opts);


	/* Build buckets. */
	build_buckets(
		&ld_sift->bucket_grps_feature_ids,
		ld_sift->bucket_grps_bucket_ids,
		sift_zero_mean_descs.size(),
		cashash_opts);
}

void 
HarrisMatching::compute_avg_descriptors(ProcessedFeatureSets const & pfs, math::Vec128f * sift_avg)
{
	math::Vec128f sift_vec_sum(0.0f);
	std::size_t num_sift_descs_total = 0;

	/* Compute sum of all SIFT descriptors. */
	for (std::size_t i = 0; i < pfs.size(); i++)
	{
		SiftDescriptors const& sift_descr = pfs[i].sift_descr;

		std::size_t num_sift_descriptors = sift_descr.size();
		num_sift_descs_total += num_sift_descriptors;

		for (std::size_t j = 0; j < num_sift_descriptors; j++)
			for (int k = 0; k < 128; k++)
				sift_vec_sum[k] += sift_descr[j][k] / 255.0f;
	}

	/* Compute average vectors for SIFT. */
	*sift_avg = sift_vec_sum / num_sift_descs_total;
}

void 
HarrisMatching::compute_zero_mean_descs(std::vector<math::Vec128f>* sift_zero_mean_descs, SiftDescriptors const & sift_descs, math::Vec128f const & sift_avg)
{
	/* Compute zero mean descriptors. */
	sift_zero_mean_descs->resize(sift_descs.size());
	for (std::size_t i = 0; i < sift_descs.size(); i++)
		for (int j = 0; j < 128; j++)
			(*sift_zero_mean_descs)[i][j] = sift_descs[i][j] / 255.0f - sift_avg[j];
}

void
HarrisMatching::build_buckets(BucketGroupsFeatures * bucket_grps_feature_ids, BucketGroupsBuckets const & bucket_grps_bucket_ids, size_t num_descs, Options const & cashash_opts)
{
	uint8_t const num_bucket_grps = this->opt_.num_bucket_groups;
	uint8_t const num_bucket_bits = this->opt_.num_bucket_bits;
	uint32_t const num_buckets_per_group = 1 << num_bucket_bits;

	bucket_grps_feature_ids->resize(num_bucket_grps);

	for (uint8_t grp_idx = 0; grp_idx < num_bucket_grps; grp_idx++)
	{
		BucketGroupFeatures &bucket_grp_features = (*bucket_grps_feature_ids)[grp_idx];
		bucket_grp_features.resize(num_buckets_per_group);
		for (size_t i = 0; i < num_descs; i++)
		{
			uint16_t bucket_id = bucket_grps_bucket_ids[grp_idx][i];
			bucket_grp_features[bucket_id].emplace_back(i);
		}
	}
}

FEA_MATCH_NAMESPACE_END
FEA_NAMESPACE_END
