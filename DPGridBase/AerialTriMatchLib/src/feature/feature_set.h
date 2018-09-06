/*
 * Copyright (C) 2015, Simon Fuhrmann
 * TU Darmstadt - Graphics, Capture and Massively Parallel Computing
 * All rights reserved.
 *
 * This software may be modified and distributed under the terms
 * of the BSD 3-Clause license. See the LICENSE.txt file for details.
 */

#ifndef FEA_FEATURE_SET_HEADER
#define FEA_FEATURE_SET_HEADER

#include <vector>

#include "math/vector.h"
#include "util/aligned_memory.h"
#include "feature/sift.h"
#include "feature/surf.h"
#include "feature/harris.h"
#include "feature/defines.h"

FEA_NAMESPACE_BEGIN

/**
 * The FeatureSet holds per-feature information for a single view, and
 * allows to transparently compute and match multiple feature types.
 */
class FeatureSet
{
public:
    /** Bitmask with feature types. */
    enum FeatureTypes
    {
		FEATURE_SIFT = 0x1,
		FEATURE_SURF = 0x2,
		FEATURE_HARRIS = 0x3,
		FEATURE_ALL = 0xFF
    };

    /** Options for feature detection and matching. */
    struct Options
    {
        Options (void);

        FeatureTypes feature_types;
        Sift::Options sift_opts;
        Surf::Options surf_opts;
		Harris::Options harris_opts;
    };

public:
    FeatureSet (void);
    explicit FeatureSet (Options const& options);
    void set_options (Options const& options);

    /** Computes the features specified in the options. */
    void compute_features (dpgrid::ByteImage::Ptr image);

    /** Clear descriptor data. */
    void clear_descriptors (void);

	/** Get the feature size*/
	std::size_t get_feature_num (void);

public:
    /** Image dimension used for feature computation. */
    int width, height;
	/** Feature scale. (used for detec image on pyrimid image and it must be set mannually )*/
	int scale;
    /** Per-feature image position. */
    std::vector<math::Vec2f> positions;
    /** Per-feature image color. */
    std::vector<math::Vec3uc> colors;
    /** The SIFT descriptors. */
    Sift::Descriptors sift_descriptors;
    /** The SURF descriptors. */
    Surf::Descriptors surf_descriptors;
	/* The Harris descriptors. */
	Harris::Descriptors harris_descriptors;

private:
    void compute_sift (dpgrid::ByteImage::ConstPtr image);
    void compute_surf (dpgrid::ByteImage::ConstPtr image);
	void compute_harris(dpgrid::ByteImage::ConstPtr image);

private:
    Options opts;
};

/** Only works for add harris feature sets to other set. */
inline void CombineFeature(FeatureSet & set1, FeatureSet const & harris_set);

typedef std::vector<FeatureSet> FeatureSets;

/* ------------------------ Implementation ------------------------ */

inline
FeatureSet::Options::Options (void)
    : feature_types(FEATURE_SIFT)
{
}

inline
FeatureSet::FeatureSet (void)
{
}

inline
FeatureSet::FeatureSet (Options const& options)
    : opts(options)
{
}

inline void
FeatureSet::set_options (Options const& options)
{
    this->opts = options;
}

inline void CombineFeature(FeatureSet & set1, FeatureSet const & harris_set)
{
	// copy position
	int oldsize = set1.positions.size();
	int harrisSize = harris_set.harris_descriptors.size();
	set1.positions.resize(oldsize + harrisSize);
	int harrisposLoc = 0 + harris_set.sift_descriptors.size() + harris_set.surf_descriptors.size();
	memcpy(&(set1.positions[oldsize]), &(harris_set.positions[harrisposLoc]), sizeof(math::Vec2f) * harrisSize);

	// copy descriptor
	set1.harris_descriptors.resize(harrisSize);
	memcpy(&set1.harris_descriptors[0], &harris_set.harris_descriptors[0], sizeof(Harris::Descriptor) * harrisSize);
}

FEA_NAMESPACE_END

#endif /* FEA_FEATURE_SET_HEADER */
