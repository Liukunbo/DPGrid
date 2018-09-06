/*
 * Copyright (C) 2015, Simon Fuhrmann
 * TU Darmstadt - Graphics, Capture and Massively Parallel Computing
 * All rights reserved.
 *
 * This software may be modified and distributed under the terms
 * of the BSD 3-Clause license. See the LICENSE.txt file for details.
 */

#include <iostream>
#include <algorithm>

#include "feature/feature_set.h"

FEA_NAMESPACE_BEGIN

namespace
{
    template <typename T>
    bool
    compare_scale (T const& descr1, T const& descr2)
    {
        return descr1.scale > descr2.scale;
    }
}  /* namespace */

void
FeatureSet::compute_features (dpgrid::ByteImage::Ptr image)
{
    this->colors.clear();
    this->positions.clear();
    this->width = image->width();
    this->height = image->height();

    /* Make sure these are in the right order. Matching relies on it. */
	switch (this->opts.feature_types)
	{
	case FEATURE_SIFT:
		this->compute_sift(image);
		break;
	case FEATURE_SURF:
		this->compute_surf(image);
		break;
	case FEATURE_HARRIS:
		this->compute_harris(image);
		break;
	case FEATURE_SIFTSURF:
		this->compute_sift(image);
		this->compute_surf(image);
		break;
	default:
		break;
	}
}

void
FeatureSet::compute_sift (dpgrid::ByteImage::ConstPtr image)
{
    /* Compute features. */
    Sift::Descriptors descr;
    {
        Sift sift(this->opts.sift_opts);
        sift.set_image(image);
        sift.process();
        descr = sift.get_descriptors();
    }

    /* Sort features by scale for low-res matching. */
    std::sort(descr.begin(), descr.end(), compare_scale<fea::Sift::Descriptor>);

    /* Prepare and copy to data structures. */
    std::size_t offset = this->positions.size();
    this->positions.resize(offset + descr.size());
    this->colors.resize(offset + descr.size());

    for (std::size_t i = 0; i < descr.size(); ++i)
    {
        Sift::Descriptor const& d = descr[i];
        this->positions[offset + i] = math::Vec2f(d.x, d.y);
        image->linear_at(d.x, d.y, this->colors[offset + i].begin());
    }

    /* Keep SIFT descriptors. */
    std::swap(descr, this->sift_descriptors);
}

void
FeatureSet::compute_surf (dpgrid::ByteImage::ConstPtr image)
{
    /* Compute features. */
    Surf::Descriptors descr;
    {
        Surf surf(this->opts.surf_opts);
        surf.set_image(image);
        surf.process();
        descr = surf.get_descriptors();
    }

    /* Sort features by scale for low-res matching. */
    std::sort(descr.begin(), descr.end(), compare_scale<fea::Surf::Descriptor>);

    /* Prepare and copy to data structures. */
    std::size_t offset = this->positions.size();
    this->positions.resize(offset + descr.size());
    this->colors.resize(offset + descr.size());

    for (std::size_t i = 0; i < descr.size(); ++i)
    {
        Surf::Descriptor const& d = descr[i];
        this->positions[offset + i] = math::Vec2f(d.x, d.y);
        image->linear_at(d.x, d.y, this->colors[offset + i].begin());
    }

    /* Keep SURF descriptors. */
    std::swap(descr, this->surf_descriptors);
}

void FeatureSet::compute_harris(dpgrid::ByteImage::ConstPtr image)
{
	/* Compute features. */
	Harris::Descriptors descr;
	{
		Harris harris(this->opts.harris_opts);
		harris.set_image(image);
		harris.process();
		descr = harris.get_descriptors();
	}

	/* Sort features by scale. */
	std::sort(descr.begin(), descr.end(), compare_scale<fea::Harris::Descriptor>);

	/* Prepare and copy to data structures. */
	std::size_t offset = this->positions.size();
	this->positions.resize(offset + descr.size());
	this->colors.resize(offset + descr.size());

	for (std::size_t i = 0; i < descr.size(); ++i)
	{
		Harris::Descriptor const& d = descr[i];
		this->positions[offset + i] = math::Vec2f(d.x, d.y);
		image->linear_at(d.x, d.y, this->colors[offset + i].begin());
	}

	/* Keep Harris descriptors. */
	std::swap(descr, this->harris_descriptors);
}

void
FeatureSet::clear_descriptors (void)
{
    this->sift_descriptors.clear();
    this->sift_descriptors.shrink_to_fit();
    this->surf_descriptors.clear();
    this->surf_descriptors.shrink_to_fit();
	this->harris_descriptors.clear();
	this->harris_descriptors.shrink_to_fit();
}

std::size_t FeatureSet::get_feature_num(void)
{
	switch (this->opts.feature_types)
	{
	case FEATURE_SIFTSURF:
		return (std::size_t)this->positions.size();
		break;
	case FEATURE_SIFT:
		return (std::size_t)this->sift_descriptors.size();
		break;
	case FEATURE_SURF:
		return (std::size_t)this->surf_descriptors.size();
		break;
	case FEATURE_HARRIS:
		return (std::size_t)this->harris_descriptors.size();
	default:
		break;
	}
	return 0;
}

FEA_NAMESPACE_END
