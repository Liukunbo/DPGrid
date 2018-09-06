#pragma once

#include <string>
#include <vector>

#if _WIN32
# ifndef MAKE_DLL
#   define FeatureMatcherAPI __declspec(dllimport)
# else
#	define FeatureMatcherAPI __declspec(dllexport)
# endif
#else
#	define FeatureMatcherAPI
#endif

class FeatureMatcherAPI FeatureMatcher
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
		/** Minimum number of matching features */
		int min_feature_matches = 24;
		/** Perform low-resolution matching to reject unlikely pairs. */
		bool use_lowres_matching = false;
		/** Matcher type. Exhaustive by default. */
		MatcherType matcher_type = MATCHER_EXHAUSTIVE;
		/** Save debug information. Default is false */
		bool debug_output = false;
		/** 
		   Do we need to reverse the y coordinates. 
		   Default is false. Any change to this will not work!
		*/
		bool reverse_y_coord = false;
		/** 
		   Do we need to normalize the coordinates. 
		   Default is false. Any change to this will not work!
		*/
		bool anti_normalize_coord = false;
	};

public:
	explicit FeatureMatcher(FeatureMatcher::Options const & options);

	~FeatureMatcher();

	void SetOptions(FeatureMatcher::Options const & options);

	int Match(std::string szFea1, std::string szFea2, std::string szMchFile);

	int Match(std::vector<std::string> szFeaList, std::string szRegionFile, std::string szMchFile, std::size_t expandSize = 0);

private:
	Options opt_;
};

