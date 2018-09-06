#pragma once

#include <string>
#include <vector>

#if _WIN32
# ifndef MAKE_DLL
#   define FeatureConnectorAPI __declspec(dllimport)
# else
#	define FeatureConnectorAPI __declspec(dllexport)
# endif
#else
#	define FeatureConnectorAPI
#endif

class FeatureConnectorAPI FeatureConnector
{
public:
	struct Option
	{
		/* The extension name of the feature file, default is ".fea" */
		std::string szFeatfileExt = ".feat";

		/* save debug information. defalut is false */
		bool debug_output = false;
	};

public:
	explicit FeatureConnector(FeatureConnector::Option const & options);
	~FeatureConnector();

	int AddMatchFile(std::string const szMatchFile);

	void RemoveAll();

	/* this function build tracks which id are mainly calculated from small image id */
	int Compute(std::string szImageListFile, std::string szFeatDir, std::string szPtsFile);

	/* if expandsize == 0, this function is recommended */
	int Compute2(std::string szImageListFile, std::string szFeatDir, std::string szPtsFile);

	/* if expandsize != 0, this function is recommended */
	int Compute3(std::string szImageListFile, std::string szFeatDir, std::string szPtsFile);

private:
	Option opt_;
	std::vector<std::string> match_file_list_;
};

