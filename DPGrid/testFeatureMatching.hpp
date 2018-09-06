#pragma once

#include <iostream>
#include <io.h>
#include <direct.h>

#include "FeatureExtractor\FeatureExtractor.h"
#include "FeatureMatcher\FeatureMatcher.h"
#include "FeatureConnector\FeatureConnector.h"

#if _DEBUG
#pragma comment(lib, "FeatureExtractorD.lib")
#pragma comment(lib, "FeatureMatcherD.lib")
#pragma comment(lib, "FeatureConnectorD.lib")
#else
#pragma comment(lib, "FeatureExtractor.lib")
#pragma comment(lib, "FeatureMatcher.lib")
#pragma comment(lib, "FeatureConnector.lib")
#endif

#include "FeatureMatcher\io_file.hpp"

int extrac_fea(std::vector<std::string> const szImageFiles, std::string szFeatDir, std::string szFeatExt = ".feat", int featureType = 1)
{
	FeatureExtractor::Options options;
	options.debug_output = false;
	options.feature_type = FeatureExtractor::FeatureTypes(featureType);
	options.feature_lvl = FeatureExtractor::FeatureLevel::LEVEL_AUTO;

	FeatureExtractor extractor(options);

	int ret = EXIT_SUCCESS;
#pragma omp parallel for schedule(dynamic)
#ifdef _MSC_VER
	for (int64_t i = 0; i < szImageFiles.size(); i++)
#else
	for (std::size_t i = 0; i < szImageFiles.size(); i++)
#endif
	{
		int loc = szImageFiles[i].rfind('\\');
		int loc2 = szImageFiles[i].rfind('.');
		std::string szFeatFile = szFeatDir + szImageFiles[i].substr(loc, loc2 - loc) + szFeatExt;

		if (_access(szFeatFile.c_str(), 0) == 0) continue;
		extractor.FeatureDetec(szImageFiles[i], szFeatFile);

#pragma omp critical
		{
			std::cout << "Detect features on image ( " << szImageFiles[i] << " ) done! " << std::endl;
		}
	}

	return ret;
}

int match(std::vector<std::string> szFeaList, std::string szRegionFile, std::string szMchFile, std::size_t expandSize = 0, int matchType = 1)
{
	if (_access(szMchFile.c_str(), 0) == 0) return EXIT_SUCCESS;

	FeatureMatcher::Options options;
	options.debug_output = false;
	options.min_feature_matches = 8;
	options.matcher_type = FeatureMatcher::MatcherType(matchType);
	FeatureMatcher matcher(options);

	return matcher.Match(szFeaList, szRegionFile, szMchFile, expandSize);
}

int match(std::vector<std::string> const & szTriidxFileList, std::string szFeaDir, std::string szMchDir, int expandSize,
	std::vector<std::string> & szMchFileList, std::string szFeatExt = ".feat", int matchType = 1)
{
#pragma omp parallel for schedule(dynamic)
#ifdef _MSC_VER
	for (int64_t i = 0; i < szTriidxFileList.size(); i++)
#else
	for (std::size_t i = 0; i < szTriidxFileList.size(); i++)
#endif
	{
		tri::triinfo tri;
		tri::read_tri_idx(szTriidxFileList[i].c_str(), tri);

		std::vector<std::string> imageFiles;
		std::vector<std::string> featFiles;
		imageFiles.resize(tri.size());
		featFiles.resize(tri.size());

		/* features */
		for (std::size_t i = 0; i < tri.size(); i++)
		{
			std::string name; name.assign(tri[i].name);
			name = name.substr(0, name.rfind('.')) + szFeatExt;

			std::string path = szFeaDir + "\\" + name;
			featFiles[i] = path;
		}

		/* do maches */
		std::string name;
		std::size_t base = szTriidxFileList[i].rfind('\\');
		name = szTriidxFileList[i].substr(base + 1, szTriidxFileList[i].length() - base);
		name = name.substr(0, name.rfind('.')) + ".mch";
		std::string mchfile = szMchDir + "\\" + name;
		if (EXIT_SUCCESS == match(featFiles, szTriidxFileList[i], mchfile, expandSize, matchType))
		{
			szMchFileList.emplace_back(mchfile);
		}
	}

	return EXIT_SUCCESS;
}

int sort_matches(std::string szImageList, std::string szFeatDir, std::vector<std::string> const & szMatchFileList, std::string szPtsFile, std::string szFeatExt = ".feat")
{
	FeatureConnector::Option options;
	options.szFeatfileExt = szFeatExt;
	options.debug_output = true;
	FeatureConnector connector(options);

	for (std::size_t i = 0; i < szMatchFileList.size(); i++)
	{
		connector.AddMatchFile(szMatchFileList[i]);
	}
	//connector.AddMatchFile(szMatchFileList[1]);

	return connector.Compute2(szImageList, szFeatDir, szPtsFile);
}

typedef struct
{
	std::string name;
	int id;
	int strip;
}ImageInfo;

bool load_list_file(std::string szImageList, std::vector<ImageInfo> & ImageinfoLists, std::string & szImageDir)
{
	FILE *fp = fopen(szImageList.c_str(), "r");
	if (fp == NULL) return false;

	char fhead[1024] = { '\0' };
	fgets(fhead, 1024, fp);
	if (strstr(fhead, "ImageList") == nullptr) { fclose(fp); return false; }

	for (int i = 0; i < 9; i++)
	{
		fgets(fhead, 1024, fp);

		if (i == 2)
		{
			char chImgDir[256];
			sscanf(fhead, "%s", chImgDir);
			szImageDir.assign(chImgDir);
		}
	}

	int nsize = 0;
	sscanf(fhead, "%d", &nsize);
	ImageinfoLists.clear();
	ImageinfoLists.resize(nsize);

	fgets(fhead, 1024, fp);
	fgets(fhead, 1024, fp);

	char fname[256] = { '\0' };
	for (int i = 0; i < nsize; i++)
	{
		ImageInfo & tinfo = ImageinfoLists[i];
		fgets(fhead, 1024, fp);
		sscanf(fhead, "%s %d %d", fname, &tinfo.strip, &tinfo.id);
		tinfo.name.assign(fname);
	}
	fclose(fp);

	return true;
}

int doTest(int argc, char* argv[])
{
	if (argc <= 1)
	{
		std::cerr << "*** Error: Parameters is error!\n";
		throw std::domain_error("*** Error: Parameters is error!\n");
		return EXIT_FAILURE;
	}
	int featType = 1;
	int matchType = 1;
	int method = 0;
	if (argc >= 4)
	{
		method = std::atoi(argv[3]);
		
	}
	if (method == 0)
	{
		std::cout << "Work use Sift Features!" << std::endl;
	}
	else
	{
		featType = 3;
		matchType = 2;
		std::cout << "Work use Harris Features!" << std::endl;
	}

	int expandSize = 100;
	std::string szImageList = "";
	szImageList.assign(argv[1]);
	std::string szImageDir = "";

	std::vector<ImageInfo> ImageinfoLists;
	if (!load_list_file(szImageList, ImageinfoLists, szImageDir))return EXIT_FAILURE;

	std::string szCurDir = szImageList.substr(0, szImageList.rfind('\\'));
	std::string szFeaDir = szCurDir + "\\feat"; if (_access(szFeaDir.c_str(), 0) != 0) _mkdir(szFeaDir.c_str());
	std::string szMchDir = szCurDir + "\\match"; if (_access(szMchDir.c_str(), 0) != 0) _mkdir(szMchDir.c_str());
	std::string szPtsFile = szCurDir + "\\pts.txt";

	std::vector<std::string> Images;
	std::vector<std::string> TriidxFileList;
	TriidxFileList.resize(ImageinfoLists.size());
	Images.resize(ImageinfoLists.size());
	for (int i = 0; i < ImageinfoLists.size(); i++)
	{
		TriidxFileList[i] = szCurDir + "\\TriMatch\\" + ImageinfoLists[i].name.substr(0, ImageinfoLists[i].name.length() - 4) + ".idx";
		Images[i] = szImageDir + "\\" + ImageinfoLists[i].name;
	}	

	/* extract features */
	std::string featExt = ".feat";
	if (argc >= 3) featExt.assign(argv[2]);
	extrac_fea(Images, szFeaDir, featExt, featType);

	/* match */
	std::vector<std::string> szMatchFileList;
	match(TriidxFileList, szFeaDir, szMchDir, expandSize, szMatchFileList, featExt, matchType);

	/* sort */
	sort_matches(szImageList, szFeaDir, szMatchFileList, szPtsFile, featExt);

	return EXIT_SUCCESS;
}

int match_harris(std::vector<std::string> szFeaList, std::string szMchFile)
{
	if (_access(szMchFile.c_str(), 0) == 0) return EXIT_SUCCESS;

	FeatureMatcher::Options options;
	options.debug_output = true;
	options.min_feature_matches = 8;
	options.matcher_type = FeatureMatcher::MatcherType::MATCHER_HARRIS;
	FeatureMatcher matcher(options);

	return matcher.Match(szFeaList, "", szMchFile, -1);
}

int doTest2()
{
	std::vector<std::string> szImgList;
	szImgList.emplace_back("F:\\Data_In\\20170328\\image\\DJI_0001.JPG");
	szImgList.emplace_back("F:\\Data_In\\20170328\\image\\DJI_0002.JPG");

	extrac_fea(szImgList, "F:\\Data_Work\\feat_harris\\feat", ".feathrs");

	std::vector<std::string> szFeatList;
	szFeatList.emplace_back("F:\\Data_Work\\feat_harris\\feat\\DJI_0001.feathrs");
	szFeatList.emplace_back("F:\\Data_Work\\feat_harris\\feat\\DJI_0002.feathrs");

	return match_harris(szFeatList, "F:\\Data_Work\\feat_harris\\match.txt");
}