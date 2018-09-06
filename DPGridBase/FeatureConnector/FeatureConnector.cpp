#include "FeatureConnector.h"

#include <algorithm>
#include <iostream>
#include <map>
#include <fstream>
#include <iomanip>

#include <io.h>
#include <direct.h>

#include "AerialTriMatchLib.h"

FeatureConnector::FeatureConnector(FeatureConnector::Option const & options):opt_(options)
{
}


FeatureConnector::~FeatureConnector()
{
}

int FeatureConnector::AddMatchFile(std::string const szMatchFile)
{
	if (! fea::fea_mch::isMatchFile(szMatchFile))
	{
		return EXIT_FAILURE;
	}
	/* If already has this file, we skip add it */
	std::vector<std::string>::iterator it;
	it = std::find(this->match_file_list_.begin(), this->match_file_list_.end(), szMatchFile);
	if (it != this->match_file_list_.end())
	{
		return EXIT_SUCCESS;
	}

	this->match_file_list_.emplace_back(szMatchFile);

	return EXIT_SUCCESS;
}

void FeatureConnector::RemoveAll()
{
	this->match_file_list_.clear();
}

#pragma region OTHER_FUNC
/*------------------------ Other Code that needed -------------------*/
typedef struct
{
	std::string name;
	int id;
	int strip;
}ImageInfo;

typedef struct
{
	int featid;
	int imgid;
	int featimgid; // ??? if it is possible, this should be replace by float pos[2]; ?
	float pos[2] = {0.0};
	char dt = 0; // dt > 0 means it will be delete
}FeaPoint;
typedef std::vector<FeaPoint> FeaPoints;

/* We maintain id_on_image-featid maps for all images(feature files), thus we can build tracks */
typedef std::map<int, int> id_id_map; // first is featid on image, second is output id
typedef std::vector<id_id_map> id_id_maps;

typedef std::map<std::string, int> name_id_map; // image(featfile) name and it's correspond id, which will be used to calculate featid
typedef std::map<int, std::pair<std::string, int> > id_name_map; // reverse to name_id_map. just used for quick searching

/* Get the title name of path */
std::string titlename(std::string const path)
{
	std::string basename;
	if (path.rfind('\\') == std::string::npos)
	{
		basename = path;		
	}
	else
	{
		int loc = path.rfind('\\');
		basename = path.substr(loc + 1, path.length() - loc);
	}

	return basename.substr(0, basename.rfind('.'));
}

/* A template finder of find a struct in vector */
template <typename T1, class T2>
class vector_finder
{
public:
	virtual void init(T1 const & value) = 0;
	
	virtual bool operator() (T2 const &value) = 0;

protected:
	T1 val_;
};

/* Find the TwoViewMatchings with same view_1_id and view_2_id */
class find_same_mch_pair : public vector_finder<std::pair<int, int>, fea::TwoViewMatching>
{
public:
	explicit find_same_mch_pair(std::pair<int, int> const & value)
	{
		val_.first = value.first;
		val_.second = value.second;
	}

	void init(std::pair<int, int> const & value) override
	{
		val_.first = value.first;
		val_.second = value.second;
	}

	bool operator() (const std::vector<fea::TwoViewMatching>::value_type &value) override
	{
		return (value.view_1_id == val_.first && value.view_2_id == val_.second);
	}
};

/* Find the matche index with same first id */
class find_matches_same_firstid : public vector_finder<int, fea::CorrespondenceIndex>
{
public:
	explicit find_matches_same_firstid(int const & value)
	{
		val_ = value;
	}

	void init(int const & value) override
	{
		val_ = value;
	}

	bool operator() (const std::vector<fea::CorrespondenceIndex>::value_type &value) override
	{
		return (value.first == val_);
	}
};

/* A pred function used to sort the matches [fea::CorrespondenceIndices] --- less */
bool vector_sort_match_index_by_first(fea::CorrespondenceIndex const & v1, fea::CorrespondenceIndex const & v2)
{
	if (v1.first == v2.first)
	{
		return v1.second < v2.second;
	}
	return v1.first < v2.first;
}

/* A pred function used to sort the pairwise matches [fea::PairwiseMatching] --- less */
bool vector_sort_pair_match_by_viewid1(fea::TwoViewMatching const & v1, fea::TwoViewMatching const & v2)
{
	if (v1.view_1_id == v2.view_1_id)
	{
		return v1.view_2_id < v2.view_2_id;
	}
	return v1.view_1_id < v2.view_1_id;
}

bool vector_sort_feat_point_by_featid(FeaPoint const &v1, FeaPoint const & v2)
{
	if (v1.featid == v2.featid)
	{
		return v1.imgid < v2.imgid;
	}
	return v1.featid < v2.featid;
}

bool vector_sort_feat_point_by_imgid(FeaPoint const &v1, FeaPoint const & v2)
{
	if (v1.imgid == v2.imgid)
	{
		return v1.featid < v2.featid;
	}
	return v1.imgid < v2.imgid;
}

bool vector_sort_feat_point_by_dt(FeaPoint const &v1, FeaPoint const & v2)
{
	if (v1.dt == v2.dt)
	{
		if (v1.imgid == v2.imgid)
		{
			return v1.featid < v2.featid;
		}
		return v1.imgid < v2.imgid;
	}
	return v1.dt < v2.dt;
}

bool load_list_file(std::string szImageList, std::vector<ImageInfo> & ImageinfoLists)
{
	FILE *fp = fopen(szImageList.c_str(), "r");
	if (fp == NULL) return false;

	char fhead[1024] = { '\0' };
	fscanf(fp, "%s\n", fhead);
	if (strstr(fhead, "ImageList") == nullptr) { fclose(fp); return false; }

	for (int i = 0; i < 8; i++)
	{
		fgets(fhead, 1024, fp);
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

/* Build feat_id_map, and generate feat points.
 * NOTE: 1. Please make sure the TwoViewMatching in PairwiseMatching is sorted by view_id_1 in Ascending
 *       //2. Please make sure feat_id_maps has called "resize", and it is recommended that featpoints call "reserve"
 */
void build_feat_id(fea::PairwiseMatching const & all_match_pair, id_name_map const & id_name_maps,
	id_id_maps & feat_id_maps, FeaPoints & featpoints)
{
	//if(feat_id_maps.capacity() <= 0)
	feat_id_maps.resize(id_name_maps.size());
	featpoints.reserve(id_name_maps.size() * 10000);

	for (fea::PairwiseMatching::const_iterator two_mch_it = all_match_pair.begin(); two_mch_it != all_match_pair.end(); two_mch_it++)
	{
		/* get features and maintain the id_id_map */
		const fea::CorrespondenceIndices & match_indices = two_mch_it->matches;

		/* get the id_id_map location by the viewid */
		id_name_map::const_iterator it;
		it = id_name_maps.find(two_mch_it->view_1_id);
		int id_id_map_loc_main = it->second.second;
		it = id_name_maps.find(two_mch_it->view_2_id);
		int id_id_map_loc_match = it->second.second;

		id_id_map & id_id_map_main = feat_id_maps[id_id_map_loc_main];
		id_id_map & id_id_map_match = feat_id_maps[id_id_map_loc_match];
		
		for (fea::CorrespondenceIndices::const_iterator match_it = match_indices.begin(); match_it != match_indices.end(); match_it++)
		{
			/* we treat the image with view_1_id as the main image.
			* A new feat id is calculated from the view_1_id and the feat id on main image.
			* When add feat point, we should first search the corresponding id_id_map.
			* If it already exists on the main image, we do not add the feat to the main image.
			*/
			FeaPoint pt;
	
			/* process the feature point on main image */		
			id_id_map::iterator id_map_it_main = id_id_map_main.find(match_it->first);
			/* if exists, we use the exist feat id */
			if (id_map_it_main != id_id_map_main.end())
			{
				pt.featid = id_map_it_main->second;
			}
			/* if not, we calculate the feat id and add it */
			else
			{
				pt.featid = two_mch_it->view_1_id * 100000 + match_it->first;
				id_id_map_main.insert(id_id_map::value_type(match_it->first, pt.featid));

				pt.imgid = two_mch_it->view_1_id;
				pt.featimgid = match_it->first;
				featpoints.emplace_back(pt);
			}

			/* procee the feature point on match image */		
			id_id_map::iterator id_map_it_match = id_id_map_match.find(match_it->second);
			/* if exists, make sure the feat id is correspond to the small main image id */
			if (id_map_it_match != id_id_map_match.end())
			{
				/* actually , we do nothing */
				/*int old_main_image_id = id_map_it2->second / 100000;
				if (old_main_image_id > two_mch_it->view_1_id)
				{
				id_map_it2->second = pt.featid;
				}*/
			}
			/* if not, we add it */
			else
			{
				id_id_map_match.insert(id_id_map::value_type(match_it->second, pt.featid));

				pt.imgid = two_mch_it->view_2_id;
				pt.featimgid = match_it->second;
				featpoints.emplace_back(pt);
			}
		}
	}
}

int feat_assign_real_coord(FeaPoints & big_feat_points, id_name_map const & id_name_maps, std::string const szFeatDir,
	std::string const szFeatExt = ".feat")
{
	int feat_file_num = id_name_maps.size();
	if (feat_file_num <= 0 || big_feat_points.size() <= 0) return EXIT_FAILURE;

	/* read feat file */
	fea::FeatureSets features;
	std::vector<int> Widths, Heights;
	std::vector<int> Scales;
	features.resize(feat_file_num);
	Widths.resize(feat_file_num); Heights.resize(feat_file_num);
	Scales.resize(feat_file_num);

#pragma omp parallel for schedule(dynamic)
	for (id_name_map::const_iterator it = id_name_maps.begin(); it != id_name_maps.end(); it++)
	{
		int index = it->second.second;
		std::string feafilepath = szFeatDir + "\\" + it->second.first + szFeatExt;

		if (EXIT_FAILURE == fea::fea_mch::readFeatureFile(feafilepath, Scales[index], Widths[index], Heights[index], features[index]))
		{
#pragma omp critical
			{
				std::cerr << "*** Error: Read feature file :" << feafilepath << " error !" << std::endl;
			}
			continue;
		}

		/* clear this to save memory */
		features[index].clear_descriptors();
	}
	Widths.clear(); Widths.shrink_to_fit();
	Heights.clear(); Heights.shrink_to_fit();

	/* get the real coord */
#pragma omp parallel for schedule(dynamic)
#ifdef _MSC_VER
	for (int64_t i = 0; i < big_feat_points.size(); i++)
#else
	for (std::size_t i = 0; i < big_feat_points.size(); i++)
#endif
	{
		FeaPoint & pt = big_feat_points[i];
		id_name_map::const_iterator it = id_name_maps.find(pt.imgid);
		if (it == id_name_maps.end())
		{
			pt.dt = 1;
			continue;
		}
		int sc = Scales[it->second.second];
		//pt.pos[0] = features[it->second.second].positions[pt.featimgid][0] * sc;
		//pt.pos[1] = features[it->second.second].positions[pt.featimgid][1] * sc;
		memcpy(pt.pos, &features[it->second.second].positions[pt.featimgid][0], sizeof(float) * 2);
		pt.pos[0] *= sc;
		pt.pos[1] *= sc;
	}

	/* sort by dt, and erase those dt > 0 */
	std::sort(big_feat_points.begin(), big_feat_points.end(), vector_sort_feat_point_by_dt);
	int remainCnt = 0;
	for (int64_t i = 0; i < big_feat_points.size(); i++)
	{
		FeaPoint const & pt = big_feat_points[i];
		if (pt.dt > 0)
		{
			break;
		}
		remainCnt++;
	}

	big_feat_points.resize(remainCnt); // erase other point

	return EXIT_SUCCESS;
}

#pragma endregion OTHER_FUNC

int FeatureConnector::Compute(std::string szImageListFile, std::string szFeatDir, std::string szPtsFile)
{
	/* 1. Load list file */
	std::vector<ImageInfo> imageinfo_lists;
	if (!load_list_file(szImageListFile, imageinfo_lists))
	{
		std::cerr << "*** Error: Cannot read list file !" << std::endl;
		return EXIT_FAILURE;
	}

	/* Check param */
	int mch_file_num = match_file_list_.size();
	if (mch_file_num <= 0)
	{
		std::cerr << "*** Error: Not enough match file!" << std::endl;
		return EXIT_FAILURE;
	}

	/* 2. Build name-id map */
	name_id_map name_id_maps;
	id_name_map id_name_maps;
	for (std::size_t i = 0; i < imageinfo_lists.size(); i++)
	{
		name_id_maps.insert(name_id_map::value_type(titlename(imageinfo_lists[i].name), imageinfo_lists[i].id));
		id_name_maps.insert(id_name_map::value_type(imageinfo_lists[i].id, std::make_pair(titlename(imageinfo_lists[i].name), i)));
	}

	fea::PairwiseMatching all_mch_pair;
	all_mch_pair.reserve(mch_file_num * 10);
	/* 3. Read match file, and do some pre-process */
	for (std::size_t i = 0; i < mch_file_num; i++)
	{
		/* 3.1. Read file */
		fea::PairwiseMatching pair_mch;
		std::vector<std::string> fea_file_list;
		if (EXIT_FAILURE == fea::fea_mch::readMatchFile(match_file_list_[i], pair_mch, &fea_file_list))
		{
			std::cout << "*** Error: Read match file " << match_file_list_[i] << "error !" << std::endl;
			continue;
		}

		/* 3.2. process the view_id */
#pragma omp parallel for
		for (fea::PairwiseMatching::iterator two_mch_it = pair_mch.begin(); two_mch_it != pair_mch.end(); )
		{
			/* correspond feature file name */
			std::string featname1 = fea_file_list[two_mch_it->view_1_id];
			std::string featname2 = fea_file_list[two_mch_it->view_2_id];

			/* find correspond id and update it */
			name_id_map::iterator id_it1 = name_id_maps.find(titlename(featname1));
			name_id_map::iterator id_it2 = name_id_maps.find(titlename(featname2));
			if (id_it1 == name_id_maps.end() || id_it2 == name_id_maps.end()) /* did not find, then we erase this match */
			{
				two_mch_it = pair_mch.erase(two_mch_it);
				continue;
			}

			/* if we find id correspond to the list, then we make sure the view_1_id point to the smaller id */
			if (id_it1->second < id_it2->second)
			{
				two_mch_it->view_1_id = id_it1->second;
				two_mch_it->view_2_id = id_it2->second;
			}
			else
			{
				two_mch_it->view_1_id = id_it2->second;
				two_mch_it->view_2_id = id_it1->second;

				/* we also need to swap the feat id */
				for (std::size_t j = 0; j < two_mch_it->matches.size(); j++)
				{
					fea::CorrespondenceIndex & cor_index = two_mch_it->matches[j];
					std::swap(cor_index.first, cor_index.second);
				}

				/* sort the feat by id */
				std::sort(two_mch_it->matches.begin(), two_mch_it->matches.end(), vector_sort_match_index_by_first);
			}
			
			/* raise the iterator */
			two_mch_it++;
		}
		/* force to clear the memory*/
		fea_file_list.swap(std::vector<std::string>());

		/* 3.3. Add to big vector */
		for (fea::PairwiseMatching::iterator two_mch_it = pair_mch.begin(); two_mch_it != pair_mch.end(); )
		{
			/* find whether it is already exist */
			std::pair<int, int> find_pair = std::make_pair(two_mch_it->view_1_id, two_mch_it->view_2_id);
			fea::PairwiseMatching::iterator big_two_mch_it = std::find_if(all_mch_pair.begin(), all_mch_pair.end(), find_same_mch_pair(find_pair));

			/* if not exist, add to end */
			if (big_two_mch_it == all_mch_pair.end())
			{
				all_mch_pair.emplace_back(*two_mch_it);
			}
			/* else, add to the exist pair */
			else
			{
				/* if the new match index's first feature id is bigger than that of the last match in big vector, then we can add the new match directly, 
				 * otherwise, we must take into consideration that those two vectors may have same first index !
				*/
				std::size_t big_size = big_two_mch_it->matches.size();
				for (fea::CorrespondenceIndices::iterator matchs_it = two_mch_it->matches.begin(); matchs_it != two_mch_it->matches.end(); )
				{
					/* add it directly */
					if (big_two_mch_it->matches[big_size - 1].first < matchs_it->first || big_two_mch_it->matches[0].first > matchs_it->first)
					{
						big_two_mch_it->matches.emplace_back(*matchs_it);
					}
					/* search if the big vector has the same first index */
					else
					{
						fea::CorrespondenceIndices::iterator big_matchs_it = std::find_if(big_two_mch_it->matches.begin(), big_two_mch_it->matches.begin() + big_size, find_matches_same_firstid(matchs_it->first));

						/* if find */
						if (big_matchs_it != (big_two_mch_it->matches.begin() + big_size))
						{
							/* we skip to add it !
							* but if the second index is not same, then this matches is not robust and it should be erased from the big vector
							*/
							if (big_matchs_it->second != matchs_it->second)
							{
								big_two_mch_it->matches.erase(big_matchs_it);
								//big_size--;
							}
						}
						/* if not find, we add it */
						else
						{
							big_two_mch_it->matches.emplace_back(*matchs_it);
						}
					}

					/* raise the iterator */
					matchs_it++;
				}

				/* after add new matches, we need sort the feat by id */
				std::sort(big_two_mch_it->matches.begin(), big_two_mch_it->matches.end(), vector_sort_match_index_by_first);
			}

			/* raise the iterator */
			two_mch_it++;
		}
		/* force to clear the memory*/
		pair_mch.swap(std::vector<fea::TwoViewMatching>());				
	}
	
	/* 4. Build featid */
	id_id_maps feat_id_maps;
	feat_id_maps.resize(imageinfo_lists.size());
	FeaPoints big_feat_points;
	big_feat_points.reserve(imageinfo_lists.size() * 100000);

	/* sort the pairwise_matching vector by viewid1 */
	std::sort(all_mch_pair.begin(), all_mch_pair.end(), vector_sort_pair_match_by_viewid1);
	build_feat_id(all_mch_pair, id_name_maps, feat_id_maps, big_feat_points);

	feat_id_maps.clear(); feat_id_maps.shrink_to_fit();

	/* 5. Read feature file, and save finnaly pts file with feat coord */
	feat_assign_real_coord(big_feat_points, id_name_maps, szFeatDir, this->opt_.szFeatfileExt);
	
	if (this->opt_.debug_output)
	{
		/* Create a temp folder under the Parent dir of szPtsFile */
		std::string szParentDir = szPtsFile.substr(0, szPtsFile.rfind('\\'));
		std::string szTempDir = szParentDir + "\\temp0";
		if (_access(szTempDir.c_str(), 0) != 0) _mkdir(szTempDir.c_str());

		/* save to different file */
		{
			std::string szDebugFile = "";

			int old_image_id = -1;
			FILE* fp = nullptr;
			for (std::size_t j = 0; j < big_feat_points.size(); j++)
			{
				if (old_image_id != big_feat_points[j].imgid)
				{
					if (nullptr != fp)
					{
						fclose(fp);
					}
					old_image_id = big_feat_points[j].imgid;
					id_name_map::iterator it = id_name_maps.find(old_image_id);
					szDebugFile = szTempDir + "\\pt_" + it->second.first + ".txt";
					fp = fopen(szDebugFile.c_str(), "w");
					if (nullptr == fp)break;
				}

				id_name_map::iterator it = id_name_maps.find(big_feat_points[j].featid / 100000);
				fprintf(fp, "%10d %8.2f %8.2f %8.2f %4d\n",
					big_feat_points[j].featid, big_feat_points[j].pos[0], big_feat_points[j].pos[1], 0.0f, it->second.second);

				if (j == big_feat_points.size() - 1)
				{
					fclose(fp);
					break;
				}
			}
		}
	}

	/*if (this->opt_.debug_output)
	{
		std::string szDebugFile = szPtsFile + ".debug.txt";

		std::ofstream out(szDebugFile);
		if (out.good())
		{
			out.flags(std::ios::left);
			out << std::setw(5) << all_mch_pair.size() << std::endl;

			for (fea::PairwiseMatching::iterator two_mch_it = all_mch_pair.begin(); two_mch_it != all_mch_pair.end(); two_mch_it++)
			{
				out << two_mch_it->view_1_id << "  " << two_mch_it->view_2_id << "   " << two_mch_it->matches.size() << std::endl;

				fea::CorrespondenceIndices  & mch_indices = two_mch_it->matches;

				for (fea::CorrespondenceIndices::iterator match_it = mch_indices.begin(); match_it != mch_indices.end(); match_it++)
				{
					out << match_it->first << "  " << match_it->second << std::endl;
				}
				out << std::endl;
			}

			out.close();
		}
	}*/

	return EXIT_SUCCESS;
}

int FeatureConnector::Compute2(std::string szImageListFile, std::string szFeatDir, std::string szPtsFile)
{
	/* 1. Load list file */
	std::vector<ImageInfo> imageinfo_lists;
	if (!load_list_file(szImageListFile, imageinfo_lists))
	{
		std::cerr << "*** Error: Cannot read list file !" << std::endl;
		return EXIT_FAILURE;
	}

	/* Check param */
	int mch_file_num = match_file_list_.size();
	if (mch_file_num <= 0)
	{
		std::cerr << "*** Error: Not enough match file!" << std::endl;
		return EXIT_FAILURE;
	}

	/* 2. Build name-id map */
	name_id_map name_id_maps;
	id_name_map id_name_maps;
	for (std::size_t i = 0; i < imageinfo_lists.size(); i++)
	{
		name_id_maps.insert(name_id_map::value_type(titlename(imageinfo_lists[i].name), imageinfo_lists[i].id));
		id_name_maps.insert(id_name_map::value_type(imageinfo_lists[i].id, std::make_pair(titlename(imageinfo_lists[i].name), i)));
	}

	/* used to store all points */
	FeaPoints big_feat_points;
	big_feat_points.reserve(mch_file_num * 100000);

	/* 3. Read match file, and do some pre-process */
#pragma omp parallel for schedule(dynamic)
#ifdef _MSC_VER
	for (int64_t i = 0; i < mch_file_num; i++)
#else
	for (std::size_t i = 0; i < mch_file_num; i++)
#endif
	{
		/* 3.1. Read match file */
		fea::PairwiseMatching pair_mch;
		std::vector<std::string> fea_file_list;
		if (EXIT_FAILURE == fea::fea_mch::readMatchFile(match_file_list_[i], pair_mch, &fea_file_list))
		{
			std::cout << "*** Error: Read match file " << match_file_list_[i] << "error !" << std::endl;
			continue;
		}
		
		/* 3.2. process the view_id */
		for (fea::PairwiseMatching::iterator two_mch_it = pair_mch.begin(); two_mch_it != pair_mch.end(); )
		{
			/* correspond feature file name */
			std::string featname1 = fea_file_list[two_mch_it->view_1_id];
			std::string featname2 = fea_file_list[two_mch_it->view_2_id];

			/* find correspond id and update it */
			name_id_map::iterator id_it1 = name_id_maps.find(titlename(featname1));
			name_id_map::iterator id_it2 = name_id_maps.find(titlename(featname2));
			if (id_it1 == name_id_maps.end() || id_it2 == name_id_maps.end()) /* did not find, then we erase this match */
			{
				two_mch_it = pair_mch.erase(two_mch_it);
				continue;
			}

			/* reset view id according to the name_id_map */
			two_mch_it->view_1_id = id_it1->second;
			two_mch_it->view_2_id = id_it2->second;		

			/* raise the iterator */
			two_mch_it++;
		}

		/* 3.4 build id_id_map and get feat */
		id_id_maps feat_id_maps;
		FeaPoints feats_match;
		id_name_map tmp_id_name_map;
		for (std::size_t j = 0; j < fea_file_list.size(); j++)
		{
			name_id_map::iterator it = name_id_maps.find(titlename(fea_file_list[j]));
			if (it == name_id_maps.end())continue;
			tmp_id_name_map.insert(id_name_map::value_type(it->second, std::make_pair(titlename(fea_file_list[j]), j)));
		}
		feat_id_maps.resize(tmp_id_name_map.size());
		feats_match.reserve(tmp_id_name_map.size() * 10000);

		/* sort the pairwise_matching vector by viewid1 */
		//std::sort(pair_mch.begin(), pair_mch.end(), vector_sort_pair_match_by_viewid1);
		build_feat_id(pair_mch, tmp_id_name_map, feat_id_maps, feats_match);

		feat_id_maps.clear(); feat_id_maps.shrink_to_fit();
		tmp_id_name_map.clear();

		/* 3.5 save the feat */		
		/* save to the big vector */
#pragma omp critical
		{
			big_feat_points.insert(big_feat_points.end(), feats_match.begin(), feats_match.end());
		}		

		/* force to clear the memory*/
		feats_match.clear(); feats_match.shrink_to_fit();
		fea_file_list.swap(std::vector<std::string>());
		pair_mch.swap(std::vector<fea::TwoViewMatching>());		
	}
	/* sort the big vector by feat image id */
	//std::sort(big_feat_points.begin(), big_feat_points.end(), vector_sort_feat_point_by_imgid);

	/* 4. Read feature file, and save finnaly pts file with feat coord */
	feat_assign_real_coord(big_feat_points, id_name_maps, szFeatDir, this->opt_.szFeatfileExt);

	if (this->opt_.debug_output)
	{
		/* Create a temp folder under the Parent dir of szPtsFile */
		std::string szParentDir = szPtsFile.substr(0, szPtsFile.rfind('\\'));
		std::string szTempDir = szParentDir + "\\temp2";
		if (_access(szTempDir.c_str(), 0) != 0) _mkdir(szTempDir.c_str());

		/* save to different file */
		{
			std::string szDebugFile = "";
			
			int old_image_id = -1;
			FILE* fp = nullptr;
			for (std::size_t j = 0; j < big_feat_points.size(); j++)
			{
				if (old_image_id != big_feat_points[j].imgid)
				{
					if (nullptr != fp)
					{
						fclose(fp);
					}
					old_image_id = big_feat_points[j].imgid;
					id_name_map::iterator it = id_name_maps.find(old_image_id);
					szDebugFile = szTempDir + "\\pt_"  + it->second.first + ".txt";
					fp = fopen(szDebugFile.c_str(), "w");
					if (nullptr == fp)break;
				}

				id_name_map::iterator it = id_name_maps.find(big_feat_points[j].featid / 100000);
				fprintf(fp, "%10d %8.2f %8.2f %8.2f %4d\n",
					big_feat_points[j].featid, big_feat_points[j].pos[0], big_feat_points[j].pos[1], 0.0f, it->second.second);

				if (j == big_feat_points.size() - 1)
				{
					fclose(fp);
					break;
				}
			}
		}
	}

	return EXIT_SUCCESS;
}

int FeatureConnector::Compute3(std::string szImageListFile, std::string szFeatDir, std::string szPtsFile)
{
	/* 1. Load list file */
	std::vector<ImageInfo> imageinfo_lists;
	if (!load_list_file(szImageListFile, imageinfo_lists))
	{
		std::cerr << "*** Error: Cannot read list file !" << std::endl;
		return EXIT_FAILURE;
	}

	/* Check param */
	int mch_file_num = match_file_list_.size();
	if (mch_file_num <= 0)
	{
		std::cerr << "*** Error: Not enough match file!" << std::endl;
		return EXIT_FAILURE;
	}

	/* 2. Build name-id map */
	name_id_map name_id_maps;
	id_name_map id_name_maps;
	for (std::size_t i = 0; i < imageinfo_lists.size(); i++)
	{
		name_id_maps.insert(name_id_map::value_type(titlename(imageinfo_lists[i].name), imageinfo_lists[i].id));
		id_name_maps.insert(id_name_map::value_type(imageinfo_lists[i].id, std::make_pair(titlename(imageinfo_lists[i].name), i)));
	}

	/* used to store all points */
	FeaPoints big_feat_points;
	big_feat_points.reserve(mch_file_num * 100000);
	id_id_maps feat_id_maps;
	feat_id_maps.resize(mch_file_num);

	/* 3. Read match file, and do some pre-process */
	for (std::size_t i = 0; i < mch_file_num; i++)
	{
		/* 3.1. Read match file */
		fea::PairwiseMatching pair_mch;
		std::vector<std::string> fea_file_list;
		if (EXIT_FAILURE == fea::fea_mch::readMatchFile(match_file_list_[i], pair_mch, &fea_file_list))
		{
			std::cout << "*** Error: Read match file " << match_file_list_[i] << "error !" << std::endl;
			continue;
		}

		/* 3.2. process the view_id */
#pragma omp parallel for
		for (fea::PairwiseMatching::iterator two_mch_it = pair_mch.begin(); two_mch_it != pair_mch.end(); )
		{
			/* correspond feature file name */
			std::string featname1 = fea_file_list[two_mch_it->view_1_id];
			std::string featname2 = fea_file_list[two_mch_it->view_2_id];

			/* find correspond id and update it */
			name_id_map::iterator id_it1 = name_id_maps.find(titlename(featname1));
			name_id_map::iterator id_it2 = name_id_maps.find(titlename(featname2));
			if (id_it1 == name_id_maps.end() || id_it2 == name_id_maps.end()) /* did not find, then we erase this match */
			{
				two_mch_it = pair_mch.erase(two_mch_it);
				continue;
			}

			/* reset view id according to the name_id_map */
			two_mch_it->view_1_id = id_it1->second;
			two_mch_it->view_2_id = id_it2->second;

			/* raise the iterator */
			two_mch_it++;
		}

		/* 3.4 build id_id_map and get feat */
		/* sort the pairwise_matching vector by viewid1 */
		//std::sort(pair_mch.begin(), pair_mch.end(), vector_sort_pair_match_by_viewid1);
		build_feat_id(pair_mch, id_name_maps, feat_id_maps, big_feat_points);	

		/* force to clear the memory*/
		fea_file_list.swap(std::vector<std::string>());
		pair_mch.swap(std::vector<fea::TwoViewMatching>());
	}
	feat_id_maps.clear(); feat_id_maps.shrink_to_fit();

	/* sort the big vector by feat image id */
	//std::sort(big_feat_points.begin(), big_feat_points.end(), vector_sort_feat_point_by_imgid);

	/* 4. Read feature file, and save finnaly pts file with feat coord */
	feat_assign_real_coord(big_feat_points, id_name_maps, szFeatDir, this->opt_.szFeatfileExt);

	if (this->opt_.debug_output)
	{
		/* Create a temp folder under the Parent dir of szPtsFile */
		std::string szParentDir = szPtsFile.substr(0, szPtsFile.rfind('\\'));
		std::string szTempDir = szParentDir + "\\temp3";
		if (_access(szTempDir.c_str(), 0) != 0) _mkdir(szTempDir.c_str());

		/* save to different file */
		{
			std::string szDebugFile = "";

			int old_image_id = -1;
			FILE* fp = nullptr;
			for (std::size_t j = 0; j < big_feat_points.size(); j++)
			{
				if (old_image_id != big_feat_points[j].imgid)
				{
					if (nullptr != fp)
					{
						fclose(fp);
					}
					old_image_id = big_feat_points[j].imgid;
					id_name_map::iterator it = id_name_maps.find(old_image_id);
					szDebugFile = szTempDir + "\\pt_" + it->second.first + ".txt";
					fp = fopen(szDebugFile.c_str(), "w");
					if (nullptr == fp)break;
				}

				id_name_map::iterator it = id_name_maps.find(big_feat_points[j].featid / 100000);
				fprintf(fp, "%10d %8.2f %8.2f %8.2f %4d\n",
					big_feat_points[j].featid, big_feat_points[j].pos[0], big_feat_points[j].pos[1], 0.0f, it->second.second);

				if (j == big_feat_points.size() - 1)
				{
					fclose(fp);
					break;
				}
			}
		}
	}

	return EXIT_SUCCESS;
}

//int FeatureConnector::Compute2(std::string szImageListFile, std::string szFeatDir, std::string szPtsFile)
//{
//	/* 1. Load list file */
//	std::vector<ImageInfo> imageinfo_lists;
//	if (!load_list_file(szImageListFile, imageinfo_lists))
//	{
//		std::cerr << "*** Error: Cannot read list file !" << std::endl;
//		return EXIT_FAILURE;
//	}
//
//	/* Check param */
//	int mch_file_num = match_file_list_.size();
//	if (mch_file_num <= 0)
//	{
//		std::cerr << "*** Error: Not enough match file!" << std::endl;
//		return EXIT_FAILURE;
//	}
//
//	/* 2. Build name-id map */
//	name_id_map name_id_maps;
//	id_name_map id_name_maps;
//	for (std::size_t i = 0; i < imageinfo_lists.size(); i++)
//	{
//		name_id_maps.insert(name_id_map::value_type(titlename(imageinfo_lists[i].name), imageinfo_lists[i].id));
//		id_name_maps.insert(id_name_map::value_type(imageinfo_lists[i].id, std::make_pair(titlename(imageinfo_lists[i].name), i)));
//	}
//
//	/* used to store all points */
//	FeaPoints big_feat_points;
//	big_feat_points.reserve(mch_file_num * 100000);
//
//	/* 3. Read match file, and do some pre-process */
//#pragma omp parallel for schedule(dynamic)
//#ifdef _MSC_VER
//	for (int64_t i = 0; i < mch_file_num; i++)
//#else
//	for (std::size_t i = 0; i < mch_file_num; i++)
//#endif
//	{
//		/* 3.1. Read match file */
//		fea::PairwiseMatching pair_mch;
//		std::vector<std::string> fea_file_list;
//		if (EXIT_FAILURE == fea::fea_mch::readMatchFile(match_file_list_[i], pair_mch, &fea_file_list))
//		{
//			std::cout << "*** Error: Read match file " << match_file_list_[i] << "error !" << std::endl;
//			continue;
//		}
//
//		id_id_maps feat_id_maps;
//		feat_id_maps.resize(fea_file_list.size());
//		FeaPoints feats_match;
//		feats_match.reserve(fea_file_list.size() * 10000);
//
//		/* 3.2. process the view_id */
//		for (fea::PairwiseMatching::iterator two_mch_it = pair_mch.begin(); two_mch_it != pair_mch.end(); )
//		{
//			/* correspond feature file name */
//			std::string featname1 = fea_file_list[two_mch_it->view_1_id];
//			std::string featname2 = fea_file_list[two_mch_it->view_2_id];
//
//			/* find correspond id and update it */
//			name_id_map::iterator id_it1 = name_id_maps.find(titlename(featname1));
//			name_id_map::iterator id_it2 = name_id_maps.find(titlename(featname2));
//			if (id_it1 == name_id_maps.end() || id_it2 == name_id_maps.end()) /* did not find, then we erase this match */
//			{
//				two_mch_it = pair_mch.erase(two_mch_it);
//				continue;
//			}
//
//			/* store the index of name_id_map, we will use it to find the right location of id_id_map */
//			std::size_t id_id_map_loc1 = two_mch_it->view_1_id;
//			std::size_t id_id_map_loc2 = two_mch_it->view_2_id;
//
//			/* reset view id according to the name_id_map */
//			two_mch_it->view_1_id = id_it1->second;
//			two_mch_it->view_2_id = id_it2->second;
//
//			/* get features and maintain the id_id_map */
//			fea::CorrespondenceIndices & match_indices = two_mch_it->matches;
//			for (fea::CorrespondenceIndices::iterator match_it = match_indices.begin(); match_it != match_indices.end(); match_it++)
//			{
//				/* we treat the image with view_1_id as the main image.
//				* A new feat id is calculated from the view_1_id and the feat id on main image.
//				* When add feat point, we should first search the corresponding id_id_map.
//				* If it already exists on the main image, we do not add the feat to the main image.
//				*/
//				FeaPoint pt;
//
//				/* process the feature point on main image */
//				id_id_map & id_id_map_main = feat_id_maps[id_id_map_loc1];
//				id_id_map::iterator id_map_it = id_id_map_main.find(match_it->first);
//				/* if exists, we use the exist feat id */
//				if (id_map_it != id_id_map_main.end())
//				{
//					pt.featid = id_map_it->second;
//				}
//				/* if not, we calculate the feat id and add it */
//				else
//				{
//					pt.featid = two_mch_it->view_1_id * 100000 + match_it->first;
//					id_id_map_main.insert(id_id_map::value_type(match_it->first, pt.featid));
//
//					pt.imgid = two_mch_it->view_1_id;
//					pt.featimgid = match_it->first;
//					//pt.pos[0] = features[id_id_map_loc1].positions[match_it->first][0] * Scales[id_id_map_loc1];
//					//pt.pos[1] = features[id_id_map_loc1].positions[match_it->first][1] * Scales[id_id_map_loc1];
//
//					feats_match.emplace_back(pt);
//				}
//
//				/* procee the feature point on match image */
//				id_id_map & id_id_map_match = feat_id_maps[id_id_map_loc2];
//				id_id_map::iterator id_map_it2 = id_id_map_match.find(match_it->second);
//				/* if exists, make sure the feat id is correspond to the small main image id */
//				if (id_map_it2 != id_id_map_match.end())
//				{
//					/* actually , we do nothing */
//					/*int old_main_image_id = id_map_it2->second / 100000;
//					if (old_main_image_id > two_mch_it->view_1_id)
//					{
//					id_map_it2->second = pt.featid;
//					}*/
//				}
//				/* if not, we add it */
//				else
//				{
//					id_id_map_match.insert(id_id_map::value_type(match_it->second, pt.featid));
//
//					pt.imgid = two_mch_it->view_2_id;
//					pt.featimgid = match_it->second;
//					//pt.pos[0] = features[id_id_map_loc2].positions[match_it->second][0] * Scales[id_id_map_loc2];
//					//pt.pos[1] = features[id_id_map_loc2].positions[match_it->second][1] * Scales[id_id_map_loc2];
//					feats_match.emplace_back(pt);
//				}
//			}
//
//			/* raise the iterator */
//			two_mch_it++;
//		}
//		/* clear to save memory */
//		feat_id_maps.clear(); feat_id_maps.shrink_to_fit();
//
//		/* 3.3 sort and save the feat */
//		/* sort the feat by image id */
//		std::sort(feats_match.begin(), feats_match.end(), vector_sort_feat_point_by_imgid);
//		/* save to the big vector */
//#pragma omp critical
//		{
//			big_feat_points.insert(big_feat_points.end(), feats_match.begin(), feats_match.end());
//		}
//
//		feats_match.clear(); feats_match.shrink_to_fit();
//
//		/* force to clear the memory*/
//		fea_file_list.swap(std::vector<std::string>());
//		pair_mch.swap(std::vector<fea::TwoViewMatching>());
//	}
//	/* sort the big vector by feat image id */
//	//std::sort(big_feat_points.begin(), big_feat_points.end(), vector_sort_feat_point_by_imgid);
//
//	/* 4. Read feature file, and save finnaly pts file with feat coord */
//	/* read feat file */
//	fea::FeatureSets features;
//	std::vector<int> Widths, Heights;
//	std::vector<int> Scales;
//	features.resize(imageinfo_lists.size());
//	Widths.resize(imageinfo_lists.size()); Heights.resize(imageinfo_lists.size());
//	Scales.resize(imageinfo_lists.size());
//#pragma omp parallel for schedule(dynamic)
//#ifdef _MSC_VER
//	for (int64_t i = 0; i < imageinfo_lists.size(); i++)
//#else
//	for (std::size_t i = 0; i < imageinfo_lists.size(); i++)
//#endif
//	{
//		std::string feafilepath = szFeatDir + "\\" + titlename(imageinfo_lists[i].name) + this->opt_.szFeatfileExt;
//
//		if (EXIT_FAILURE == fea::fea_mch::readFeatureFile(feafilepath, Scales[i], Widths[i], Heights[i], features[i]))
//		{
//#pragma omp critical
//			{
//				std::cerr << "*** Error: Read feature file :" << feafilepath << " error !" << std::endl;
//			}
//			continue;
//		}
//
//		/* clear this to save memory */
//		features[i].clear_descriptors();
//	}
//
//	/* get the real coord */
//#pragma omp parallel for schedule(dynamic)
//#ifdef _MSC_VER
//	for (int64_t i = 0; i < big_feat_points.size(); i++)
//#else
//	for (std::size_t i = 0; i < big_feat_points.size(); i++)
//#endif
//	{
//		FeaPoint & pt = big_feat_points[i];
//		id_name_map::iterator it = id_name_maps.find(pt.imgid);
//		if (it == id_name_maps.end())
//		{
//			pt.dt = 1;
//			continue;
//		}
//		memcpy(pt.pos, &features[it->second.second].positions[pt.featimgid][0], sizeof(float) * 2);
//	}
//	/* sort by dt, and erase those dt > 0 */
//	std::sort(big_feat_points.begin(), big_feat_points.end(), vector_sort_feat_point_by_dt);
//	int remainCnt = 0;
//	for (int64_t i = 0; i < big_feat_points.size(); i++)
//	{
//		FeaPoint const & pt = big_feat_points[i];
//		if (pt.dt > 0)
//		{
//			break;
//		}
//		remainCnt++;
//	}
//	big_feat_points.resize(remainCnt); // erase other point
//
//	if (this->opt_.debug_output)
//	{
//		/* Create a temp folder under the Parent dir of szPtsFile */
//		std::string szParentDir = szPtsFile.substr(0, szPtsFile.rfind('\\'));
//		std::string szTempDir = szParentDir + "\\temp";
//		if (_access(szTempDir.c_str(), 0) != 0) _mkdir(szTempDir.c_str());
//
//		/* save to different file */
//		{
//			std::string szDebugFile = "";
//
//			int old_image_id = -1;
//			FILE* fp = nullptr;
//			for (std::size_t j = 0; j < big_feat_points.size(); j++)
//			{
//				if (old_image_id != big_feat_points[j].imgid)
//				{
//					if (nullptr != fp)
//					{
//						fclose(fp);
//					}
//					old_image_id = big_feat_points[j].imgid;
//					id_name_map::iterator it = id_name_maps.find(old_image_id);
//					szDebugFile = szTempDir + "\\pt_" + it->second.first + ".txt";
//					fp = fopen(szDebugFile.c_str(), "w");
//					if (nullptr == fp)break;
//				}
//
//				id_name_map::iterator it = id_name_maps.find(big_feat_points[j].featid / 100000);
//				fprintf(fp, "%10d %8.2f %8.2f %8.2f %4d\n",
//					big_feat_points[j].featid, big_feat_points[j].pos[0], big_feat_points[j].pos[1], 0.0f, it->second.second);
//
//				if (j == big_feat_points.size() - 1)
//				{
//					fclose(fp);
//					break;
//				}
//			}
//		}
//	}
//
//	return EXIT_SUCCESS;
//}
