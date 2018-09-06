#ifndef DPGRID_COMMON_HEADER
#define DPGRID_COMMON_HEADER

#include <vector>
#include <fstream>
#include <iostream>

#include "feature\defines.h"
#include "feature\feature_set.h"

FEA_NAMESPACE_BEGIN

struct Correspondence2D2D;
typedef std::vector<Correspondence2D2D> Correspondences2D2D;

struct Correspondence2D3D;
typedef std::vector<Correspondence2D3D> Correspondences2D3D;

/** The IDs of a matching feature pair in two images. */
typedef std::pair<int, int> CorrespondenceIndex;
/** A list of all matching feature pairs in two images. */
typedef std::vector<CorrespondenceIndex> CorrespondenceIndices;

/**
* Two image coordinates which correspond to each other in terms of observing
* the same point in the scene.
*/
struct Correspondence2D2D
{
	double p1[2];
	double p2[2];
};

/**
* A 3D point and an image coordinate which correspond to each other in terms
* of the image observing this 3D point in the scene.
*/
struct Correspondence2D3D
{
	double p3d[3];
	double p2d[2];
};

/** The matching result between two views. */
struct TwoViewMatching
{
	bool operator< (TwoViewMatching const& rhs) const;

	int view_1_id;
	int view_2_id;
	CorrespondenceIndices matches;
};

/** The matching result between several pairs of views. */
typedef std::vector<TwoViewMatching> PairwiseMatching;

/* ------------------------ Implementation ------------------------ */

inline bool
TwoViewMatching::operator< (TwoViewMatching const& rhs) const
{
	return this->view_1_id == rhs.view_1_id
		? this->view_2_id < rhs.view_2_id
		: this->view_1_id < rhs.view_1_id;
}

FEA_MATCH_NAMESPACE_BEGIN

#ifndef FeatureFile_Header
#define FeatureFile_Header "Feature_Detect_File\n"
#define FeatureFile_Header_Len 20
#endif // !FeatureFile_Header

#ifndef MatchFile_Header
#define MatchFile_Header "Feature_Match_File\n"
#define MatchFile_Header_LEN 19
#endif // !MatchFile_Header

/* save feature file */
inline
int saveFeatureFile(std::string szFeaturePath, int Scale, int FeaW, int FeaH, fea::FeatureSet const& features);

inline
int saveFeatureFileTxt(std::string szFeaturePath, int Scale, int FeaW, int FeaH, fea::FeatureSet const& features);

/* read feature file */
inline
bool isFeatureFile(std::string szFeaturePath);

inline
int readFeatureFile(std::string szFeaturePath, int & Scale, int & FeaW, int & FeaH, fea::FeatureSet & features);

/* save match file */
inline
int saveMatchFile(std::string szMatchFile, fea::PairwiseMatching const & matching, std::vector<std::string> const & feat_file_list);

inline
int saveMatchFileTxt(std::string szMatchFile, fea::PairwiseMatching const & matching, std::vector<std::string> const & feat_file_list, fea::FeatureSets * features = nullptr);

/* read match file */
inline
bool isMatchFile(std::string szMatchFile);

inline
int readMatchFile(std::string szMatchFile, fea::PairwiseMatching & matching, std::vector<std::string> * feat_file_list = nullptr);

/* ------------------------ Implementation ------------------------ */

inline
std::string basename(std::string const & path)
{
	int len = path.length();
	int last_backslashed = path.rfind('\\');

	return path.substr(last_backslashed + 1, len - last_backslashed);
}

template <class T>
void ClearVector(std::vector<T> &vec)
{
	vec.clear();
	vec.shrink_to_fit();
}

inline
int saveFeatureFile(std::string szFeaturePath, int Scale, int FeaW, int FeaH, fea::FeatureSet const& features)
{
	std::ofstream out(szFeaturePath, std::ios::binary);
	if (!out.good()) { std::cerr << "*** Error: Can't save file: " << szFeaturePath << std::endl; return EXIT_FAILURE; }

	/* write header */
	out.write(FeatureFile_Header, FeatureFile_Header_Len);

	/* Write scale and scaled image size */
	int32_t sc = Scale;	out.write(reinterpret_cast<char const*>(&sc), sizeof(int32_t));
	int32_t sw = FeaW;	out.write(reinterpret_cast<char const*>(&sw), sizeof(int32_t));
	int32_t sh = FeaH;	out.write(reinterpret_cast<char const*>(&sh), sizeof(int32_t));

	/* Write number of features. */
	int32_t num_positions = static_cast<int32_t>(features.positions.size());
	out.write(reinterpret_cast<char const*>(&num_positions), sizeof(int32_t));

	/* Write number of sift features. */
	int32_t num_sifts = static_cast<int32_t>(features.sift_descriptors.size());
	out.write(reinterpret_cast<char const*>(&num_sifts), sizeof(int32_t));

	/* Write number of surf features. */
	int32_t num_surfs = static_cast<int32_t>(features.surf_descriptors.size());
	out.write(reinterpret_cast<char const*>(&num_surfs), sizeof(int32_t));

	/* Write number of harris features. */
	int32_t num_harris = static_cast<int32_t>(features.harris_descriptors.size());
	out.write(reinterpret_cast<char const*>(&num_harris), sizeof(int32_t));

	/* Write positions and colors . */
	for (std::size_t i = 0; i < num_positions; i++)
	{
		math::Vec2f const& pos = features.positions[i];
		out.write(reinterpret_cast<char const*>(&pos), sizeof(math::Vec2f));

		math::Vec3uc const& col = features.colors[i];
		out.write(reinterpret_cast<char const*>(&col), sizeof(math::Vec3uc));
	}

	/* Write sift descriptors. */
	for (std::size_t i = 0; i < num_sifts; i++)
	{
		fea::Sift::Descriptor const& des = features.sift_descriptors[i];
		out.write(reinterpret_cast<char const*>(&des), sizeof(fea::Sift::Descriptor));
	}

	/* Write surf descriptors. */
	for (std::size_t i = 0; i < num_surfs; i++)
	{
		fea::Surf::Descriptor const& des = features.surf_descriptors[i];
		out.write(reinterpret_cast<char const*>(&des), sizeof(fea::Surf::Descriptor));
	}

	/* Write harris descriptors. */
	for (std::size_t i = 0; i < num_harris; i++)
	{
		fea::Harris::Descriptor const& des = features.harris_descriptors[i];
		out.write(reinterpret_cast<char const*>(&des), sizeof(fea::Harris::Descriptor));
	}

	/* write feature name */
	std::string feature_name = basename(szFeaturePath);
	int32_t len = feature_name.length();
	out.write(reinterpret_cast<char const*>(&len), sizeof(int32_t));
	out.write(feature_name.c_str(), len);

	out.close(); std::cout << "\rSave file: " << szFeaturePath << " Success" << std::endl;

	return EXIT_SUCCESS;
}

inline
int saveFeatureFileTxt(std::string szFeaturePath, int Scale, int FeaW, int FeaH, fea::FeatureSet const & features)
{
	std::ofstream out(szFeaturePath);
	if (!out.good()) { std::cerr << "*** Error: Can't save file: " << szFeaturePath << std::endl; return EXIT_FAILURE; }

	/* write header */
	out << FeatureFile_Header;

	/* Write scale and scaled image size */
	out.flags(std::ios::right);
	int32_t sc = Scale; out << std::setw(5) << sc; 
	int32_t sw = FeaW; out << std::setw(5) << sw; 
	int32_t sh = FeaH; out << std::setw(5) << sh << std::endl; 

	/* Write number of features. */
	int32_t num_positions = static_cast<int32_t>(features.positions.size());
	out << std::setw(7) << num_positions; 

	/* Write number of sift features. */
	int32_t num_sifts = static_cast<int32_t>(features.sift_descriptors.size());
	out << std::setw(7) << num_sifts; 

	/* Write number of surf features. */
	int32_t num_surfs = static_cast<int32_t>(features.surf_descriptors.size());
	out << std::setw(7) << num_surfs << std::endl; 

	/* Write positions and colors . */
	float const fwidth = static_cast<float>(features.width);
	float const fheight = static_cast<float>(features.height);
	float const fnorm = std::max(fwidth, fheight);
	double x, y;
	for (std::size_t i = 0; i < num_positions; i++)
	{
		math::Vec2f const& pos = features.positions[i];
		
		x = pos[0];
		y = sh - 1 - pos[1];
		out << std::setw(7) << x << "  " << y;

		math::Vec3uc const& col = features.colors[i];
		out << std::setw(5) << (int)col[0] << "  " << (int)col[1] << "  " << (int)col[2] << std::endl; 
	}

	/* write feature name */
	std::string feature_name = basename(szFeaturePath);
	feature_name = feature_name + "\n";
	out << feature_name;

	out.close(); std::cout << "\rSave file: " << szFeaturePath << " Success" << std::endl;

	return EXIT_SUCCESS;
}

//inline 
//int saveHarrisFeature(std::string szFeaturePath, int Scale, int FeaW, int FeaH, fea::FeatureSet const& features)
//{
//	std::ofstream out(szFeaturePath, std::ios::binary);
//	if (!out.good()) { std::cerr << "*** Error: Can't save file: " << szFeaturePath << std::endl; return EXIT_FAILURE; }
//
//	/* write header */
//	out.write(FeatureFile_Header, FeatureFile_Header_Len);
//
//	/* Write scale and scaled image size */
//	int32_t sc = Scale;	out.write(reinterpret_cast<char const*>(&sc), sizeof(int32_t));
//	int32_t sw = FeaW;	out.write(reinterpret_cast<char const*>(&sw), sizeof(int32_t));
//	int32_t sh = FeaH;	out.write(reinterpret_cast<char const*>(&sh), sizeof(int32_t));
//
//	/* Write number of features. */
//	int32_t num_positions = static_cast<int32_t>(features.positions.size());
//	out.write(reinterpret_cast<char const*>(&num_positions), sizeof(int32_t));
//
//	/* Write number of sift features. */
//	int32_t num_sifts = static_cast<int32_t>(features.sift_descriptors.size());
//	out.write(reinterpret_cast<char const*>(&num_sifts), sizeof(int32_t));
//
//	/* Write number of surf features. */
//	int32_t num_harris = static_cast<int32_t>(features.harris_descriptors.size());
//	out.write(reinterpret_cast<char const*>(&num_harris), sizeof(int32_t));
//
//	/* Write positions and colors . */
//	for (std::size_t i = 0; i < num_positions; i++)
//	{
//		math::Vec2f const& pos = features.positions[i];
//		out.write(reinterpret_cast<char const*>(&pos), sizeof(math::Vec2f));
//
//		math::Vec3uc const& col = features.colors[i];
//		out.write(reinterpret_cast<char const*>(&col), sizeof(math::Vec3uc));
//	}
//
//	/* Write sift descriptors. */
//	for (std::size_t i = 0; i < num_sifts; i++)
//	{
//		fea::Sift::Descriptor const& des = features.sift_descriptors[i];
//		out.write(reinterpret_cast<char const*>(&des), sizeof(fea::Sift::Descriptor));
//	}
//
//	/* Write harris descriptors. */
//	for (std::size_t i = 0; i < num_harris; i++)
//	{
//		fea::Harris::Descriptor const& des = features.harris_descriptors[i];
//		out.write(reinterpret_cast<char const*>(&des), sizeof(fea::Harris::Descriptor));
//	}
//
//	/* write feature name */
//	std::string feature_name = basename(szFeaturePath);
//	int32_t len = feature_name.length();
//	out.write(reinterpret_cast<char const*>(&len), sizeof(int32_t));
//	out.write(feature_name.c_str(), len);
//
//	out.close(); std::cout << "\rSave file: " << szFeaturePath << " Success" << std::endl;
//
//	return EXIT_SUCCESS;
//}

inline 
int saveHarrisFeatureTxt(std::string szFeaturePath, int Scale, int FeaW, int FeaH, fea::FeatureSet const& features)
{
	std::ofstream out(szFeaturePath);
	if (!out.good()) { std::cerr << "*** Error: Can't save file: " << szFeaturePath << std::endl; return EXIT_FAILURE; }

	/* write header */
	out << FeatureFile_Header;

	/* Write scale and scaled image size */
	out.flags(std::ios::right);
	int32_t sc = Scale; out << std::setw(5) << sc;
	int32_t sw = FeaW; out << std::setw(5) << sw;
	int32_t sh = FeaH; out << std::setw(5) << sh << std::endl;

	/* Write number of features. */
	int32_t num_positions = static_cast<int32_t>(features.positions.size());
	out << std::setw(7) << num_positions;

	/* Write number of sift features. */
	int32_t num_sifts = static_cast<int32_t>(features.sift_descriptors.size());
	out << std::setw(7) << num_sifts;

	/* Write number of surf features. */
	int32_t num_harris = static_cast<int32_t>(features.harris_descriptors.size());
	out << std::setw(7) << num_harris << std::endl;

	/* Write positions and colors . */
	double x, y;
	for (std::size_t i = 0; i < num_harris; i++)
	{
		out << i << "  ";
		math::Vec2f const& pos = features.positions[i + num_sifts];

		x = pos[0];
		y = sh * sc - 1 - pos[1];
		out << 0.0 << "  " << 0.0 << "  " << 0.0 << "  ";
		out << std::setw(7) << x << "  " << y;

		math::Vec3uc const& col = features.colors[i + num_sifts];
		out << std::setw(5) << (int)col[0] << "  " << (int)col[1] << "  " << (int)col[2] << std::endl;
	}

	/* write feature name */
	std::string feature_name = basename(szFeaturePath);
	feature_name = feature_name + "\n";
	out << feature_name;

	out.close(); std::cout << "\rSave file: " << szFeaturePath << " Success" << std::endl;

	return EXIT_SUCCESS;
}

inline bool isFeatureFile(std::string szFeaturePath)
{
	std::ifstream in(szFeaturePath, std::ios::binary);
	if (!in.good()) return false;

	char fhead[32] = { '\0' };
	in.read(fhead, FeatureFile_Header_Len);
	if (strcmp(fhead, FeatureFile_Header) != 0) { std::cerr << "*** Error: Error fileformat" << std::endl; in.close(); return false; };

	in.close();
	return true;
}

inline
int readFeatureFile(std::string szFeaturePath, int & Scale, int & FeaW, int & FeaH, fea::FeatureSet & features)
{
	std::ifstream in(szFeaturePath, std::ios::binary);
	if (!in.good()) { std::cerr << "*** Error: Can't open file: " << szFeaturePath << std::endl; return EXIT_FAILURE; }
	std::cout << "Open file: " << szFeaturePath << " Success" << std::endl;

	char fhead[32] = { '\0' };
	in.read(fhead, FeatureFile_Header_Len);
	if (strcmp(fhead, FeatureFile_Header) != 0) { std::cerr << "*** Error: Error fileformat" << std::endl; in.close(); return EXIT_FAILURE; };

	/* Read scale and scaled image size */
	int32_t sc = 0, sw = 0, sh = 0;
	in.read(reinterpret_cast<char *>(&sc), sizeof(int32_t)); Scale = sc;
	in.read(reinterpret_cast<char *>(&sw), sizeof(int32_t)); FeaW = sw;
	in.read(reinterpret_cast<char *>(&sh), sizeof(int32_t)); FeaH = sh;

	features.width = sw;
	features.height = sh;
	features.scale = sc;

	/* Read number of features. */
	int32_t num_positions = 0;
	in.read(reinterpret_cast<char *>(&num_positions), sizeof(int32_t));
	features.positions.clear();	features.positions.resize(num_positions);
	features.colors.clear(); features.colors.resize(num_positions);

	/* Read number of sift features. */
	int32_t num_sifts = 0;
	in.read(reinterpret_cast<char *>(&num_sifts), sizeof(int32_t));
	features.sift_descriptors.clear();	features.sift_descriptors.resize(num_sifts);

	/* Read number of surf features. */
	int32_t num_surfs = 0;
	in.read(reinterpret_cast<char *>(&num_surfs), sizeof(int32_t));
	features.surf_descriptors.clear();	features.surf_descriptors.resize(num_surfs);

	/* Read number of harris features. */
	int32_t num_harris = 0;
	in.read(reinterpret_cast<char *>(&num_harris), sizeof(int32_t));
	features.harris_descriptors.clear();	features.harris_descriptors.resize(num_harris);

	/* Read positions and colors . */
	for (std::size_t i = 0; i < num_positions; i++)
	{
		math::Vec2f& pos = features.positions[i];
		in.read(reinterpret_cast<char *>(&pos), sizeof(math::Vec2f));

		math::Vec3uc& col = features.colors[i];
		in.read(reinterpret_cast<char *>(&col), sizeof(math::Vec3uc));
	}

	/* Read sift descriptors. */
	for (std::size_t i = 0; i < num_sifts; i++)
	{
		fea::Sift::Descriptor& des = features.sift_descriptors[i];
		in.read(reinterpret_cast<char *>(&des), sizeof(fea::Sift::Descriptor));
	}

	/* Read surf descriptors. */
	for (std::size_t i = 0; i < num_surfs; i++)
	{
		fea::Surf::Descriptor& des = features.surf_descriptors[i];
		in.read(reinterpret_cast<char *>(&des), sizeof(fea::Surf::Descriptor));
	}

	/* Read harris descriptors. */
	for (std::size_t i = 0; i < num_harris; i++)
	{
		fea::Harris::Descriptor& des = features.harris_descriptors[i];
		in.read(reinterpret_cast<char *>(&des), sizeof(fea::Harris::Descriptor));
	}

	/* Read feature name */
	int32_t len = 0;
	in.read(reinterpret_cast<char *>(&len), sizeof(int32_t));
	char chFeaName[128] = { '\0' };
	in.read(chFeaName, len);
	std::string thisname = basename(szFeaturePath);

	if (strcmp(chFeaName, thisname.c_str()) != 0)
	{
		std::cerr << "*** Error: File name do not match! " << std::endl; in.close(); return EXIT_FAILURE;
	}
	

	in.close();

	/*fea::FeatureSet::Options option;
	if (num_sifts > 0 && num_surfs > 0)
	{
		option.feature_types = fea::FeatureSet::FEATURE_ALL;
	}
	else
	{
		if (num_sifts > 0)
		{
			option.feature_types = fea::FeatureSet::FEATURE_SIFT;
		}
		else
		{
			option.feature_types = fea::FeatureSet::FEATURE_SURF;
		}

		features.set_options(option);
	}*/
	return EXIT_SUCCESS;
}

//inline 
//int readHarrisFeature(std::string szFeaturePath, int & Scale, int & FeaW, int & FeaH, fea::FeatureSet & features)
//{
//	std::ifstream in(szFeaturePath, std::ios::binary);
//	if (!in.good()) { std::cerr << "*** Error: Can't open file: " << szFeaturePath << std::endl; return EXIT_FAILURE; }
//	std::cout << "Open file: " << szFeaturePath << " Success" << std::endl;
//
//	char fhead[32] = { '\0' };
//	in.read(fhead, FeatureFile_Header_Len);
//	if (strcmp(fhead, FeatureFile_Header) != 0) { std::cerr << "*** Error: Error fileformat" << std::endl; in.close(); return EXIT_FAILURE; };
//
//	/* Read scale and scaled image size */
//	int32_t sc = 0, sw = 0, sh = 0;
//	in.read(reinterpret_cast<char *>(&sc), sizeof(int32_t)); Scale = sc;
//	in.read(reinterpret_cast<char *>(&sw), sizeof(int32_t)); FeaW = sw;
//	in.read(reinterpret_cast<char *>(&sh), sizeof(int32_t)); FeaH = sh;
//
//	features.width = sw;
//	features.height = sh;
//	features.scale = sc;
//	features.surf_descriptors.clear();	features.surf_descriptors.resize(0);
//	
//	/* Read number of features. */
//	int32_t num_positions = 0;
//	in.read(reinterpret_cast<char *>(&num_positions), sizeof(int32_t));
//	features.positions.clear();	features.positions.resize(num_positions);
//	features.colors.clear(); features.colors.resize(num_positions);
//
//	/* Read number of sift features. */
//	int32_t num_sifts = 0;
//	in.read(reinterpret_cast<char *>(&num_sifts), sizeof(int32_t));
//	features.sift_descriptors.clear();	features.sift_descriptors.resize(num_sifts);
//
//	/* Read number of surf features. */
//	int32_t num_harris = 0;
//	in.read(reinterpret_cast<char *>(&num_harris), sizeof(int32_t));
//	features.harris_descriptors.clear();	features.harris_descriptors.resize(num_harris);
//
//	/* Read positions and colors . */
//	for (std::size_t i = 0; i < num_positions; i++)
//	{
//		math::Vec2f& pos = features.positions[i];
//		in.read(reinterpret_cast<char *>(&pos), sizeof(math::Vec2f));
//
//		math::Vec3uc& col = features.colors[i];
//		in.read(reinterpret_cast<char *>(&col), sizeof(math::Vec3uc));
//	}
//
//	/* Read sift descriptors. */
//	for (std::size_t i = 0; i < num_sifts; i++)
//	{
//		fea::Sift::Descriptor& des = features.sift_descriptors[i];
//		in.read(reinterpret_cast<char *>(&des), sizeof(fea::Sift::Descriptor));
//	}
//
//	/* Read harris descriptors. */
//	for (std::size_t i = 0; i < num_harris; i++)
//	{
//		fea::Harris::Descriptor& des = features.harris_descriptors[i];
//		in.read(reinterpret_cast<char *>(&des), sizeof(fea::Harris::Descriptor));
//	}
//
//	/* Read feature name */
//	int32_t len = 0;
//	in.read(reinterpret_cast<char *>(&len), sizeof(int32_t));
//	char chFeaName[128] = { '\0' };
//	in.read(chFeaName, len);
//	std::string thisname = basename(szFeaturePath);
//
//	if (strcmp(chFeaName, thisname.c_str()) != 0)
//	{
//		std::cerr << "*** Error: File name do not match! " << std::endl; in.close(); return EXIT_FAILURE;
//	}
//
//
//	in.close();
//
//	/*fea::FeatureSet::Options option;
//	if (num_sifts > 0 && num_surfs > 0)
//	{
//	option.feature_types = fea::FeatureSet::FEATURE_ALL;
//	}
//	else
//	{
//	if (num_sifts > 0)
//	{
//	option.feature_types = fea::FeatureSet::FEATURE_SIFT;
//	}
//	else
//	{
//	option.feature_types = fea::FeatureSet::FEATURE_SURF;
//	}
//
//	features.set_options(option);
//	}*/
//	return EXIT_SUCCESS;
//}

inline 
int saveMatchFile(std::string szMatchFile, fea::PairwiseMatching const & matching, std::vector<std::string> const & feat_file_list)
{
	std::ofstream out(szMatchFile, std::ios::binary);
	if (!out.good()) { std::cerr << "*** Error: Can't save file: " << szMatchFile << std::endl; return EXIT_FAILURE; }

	/* write header */
	out.write(MatchFile_Header, MatchFile_Header_LEN);

	/* Write number of matching pairs. */
	int32_t num_pairs = static_cast<int32_t>(matching.size());
	out.write(reinterpret_cast<char const*>(&num_pairs), sizeof(int32_t));

	/* Write per-matching pair data. */
	for (std::size_t i = 0; i < num_pairs; ++i)
	{
		fea::TwoViewMatching const& tvr = matching[i];
		int32_t id1 = static_cast<int32_t>(tvr.view_1_id);
		int32_t id2 = static_cast<int32_t>(tvr.view_2_id);
		int32_t num_matches = static_cast<int32_t>(tvr.matches.size());
		out.write(reinterpret_cast<char const*>(&id1), sizeof(int32_t));
		out.write(reinterpret_cast<char const*>(&id2), sizeof(int32_t));
		out.write(reinterpret_cast<char const*>(&num_matches), sizeof(int32_t));

		for (std::size_t j = 0; j < num_matches; ++j)
		{
			fea::CorrespondenceIndex const& c = tvr.matches[j];
			int32_t i1 = static_cast<int32_t>(c.first);
			int32_t i2 = static_cast<int32_t>(c.second);
			out.write(reinterpret_cast<char const*>(&i1), sizeof(int32_t));
			out.write(reinterpret_cast<char const*>(&i2), sizeof(int32_t));
		}
	}

	/* Write feat file name */
	int32_t file_size = feat_file_list.size();
	out.write(reinterpret_cast<char const*>(&file_size), sizeof(int32_t));
	for (std::size_t i = 0; i < feat_file_list.size(); i++)
	{
		std::string feature_name = basename(feat_file_list[i]);
		int32_t len = feature_name.length();
		out.write(reinterpret_cast<char const*>(&len), sizeof(int32_t));
		out.write(feature_name.c_str(), len);
	}

	out.close(); std::cout << "\rSave file: " << szMatchFile << " Success.\n" << std::endl;

	return EXIT_SUCCESS;
}

inline 
int saveMatchFileTxt(std::string szMatchFile, fea::PairwiseMatching const & matching, std::vector<std::string> const & feat_file_list, fea::FeatureSets * features)
{
	std::ofstream out(szMatchFile);
	if (!out.good()) { std::cerr << "*** Error: Can't save file: " << szMatchFile << std::endl; return EXIT_FAILURE; }

	/* write header */
	out << MatchFile_Header; 

	out.flags(std::ios::left);
	/* Write number of matching pairs. */
	int32_t num_pairs = static_cast<int32_t>(matching.size());
	out << std::setw(2) << num_pairs << std::endl; 

	/* Write match pair name */
	for (std::size_t i = 0; i < num_pairs; ++i)
	{
		fea::TwoViewMatching const& tvr = matching[i];
		int32_t id1 = static_cast<int32_t>(tvr.view_1_id);
		int32_t id2 = static_cast<int32_t>(tvr.view_2_id);
		out << basename(feat_file_list[id1]) << "-----" << basename(feat_file_list[id2]) << std::endl;
	}

	/* Write per-matching pair data. */
	for (std::size_t i = 0; i < num_pairs; ++i)
	{
		fea::TwoViewMatching const& tvr = matching[i];
		int32_t id1 = static_cast<int32_t>(tvr.view_1_id);
		int32_t id2 = static_cast<int32_t>(tvr.view_2_id);
		int32_t num_matches = static_cast<int32_t>(tvr.matches.size());
		out << std::setw(6) << id1 << "  " << id2 << "  " << num_matches << std::endl;

		if (nullptr != features)
		{
			fea::FeatureSet const & feature1 = features->at(id1);
			int hei1 = feature1.height * feature1.scale;
			fea::FeatureSet const & feature2 = features->at(id2);
			int hei2 = feature2.height * feature2.scale;

			for (std::size_t j = 0; j < num_matches; ++j)
			{
				fea::CorrespondenceIndex const& c = tvr.matches[j];		
				int32_t i1 = static_cast<int32_t>(c.first);
				int32_t i2 = static_cast<int32_t>(c.second);

				out << std::setw(6) << j << "  ";
				out << std::setw(7) << feature1.positions[i1][0] * feature1.scale << "  " << (hei1 - 1) - feature1.positions[i1][1] * feature1.scale << "  "
					<< feature2.positions[i2][0] * feature2.scale << "  " << (hei2 - 1) -  feature2.positions[i2][1] * feature2.scale << "  ";
				out << std::setw(6) << i1 << "  " << i2 << "  "<< 0.0 <<std::endl;
			}
		}
		else
		{
			for (std::size_t j = 0; j < num_matches; ++j)
			{
				fea::CorrespondenceIndex const& c = tvr.matches[j];
				int32_t i1 = static_cast<int32_t>(c.first);
				int32_t i2 = static_cast<int32_t>(c.second);
				out << std::setw(6) << j << "  " << i1 << "  " << i2 << std::endl;
			}
		}
		
	}
	out.close(); std::cout << "\rSave file: " << szMatchFile << " Success" << std::endl;

	return EXIT_SUCCESS;
}

inline 
int saveHarrisMatchTxt(std::string szMatchFile, fea::PairwiseMatching const & matching, std::vector<std::string> const & feat_file_list, fea::FeatureSets * features)
{
	std::ofstream out(szMatchFile);
	if (!out.good()) { std::cerr << "*** Error: Can't save file: " << szMatchFile << std::endl; return EXIT_FAILURE; }

	/* write header */
	out << MatchFile_Header;

	out.flags(std::ios::left);
	/* Write number of matching pairs. */
	int32_t num_pairs = static_cast<int32_t>(matching.size());
	out << std::setw(2) << num_pairs << std::endl;

	/* Write match pair name */
	for (std::size_t i = 0; i < num_pairs; ++i)
	{
		fea::TwoViewMatching const& tvr = matching[i];
		int32_t id1 = static_cast<int32_t>(tvr.view_1_id);
		int32_t id2 = static_cast<int32_t>(tvr.view_2_id);
		out << basename(feat_file_list[id1]) << "-----" << basename(feat_file_list[id2]) << std::endl;
	}

	/* Write per-matching pair data. */
	for (std::size_t i = 0; i < num_pairs; ++i)
	{
		fea::TwoViewMatching const& tvr = matching[i];
		int32_t id1 = static_cast<int32_t>(tvr.view_1_id);
		int32_t id2 = static_cast<int32_t>(tvr.view_2_id);
		int32_t num_matches = static_cast<int32_t>(tvr.matches.size());
		out << std::setw(6) << id1 << "  " << id2 << "  " << num_matches << std::endl;

		if (nullptr != features)
		{
			for (std::size_t j = 0; j < num_matches; ++j)
			{
				fea::CorrespondenceIndex const& c = tvr.matches[j];
				fea::FeatureSet const & feature1 = features->at(id1);
				int sift1 = feature1.sift_descriptors.size();
				int hei1 = feature1.height * feature1.scale;
				fea::FeatureSet const & feature2 = features->at(id2);
				int sift2 = feature2.sift_descriptors.size();
				int hei2 = feature2.height * feature2.scale;
				int32_t i1 = static_cast<int32_t>(c.first);
				int32_t i2 = static_cast<int32_t>(c.second);
				out << std::setw(6) << j << "  ";
				out << std::setw(7) << feature1.positions[i1][0] << "  " << (hei1 - 1 - feature1.positions[i1][1])  << "  "
					<< feature2.positions[i2][0] << "  " << (hei2 - 1 - feature2.positions[i2][1]) << "  ";
				out << std::setw(6) << i1 << "  " << i2 << std::endl;
			}
		}
		else
		{
			for (std::size_t j = 0; j < num_matches; ++j)
			{
				fea::CorrespondenceIndex const& c = tvr.matches[j];
				int32_t i1 = static_cast<int32_t>(c.first);
				int32_t i2 = static_cast<int32_t>(c.second);
				out << std::setw(6) << j << "  " << i1 << "  " << i2 << std::endl;
			}
		}

	}
	out.close(); std::cout << "\rSave file: " << szMatchFile << " Success" << std::endl;

	return EXIT_SUCCESS;
}

inline
bool isMatchFile(std::string szMatchFile)
{
	std::ifstream in(szMatchFile, std::ios::binary);
	if (!in.good()) return false;

	char fhead[32] = { '\0' };
	in.read(fhead, MatchFile_Header_LEN);
	if (strcmp(fhead, MatchFile_Header) != 0) { std::cerr << "*** Error: Error fileformat" << std::endl; in.close(); return false; };

	in.close();

	return true;
}

inline
int readMatchFile(std::string szMatchFile, fea::PairwiseMatching & matching, std::vector<std::string> * feat_file_list)
{
	std::ifstream in(szMatchFile, std::ios::binary);
	if (!in.good()) { std::cerr << "*** Error: Can't open file: " << szMatchFile << std::endl; return EXIT_FAILURE; }
	std::cout << "Open file: " << szMatchFile << " Success" << std::endl;

	char fhead[32] = { '\0' };
	in.read(fhead, MatchFile_Header_LEN);
	if (strcmp(fhead, MatchFile_Header) != 0) { std::cerr << "*** Error: Error fileformat" << std::endl; in.close(); return EXIT_FAILURE; };

	/* Read number of matching pairs. */
	int32_t num_pairs = 0;	in.read(reinterpret_cast<char *>(&num_pairs), sizeof(int32_t));
	if (num_pairs <= 0) { std::cerr << "*** Error: No matching pairs !" << std::endl; in.close(); return EXIT_FAILURE; };

	/* Resize the space of matching */
	matching.resize(num_pairs);

	/* Read per-matching pair data. */
	for (std::size_t i = 0; i < num_pairs; ++i)
	{
		/* Read matching data. */
		int32_t id1 = 0, id2 = 0, num_matches = 0;
		in.read(reinterpret_cast<char *>(&id1), sizeof(int32_t)); matching[i].view_1_id = id1;
		in.read(reinterpret_cast<char *>(&id2), sizeof(int32_t)); matching[i].view_2_id = id2;
		in.read(reinterpret_cast<char *>(&num_matches), sizeof(int32_t)); matching[i].matches.clear(); matching[i].matches.resize(num_matches);

		for (std::size_t j = 0; j < num_matches; ++j)
		{
			int32_t i1 = 0, i2 = 0;
			in.read(reinterpret_cast<char *>(&i1), sizeof(int32_t)); matching[i].matches[j].first = i1;
			in.read(reinterpret_cast<char *>(&i2), sizeof(int32_t)); matching[i].matches[j].second = i2;
		}
	}

	if (feat_file_list != nullptr)
	{
		int32_t file_size = 0;
		in.read(reinterpret_cast<char *>(&file_size), sizeof(int32_t));
	
		int32_t len = 0;
		char chFeaName[128] = { '\0' };
		std::string strstr;
		for (std::size_t i = 0; i < file_size; i++)
		{		
			in.read(reinterpret_cast<char *>(&len), sizeof(int32_t));		
			in.read(chFeaName, len); 
			chFeaName[len] = '\0';
			strstr.assign(chFeaName);

			feat_file_list->emplace_back(chFeaName);
		}
	}
	in.close();

	return EXIT_SUCCESS;
}

FEA_MATCH_NAMESPACE_END

FEA_NAMESPACE_END

#endif // DPGRID_COMMON_HEADER


