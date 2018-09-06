
#include <iostream>
#include "AerialTriMatchLib.h"
#if _DEBUG
#pragma comment(lib, "AerialTriMatchLibD.lib")
#else
#pragma comment(lib, "AerialTriMatchLib.lib")
#endif

#include "image\image_drawing.h"

struct tfw
{
	double x0, y0;
	double intervalX, intervalY;
	double rotateX, rotateY;
};

int is_regular_file(std::string const & file)
{
	struct stat st;
	stat(file.c_str(), &st);
	if (st.st_mode == 33206)
	{
		return 1;
	}
	return 0;
}

int is_directory(std::string const & file)
{
	struct stat st;
	stat(file.c_str(), &st);
	if (st.st_mode == 16895)
	{
		return 1;
	}
	return 0;
}

int read_tfw_file(std::string const & tifffile, tfw & para)
{
	if (!is_regular_file(tifffile))
	{
		return 0;
	}

	std::string tfwfile = tifffile.substr(0, tifffile.rfind('.')) + ".tfw";
	FILE* fp = fopen(tfwfile.c_str(), "r");
	if (nullptr == fp)
	{
		return 0;
	}

	fscanf(fp, "%lf", &para.intervalX);
	fscanf(fp, "%lf", &para.rotateX);
	fscanf(fp, "%lf", &para.rotateY);
	fscanf(fp, "%lf", &para.intervalY);
	fscanf(fp, "%lf", &para.x0);
	fscanf(fp, "%lf", &para.y0);

	fclose(fp);
	return 1;
}

typedef std::pair<math::Vec2i, math::Vec2i> Correspondense;
typedef std::vector<Correspondense> Correspondenses;

template<typename T>
float ncc_harris(T* const des1, T* const des2, int size)
{
	double w1 = 0.0f, w2 = 0.0f, w3 = 0.0f;
	const T* ds1 = des1;
	const T* ds2 = des2;
	for (int j = 0; j < size; j++)
	{
		w1 += (*ds1 * *ds1);
		w2 += (*ds2 * *ds2);
		w3 += (*ds1 * *ds2);
		ds1++;
		ds2++;
	}

	double ncc_val = w3 / (std::sqrt(w1) * std::sqrt(w2));
	return float(ncc_val);
}

int match_dom(dpgrid::ByteImage::Ptr const & dom1, tfw const & tfwpara1, dpgrid::ByteImage::Ptr const & dom2, tfw const & tfwpara2, Correspondenses & match_corres,
	int search_range = 1, int window_size = 7)
{
	if (window_size > 31) return 0;

	int interval = 10;
	int border_threshold = window_size / 2 + search_range;
	int wid1 = dom1->width(); int hei1 = dom1->height();
	int wid2 = dom2->width(); int hei2 = dom2->height();
	int sub_window = window_size / 2;

	int match_type = 0;
	if (match_corres.size() > 0)
	{
		match_type = 1;
	}

	if (match_type == 0) // find match for every point
	{
#pragma omp parallel for schedule(dynamic)
#if _MSC_VER
		for (int64_t j = border_threshold; j < hei1 - border_threshold; j+= interval)
#else
		for (std::size_t j = border_threshold; j < hei1 - border_threshold; j+= interval)
#endif
		{
			for (std::size_t i = border_threshold; i < wid1 - border_threshold; i+= interval)
			{
				uint8_t data1[1000];
				uint8_t data2[1000];
				float ncc[81] = { 0.0 };

				// get data 1
				uint8_t* src1 = dom1->begin();
				for (int r = -sub_window; r <= sub_window; r++)
				{
					memcpy(data1 + window_size * (r + sub_window), src1 + (j + r) * wid1 + (i - sub_window), sizeof(uint8_t) * window_size);
				}

				int x = (int)(((i * tfwpara1.intervalX + tfwpara1.x0) - tfwpara2.x0) / tfwpara2.intervalX);
				int y = (int)(((j * tfwpara1.intervalY + tfwpara1.y0) - tfwpara2.y0) / tfwpara2.intervalY);

				if (x < border_threshold || x >= (wid2 - border_threshold) || y < border_threshold || y >= (hei2 - border_threshold))continue;
				for (int jj = -search_range; jj <= search_range; jj++)
				{
					for (int ii = -search_range; ii <= search_range; ii++)
					{
						// get data 2
						for (int r = -sub_window; r <= sub_window; r++)
						{
							memcpy(data2 + window_size * (r + sub_window), dom2->begin() + (y + jj + r) * wid2 + (x + ii - sub_window), sizeof(uint8_t) * window_size);
						}

						// ncc
						ncc[(jj + search_range) * search_range + ii + search_range] = ncc_harris(data1, data2, window_size);
					}
				}

				int search_wid = search_range * 2 + 1;
				float* max_loc = std::max_element(ncc, ncc + search_wid * search_wid);
				int loc_x = (max_loc - ncc) % search_wid - search_range + x;
				int loc_y = (max_loc - ncc) / search_wid - search_range + y;

				if (*max_loc > 0.9)
				{
					match_corres.emplace_back(math::Vec2i(i, j), math::Vec2i(loc_x, loc_y));
				}
			}
		}
	}
	else // refine match
	{
#pragma omp parallel for schedule(dynamic)
#if _MSC_VER
		for(int64_t i = 0; i < match_corres.size(); i++)
#else
		for (std::size_t i = 0; i < match_corres.size(); i++)
#endif
		{
			Correspondense & cor = match_corres[i];
			if (cor.first[0] < 0 || cor.first[1] < 0 || cor.second[0] < 0 || cor.second[1] < 0)
			{
				continue;
			}

			uint8_t data1[1000];
			uint8_t data2[1000];
			float ncc[81] = { 0.0 };

			int x1 = cor.first[0]; int y1 = cor.first[1];
			int x2 = cor.second[0]; int y2 = cor.second[1];

			// get data 1
			for (int r = -sub_window; r <= sub_window; r++)
			{
				memcpy(data1 + window_size * (r + sub_window), dom1->begin() + (y1 + r) * wid1 + (x1 - sub_window), sizeof(uint8_t) * window_size);
			}

			for (int jj = -search_range; jj <= search_range; jj++)
			{
				if (jj + y2 < 0 || jj + y2 >= hei2) continue;

				for (int ii = -search_range; ii <= search_range; ii++)
				{
					if (ii + x2 < 0 || ii + x2 >= wid2) continue;

					// get data 2
					for (int r = -sub_window; r <= sub_window; r++)
					{
						memcpy(data2 + window_size * (r + sub_window), dom2->begin() + (y2 + jj + r) * wid2 + (x2 + ii - sub_window), sizeof(uint8_t) * window_size);
					}

					// ncc
					ncc[(jj + search_range) * search_range + ii + search_range] = ncc_harris(data1, data2, window_size);
				}
			}

			int search_wid = search_range * 2 + 1;
			float* max_loc = std::max_element(ncc, ncc + search_wid * search_wid);
			int loc_x = (max_loc - ncc) % search_wid - search_range + x2;
			int loc_y = (max_loc - ncc) / search_wid - search_range + y2;

			if (*max_loc > 0.95)
			{
				cor.second[0] = loc_x; cor.second[1] = loc_y;
			}
			else
			{
				cor.first[0] = -1; cor.first[1] = -1;
				cor.second[0] = -1; cor.second[1] = -1;
			}
		}
	}

	int succ_cnt = 0;
	for (int64_t i = 0; i < match_corres.size(); i++)
	{
		Correspondense const & cor = match_corres[i];
		if (cor.first[0] < 0 || cor.first[1] < 0 || cor.second[0] < 0 || cor.second[1] < 0)
		{
			continue;
		}
		succ_cnt++;
	}
	return succ_cnt;
}

int save_result_file(std::string const & dom_file1, std::string const & dom_file2, Correspondenses match_corres)
{
	int succ_cnt = 0;
	for (int64_t i = 0; i < match_corres.size(); i++)
	{
		Correspondense const & cor = match_corres[i];
		if (cor.first[0] < 0 || cor.first[1] < 0 || cor.second[0] < 0 || cor.second[1] < 0)
		{
			continue;
		}
		succ_cnt++;
	}
	//if (succ_cnt <= 0) return 0;

	std::string result = dom_file1.substr(0, dom_file1.rfind('.')) + "_" + dom_file2.substr(dom_file2.rfind('\\') + 1, dom_file2.length() - dom_file2.rfind('\\') - 5) + ".txt";

	FILE* fp = fopen(result.c_str(), "w+");
	if (fp == nullptr) return 0;

	fprintf(fp, "%s\n", dom_file1.c_str());
	fprintf(fp, "%s\n", dom_file2.c_str());
	fprintf(fp, "%d\n", succ_cnt);

	for (int64_t i = 0; i < match_corres.size(); i++)
	{
		Correspondense const & cor = match_corres[i];
		if (cor.first[0] < 0 || cor.first[1] < 0 || cor.second[0] < 0 || cor.second[1] < 0)
		{
			continue;
		}
		fprintf(fp, "%5d  %8.3lf %8.3lf %8.3lf %8.3lf \n",
			i, (double)cor.first[0], (double)cor.first[1], (double)cor.second[0], (double)cor.second[1]);
	}

	fclose(fp);

	return 1;
}

int match_dom_dense(std::string const & dom_file1, std::string const & dom_file2, int approximate_misalignment = 2)
{
	/* Read Dom file and para. */
	std::cout << "Reading DOMs...\n" << std::endl;
	tfw tfwpara1, tfwpara2;
	if (!read_tfw_file(dom_file1, tfwpara1) || !read_tfw_file(dom_file2, tfwpara2))
	{
		std::cout << "Input images are not dom file.." << std::endl;
		return 0;
	}

	dpgrid::ByteImage::Ptr dom1, dom2;
	try 
	{ 
		dom1 = dpgrid::image::load_file(dom_file1);
		dom2 = dpgrid::image::load_file(dom_file2);
	}
	catch (const std::exception&) { std::cerr << "*** Error: Can't DOM file !" << std::endl;	return EXIT_FAILURE; }
	if (dom1->channels() == 3) { dom1 = dpgrid::image::desaturate<uint8_t>(dom1, dpgrid::image::DESATURATE_LUMINANCE);}
	if (dom2->channels() == 3) { dom2 = dpgrid::image::desaturate<uint8_t>(dom2, dpgrid::image::DESATURATE_LUMINANCE); }

	/* do matching pixel by pixel */
	int downscale_times = 0;
	downscale_times = approximate_misalignment / 2;

	// downscale image
	std::cout << "Downsampling DOMs...\n" << std::endl;
	dpgrid::ByteImage::Ptr scale_image1[10];
	dpgrid::ByteImage::Ptr scale_image2[10];
	scale_image1[0] = dom1;
	scale_image2[0] = dom2;
	for (std::size_t i = 1; i <= downscale_times; i++)
	{
		scale_image1[i] = dpgrid::image::rescale_half_size<uint8_t>(scale_image1[i - 1]);
		scale_image2[i] = dpgrid::image::rescale_half_size<uint8_t>(scale_image2[i - 1]);
	}

	// match
	Correspondenses match_corres; match_corres.empty();
	std::cout << "Matching DOMs...\n" << std::endl;
	for (int i = downscale_times; i >= 0; i--)
	{
		int ratio = std::pow(2, i);
		tfw tfpara1 = tfwpara1; tfpara1.intervalX *= ratio; tfpara1.intervalY *= ratio;
		tfw tfpara2 = tfwpara2; tfpara2.intervalX *= ratio; tfpara2.intervalY *= ratio;
		int search_range = approximate_misalignment / ratio;

		for (int64_t i = 0; i < match_corres.size(); i++)
		{
			Correspondense  & cor = match_corres[i];
			if (cor.first[0] < 0 || cor.first[1] < 0 || cor.second[0] < 0 || cor.second[1] < 0)
			{
				continue;
			}
			cor.first[0] *= 2; cor.first[1] *= 2;
			cor.second[0] *= 2; cor.second[1] *= 2;
		}

		int succ_num = match_dom(scale_image1[i], tfpara1, scale_image2[i], tfpara2, match_corres, search_range);
		if (succ_num <= 0) break;		
		std::cout << "In the " << downscale_times - i + 1 << "-th times interation, find " << succ_num << " match points!\n" << std::endl;
	}

	std::cout << "Saving results...\n"<< std::endl;
	for (int64_t i = 0; i < match_corres.size(); i++)
	{
		Correspondense  & cor = match_corres[i];
		if (cor.first[0] < 0 || cor.first[1] < 0 || cor.second[0] < 0 || cor.second[1] < 0)
		{
			continue;
		}
		cor.first[1] = dom1->height() - 1 - cor.first[1];
		cor.second[1] = dom2->height() - 1 - cor.second[1];
	}

	save_result_file(dom_file1, dom_file2, match_corres);
}

int main(int argc, char* argv[])
{
	/*dpgrid::ByteImage::Ptr image;
	try
	{
		image = dpgrid::image::load_file(argv[1]);
	}
	catch (const std::exception&) { std::cerr << "*** Error: Can't DOM file !" << std::endl;	return EXIT_FAILURE; }

	uint8_t cr[3] = { 255,0,0 };
	
	dpgrid::image::draw_line<uint8_t>(*image, 0, 0, 100, 50, cr);

	std::string sav; sav.assign(argv[2]);
	dpgrid::image::save_file(image, sav);*/

	/*Correspondenses match_corres; match_corres.empty();
	save_result_file(argv[1], argv[2], match_corres);*/

	match_dom_dense(argv[1], argv[2], 2);
	return 1;
}