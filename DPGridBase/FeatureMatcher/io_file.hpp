#pragma once

#include "vector\tri_define.hpp"
using namespace tri;

namespace tri
{
    bool read_image_list_txt(const char* szImageListTxt, imginfo_list &imginfolist, char* imgDir = nullptr)
    {
        //! open the file
        FILE* fp = fopen(szImageListTxt, "r");
        if(nullptr == fp)
            return false;
        
        //! clear the memory
        imginfolist.clear();

        int i = 0;
        char chLine[1024];

        if (imgDir != nullptr)
        {
            for (i = 0; i < 9; i++)
            {
                fgets(chLine, 1024, fp);
                if (i == 3)
                {
                    strcpy(imgDir, chLine);
                    if (imgDir[strlen(imgDir) - 1] == '\n')
                    {
                        imgDir[strlen(imgDir) - 1] = '\0';
                    }
                }
            }
        }
        else
        {
            for (i = 0; i < 9; i++) fgets(chLine, 1024, fp);
        }

        //! read the num of image info
        int num = 0;
        fgets(chLine, 1024, fp);
        sscanf(chLine, "%d\n", &num);
        if (num <= 0)
        {
            fclose(fp);
            return false;
        }

        imginfo infoTemp;
		infoTemp.dt = 5;
        //! read other lines
        fgets(chLine, 1024, fp);
        fgets(chLine, 1024, fp);
        for (i = 0; i < num; i++)
        {
            fgets(chLine, 1024, fp);
            sscanf(chLine, "%s %d %d\n", infoTemp.name, &infoTemp.stripId, &infoTemp.id);
            imginfolist.push_back(infoTemp);
        }

        fclose(fp);

        return true;
    }

	bool write_image_list_txt(const char* szImageListTxt, imginfo_list &imginfolist, char* imgDir = nullptr)
	{
		//! open the file
		FILE* fp = fopen(szImageListTxt, "w");
		if (nullptr == fp)
			return false;

		int i = 0;
		int num = 0;
		num = imginfolist.size();

		fprintf(fp, "ImageListV1.0\n\n");
		if (imgDir == nullptr)
		{
			fprintf(fp, "ImageDir:\n%s\n\n", "C:\\Images");
		}
		else
		{
			fprintf(fp, "ImageDir:\n%s\n\n", imgDir);
		}
		fprintf(fp, "StripNumber:\n 1\n\n");
		fprintf(fp, "ImageNumber:\n %5d\n\n", num);
		fprintf(fp, "ImageName		StripID		ImageID\n");

		imginfo infoTemp;
		int ireal = 0;
		for (i = 0; i < num; i++)
		{
			infoTemp = imginfolist[i];
			if (infoTemp.dt == 0)continue;

			fprintf(fp, "%s %5d \t %4d\n", infoTemp.name, infoTemp.stripId, infoTemp.id);
			ireal++;
		}

		fseek(fp, 0L, SEEK_SET);
		fprintf(fp, "ImageListV1.0\n\n");
		if (imgDir == nullptr)
		{
			fprintf(fp, "ImageDir:\n%s\n\n", "C:\\Images");
		}
		else
		{
			fprintf(fp, "ImageDir:\n%s\n\n", imgDir);
		}
		fprintf(fp, "StripNumber:\n 1\n\n");
		fprintf(fp, "ImageNumber:\n %5d", ireal);

		fclose(fp);

		return true;
	}

    bool read_image_exts(const char* szExtsTxt, imginfo_list &imginfolist)
    {
        //! open the file
        FILE* fp = fopen(szExtsTxt, "r");
        if(nullptr == fp)
            return false;
        
        //! clear the memory
        imginfolist.clear();

        //! exts number
        int N;
        fscanf(fp, "%d\n", &N);
        if (0 >= N)
        {
            fclose(fp);
            return false;
        }

        //! read other lines
        int i = 0;
        imginfo infoTemp;
        for (i = 0; i < N; i++)
        {
            fscanf(fp, "%s cmr %lf %lf %lf %lf %lf %lf\n",
                infoTemp.name,
                &infoTemp.Xs, &infoTemp.Ys, &infoTemp.Zs, &infoTemp.fai, &infoTemp.wmg, &infoTemp.kav);

            imginfolist.push_back(infoTemp);
        }

        fclose(fp);

        return true;
    }

	bool write_image_exts(const char* szExtsTxt, imginfo_list &imginfolist)
	{
		//! open the file
		FILE* fp = fopen(szExtsTxt, "w");
		if (nullptr == fp)
			return false;

		//! exts number
		int N = imginfolist.size();
		fprintf(fp, "%d\n", N);

		//! read other lines
		int i = 0;
		imginfo infoTemp;
		int ireal = 0;
		for (i = 0; i < N; i++)
		{
			infoTemp = imginfolist[i];
			if (infoTemp.dt == 0)continue;

			fprintf(fp, "%s cmr %14.6lf %14.6lf %14.6lf %10.6lf %10.6lf %10.6lf\n",
				infoTemp.name,
				infoTemp.Xs, infoTemp.Ys, infoTemp.Zs, infoTemp.fai, infoTemp.wmg, infoTemp.kav);
			ireal++;
		}

		fseek(fp, 0L, SEEK_SET);
		fprintf(fp, "%d\n", ireal);

		fclose(fp);

		return true;
	}

    bool read_camera_cmr(const char* szCmrCmr, caminfo_list &camlist)
    {
        //! open the file
        FILE* fp = fopen(szCmrCmr, "r");
        if(nullptr == fp)
            return false;
        
        //! clear the memory
        camlist.clear();

        //! read info line
        char chLine[1024];
        while (1)
        {
            fgets(chLine, 1024, fp);
            if (chLine[0] != '$')
            {
                break;
            }
        }

        //! read cmr number
        int num;
        sscanf(chLine, "%d\n", &num);

        //! read other lines
        caminfo camTemp;

        int i = 0;
        for (i = 0; i < num; i++)
        {
            //! 读点的space coordinate
            fgets(chLine, 1024, fp);
            sscanf(chLine, "%d %lf %lf %lf %lf %lf %lf %lf %lf %lf %lf %lf %lf\n",
                &camTemp.id,
                &camTemp.x0, &camTemp.y0, &camTemp.f,
                &camTemp.format_X, &camTemp.format_Y,
                &camTemp.pixelSize,
                &camTemp.k0, &camTemp.k1, &camTemp.k2, &camTemp.k3, &camTemp.p1, &camTemp.p2);

            camlist.push_back(camTemp);
        }

        fclose(fp);

        return true;
    }

	bool write_camera_cmr(const char* szCmrCmr, caminfo_list &camlist)
	{
		//! open the file
		FILE* fp = fopen(szCmrCmr, "w");
		if (nullptr == fp)
			return false;

		//! read cmr number
		int num = camlist.size();
		fprintf(fp, "%d\n", num);

		//! read other lines
		caminfo camTemp;

		int i = 0;
		for (i = 0; i < num; i++)
		{
			//! 读点的space coordinate
			camTemp = camlist[i];

			fprintf(fp, "%-5d %8.6lf %8.6lf %8.4lf %8.4lf %8.4lf %8.6lf %10.6e %10.6e %10.6e %10.6e %10.6e %10.6e\n",
				camTemp.id,
				camTemp.x0, camTemp.y0, camTemp.f,
				camTemp.format_X, camTemp.format_Y,
				camTemp.pixelSize,
				camTemp.k0, camTemp.k1, camTemp.k2, camTemp.k3, camTemp.p1, camTemp.p2);
		}

		fclose(fp);
		return true;
	}

    bool read_tri_idx(const char* szTriIdx, triinfo &triinfo)
    {
        //! open the file
        FILE* fp = fopen(szTriIdx, "r");
        if(nullptr == fp)
            return false;
        
        //! clear the memory
        triinfo.clear();

		char chline[1024];

        int iNum = 0;	
		fgets(chline, 1024, fp);
        sscanf(chline, "Images:%d", &iNum);

        int i = 0, j = 0;
		
        tri_contour conTemp;
        for(i = 0; i < iNum; i++)
        {
            conTemp.contour.clear();

			fgets(chline, 1024, fp);
			fgets(chline, 1024, fp);
           int ret = sscanf(chline, "Name:%s	nPts:%d  ID:%d",
            conTemp.name, &conTemp.vertex_num, &conTemp.id);
		   if (ret == 2) conTemp.id = i;
            vertexd vert;
			for (j = 0; j < conTemp.vertex_num; j++)
            {
				fgets(chline, 1024, fp);
                sscanf(chline, "%lf %lf\n", &vert.x, &vert.y);
                conTemp.contour.push_back(vert);
            }

            triinfo.push_back(conTemp);
        }

        fclose(fp);

        return true;
    }

	bool write_lst_lst(const char* szLstLst, const char* szOriginImageDir, imginfo_list &imginfolist, caminfo &cam)
	{
		FILE* fp = fopen(szLstLst, "w+");
		if (nullptr == fp)
			return false;

		fprintf(fp, "%12.6lf %12.6lf %12.6lf %12.6lf %5d %5d\n",
			cam.x0, cam.y0, cam.f, cam.pixelSize, (int)(cam.format_X / cam.pixelSize), (int)(cam.format_Y / cam.pixelSize));
		fprintf(fp, "%.6e %.6e %.6e %.6e\n",
			cam.k1, cam.k2, cam.p1, cam.p2);

		int iImgNum = imginfolist.size();
		fprintf(fp, "%d\n", iImgNum);

		int flag = 2;
		if (imginfolist.at(0).Xs == NO_USEAGE)
		{
			flag = 1;
		}
		fprintf(fp, "%d\n", flag);

		int i = 0;
		int uid = 0;
		if (flag == 1)
		{
			for (i = 0; i < imginfolist.size(); i++)
			{
				if (imginfolist.at(i).dt == 0)
				{
					continue;
				}

				fprintf(fp, "%s\\%s %4d %12.6lf\n",
					szOriginImageDir, imginfolist.at(i).name, uid, cam.f);
				uid++;
			}				
		}
		else if (flag == 2)
		{
			for (i = 0; i < imginfolist.size(); i++)
			{
				if (imginfolist.at(i).dt == 0)
				{
					continue;
				}

				fprintf(fp, "%s\\%s %4d %12.6lf %12.6lf %12.6lf %12.6lf %12.6lf %12.6lf\n",
					szOriginImageDir, imginfolist.at(i).name, uid, 
					imginfolist.at(i).Xs, imginfolist.at(i).Ys, imginfolist.at(i).Zs,
					imginfolist.at(i).fai, imginfolist.at(i).wmg, imginfolist.at(i).kav);
				uid++;
			}
		}

		fseek(fp, 0L, SEEK_SET);
		fprintf(fp, "%12.6lf %12.6lf %12.6lf %12.6lf %5d %5d\n",
			cam.x0, cam.y0, cam.f, cam.pixelSize, (int)(cam.format_X / cam.pixelSize), (int)(cam.format_Y / cam.pixelSize));
		fprintf(fp, "%.6e %.6e %.6e %.6e\n",
			cam.k1, cam.k2, cam.p1, cam.p2);

		fprintf(fp, "%d\n", uid);

		fclose(fp);
		
		return true;
	}
}