#include "harris.h"

#include "image/image_tools.h"
#include "image/image_io.h"
#include "util/timer.h"

FEA_NAMESPACE_BEGIN

/* wallis filter */
namespace wallis
{
	typedef unsigned char BYTE;
	typedef unsigned int WORD;
	typedef unsigned int UINT;
	typedef int BOOL;

#ifndef TRUE
#define TRUE 1
#endif

#ifndef FALSE
#define FALSE 0
#endif

#ifndef _BND2RGB
#define _BND2RGB
	static inline void Bnd2RGB(BYTE *pBuf, int cols, int band)
	{
		if (cols <= 0 || band <= 3) return;
#ifdef _X64
		BYTE *pRGB = pBuf; int i;
		for (i = 0; i < cols; i++, pRGB += 3, pBuf += band) { *((WORD*)pRGB) = *((WORD*)pBuf); *(pRGB + 2) = *(pBuf + 2); }
#else
		_asm {
			mov edi, pBuf
			mov esi, pBuf
			mov edx, cols
			_loop :
			dec edx
				mov eax, dword ptr[esi]
				mov byte ptr[edi], al
				shr eax, 08h
				inc edi
				mov word ptr[edi], ax
				add esi, band
				inc edi
				inc edi
				cmp edx, 0h
				jne _loop
		}
#endif
	}
#endif

	class CWuWlsImage
	{
	public:
		virtual int  GetRows() = 0;
		virtual int  GetCols() = 0;
		virtual int  GetPixelBytes() = 0;
		virtual BOOL Read(BYTE *pBuf, int rowIdx) = 0;

		virtual void SetRows(int nRows) {};
		virtual void SetCols(int nCols) {};
		virtual void SetPixelBytes(int nPxlBytes) {};
	};

	class CWuWlsFile
	{
	public:
		CWuWlsFile() 
		{
			m_pGridR0 = NULL; m_pGridR1 = NULL;
			m_gridRow = m_gridCol = 0;

			m_meanValue = 137;
			m_sigmaValue = 190;
			m_CValue = 0.8;
			m_BValue = 0.9;
			m_gridWZ = 16; 
			m_filterWZ = 32;
		};

		CWuWlsFile(double const para[6])
		{
			CWuWlsFile();

			m_meanValue = para[0];
			m_sigmaValue = para[1];
			m_CValue = para[2];
			m_BValue = para[3];
			m_gridWZ = para[4];
			m_filterWZ = para[5];
		}

		virtual ~CWuWlsFile() { Reset(); };
		virtual void    Set(double const para[6])
		{
			m_meanValue = para[0];
			m_sigmaValue = para[1];
			m_CValue = para[2];
			m_BValue = para[3];
			m_gridWZ = para[4];
			m_filterWZ = para[5];
		};

		virtual void    Reset() 
		{
			if (m_pGridR0) delete m_pGridR0; m_pGridR0 = NULL;
			if (m_pGridR1) delete m_pGridR1; m_pGridR1 = NULL;
			m_gridRow = m_gridCol = 0;     
		};

		BOOL CalWlsPar(CWuWlsImage *pImgFile) 
		{
			int colsS = pImgFile->GetCols(), rowsS = pImgFile->GetRows(), pxlBytes = pImgFile->GetPixelBytes();
			int lineSize = colsS*pxlBytes; int cw, rw, br, er, bc, ec, sum;

			m_gridRow = rowsS / m_gridWZ + 1;
			m_gridCol = colsS / m_gridWZ + 1;
			m_pGridR0 = new float[(m_gridRow + 1)*(m_gridCol + 1)];
			m_pGridR1 = new float[(m_gridRow + 1)*(m_gridCol + 1)];

			ProgBegin(m_gridRow); int cancel = 0; PrintMsg("Calculate WallisFilter Parameter ...");

			float meanV = m_meanValue, sigmaV = m_sigmaValue, cV = m_CValue, bV = m_BValue;
			float r0, r1, mean0, sigma0, mean, sigma;
			BYTE *imgBuf = new BYTE[lineSize*(1 + m_filterWZ) + 64];
			BYTE *pRow = new BYTE[lineSize];
			BYTE *pR, *pC, *pD, *pS;
			for (int vv, i = 0; i < m_gridRow; i++, ProgStep(cancel)) {
				br = i*m_gridWZ; er = br + m_filterWZ; if (cancel) break;
				if (er > rowsS) { er = rowsS; br = er - m_filterWZ; }
				for (int r = br; r < er; r++) {
					pImgFile->Read(pRow, r);
					if (pxlBytes > 1) {
						Bnd2RGB(pRow, colsS, pxlBytes);
						for (pD = pRow, pS = pRow, vv = 0; vv < colsS; vv++, pD++, pS += pxlBytes) { *pD = BYTE((UINT(*pS) + *(pS + 1) + *(pS + 2)) / 3); }
					}
					memcpy(imgBuf + (r - br)*colsS, pRow, colsS);
				}
				for (int j = 0; j < m_gridCol; j++) {
					bc = j*m_gridWZ; ec = bc + m_filterWZ;
					if (ec > colsS) { ec = colsS; bc = ec - m_filterWZ; }

					mean0 = sigma0 = 0; sum = m_filterWZ*m_filterWZ;
					for (rw = br, pR = imgBuf; rw < er; rw++, pR += colsS) {
						for (cw = bc, pC = pR + bc; cw < ec; cw++, pC++) {
							mean0 += *pC;
							sigma0 += *pC * *pC;
						}
					}
					mean = mean0 / sum;
					if (sigma0 / sum <= mean*mean) {
						r1 = 1.0;
						r0 = bV*meanV + (1.0f - bV - r1)*mean;
					}
					else {
						sigma = (float)sqrt(sigma0 / sum - mean*mean);
						r1 = cV*sigmaV / (cV*sigma + (1.0f - cV)*sigmaV);
						r0 = bV*meanV + (1.0f - bV - r1)*mean;
					}
					if (mean < 1) r1 = r0 = 0;
					*(m_pGridR0 + i*m_gridCol + j) = r0;
					*(m_pGridR1 + i*m_gridCol + j) = r1;
				}
			}
			delete pRow;
			delete imgBuf;
			ProgEnd();

			return TRUE;
		};
		inline void WlsFlt(BYTE *pS, int r, int cols)
		{
			float val, r0, r1;
			for (int c = 0; c < cols; c++, pS++) {
				InterplotWallisParameter(m_pGridR0, m_pGridR1, m_gridRow, m_gridCol, m_gridWZ, &r0, &r1, c - m_filterWZ / 2, r - m_filterWZ / 2);

				val = *pS * r1 + r0;
				if (val > 255) val = 255.f;
				else { if (val < 0) val = 0.f; }

				*pS = BYTE(val);
			}
		};


		float   m_meanValue, m_sigmaValue, m_CValue, m_BValue;
		float   m_rmean, m_rsigma;
		int     m_filterWZ, m_gridWZ;
		int     m_gridRow, m_gridCol;
		float   *m_pGridR0, *m_pGridR1;

	private:
		void     InterplotWallisParameter(float *pR0, float *pR1,
			int gridRow, int gridCol, int grdSz,
			float *r0, float *r1, int x, int y) 
		{
			if (x < 0) x = 0; if (y < 0) y = 0;
			float Z00, Z10, Z01, Z11, *pC;
			int   grid_r = y / grdSz, grid_c = x / grdSz;
			float dx = float(x - grid_c*grdSz) / grdSz;
			float dy = float(y - grid_r*grdSz) / grdSz;
			if (grid_r >= gridRow - 1) { grid_r = gridRow - 2; dy = 0.999f; }
			if (grid_c >= gridCol - 1) { grid_c = gridCol - 2; dx = 0.999f; }

			pC = pR0 + grid_r*gridCol + grid_c;
			Z00 = *(pC + 0);       Z10 = *(pC + 1);
			Z01 = *(pC + gridCol); Z11 = *(pC + gridCol + 1);
			if (Z00 == 0 || Z10 == 0 || Z01 == 0 || Z11 == 0)  *r0 = 0;
			else *r0 = (1 - dx)*(1 - dy)*Z00 + dx*(1 - dy)*Z10 + (1 - dx)*dy*Z01 + dx*dy*Z11;

			pC = pR1 + grid_r*gridCol + grid_c;
			Z00 = *(pC + 0);       Z10 = *(pC + 1);
			Z01 = *(pC + gridCol); Z11 = *(pC + gridCol + 1);
			if (Z00 == 0 || Z10 == 0 || Z01 == 0 || Z11 == 0)  *r1 = 0;
			else *r1 = (1 - dx)*(1 - dy)*Z00 + dx*(1 - dy)*Z10 + (1 - dx)*dy*Z01 + dx*dy*Z11;
		};

	protected:
		virtual void ProgBegin(int range) {};//if ( ::IsWindow(m_hWndRec) )::SendMessage( m_hWndRec,m_msgID,PROG_START,range );          };
		virtual void ProgStep(int& cancel) {};//if ( ::IsWindow(m_hWndRec) )::SendMessage( m_hWndRec,m_msgID,PROG_STEP ,LONG(&cancel) );  };
		virtual void ProgEnd() {};//if ( ::IsWindow(m_hWndRec) )::SendMessage( m_hWndRec,m_msgID,PROG_OVER ,0 );              };
		virtual void PrintMsg(const char* lpstrMsg) {};//if ( ::IsWindow(m_hWndRec) )::SendMessage( m_hWndRec,m_msgID,PROG_MSG  ,UINT(lpstrMsg) ); };
	};


	inline BOOL WallisFlt(BYTE *pImg, int colsS, int rowsS, double const * para = nullptr)
	{
		class CWuVZImageWls : public CWuWlsImage
		{
		public:
			CWuVZImageWls() { m_pImg = NULL; };
			virtual ~CWuVZImageWls() { m_pImg = NULL; };
			void    Attach(BYTE *pImg, int cols, int rows) { m_pImg = pImg; m_nCols = cols; m_nRows = rows; };
			BOOL    Read(BYTE *pBuf, int rowIdx) {
				if (rowIdx < 0 || rowIdx >= m_nRows) memset(pBuf, 0, m_nCols);
				else memcpy(pBuf, m_pImg + rowIdx*m_nCols, m_nCols);
				return m_nCols;
			};
			int  GetRows() { return m_nRows; };
			int  GetCols() { return m_nCols; };
			int  GetPixelBytes() { return 1; };
		protected:
			BYTE *m_pImg;
			int m_nCols, m_nRows;
		}wlsImage; wlsImage.Attach(pImg, colsS, rowsS);
		CWuWlsFile wlsFlt; 
		if (para != nullptr) wlsFlt.Set(para);
		if (!wlsFlt.CalWlsPar(&wlsImage)) return FALSE;
		for (int r = 0; r < rowsS; r++) wlsFlt.WlsFlt(pImg + r*colsS, r, colsS);
		return TRUE;
	}

}

#define EDGE_SZ     3

// 高斯模板(0.7)
static float gaussTemplate07[9] = {
	0.04397081413862f, 0.12171198028232f, 0.04387081413862f,
	0.12171198028232f, 0.33766882231624f, 0.12171198028232f,
	0.04387081413862f, 0.12171198028232f, 0.04387081413862f
};

Harris::Harris(Options const & options)
{
	this->opt_ = options;
}

void Harris::set_image(dpgrid::ByteImage::ConstPtr img)
{
	if (img->channels() != 1 && img->channels() != 3)
		throw std::invalid_argument("Gray or color image expected");

	//this->orig_ = img;
	if (img->channels() == 3)
	{
		this->orig_ = dpgrid::image::desaturate<uint8_t>
			(img, dpgrid::image::DESATURATE_LUMINANCE);
	}
	else
	{
		this->orig_ = img->duplicate();
	}
}

void Harris::process(void)
{
	util::ClockTimer timer;

	if (this->opt_.do_wallis_filte)
	{
		double para[6];

		para[0] = this->opt_.meanValue;
		para[1] = this->opt_.sigmaValue;
		para[2] = this->opt_.CValue;
		para[3] = this->opt_.BValue;
		para[4] = this->opt_.gridWZ;
		para[5] = this->opt_.filterWZ;

		if (!wallis::WallisFlt(&this->orig_->at(0), this->orig_->width(), this->orig_->height(), para))
		{
			if (this->opt_.debug_output)
			{
				std::cout << "Harris: Warning: "
					<< " Cannot do Wallis Filter and We will Skip it !" << std::endl;
			}	
		}

		if (this->opt_.debug_output)
		{
			std::cout << "Harris: Doing Wallis Filter took "
				<< timer.get_elapsed() << " ms." << std::endl;
		}
	}

	if (!Harris_ExtrFeat(&this->orig_->at(0), this->orig_->width(), this->orig_->height(), this->opt_.window_size, this->opt_.overlap_size))
	{
		throw std::domain_error("Harris Feature Extract Failure");
	}

	if (this->opt_.debug_output)
	{
		std::cout << "Harris: Get  " << descriptors_.size()
			<< " keypoints, took " << timer.get_elapsed() << " ms." << std::endl;
	}
}

bool Harris::Harris_ExtrFeatPt(unsigned char * pImg, float * pIx2, float * pIy2, float * pIxy, float * pM, int grdSz, float * px, float * py, float * xvalue, float mxv)
{
	// 计算灰度梯度矩阵
	int cc, rr, c, r, cur, bk; float dx, dy, *xx, mmx = 0;
	unsigned char *pL = nullptr, *pC = nullptr, *pN = nullptr;
	for (r = 1; r<grdSz - 1; r++) {
		for (c = 1; c<grdSz - 1; c++) {
			cur = r*grdSz + c;
			pC = pImg + cur; pL = pC - grdSz; pN = pC + grdSz;
			dx = (*(pL - 1) + *(pC - 1) + *(pN - 1) - *(pL + 1) - *(pC + 1) - *(pN + 1)) / 9.f;
			dy = (*(pL - 1) + *(pL)+*(pL + 1) - *(pN - 1) - *(pN)-*(pN + 1)) / 9.f;

			pIx2[cur] = dx*dx;
			pIy2[cur] = dy*dy;
			pIxy[cur] = dx*dy;
		}
	}
	// Gauss滤波
	float* pt = gaussTemplate07;
	memcpy(pM, pIx2, sizeof(float)*grdSz*grdSz);
	for (r = 1; r<grdSz - 1; r++) {
		for (c = 1; c<grdSz - 1; c++) {
			cur = r*grdSz + c; xx = pM + cur;
			pIx2[cur] = (*(pt + 0))*(*(xx - grdSz - 1)) + (*(pt + 1))*(*(xx - grdSz)) + (*(pt + 2))*(*(xx - grdSz + 1)) +
				(*(pt + 3))*(*(xx - 1)) + (*(pt + 4))*(*(xx)) + (*(pt + 5))*(*(xx + 1)) +
				(*(pt + 6))*(*(xx + grdSz - 1)) + (*(pt + 7))*(*(xx + grdSz)) + (*(pt + 8))*(*(xx + grdSz + 1));
		}
	}
	memcpy(pM, pIy2, sizeof(float)*grdSz*grdSz);
	for (r = 1; r<grdSz - 1; r++) {
		for (c = 1; c<grdSz - 1; c++) {
			cur = r*grdSz + c; xx = pM + cur;
			pIy2[cur] = (*(pt + 0))*(*(xx - grdSz - 1)) + (*(pt + 1))*(*(xx - grdSz)) + (*(pt + 2))*(*(xx - grdSz + 1)) +
				(*(pt + 3))*(*(xx - 1)) + (*(pt + 4))*(*(xx)) + (*(pt + 5))*(*(xx + 1)) +
				(*(pt + 6))*(*(xx + grdSz - 1)) + (*(pt + 7))*(*(xx + grdSz)) + (*(pt + 8))*(*(xx + grdSz + 1));
		}
	}
	memcpy(pM, pIxy, sizeof(float)*grdSz*grdSz);
	for (r = 1; r<grdSz - 1; r++) {
		for (c = 1; c<grdSz - 1; c++) {
			cur = r*grdSz + c; xx = pM + cur;
			pIxy[cur] = (*(pt + 0))*(*(xx - grdSz - 1)) + (*(pt + 1))*(*(xx - grdSz)) + (*(pt + 2))*(*(xx - grdSz + 1)) +
				(*(pt + 3))*(*(xx - 1)) + (*(pt + 4))*(*(xx)) + (*(pt + 5))*(*(xx + 1)) +
				(*(pt + 6))*(*(xx + grdSz - 1)) + (*(pt + 7))*(*(xx + grdSz)) + (*(pt + 8))*(*(xx + grdSz + 1));
		}
	}
	//计算兴趣值 
	memset(pM, 0, sizeof(float)*grdSz*grdSz);
	for (mmx = 0, r = 1; r<grdSz - 1; r++) {
		for (c = 1; c<grdSz - 1; c++) {
			cur = r*grdSz + c;
			pM[cur] = (pIx2[cur] * pIy2[cur] - pIxy[cur] * pIxy[cur]) - 0.04f*(pIx2[cur] + pIy2[cur])*(pIx2[cur] + pIy2[cur]);
			if (mmx<pM[cur]) {
				bk = 0; pC = pImg + cur;
				if (c>2 && r>2 && c<grdSz - 4 && r<grdSz - 4 && grdSz>9) {
					for (rr = r - 3; rr<r + 4; rr++) {
						for (cc = c - 3; cc <= c + 3; cc++) {
							if (pImg[rr*grdSz + cc] == 0) bk++;
						}
					}
					if (bk<3) bk = 0;
				}
				else {
					if (*(pC - 1) == 0) bk++; if (*(pC + 1) == 0) bk++;
					if (*(pC - grdSz) == 0) bk++; if (*(pC - grdSz - 1) == 0) bk++; if (*(pC - grdSz + 1) == 0) bk++;
					if (*(pC + grdSz) == 0) bk++; if (*(pC + grdSz - 1) == 0) bk++; if (*(pC + grdSz + 1) == 0) bk++;
				}
				if (!bk) {
					mmx = pM[cur]; *px = (float)c; *py = (float)r;
				}
			}
		}
	}

	*xvalue = mmx;
	return (mmx>mxv);
}

bool get_feature_data(unsigned char * pImg, int cols, int rows, Harris::Descriptor &des)
{
	int data_size = 11;
	int sub_window = data_size / 2;

	if (des.y - sub_window < 0 || des.y + sub_window >= rows
		|| des.x - sub_window < 0 || des.y + sub_window >= cols)
	{
		return false;
	}

	int r;
	int rr, sc;
	sc = des.x - sub_window;
	uint8_t* pDes = nullptr;
	for (r = -sub_window; r <= sub_window; r++)
	{
		rr = des.y + r;
		pDes = &(des.data((r + sub_window) * data_size));
		memcpy(pDes, pImg + rr * cols + sc, sizeof(unsigned char) * 11);
	}

	return true;
}

bool Harris::Harris_ExtrFeat(unsigned char * pImg, int cols, int rows, int grdSz0, int ovlp, int sC, int eC, int sR, int eR, float mxv)
{
	int r, c, br, bc, rw, gzc, gzr, grdSz = grdSz0; float x, y;
	int gc = cols / grdSz + 1, gr = rows / grdSz + 1;
	unsigned char *pR = nullptr, *pW = new unsigned char[(grdSz0 + 2 + ovlp)*(grdSz0 + 2 + ovlp) + 8];
	float *pH = new float[(grdSz0 + 2 + ovlp)*(grdSz0 + 2 + ovlp) * 4 + 8];
	
	descriptors_.clear();
	keypoints_.clear();

	if (eC == -1) eC = cols; if (eR == -1) eR = rows;
	for (r = 0; r<gr; r++) 
	{
		br = r*grdSz0 - 1 + EDGE_SZ; gzr = grdSz0 + 2; if (br + grdSz0 + EDGE_SZ + ovlp>rows) { gzr = rows - br; br = rows - gzr - EDGE_SZ; }
		for (c = 0; c<gc; c++)
		{
			bc = c*grdSz0 - 2 + EDGE_SZ; gzc = grdSz0 + 2; if (bc + grdSz0 + EDGE_SZ + ovlp>cols) { gzc = cols - bc; bc = cols - gzc - EDGE_SZ; }
			if (br >= sR&&bc >= sC&&br<eR&&bc<eC)
			{
				grdSz = std::min(gzr, gzc); if (grdSz<9 || grdSz - 2>grdSz0) continue;
				memset(pH, 0, sizeof(float)*grdSz*grdSz * 4);
				for (pR = pImg + br*cols + bc, rw = 0; rw<grdSz + ovlp; rw++, pR += cols)
					memcpy(pW + rw*(grdSz + ovlp), pR, (grdSz + ovlp));
				float xvalue = 0.0f;
				if (Harris_ExtrFeatPt(pW, pH, pH + (grdSz + ovlp)*(grdSz + ovlp), pH + (grdSz + ovlp)*(grdSz + ovlp) * 2, pH + (grdSz + ovlp)*(grdSz + ovlp) * 3, grdSz + ovlp, &x, &y, &xvalue, mxv))
				{
					Descriptor des;
					des.x = short(bc + x); des.y = short(br + y);
					des.scale = xvalue;
					if (get_feature_data(pImg, cols, rows, des))
					{
						descriptors_.push_back(des);
						Keypoint pt;
						pt.x = short(bc + x); pt.y = short(br + y);
						keypoints_.push_back(pt);
					}
				}
			}
		}
	}

	delete[]pW; pW = nullptr;
	delete[]pH; pH = nullptr;
	
	return descriptors_.size() > 0;
}


FEA_NAMESPACE_END