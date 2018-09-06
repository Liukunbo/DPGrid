/*----------------------------------------------------------------------+
|		WuHarris											 		    |
|       Author:     DuanYanSong  2013/05/26		  						|
|            Ver 1.0												    |
|       Copyright (c) 2013,WHU RSGIS DPGrid Group All rights reserved.  |
|		ysduan@whu.edu.cn; ysduan@sohu.com                              |
+----------------------------------------------------------------------*/
#ifndef WUHARRIS_H_DUANYANSONG_2013_05_26_08_59_324983098759874325987432098543543
#define WUHARRIS_H_DUANYANSONG_2013_05_26_08_59_324983098759874325987432098543543

// 1.
//    Ix = (p[c-1])-(p[c+1]) // dX
//    Iy = (p[r-1])-(p[r+1]) // dY
// 2.
//       Ix*Ix   Ix*Iy
// M = [               ]
//       Ix*Iy   Iy*Iy
// 3.
//                       (x*x+y*y)
// Gauss = exp( (-1) * --------------- , sigm = 0.7
//                      2* sigm*sigm
// 4.
//  cim = det(M) - r(trace(M)        
//   = ((Ix*Ix)*(Iy*Iy) - (Ix*Iy)*(Ix*Iy)) - 0.04*(Ix*Ix+Iy*Iy)*(Ix*Ix+Iy*Iy);
//                   

#ifndef _SPT
#define _SPT
//!二维点数据定义（short型）
typedef struct tagSPT
{
	short x,y;
}SPT;
#else 
#pragma message("WuHarris.hpp, Warning: SPT alread define, be sure it was define as: struct tagSPT{ short x,y; }. \
               \nWuHarris.hpp, 警告:类型 SPT 已经定义过,请确保其定义为: struct tagSPT{ short x,y; }. ") 
#endif

/*
#ifdef _DEBUG
#include "stdio.h"
#include "SpBmpFile.hpp"
inline void Save2BMP( const char *lpstrPN,float *pImg,int cols,int rows)
{
    BYTE *pI = new BYTE[cols*rows];
    int i; float mx=-99999999.f,mi=99999999.f;
    for( i=0;i<cols*rows;i++ ){
        if ( mi>pImg[i] ) mi = pImg[i];
        if ( mx<pImg[i] ) mx = pImg[i];
    }
    for( i=0;i<cols*rows;i++ ){
        pI[i]= int((pImg[i]-mi)/(mx-mi)*255);
    }
    Save2Bmp( lpstrPN,pI,cols,rows,1 );
    delete pI;
}
#endif
*/

inline void ReleaseHarrisPts(SPT *pPts){ if (pPts) delete pPts;  }
inline bool Harris_ExtrFeatPt(BYTE *pImg,float *pIx2,float *pIy2,float *pIxy,float *pM,int grdSz,float *px,float *py )
{
    // 高斯模板(0.7)
    static float gaussTemplate07[9] = {
        0.04397081413862f, 0.12171198028232f, 0.04387081413862f,
        0.12171198028232f, 0.33766882231624f, 0.12171198028232f,
        0.04387081413862f, 0.12171198028232f, 0.04387081413862f
    };
    // 计算灰度梯度矩阵
    int c,r,cur;float dx,dy,*xx,mmx=0; BYTE *pL,*pC,*pN;
    for( r=1;r<grdSz-1;r++ ){
        for( c=1;c<grdSz-1;c++ ){
            cur = r*grdSz+c;
            pC = pImg+cur; pL = pC-grdSz; pN = pC+grdSz;
            dx = (*(pL-1)+*(pC-1)+*(pN-1)-*(pL+1)-*(pC+1)-*(pN+1))/9.f;
            dy = (*(pL-1)+*(pL  )+*(pL+1)-*(pN-1)-*(pN  )-*(pN+1))/9.f;
            
            pIx2[cur] = dx*dx;
            pIy2[cur] = dy*dy;
            pIxy[cur] = dx*dy;
        }        
    }
    // Gauss滤波
    float* pt = gaussTemplate07;
    memcpy( pM,pIx2,sizeof(float)*grdSz*grdSz );
    for ( r=1;r<grdSz-1;r++ ){
        for ( c=1;c<grdSz-1;c++ ){
            cur = r*grdSz+c; xx = pM+cur;
            pIx2[cur] = (*(pt+0))*(*(xx-grdSz-1)) + (*(pt+1))*(*(xx-grdSz)) + (*(pt+2))*(*(xx-grdSz+1)) +
                        (*(pt+3))*(*(xx      -1)) + (*(pt+4))*(*(xx      )) + (*(pt+5))*(*(xx+1      )) + 
                        (*(pt+6))*(*(xx+grdSz-1)) + (*(pt+7))*(*(xx+grdSz)) + (*(pt+8))*(*(xx+grdSz+1));            
        }
	}
    memcpy( pM,pIy2,sizeof(float)*grdSz*grdSz );
    for ( r=1;r<grdSz-1;r++ ){
        for ( c=1;c<grdSz-1;c++ ){
            cur = r*grdSz+c; xx = pM+cur;
            pIy2[cur] = (*(pt+0))*(*(xx-grdSz-1)) + (*(pt+1))*(*(xx-grdSz)) + (*(pt+2))*(*(xx-grdSz+1)) +
                (*(pt+3))*(*(xx      -1)) + (*(pt+4))*(*(xx      )) + (*(pt+5))*(*(xx+1      )) + 
                (*(pt+6))*(*(xx+grdSz-1)) + (*(pt+7))*(*(xx+grdSz)) + (*(pt+8))*(*(xx+grdSz+1));            
        }
	}
    memcpy( pM,pIxy,sizeof(float)*grdSz*grdSz );
    for ( r=1;r<grdSz-1;r++ ){
        for ( c=1;c<grdSz-1;c++ ){
            cur = r*grdSz+c; xx = pM+cur;
            pIxy[cur] = (*(pt+0))*(*(xx-grdSz-1)) + (*(pt+1))*(*(xx-grdSz)) + (*(pt+2))*(*(xx-grdSz+1)) +
                (*(pt+3))*(*(xx      -1)) + (*(pt+4))*(*(xx      )) + (*(pt+5))*(*(xx+1      )) + 
                (*(pt+6))*(*(xx+grdSz-1)) + (*(pt+7))*(*(xx+grdSz)) + (*(pt+8))*(*(xx+grdSz+1));            
        }
	}
    //计算兴趣值 
    memset( pM,0,sizeof(float)*grdSz*grdSz );
    for ( mmx=0,r=1;r<grdSz-1;r++ ){
        for ( c=1;c<grdSz-1;c++ ){
            cur = r*grdSz+c;
            pM[cur]= (pIx2[cur]*pIy2[cur]-pIxy[cur]*pIxy[cur]) - 0.04f*(pIx2[cur]+pIy2[cur])*(pIx2[cur]+pIy2[cur]);
            if ( mmx<pM[cur] ){ mmx=pM[cur]; *px=(float)c; *py=(float)r; }
        }
    }
    return (mmx>3);
}

inline bool Harris_ExtrFeat(unsigned char *pImg,int cols,int rows,int grdSz,int ovlp,SPT **pPts,int *ptSum,int sC=0,int eC=-1,int sR=0,int eR=-1)
{
    int r,c,br,bc,rw; float x,y; 
    int gc = cols/grdSz+1, gr = rows/grdSz+1; if(grdSz<5){ grdSz=3; ovlp=2; }    
    BYTE *pR,*pW = new BYTE[ (grdSz+ovlp)*(grdSz+ovlp) +8];
    float *pH = new float[ (grdSz+ovlp)*(grdSz+ovlp)*4 +8];
    SPT *pFts = *pPts = new SPT[gc*gr +8]; 

    if ( eC==-1 ) eC = cols; if (eR==-1) eR = rows;
    for( r=0; r<gr; r++ ){
        br = r*grdSz; if ( br+grdSz+ovlp>rows ){  br = rows-grdSz-ovlp; }
        for (c=0; c<gc; c++ ){
            bc = c*grdSz; if ( bc+grdSz+ovlp>cols ){ bc = cols-grdSz-ovlp;  }
            if ( br>sR&&bc>sC&&br<eR&&bc<eC ){
                
                memset( pH,0,sizeof(float)*grdSz*grdSz*4 );
                for ( pR=pImg+br*cols+bc,rw=0;rw<grdSz+ovlp;rw++,pR+=cols ) memcpy( pW+rw*(grdSz+ovlp),pR,(grdSz+ovlp) );                        
                if ( Harris_ExtrFeatPt(pW,pH,pH+(grdSz+ovlp)*(grdSz+ovlp),pH+(grdSz+ovlp)*(grdSz+ovlp)*2,pH+(grdSz+ovlp)*(grdSz+ovlp)*3,grdSz+ovlp,&x,&y) ){
                    pFts->x = short( bc+x );
                    pFts->y = short( br+y );
                    pFts++;                    
                }
            }
        }
    }

	delete pW; delete pH; *ptSum = pFts-(*pPts);
    if ( pFts==*pPts ){ delete *pPts; *pPts = NULL; }
    return (*ptSum>0);
}

#endif
