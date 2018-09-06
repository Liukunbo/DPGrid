/*----------------------------------------------------------------------+
|     SpWallisFlt                                                       |
|       Author:     DuanYanSong  2009/08/18                             |
|            Ver 1.0                                                    |
|       Copyright (c) 2009, Supresoft Corporation                       |
|          All rights reserved.                                         |
|       http://www.supresoft.com.cn                                     |
|     ysduan@supresoft.com.cn;ysduan@sohu.com                           |
+----------------------------------------------------------------------*/
#ifndef SPWALLISFLT_H_DUANYANSONG_2009_08_18_16_50_245678798097654
#define SPWALLISFLT_H_DUANYANSONG_2009_08_18_16_50_245678798097654

#include "stdio.h"


#ifndef _BND2RGB
#define _BND2RGB
static inline void Bnd2RGB(BYTE *pBuf, int cols, int band)
{
	if (cols <= 0 || band <= 3) return;
#ifdef _X64
	BYTE *pRGB = pBuf; int i;
	for (i = 0; i<cols; i++, pRGB += 3, pBuf += band){ *((WORD*)pRGB) = *((WORD*)pBuf); *(pRGB + 2) = *(pBuf + 2); }
#else
	_asm{
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

class CSpWlsImage
{
public:
    virtual int  GetRows()=0;
    virtual int  GetCols()=0;
    virtual int  GetPixelBytes()=0;
    virtual BOOL Read ( BYTE *pBuf,int rowIdx )=0;
    
    virtual void SetRows(int nRows){};
    virtual void SetCols(int nCols){};
    virtual void SetPixelBytes( int nPxlBytes ){};
    virtual BOOL Write( BYTE *pBuf,int rowIdx ){ return FALSE; };
};

class CSpWlsFile 
{ 
public:
    CSpWlsFile(){
        m_pGridR0 = NULL; m_pGridR1 = NULL;  
        m_gridRow = m_gridCol = 0;
        
        char iniFN[256]; GetModuleFileName( NULL,iniFN,sizeof(iniFN)); strcpy( strrchr(iniFN,'.'),".ini" );
        m_meanValue   = float( atof( GetPrivateProfileVal("Wallis","meanV" ,"137" ,iniFN) ) );
        m_sigmaValue  = float( atof( GetPrivateProfileVal("Wallis","sigmaV","190"  ,iniFN) ) );
        m_CValue      = float( atof( GetPrivateProfileVal("Wallis","C_Val" ,"0.8" ,iniFN) ) );
        m_BValue      = float( atof( GetPrivateProfileVal("Wallis","B_Val" ,"0.9",iniFN) ) );
        m_gridWZ      = atoi( GetPrivateProfileVal("Wallis","Grid"  ,"16",iniFN ) );
        m_filterWZ    = atoi( GetPrivateProfileVal("Wallis","Filter","31",iniFN ) );
    };
    virtual ~CSpWlsFile(){ Reset(); };
    virtual void    Reset(){
        if (m_pGridR0) delete m_pGridR0; m_pGridR0 = NULL;
        if (m_pGridR1) delete m_pGridR1; m_pGridR1 = NULL;
        m_gridRow = m_gridCol = 0;
        
        char iniFN[256]; GetModuleFileName( NULL,iniFN,sizeof(iniFN)); strcpy( strrchr(iniFN,'.'),".ini" );
        m_meanValue   = float( atof( GetPrivateProfileVal("Wallis","meanV" ,"137" ,iniFN) ) );
        m_sigmaValue  = float( atof( GetPrivateProfileVal("Wallis","sigmaV","190"  ,iniFN) ) );
        m_CValue      = float( atof( GetPrivateProfileVal("Wallis","C_Val" ,"0.8" ,iniFN) ) );
        m_BValue      = float( atof( GetPrivateProfileVal("Wallis","B_Val" ,"0.9",iniFN) ) );
        m_gridWZ      = atoi( GetPrivateProfileVal("Wallis","Grid"  ,"16",iniFN ) );
        m_filterWZ    = atoi( GetPrivateProfileVal("Wallis","Filter","31",iniFN ) );        
    };

    BOOL CalWlsPar( CSpWlsImage *pImgFile ){
        
        int colsS=pImgFile->GetCols(),rowsS=pImgFile->GetRows(),pxlBytes=pImgFile->GetPixelBytes();
        int lineSize = colsS*pxlBytes; int cw,rw,br,er,bc,ec,sum;
        
        m_gridRow = rowsS/m_gridWZ +1;
        m_gridCol = colsS/m_gridWZ +1;                        
        m_pGridR0 = new float[ (m_gridRow+1)*(m_gridCol+1) ];
        m_pGridR1 = new float[ (m_gridRow+1)*(m_gridCol+1) ];
        
        ProgBegin( m_gridRow );int cancel=0; PrintMsg( "Calculate WallisFilter Parameter ..." );

        float meanV=m_meanValue,sigmaV=m_sigmaValue,cV=m_CValue,bV=m_BValue;
        float r0,r1,mean0,sigma0,mean,sigma;
        BYTE *imgBuf = new BYTE[ lineSize*(1+m_filterWZ) +64];
        BYTE *pRow = new BYTE[ lineSize ];
        BYTE *pR,*pC,*pD,*pS;
        for(int vv,i=0; i<m_gridRow; i++,ProgStep(cancel) ){
            br = i*m_gridWZ; er = br+m_filterWZ; if(cancel) break;
            if ( er>rowsS ){ er = rowsS; br = er-m_filterWZ; }
            for ( int r=br;r<er;r++ ){
                pImgFile->Read( pRow,r );
                if ( pxlBytes>1 ){
                    Bnd2RGB( pRow,colsS,pxlBytes );
                    for ( pD=pRow,pS=pRow, vv=0;vv<colsS;vv++,pD++,pS+=pxlBytes){ *pD = BYTE( ( UINT(*pS)+*(pS+1)+*(pS+2) )/3 ); } 
                }
                memcpy( imgBuf+(r-br)*colsS,pRow,colsS );
            }
            for (int j=0; j<m_gridCol; j++){
                bc = j*m_gridWZ; ec = bc+m_filterWZ;
                if ( ec>colsS ){ ec = colsS; bc = ec-m_filterWZ;  }
                
                mean0=sigma0=0; sum = m_filterWZ*m_filterWZ;
                for ( rw=br,pR=imgBuf;rw<er;rw++,pR+=colsS ){
                    for( cw=bc,pC=pR+bc;cw<ec;cw++,pC++ ){
                        mean0  += *pC;
                        sigma0 += *pC * *pC;
                    }
                }
                mean = mean0/sum;
                if ( sigma0/sum<=mean*mean ){
                    r1 = 1.0;
                    r0 = bV*meanV+(1.0f-bV-r1)*mean;
                }else{
                    sigma = (float)sqrt( sigma0/sum - mean*mean  );    
                    r1 = cV*sigmaV/(cV*sigma+(1.0f-cV)*sigmaV);
                    r0 = bV*meanV+(1.0f-bV-r1)*mean;
                }
                
                *(m_pGridR0+i*m_gridCol+j) = r0;
                *(m_pGridR1+i*m_gridCol+j) = r1;                                
            }
        }
        delete pRow;         
        delete imgBuf;
        ProgEnd();
        
        return TRUE;
    };
    inline void WlsFlt(BYTE *pS,int r,int cols){
        float val,r0,r1;
        for ( int c=0;c<cols;c++,pS++ ){
            InterplotWallisParameter( m_pGridR0,m_pGridR1,m_gridRow,m_gridCol,m_gridWZ,&r0,&r1,c-m_filterWZ/2,r-m_filterWZ/2 );
            
            val = *pS * r1 + r0;
            if (val>255) val = 255.f;
            else{ if (val<0) val = 0.f; }
            
            *pS = BYTE( val );
        }
    };

    BOOL Load4File( LPCSTR lpstrPathName ){
        DWORD rw; char str[128]; memset(str,0,sizeof(str)); Reset(); 
        HANDLE hFile = ::CreateFile( lpstrPathName,GENERIC_READ,FILE_SHARE_READ,NULL,OPEN_EXISTING,0,NULL ); 
        if ( hFile == INVALID_HANDLE_VALUE ) return FALSE;   
        ::ReadFile( hFile,str,128,&rw,NULL ); 
        if ( strncmp(str,"SPF_WLS_VER:1.0",15)==0 ){ 
            ::ReadFile( hFile,str,128,&rw,NULL ); 
            sscanf( str,"%f%f%f%f%f%f%d%d%d%d",&m_meanValue,&m_sigmaValue,&m_CValue,&m_BValue,&m_rmean,&m_rsigma,
                &m_filterWZ,&m_gridWZ,&m_gridRow,&m_gridCol );
            m_pGridR0 = new float[ m_gridRow*m_gridCol +8];
            m_pGridR1 = new float[ m_gridRow*m_gridCol +8];
            ::ReadFile( hFile,m_pGridR0,m_gridRow*m_gridCol*sizeof(float),&rw,NULL );
            ::ReadFile( hFile,m_pGridR1,m_gridRow*m_gridCol*sizeof(float),&rw,NULL );
        }
        ::CloseHandle( hFile ); 
        return (m_gridRow*m_gridCol)>0;
    }
    BOOL Save2File( LPCSTR lpstrPathName ){
        DWORD rw; char str[128]; memset(str,0,sizeof(str));
        strcpy( str,"SPF_WLS_VER:1.0_WallisFilter_Parameter:mean,sigma,C,B,rmean,rsigma,nFilterWZ,nGridWZ,nCol,nRow,pfGridR0,pfGridR1" );
        HANDLE hFile = ::CreateFile( lpstrPathName,GENERIC_READ|GENERIC_WRITE,FILE_SHARE_READ,NULL,CREATE_ALWAYS,0,NULL ); 
        if ( hFile == INVALID_HANDLE_VALUE ) return FALSE;   
        ::WriteFile( hFile,str,128,&rw,NULL ); 
        memset(str,0,sizeof(str));
        sprintf( str,"%f %f %f %f %f %f %d %d %d %d",m_meanValue,m_sigmaValue,m_CValue,m_BValue,m_rmean,m_rsigma,
            m_filterWZ,m_gridWZ,m_gridRow,m_gridCol );
        ::WriteFile( hFile,str,128,&rw,NULL );
        ::WriteFile( hFile,m_pGridR0,m_gridRow*m_gridCol*sizeof(float),&rw,NULL );
        ::WriteFile( hFile,m_pGridR1,m_gridRow*m_gridCol*sizeof(float),&rw,NULL );
        ::FlushFileBuffers( hFile );
        ::CloseHandle( hFile );         
        return TRUE;
    }
    
    float   m_meanValue,m_sigmaValue,m_CValue,m_BValue;
    float   m_rmean,m_rsigma;
    int     m_filterWZ,m_gridWZ;
    int     m_gridRow,m_gridCol;
    float   *m_pGridR0,*m_pGridR1;

private:
    LPCSTR GetPrivateProfileVal(LPCSTR lpAppName,LPCSTR lpKeyName,LPCSTR lpDefault,LPCSTR lpFileName){
        static char strVal[64]; strVal[0]=0;
        ::GetPrivateProfileString( lpAppName,lpKeyName,lpDefault,strVal,64,lpFileName );
        ::WritePrivateProfileString( lpAppName,lpKeyName,strVal,lpFileName );
        return strVal;
    };
    void     InterplotWallisParameter( float *pR0,float *pR1, 
                                       int gridRow, int gridCol,int grdSz,
                                       float *r0,float *r1, int x, int y ){
        if ( x<0 ) x = 0; if ( y<0 ) y = 0;    
        float Z00,Z10,Z01,Z11,*pC;
        int   grid_r = y/grdSz,grid_c = x/grdSz;
        float dx = float( x-grid_c*grdSz )/grdSz;
        float dy = float( y-grid_r*grdSz )/grdSz;
        if (grid_r>=gridRow-1){ grid_r = gridRow-2; dy=0.999f; }
        if (grid_c>=gridCol-1){ grid_c = gridCol-2; dx=0.999f; }
        
        pC  = pR0+grid_r*gridCol+grid_c;
        Z00 = *(pC+0);         Z10 = *(pC+1);
        Z01 = *(pC+gridCol); Z11 = *(pC+gridCol+1);
        *r0 = (1-dx)*(1-dy)*Z00+dx*(1-dy)*Z10+(1-dx)*dy*Z01+dx*dy*Z11;    

        pC  = pR1+grid_r*gridCol+grid_c;
        Z00 = *(pC+0);         Z10 = *(pC+1);
        Z01 = *(pC+gridCol); Z11 = *(pC+gridCol+1);
        *r1 = (1-dx)*(1-dy)*Z00+dx*(1-dy)*Z10+(1-dx)*dy*Z01+dx*dy*Z11;    
    };

public:
    enum OUTMSG{
         PROG_MSG   =   10,
         PROG_START =   11,
         PROG_STEP  =   12,
         PROG_OVER  =   13,
    };
    void SetRevMsgWnd( HWND hWnd,UINT msgID ){   m_hWndRec=hWnd; m_msgID=msgID; };
protected:
    virtual void ProgBegin(int range)       {if ( ::IsWindow(m_hWndRec) )::SendMessage( m_hWndRec,m_msgID,PROG_START,range );          };
    virtual void ProgStep(int& cancel)      {if ( ::IsWindow(m_hWndRec) )::SendMessage( m_hWndRec,m_msgID,PROG_STEP ,LONG(&cancel) );  };
    virtual void ProgEnd()                  {if ( ::IsWindow(m_hWndRec) )::SendMessage( m_hWndRec,m_msgID,PROG_OVER ,0 );              };
    virtual void PrintMsg(LPCSTR lpstrMsg ) {if ( ::IsWindow(m_hWndRec) )::SendMessage( m_hWndRec,m_msgID,PROG_MSG  ,UINT(lpstrMsg) ); };
private:
    HWND            m_hWndRec;
    UINT            m_msgID;
};

inline BOOL WallisFlt(CSpWlsImage *pImgFileS,CSpWlsImage *pImgFileD,HWND hWnd=NULL,UINT msgID=0)
{
    int colsS=pImgFileS->GetCols(),rowsS=pImgFileS->GetRows(),pxlBytes=pImgFileS->GetPixelBytes();
    pImgFileD->SetCols( colsS ); pImgFileD->SetRows( rowsS ); pImgFileD->SetPixelBytes( 1 );
    CSpWlsFile wlsFlt; wlsFlt.SetRevMsgWnd(hWnd,msgID);
    if ( !wlsFlt.CalWlsPar(pImgFileS) ) return FALSE;
    
    if ( ::IsWindow(hWnd) ) ::SendMessage( hWnd,msgID,CSpWlsFile::PROG_START,rowsS/512 );

    int lineSize = colsS*pxlBytes,vv; 
    BYTE *pRow = new BYTE[ lineSize ],*pD,*pS;
    for ( int r=0;r<rowsS;r++ ){
        pImgFileS->Read( pRow,r );
        Bnd2RGB( pRow,colsS,pxlBytes ); 
        if ( pxlBytes>1 ){ for ( pD=pRow,pS=pRow, vv=0;vv<colsS;vv++,pD++,pS+=3){ *pD = BYTE( ( UINT(*pS)+*(pS+1)+*(pS+2) )/3 ); } }

        wlsFlt.WlsFlt( pRow,r,colsS );        
        pImgFileD->Write( pRow,r );

        if ( ::IsWindow(hWnd) && (r%512)==0 ) ::SendMessage( hWnd,msgID,CSpWlsFile::PROG_STEP ,0 ); 
    }
    delete pRow;

    if ( ::IsWindow(hWnd) ) ::SendMessage( hWnd,msgID,CSpWlsFile::PROG_OVER ,0 );
    return TRUE;
}

inline BOOL WallisFlt(BYTE *pImg,int colsS,int rowsS)
{
    class CSpVZImageWls : public CSpWlsImage
    {
    public:
        CSpVZImageWls(){m_pImg=NULL;};
        virtual ~CSpVZImageWls(){m_pImg=NULL;};    
        void    Attach( BYTE *pImg,int cols,int rows ){ m_pImg = pImg; m_nCols=cols; m_nRows=rows; };
        BOOL    Read ( BYTE *pBuf,int rowIdx ){
            if ( rowIdx<0 || rowIdx>=m_nRows ) memset( pBuf,0,m_nCols );
            else memcpy( pBuf,m_pImg+rowIdx*m_nCols,m_nCols );
            return m_nCols;
        };
        int  GetRows(){ return m_nRows; };
        int  GetCols(){ return m_nCols; };
        int  GetPixelBytes(){ return 1; };                
    protected:
        BYTE *m_pImg;
        int m_nCols,m_nRows;
    }wlsImage; wlsImage.Attach( pImg,colsS,rowsS );
    CSpWlsFile wlsFlt;  if ( !wlsFlt.CalWlsPar(&wlsImage) ) return FALSE;    
    for ( int r=0;r<rowsS;r++ ) wlsFlt.WlsFlt( pImg+r*colsS,r,colsS );                    
    return TRUE;
}

inline BOOL WallisFlt(int rowsS,int colsS,BYTE* pImg,
                      float meanV,float sigmaV,float cV,float bV,
                      int grdSz,int fltSz,HWND hWnd=NULL,UINT msgID=0)
{
    class CSpVZImageWls : public CSpWlsImage
    {
    public:
        CSpVZImageWls(){m_pImg=NULL;};
        virtual ~CSpVZImageWls(){m_pImg=NULL;};    
        void    Attach( BYTE *pImg,int cols,int rows ){ m_pImg = pImg; m_nCols=cols; m_nRows=rows; };
        BOOL    Read ( BYTE *pBuf,int rowIdx ){
            if ( rowIdx<0 || rowIdx>=m_nRows ) memset( pBuf,0,m_nCols );
            else memcpy( pBuf,m_pImg+rowIdx*m_nCols,m_nCols );
            return m_nCols;
        };
        int  GetRows(){ return m_nRows; };
        int  GetCols(){ return m_nCols; };
        int  GetPixelBytes(){ return 1; };                
    protected:
        BYTE *m_pImg;
        int m_nCols,m_nRows;
    }wlsImage; wlsImage.Attach( pImg,colsS,rowsS );

    CSpWlsFile wlsFlt; wlsFlt.SetRevMsgWnd(hWnd,msgID);
    wlsFlt.m_meanValue   = meanV;
    wlsFlt.m_sigmaValue  = sigmaV;
    wlsFlt.m_CValue      = cV;
    wlsFlt.m_BValue      = bV;
    wlsFlt.m_gridWZ      = min(colsS,rowsS)<512?8:grdSz;
    wlsFlt.m_filterWZ    = fltSz;
    if ( !wlsFlt.CalWlsPar(&wlsImage) ) return FALSE;
    
    if ( ::IsWindow(hWnd) ) ::SendMessage( hWnd,msgID,CSpWlsFile::PROG_START,rowsS/512 );
    for ( int r=0;r<rowsS;r++ ){
        wlsFlt.WlsFlt( pImg+r*colsS,r,colsS );        
        if ( ::IsWindow(hWnd) && (r%512)==0 ) ::SendMessage( hWnd,msgID,CSpWlsFile::PROG_STEP ,0 ); 
    }
    if ( ::IsWindow(hWnd) ) ::SendMessage( hWnd,msgID,CSpWlsFile::PROG_OVER ,0 );
    return TRUE;
}

#endif

