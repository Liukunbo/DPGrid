
// AerialTriMatch.h : PROJECT_NAME Ӧ�ó������ͷ�ļ�
//

#pragma once

#ifndef __AFXWIN_H__
	#error "�ڰ������ļ�֮ǰ������stdafx.h�������� PCH �ļ�"
#endif

#include "resource.h"		// ������


// CAerialTriMatchApp: 
// �йش����ʵ�֣������ AerialTriMatch.cpp
//

class CAerialTriMatchApp : public CWinApp
{
public:
	CAerialTriMatchApp();

// ��д
public:
	virtual BOOL InitInstance();

// ʵ��

	DECLARE_MESSAGE_MAP()
};

extern CAerialTriMatchApp theApp;