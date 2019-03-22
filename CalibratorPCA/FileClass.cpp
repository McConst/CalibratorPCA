#include "pch.h"
#include "FileClass.h"
#include <iostream>
#include <string>
#include <Windows.h>
#include <Eigen/Dense>//���������� Eigen

using namespace std;
using namespace Eigen;

FileClass::FileClass()
{
}


FileClass::~FileClass()
{
}

LPWSTR FileClass::StringToW_Char(const std::string &Str)
//������� ������������ ��������� �������� � ���� w_char
{
	//����������� ������ FileName �� ���� string � ���� widechar ����� WinAPI
	int bufferlen = ::MultiByteToWideChar(CP_ACP, 0, Str.c_str(), Str.size(), NULL, 0);//�������� ����� ������
	if (bufferlen == 0)
	{
		// ���-�� ������, ��������� GetLastError() and log.
		return 0;//������. ������� �� �������� ����� ��� ������. ���������� ������� ���������
	}

	// ������� �������� ������, ������� ���������� ����� ������� ����� delete ��� �������
	LPWSTR widestr = new WCHAR[bufferlen + 1];
	::MultiByteToWideChar(CP_ACP, 0, Str.c_str(), Str.size(), widestr, bufferlen);//��������� �������� ������ ������� ���� w_char

	// ��������� ���� � �����, ����� �������� ����-������
	widestr[bufferlen] = 0;
	return widestr;
}

bool FileClass::OpenForSave(const string FullFileName) 
//�������� ����� �� ������� ����� ����� ��� ����������
//���� �������� ������ �������, ������������ True
{
	FileName = FullFileName;

	LPWSTR widestrFileName = StringToW_Char(FileName);
	hFile;//��������� ���������� API
	hFile = CreateFile(widestrFileName, GENERIC_WRITE, 0, NULL, CREATE_ALWAYS,
		FILE_ATTRIBUTE_NORMAL, 0);
	delete[] widestrFileName;

	return (bool) hFile;
}

bool FileClass::SaveObject(int INT)
//���������������� ������� ������ � ���� �������� ���� int
{
	int DataSize = sizeof(INT);//������ ���� ������������� �� ������ �� �����
	DWORD nBytesSave;//������, ���������� ������� �� �����
	bool res = WriteFile(hFile, &INT, DataSize, &nBytesSave, NULL);
	return (nBytesSave == DataSize && res);
}

bool FileClass::SaveObject(const MatrixXd &X)
//���������������� ������� ������ � ���� �������� ���� int
{
	int N = X.rows();
	int M = X.cols();
	int DataSize =N * M * sizeof(double);//������ ���� ������������� �� ������ �� �����
	//��������� ������������ ������, � ������� ����������� ����� ��������� �������� �������
	double* Arr = new double[N*M];//������� ������ ��� ������ ������������� ��������
	
	for (int i = 0; i < M; ++i)
	{
		{
			for (int j = 0; j < N; ++j)
				Arr[i*N + j] = X(j, i);
		}
	}
		DWORD nBytesSave;//������, ���������� ������� �� �����
	bool res = WriteFile(hFile, Arr, DataSize, &nBytesSave, NULL);
	delete[] Arr;
	return (nBytesSave == DataSize && res);
}