#include "pch.h"
#include "FileClass.h"
#include "MyFuncs.h"
#include "Calibr.h"//���������� �� ���������� StructPLS ��� ���������� ���������� �������
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

	hFile = CreateFile(widestrFileName, GENERIC_WRITE, 0, NULL, CREATE_ALWAYS, FILE_ATTRIBUTE_NORMAL, 0);
	delete[] widestrFileName;
	bool res;
	res=((hFile==0) || (hFile == INVALID_HANDLE_VALUE)) ?  false: true;
	
	return res;
}

bool FileClass::OpenForRead(const string FullFileName)
//��������� ���� ��� ������
{
	FileName = FullFileName;

	LPWSTR widestrFileName = StringToW_Char(FileName);

	hFile = CreateFile(widestrFileName, GENERIC_READ, 0, 0, OPEN_EXISTING, FILE_ATTRIBUTE_NORMAL, 0);
	delete[] widestrFileName;
	bool res;
	res = ((hFile == 0) || (hFile == INVALID_HANDLE_VALUE)) ? false : true;

	return res;
}


bool FileClass::SaveObject(Index INT)
//���������������� ������� ������ � ���� �������� ���� int
{
	int DataSize = sizeof(INT);//������ ���� ������������� �� ������ �� �����
	DWORD nBytesSave;//������, ���������� ������� �� �����
	bool res = WriteFile(hFile, &INT, DataSize, &nBytesSave, NULL);
	return (nBytesSave == DataSize && res);
}

bool FileClass::SaveObject(int INT)
//���������������� ������� ������ � ���� �������� ���� int
{
	int DataSize = sizeof(INT);//������ ���� ������������� �� ������ �� �����
	DWORD nBytesSave;//������, ���������� ������� �� �����
	bool res = WriteFile(hFile, &INT, DataSize, &nBytesSave, NULL);
	return (nBytesSave == DataSize && res);
}

bool FileClass::SaveObject(char CHAR)
//���������� � ���� ������� ���� double
{
	int DataSize = sizeof(CHAR);//������ ���� ������������� �� ������ �� �����
	DWORD nBytesSave;//������, ���������� ������� �� �����
	bool res = WriteFile(hFile, &CHAR, DataSize, &nBytesSave, NULL);
	return (nBytesSave == DataSize && res);
}

bool FileClass::SaveObject(char *CHAR, int charSize)
//���������� � ���� ����-������ *CHAR ������ charSize
{
	DWORD nBytesSave;//������, ���������� ������� �� �����
	bool res = WriteFile(hFile, CHAR, charSize, &nBytesSave, NULL);
	return (nBytesSave == charSize && res);
}

bool FileClass::SaveObject(double DOUBLE)
//���������� � ���� ������� ���� double
{
	int DataSize = sizeof(DOUBLE);//������ ���� ������������� �� ������ �� �����
	DWORD nBytesSave;//������, ���������� ������� �� �����
	bool res = WriteFile(hFile, &DOUBLE, DataSize, &nBytesSave, NULL);
	return (nBytesSave == DataSize && res);
}

bool FileClass::SaveObject(const MatrixXd &X)
//���������������� ������� ������ � ���� �������� ���� int
{
	Index N = X.rows();
	Index M = X.cols();
	DWORD DataSize =N * M * sizeof(double);//������ ���� ������������� �� ������ �� �����
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

bool FileClass::SaveObject(const Calibr::StructPLS &X)
/*���������� � ���� ��������� � ����������� ������� �� PLS
��������� ������ � �����
int N				4 ����		���������� ��������
int M				4 ����		���������� �������
int A				4 ����		����� ��
VectorXd Xmean		M * 8 ����	������ - ������ ������� �������� �������
double Ymean		8 ����		������� ������(������������ / LE) � �������
MatrixXd T			N*A * 8 ����	������� ������
MatrixXd P			M*A * 8 ����	������� ��������
VectorXd Q			A * 8 ����	������ �������� ��� ������� ������� Y
MatrixXd W			M*A * 8 ����	������� ���������� �������� W
MatrixXd E			N*M * 8 ����	������� ������(��������) �� ������� X
VectorXd F			N * 8 ����	������ �������� ��� Y
VectorXd B			A * 8 ����	������ ������������� ������������� B
*/
{
	bool res{ true };
	res=res && SaveObject(X.N);
	res=res && SaveObject(X.M);
	res=res && SaveObject(X.A);
	res=res && SaveObject(X.Xmean);
	res=res && SaveObject(X.Ymean);
	res=res && SaveObject(X.T);
	res=res && SaveObject(X.P);
	res=res && SaveObject(X.Q);
	res=res && SaveObject(X.W);
	res=res && SaveObject(X.E);
	res=res && SaveObject(X.F);
	res=res && SaveObject(X.B);
	return res;
}

bool FileClass::LoadObject(Index &INT)
//���������������� ������� ������ � ���� �������� ���� int
{
	int DataSize = sizeof(INT);//������ ���� ������������� �� ������ �� �����
	//int buffer;
	DWORD nBytesSave;//������, ���������� ������� �� �����
	bool res = ReadFile(hFile, &INT, DataSize, &nBytesSave, NULL);
	//INT = buffer;
	return (nBytesSave == DataSize && res);
}

bool FileClass::LoadObject(int &INT)
//���������������� ������� ������ � ���� �������� ���� int
{
	int DataSize = sizeof(INT);//������ ���� ������������� �� ������ �� �����
	//int buffer;
	DWORD nBytesSave;//������, ���������� ������� �� �����
	bool res = ReadFile(hFile, &INT, DataSize, &nBytesSave, NULL);
	//INT = buffer;
	return (nBytesSave == DataSize && res);
}

bool FileClass::LoadObject(char &CHAR)
//���������������� ������� ������ � ���� �������� ���� int
{
	int DataSize = sizeof(char);//������ ���� ������������� �� ������ �� �����
	DWORD nBytesSave;//������, ���������� ������� �� �����
	bool res = ReadFile(hFile, &CHAR, DataSize, &nBytesSave, NULL);
	return (nBytesSave == DataSize && res);
}

bool FileClass::LoadObject(char *CHAR, int charSize)
//������ �� ����� ����-������ ������ charSize ����
{
	DWORD nBytesSave;//������, ���������� ������� �� �����
	bool res = ReadFile(hFile, CHAR, charSize, &nBytesSave, NULL);
	return (nBytesSave == charSize && res);
}


bool FileClass::LoadObject(MatrixXd &X, const Index N, const Index M)
//���������������� ������� ������ ������� �� ����� �������� ���� double
{

	DWORD DataSize = N * M * sizeof(double);//������ ���� ������������� �� ������ �� �����
	X.resize(N, M);
	//��������� ������������ ������, � ������� ����������� ����� ��������� �������� �������
	double* Arr = new double[N*M];//������� ������ ��� ������ ������������� ��������
	DWORD nBytesSave;//������, ���������� ������� �� �����
	bool res = ReadFile(hFile, Arr, DataSize, &nBytesSave, NULL);
	//��������� ������� ������� �� �������
	for (int i = 0; i < M; ++i)
	{
		{
			for (int j = 0; j < N; ++j)
				X(j, i)=Arr[i*N+j];
		}
	}
		
	delete[] Arr;
	return (nBytesSave == DataSize && res);
}

bool FileClass::LoadObject(double &DOUBLE)
//���������������� ������� ������ � ���� �������� ���� int
{
	int DataSize = sizeof(DOUBLE);//������ ���� ������������� �� ������ �� �����
	DWORD nBytesSave;//������, ���������� ������� �� �����
	bool res = ReadFile(hFile, &DOUBLE, DataSize, &nBytesSave, NULL);
	return (nBytesSave == DataSize && res);
}


bool  FileClass::LoadObject(Calibr::StructPLS &X)
/*�������� ������ PLS ���������� � ���������
��������� ������ � �����
int N				4 ����		���������� ��������
int M				4 ����		���������� �������
int A				4 ����		����� ��
VectorXd Xmean		M * 8 ����	������ - ������ ������� �������� �������
double Ymean		8 ����		������� ������(������������ / LE) � �������
MatrixXd T			N*A * 8 ����	������� ������
MatrixXd P			M*A * 8 ����	������� ��������
VectorXd Q			A * 8 ����	������ �������� ��� ������� ������� Y
MatrixXd W			M * A * 8 ����	������� ���������� �������� W
MatrixXd E			N * M * 8 ����	������� ������(��������) �� ������� X
VectorXd F			N * 8 ����	������ �������� ��� Y
VectorXd B			A * 8 ����	������ ������������� ������������� B
*/
{
	bool res{ true };
	res = res && LoadObject(X.N);//������ �������� ���������� ��������
	res = res && LoadObject(X.M);//������ �������� ���������� �������
	res = res && LoadObject(X.A);//������ ���������� ������� ���������
	MatrixXd Temp;//������� ������������ ������� ��� ������ ������ �� �����
	res = res && LoadObject(Temp, 1, X.M);//������ ������ 
	X.Xmean = Temp;
	res = res && LoadObject(X.Ymean);//������ ������� ������
	res = res && LoadObject(X.T, X.N, X.A);//������ ������� ������ T
	res = res && LoadObject(X.P, X.M, X.A);//������ ������� �������� P
	res = res && LoadObject(Temp, X.A, 1);//������ ������ ���������� �������� Q �� ��������� �������
	X.Q = Temp;
	res = res && LoadObject(X.W, X.M, X.A);//������ ������� ���������� �������� W
	res = res && LoadObject(X.E, X.N, X.M);//������ ������� ������ E
	res = res && LoadObject(Temp, X.N,1);//������ ������ �������� Y �� ��������� �������
	X.F = Temp;
	res = res && LoadObject(Temp, X.A, 1);//������ ������ ������������� B �� ��������� �������
	X.B = Temp;
	return res;
}

void FileClass::FindFileList(const string &FullFileMask, std::vector<std::string> &Result)
//������� ������ ������ � �������� �� �������
{
	Result.clear();
	LPWSTR widestrFileName = StringToW_Char(FullFileMask);
	WIN32_FIND_DATA FindFileData;//�����, �������� ����� ������������ ���������� ������ �� WinAPI
	HANDLE hFind;//����� �������� ����������������� ������
		
	hFind = FindFirstFile(widestrFileName, &FindFileData);
	
	if (hFind == INVALID_HANDLE_VALUE)
	//����� ����� �� ������
	{
		printf("������ ��� ������ ������� �����. �������� ������: (%d)\n", GetLastError());
		return;
	}

	//��������� ����������, ���� ������ � �����
	do
	{
		Result.push_back(W_charToString(FindFileData.cFileName));
	} while (FindNextFile(hFind, &FindFileData));

	FindClose(hFind);
	delete[] widestrFileName;

}