#pragma once
#include <string>
#include <Windows.h>
#include "Calibr.h"//���������� �� ���������� StructPLS ��� ���������� ���������� �������
#include <Eigen/Dense>//���������� Eigen
using namespace std;
using namespace Eigen;

class FileClass
{
	HANDLE hFile;//����� ��������� �����
	string FileName;//������ ��� �����, � ������� ���������� ������

	LPWSTR StringToW_Char(const std::string &Str);//����� �������� ������ ���� string � w_char

public:
	FileClass();
	~FileClass();

	bool OpenForSave(const string FullFileName);//������� ���� ��� ������
	bool SaveObject(int INT);//���������� � ���� ������� ���� int
	bool SaveObject(double DOUBLE);//���������� � ���� ������� ���� double
	bool SaveObject(const MatrixXd &X);//���������� ������������ ������� double � ����
	bool SaveObject(const Calibr::StructPLS X);

	void Close() //������� ����
	{
		CloseHandle(hFile);
	}
};

