#pragma once
#include <string>
#include <Windows.h>
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
	bool SaveObject(const MatrixXd &X);//���������� ������������ ������� double � ����

	void Close() //������� ����
	{
		CloseHandle(hFile);
	}
};

