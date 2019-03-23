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
	bool OpenForRead(const string FullFileName);//������� ���� ��� ������
	bool SaveObject(Index INT);//���������� � ���� ������� ���� int
	bool SaveObject(double DOUBLE);//���������� � ���� ������� ���� double
	bool SaveObject(const MatrixXd &X);//���������� ������������ ������� double � ����
	bool SaveObject(const Calibr::StructPLS &X);//������ ��������� StructPLS � �������� ����
	bool LoadObject(Calibr::StructPLS &X);//�������� ������ ��������� StructPLS �� ��������� �����
	bool LoadObject(Index &N);//������ ���������� � ���������� int, �������� ����������
	bool LoadObject(MatrixXd &X, const Index N, const Index M);//������ �� ����� �������� ������� NxM
	bool LoadObject(double &DOUBLE);//������ �������� � ���������� ���� double
	

	void Close() //������� ����
	{
		CloseHandle(hFile);
	}
};

