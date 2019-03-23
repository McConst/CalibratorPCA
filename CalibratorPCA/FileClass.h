#pragma once
#include <string>
#include <Windows.h>
#include "Calibr.h"//Библиотека со структурой StructPLS для сохранения параметров расчета
#include <Eigen/Dense>//Библиотека Eigen
using namespace std;
using namespace Eigen;

class FileClass
{
	HANDLE hFile;//Хэндл открытого файла
	string FileName;//Полное имя файла, с которым происходит работа

	LPWSTR StringToW_Char(const std::string &Str);//Метод перевода строки типа string в w_char

public:
	FileClass();
	~FileClass();

	bool OpenForSave(const string FullFileName);//Открыть файл для записи
	bool SaveObject(int INT);//Сохранение в файл объекта типа int
	bool SaveObject(double DOUBLE);//Сохранение в файл объекта типа double
	bool SaveObject(const MatrixXd &X);//Сохранение динамической матрицы double в файл
	bool SaveObject(const Calibr::StructPLS X);

	void Close() //Закрыть файл
	{
		CloseHandle(hFile);
	}
};

