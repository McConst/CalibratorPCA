#pragma once
#include <string>
#include <Windows.h>
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
	bool SaveObject(const MatrixXd &X);//Сохранение динамической матрицы double в файл

	void Close() //Закрыть файл
	{
		CloseHandle(hFile);
	}
};

