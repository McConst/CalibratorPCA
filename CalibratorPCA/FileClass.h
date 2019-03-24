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
	bool OpenForRead(const string FullFileName);//Открыть файл для чтения
	bool SaveObject(Index INT);//Сохранение в файл объекта типа Index
	bool SaveObject(int INT);//Сохранение в файл объекта типа int
	bool SaveObject(double DOUBLE);//Сохранение в файл объекта типа double
	bool SaveObject(char CHAR);//Сохранение байта в файл
	bool SaveObject(const MatrixXd &X);//Сохранение динамической матрицы double в файл
	bool SaveObject(const Calibr::StructPLS &X);//Запись структуры StructPLS в открытый файл
	bool LoadObject(Calibr::StructPLS &X);//Загрузка данных структуры StructPLS из открытого файла
	bool LoadObject(Index &N);//Чтение результата в переменную Index, заданную параметром
	bool LoadObject(int &N);//Чтение результата в переменную int, заданную параметром
	bool LoadObject(MatrixXd &X, const Index N, const Index M);//Чтение из файла значений матрицы NxM
	bool LoadObject(double &DOUBLE);//Читаем значение в переменную типа double
	bool LoadObject(char &CHAR);//Чтение байта из файла
	
	void FindFileList(const string &FullFileMask, std::vector<std::string> &Result);//Функция поиска файлов в каталоге по шаблону

	void Close() //Закрыть файл
	{
		CloseHandle(hFile);
	}
};

