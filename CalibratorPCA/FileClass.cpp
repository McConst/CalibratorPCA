#include "pch.h"
#include "FileClass.h"
#include "Calibr.h"//Библиотека со структурой StructPLS для сохранения параметров расчета
#include <iostream>
#include <string>
#include <Windows.h>
#include <Eigen/Dense>//Библиотека Eigen

using namespace std;
using namespace Eigen;

FileClass::FileClass()
{
}


FileClass::~FileClass()
{
}

LPWSTR FileClass::StringToW_Char(const std::string &Str)
//Функция конвертирует строковое значение к типу w_char
{
	//Преобразуем строку FileName от типа string к типу widechar через WinAPI
	int bufferlen = ::MultiByteToWideChar(CP_ACP, 0, Str.c_str(), Str.size(), NULL, 0);//Получаем длину буфера
	if (bufferlen == 0)
	{
		// Где-то ошибка, проверьте GetLastError() and log.
		return 0;//Ошибка. Система не выделила буфер под строку. Возвращаем нулевой указатель
	}

	// Создаем буферный массив, который необходимо будет удалить через delete вне функции
	LPWSTR widestr = new WCHAR[bufferlen + 1];
	::MultiByteToWideChar(CP_ACP, 0, Str.c_str(), Str.size(), widestr, bufferlen);//Заполняем буферный массив строкой типа w_char

	// Добавляем нуль в конце, чтобы получить нуль-строку
	widestr[bufferlen] = 0;
	return widestr;
}

bool FileClass::OpenForSave(const string FullFileName) 
//Открытие файла по полному имени файла для перезаписи
//Если операция прошла успешно, возвращается True
{
	FileName = FullFileName;

	LPWSTR widestrFileName = StringToW_Char(FileName);

	hFile = CreateFile(widestrFileName, GENERIC_WRITE, 0, NULL, CREATE_ALWAYS, FILE_ATTRIBUTE_NORMAL, 0);
	delete[] widestrFileName;
	bool res;
	res=((hFile==0) || (hFile == INVALID_HANDLE_VALUE)) ?  false: true;
	
	return res;
}

bool FileClass::SaveObject(int INT)
//Переопределенная функция записи в файл значения типа int
{
	int DataSize = sizeof(INT);//Размер байт запрашиваемый на чтение из файла
	DWORD nBytesSave;//Размер, полученный чтением из файла
	bool res = WriteFile(hFile, &INT, DataSize, &nBytesSave, NULL);
	return (nBytesSave == DataSize && res);
}

bool FileClass::SaveObject(double DOUBLE)
//Сохранение в файл объекта типа double
{
	int DataSize = sizeof(DOUBLE);//Размер байт запрашиваемый на чтение из файла
	DWORD nBytesSave;//Размер, полученный чтением из файла
	bool res = WriteFile(hFile, &DOUBLE, DataSize, &nBytesSave, NULL);
	return (nBytesSave == DataSize && res);
}

bool FileClass::SaveObject(const MatrixXd &X)
//Переопределенная функция записи в файл значения типа int
{
	int N = X.rows();
	int M = X.cols();
	int DataSize =N * M * sizeof(double);//Размер байт запрашиваемый на чтение из файла
	//Объявляем динамический массив, в котором поколоночно будут храниться значения матрицы
	double* Arr = new double[N*M];//Создаем массив для чтения целочисленных спектров
	
	for (int i = 0; i < M; ++i)
	{
		{
			for (int j = 0; j < N; ++j)
				Arr[i*N + j] = X(j, i);
		}
	}
		DWORD nBytesSave;//Размер, полученный чтением из файла
	bool res = WriteFile(hFile, Arr, DataSize, &nBytesSave, NULL);
	delete[] Arr;
	return (nBytesSave == DataSize && res);
}

bool FileClass::SaveObject(Calibr::StructPLS X)
/*Сохранение в файл структуры с параметрами расчета по PLS
Структура данных в файле
int N				4 байт		количество спектров
int M				4 байт		количество каналов
int A				4 байт		число ГК
VectorXd Xmean		M * 8 байт	вектор - строка средних значений спектра
double Ymean		8 байт		средний отклик(концентрация / LE) в спектре
MatrixXd T			N*A * 8 байт	матрица счетов
MatrixXd P			N*M * 8 байт	матрица нагрузки
VectorXd Q			A * 8 байт	вектор нагрузок для вектора отклика Y
VectorXd W			M*A * 8 байт	матрица взвешенных нагрузок W
MatrixXd E			N*M * 8 байт	матрица ошибок(остатков) от матрицы X
VectorXd F			N * 8 байт	вектор остатков для Y
VectorXd B			A * 8 байт	вектор регрессионных коэффициентов B
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