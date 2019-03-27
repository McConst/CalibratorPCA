#include "pch.h"
#include "FileClass.h"
#include "MyFuncs.h"
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

bool FileClass::OpenForRead(const string FullFileName)
//Открываем файл для чтения
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
//Переопределенная функция записи в файл значения типа int
{
	int DataSize = sizeof(INT);//Размер байт запрашиваемый на чтение из файла
	DWORD nBytesSave;//Размер, полученный чтением из файла
	bool res = WriteFile(hFile, &INT, DataSize, &nBytesSave, NULL);
	return (nBytesSave == DataSize && res);
}

bool FileClass::SaveObject(int INT)
//Переопределенная функция записи в файл значения типа int
{
	int DataSize = sizeof(INT);//Размер байт запрашиваемый на чтение из файла
	DWORD nBytesSave;//Размер, полученный чтением из файла
	bool res = WriteFile(hFile, &INT, DataSize, &nBytesSave, NULL);
	return (nBytesSave == DataSize && res);
}

bool FileClass::SaveObject(char CHAR)
//Сохранение в файл объекта типа double
{
	int DataSize = sizeof(CHAR);//Размер байт запрашиваемый на чтение из файла
	DWORD nBytesSave;//Размер, полученный чтением из файла
	bool res = WriteFile(hFile, &CHAR, DataSize, &nBytesSave, NULL);
	return (nBytesSave == DataSize && res);
}

bool FileClass::SaveObject(char *CHAR, int charSize)
//Сохранение в файл нуль-строки *CHAR длиной charSize
{
	DWORD nBytesSave;//Размер, полученный чтением из файла
	bool res = WriteFile(hFile, CHAR, charSize, &nBytesSave, NULL);
	return (nBytesSave == charSize && res);
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
	Index N = X.rows();
	Index M = X.cols();
	DWORD DataSize =N * M * sizeof(double);//Размер байт запрашиваемый на чтение из файла
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

bool FileClass::SaveObject(const Calibr::StructPLS &X)
/*Сохранение в файл структуры с параметрами расчета по PLS
Структура данных в файле
int N				4 байт		количество спектров
int M				4 байт		количество каналов
int A				4 байт		число ГК
VectorXd Xmean		M * 8 байт	вектор - строка средних значений спектра
double Ymean		8 байт		средний отклик(концентрация / LE) в спектре
MatrixXd T			N*A * 8 байт	матрица счетов
MatrixXd P			M*A * 8 байт	матрица нагрузки
VectorXd Q			A * 8 байт	вектор нагрузок для вектора отклика Y
MatrixXd W			M*A * 8 байт	матрица взвешенных нагрузок W
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

bool FileClass::LoadObject(Index &INT)
//Переопределенная функция записи в файл значения типа int
{
	int DataSize = sizeof(INT);//Размер байт запрашиваемый на чтение из файла
	//int buffer;
	DWORD nBytesSave;//Размер, полученный чтением из файла
	bool res = ReadFile(hFile, &INT, DataSize, &nBytesSave, NULL);
	//INT = buffer;
	return (nBytesSave == DataSize && res);
}

bool FileClass::LoadObject(int &INT)
//Переопределенная функция записи в файл значения типа int
{
	int DataSize = sizeof(INT);//Размер байт запрашиваемый на чтение из файла
	//int buffer;
	DWORD nBytesSave;//Размер, полученный чтением из файла
	bool res = ReadFile(hFile, &INT, DataSize, &nBytesSave, NULL);
	//INT = buffer;
	return (nBytesSave == DataSize && res);
}

bool FileClass::LoadObject(char &CHAR)
//Переопределенная функция записи в файл значения типа int
{
	int DataSize = sizeof(char);//Размер байт запрашиваемый на чтение из файла
	DWORD nBytesSave;//Размер, полученный чтением из файла
	bool res = ReadFile(hFile, &CHAR, DataSize, &nBytesSave, NULL);
	return (nBytesSave == DataSize && res);
}

bool FileClass::LoadObject(char *CHAR, int charSize)
//чтение из файла нуль-строки длиной charSize байт
{
	DWORD nBytesSave;//Размер, полученный чтением из файла
	bool res = ReadFile(hFile, CHAR, charSize, &nBytesSave, NULL);
	return (nBytesSave == charSize && res);
}


bool FileClass::LoadObject(MatrixXd &X, const Index N, const Index M)
//Переопределенная функция чтения матрицы из файла значения типа double
{

	DWORD DataSize = N * M * sizeof(double);//Размер байт запрашиваемый на чтение из файла
	X.resize(N, M);
	//Объявляем динамический массив, в котором поколоночно будут храниться значения матрицы
	double* Arr = new double[N*M];//Создаем массив для чтения целочисленных спектров
	DWORD nBytesSave;//Размер, полученный чтением из файла
	bool res = ReadFile(hFile, Arr, DataSize, &nBytesSave, NULL);
	//Наполняем матрицу данными из массива
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
//Переопределенная функция записи в файл значения типа int
{
	int DataSize = sizeof(DOUBLE);//Размер байт запрашиваемый на чтение из файла
	DWORD nBytesSave;//Размер, полученный чтением из файла
	bool res = ReadFile(hFile, &DOUBLE, DataSize, &nBytesSave, NULL);
	return (nBytesSave == DataSize && res);
}


bool  FileClass::LoadObject(Calibr::StructPLS &X)
/*Загрузка данных PLS калибровки в структуру
Структура данных в файле
int N				4 байт		количество спектров
int M				4 байт		количество каналов
int A				4 байт		число ГК
VectorXd Xmean		M * 8 байт	вектор - строка средних значений спектра
double Ymean		8 байт		средний отклик(концентрация / LE) в спектре
MatrixXd T			N*A * 8 байт	матрица счетов
MatrixXd P			M*A * 8 байт	матрица нагрузки
VectorXd Q			A * 8 байт	вектор нагрузок для вектора отклика Y
MatrixXd W			M * A * 8 байт	матрица взвешенных нагрузок W
MatrixXd E			N * M * 8 байт	матрица ошибок(остатков) от матрицы X
VectorXd F			N * 8 байт	вектор остатков для Y
VectorXd B			A * 8 байт	вектор регрессионных коэффициентов B
*/
{
	bool res{ true };
	res = res && LoadObject(X.N);//Читаем значение количества спектров
	res = res && LoadObject(X.M);//Читаем значение количества каналов
	res = res && LoadObject(X.A);//Читаем количество главных компонент
	MatrixXd Temp;//Создаем динамическую матрицу для чтения данных из файла
	res = res && LoadObject(Temp, 1, X.M);//Читаем вектор 
	X.Xmean = Temp;
	res = res && LoadObject(X.Ymean);//Читаем средний сигнал
	res = res && LoadObject(X.T, X.N, X.A);//Читаем матрицу счетов T
	res = res && LoadObject(X.P, X.M, X.A);//Читаем матрицу нагрузок P
	res = res && LoadObject(Temp, X.A, 1);//Читаем вектор химических нагрузок Q во временную матрицу
	X.Q = Temp;
	res = res && LoadObject(X.W, X.M, X.A);//Читаем матрицу взвешенных нагрузок W
	res = res && LoadObject(X.E, X.N, X.M);//Читаем матрицу ошибок E
	res = res && LoadObject(Temp, X.N,1);//Читаем вектор остатков Y во временную матрицу
	X.F = Temp;
	res = res && LoadObject(Temp, X.A, 1);//Читаем вектор коэффициентов B во временную матрицу
	X.B = Temp;
	return res;
}

void FileClass::FindFileList(const string &FullFileMask, std::vector<std::string> &Result)
//Функция поиска файлов в каталоге по шаблону
{
	Result.clear();
	LPWSTR widestrFileName = StringToW_Char(FullFileMask);
	WIN32_FIND_DATA FindFileData;//Буфер, которому будут возвращаться результаты поиска из WinAPI
	HANDLE hFind;//Хэндл операции последовательного поиска
		
	hFind = FindFirstFile(widestrFileName, &FindFileData);
	
	if (hFind == INVALID_HANDLE_VALUE)
	//Поиск файла не удался
	{
		printf("Ошибка при поиске первого файла. Описание ошибки: (%d)\n", GetLastError());
		return;
	}

	//Результат существует, ищем дальше в цикле
	do
	{
		Result.push_back(W_charToString(FindFileData.cFileName));
	} while (FindNextFile(hFind, &FindFileData));

	FindClose(hFind);
	delete[] widestrFileName;

}