// CalibratorPCA.cpp : Этот файл содержит функцию "main". Здесь начинается и заканчивается выполнение программы.
//

#define _CRT_SECURE_NO_WARNINGS//Отключаем предупреждения CRT в компиляторе

#include "pch.h"
#include <iostream>
#include <string>
#include <locale.h> /* Для русского языка */
#include "const.h"
#include "Calibr.h"
#include "Split.h"
#include "Eigen/Dense"

using namespace Eigen;

int main(int argc, char* argv[])
{
	//В argv находится путь, который будет использоваться для калибровки. Извлекаем путь для инициализации класса Calibr
	std::string Path; //Переменная для получения пути
	std::string::size_type pos; //Позиция поиска в тексте

	setlocale(LC_ALL, "Rus");//Установка русской кодировки
	//setlocale(LC_NUMERIC, "C");//Разделителем дробной части делаем точку

	Path = argv[0];
	pos = Path.rfind("\\") + 1;
	Path = Path.substr(0, pos); //Получили каталог, в котором находятся файлы для расчета калибровки

	Calibr Calc(Path.append(InitFile));//Инициализируем новый объект Calc данными из файла настроек
	
	Calc.LoadInitDataForCalibrat();//Подгружаем массив спектров или спектральные файлы, присваиваем аттестованные значения
	
	//Блок выбора способа инициализации
	if (Calc.LEInitialType == "BckSctrng")
	{
		//Читаем спектральные данные и аттестованных значения на них из файлов Spectra.dat и Y.txt
		Calc.InitLE_NonCoherentBackScatter();
	}
	else if (Calc.LEInitialType == "Summ")
	{
		Calc.InitLE_SetSumX();
	}
	else if (Calc.LEInitialType == "Max")
	{
		Calc.InitLE_SetMaxX();
	}
	else if (Calc.LEInitialType == "Mean")
	{
		Calc.InitLE_SetMeanX();
	}

	if (Calc.CalcParametersFile != "")
	{
		Calc.LoadResultsPLS(Calc.CalcParametersFile);//Подгружаем файл с промежуточными параметрами калибровки
		//Calc.MatrixToClipboard(Calc.LEcoeff.T);//Копируем матрицу счетов в буфер обмена для анализа в Excel
		//Calc.VectorToClipboard(Calc.LEcoeff.B);//Копируем вектор коэфф. в буфер

	}
	else
	//Так как файл с промежуточными расчетами отсутствует, выполняем расчет начальных счетов и коэффициентов для LE
	//по LE, инициализированным из файлов
	{
		if (!_stricmp(Calc.CalibrMethod, "PLS"))
		{
			Calc.DecomposePLS(Calc.Spectra, Calc.LE, Calc.LEcoeff);//Калибровка по LE для получения начальних значений LECoeff
		}
	}

		//Блок выбора способа калибровки
	if (!_stricmp(Calc.CalibrMethod, "PLS"))
	{
		if (!Calc.FinalPLS)
			//Расчеты ещё не выполнены
		{
			Calc.MainCalibrationPLS(1);
		}
		else
		{
			std::cout << "Расчеты градуировочных коэффициентов в файле\n" << Calc.CalcParametersFile << "\nзавершены\n";
		}
	}
	else if (Calc.CalibrMethod == "PCR")
	{

	}

	//system("Pause");
	
}

