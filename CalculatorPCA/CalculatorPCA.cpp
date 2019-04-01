#include "pch.h"
#include <iostream>
#include <string>
#include <locale.h> /* Для русского языка */
#include "const.h"
#include "Calibr.h"
#include "Split.h"
#include "Eigen/Dense"

using namespace std;
using namespace Eigen;

int main(int argc, char* argv[])
{
	//В argv находится путь, который будет использоваться для подключения файла инициализации для класса Calibr
	
	setlocale(LC_ALL, "Rus");//Установка русской кодировки
	//setlocale(LC_NUMERIC, "C");//Разделителем дробной части делаем точку

	string Path = argv[0];
	string::size_type pos = Path.rfind("\\") + 1;
	Path = Path.substr(0, pos); //Получили каталог, в котором находится файл инициализации
	Calibr Calc (Path.append(CalculatorInitFile));//Инициализируем новый объект Calc данными из файла .ini
	Calc.LoadElvaXSpectra(Calc.WorkingPath, Calc.Spectra);//Читаем спектры из рабочего каталога в матрицу Spectra
	Calc.LoadResultsPLS(Calc.CalcParametersFile);//Подгружаем градуировочные коэффициенты PLS
	//Calc.LoadInitDataForCalibrat();//Подгружаем начальный спектр для проверки в Spectra

	//выполняем расчеты
	cout << "Концентрации:\n" << Calc.SpectraPredictPLS(Calc.Spectra, Calc.LEcoeff, Calc.XNormcoeff);
}