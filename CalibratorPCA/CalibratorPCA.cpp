// CalibratorPCA.cpp : Этот файл содержит функцию "main". Здесь начинается и заканчивается выполнение программы.
//

#include "pch.h"
#include <iostream>
#include <cstring>
#include <locale.h> /* Для русского языка */
#include "const.h"
#include "Calibr.h"
#include "Split.h"


int main(int argc, char* argv[])
{
	//В argv находится путь, который будет использоваться для калибровки. Извлекаем путь для инициализации класса Calibr
	std::string Path; //Переменная для получения пути
	std::string::size_type pos; //Позиция поиска в тексте

	setlocale(LC_ALL, "Rus");//Установка русской кодировки
	//setlocale(LC_NUMERIC, "C");//Разделителем дробной части делаем точку

	Path = argv[0];
	pos = Path.rfind("\\")+1;
	Path = Path.substr(0, pos); //Получили каталог, в котором находятся файлы для расчета калибровки
	Calibr Calc(Path.append(InitFile));//Инициализируем новый объект Calc данными из файла
	Calc.LoadSpectra();//Читаем спектральные данные и информацию для калибровки++
	Calc.SetMaxXtoLE();
	RowVectorXd Test;
	Calc.LoadElvaXSpectrum("1K0131701BK#2.evt", Test);
	std::cout << Test;
	//Методы для копирования матриц в буфер обмена и последующей проверки в Unscrambler или Excel
	//Calc.MatrixToClipboard(Calc.Spectra);
	//Calc.VectorToClipboard(Calc.LE);

	//Calc.DecomposePLS(Calc.Spectra, Calc.LE, Calc.LEcoeff);//Калибровка по LE для получения начальних значений LECoeff
	Calc.LoadResultsPLS("test.dat");

	Calc.MainCalibrationPLS();
	//Calc.SaveResultsPLS("test.dat");//Сохраняем информацию о параметрах разложения LE
	//system("Pause");
	
}

