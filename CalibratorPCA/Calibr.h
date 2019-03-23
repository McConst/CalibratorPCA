#pragma once

#include <string>
#include <vector>
#include "const.h"
#include "Eigen/Dense"

using namespace Eigen;




class Calibr
{
	//Свойства класса
	std::string InitFile;// Полный путь к файлу инициализации
	std::string SpectraPath;//Путь к файлам спектров и градуировки;
	int SpectraCount;//Количество спектров с аттестованными значениями концентрации
	

	//Методы
	int LoadMatrixLong(const std::string &FileName, MatrixXd &X);//Функция заполнения двумерного массива данными из файла типа long (спектральными интенсивностями)
	int LoadMatrixDouble(const std::string &FileName);
	int LoadVectorLong(const std::string &FileName, VectorXd &X);//Функция инициализация вектора данными из файла

public:
	//Свойства класса

	struct StructPLS//Структура, в которой хранятся результаты PLS-разложения
	{
		int N;//Количество спектров/образцов
		int M;//Количество каналов
		int A;//Количество ГК
		RowVectorXd Xmean;//Средний спектр
		double Ymean;//Средний отклик
		Matrix <double, Dynamic, Dynamic> T;//матрица счетов размерностью NxA (N-кол-во образцов, A - кол-во ГК)
		Matrix <double, Dynamic, Dynamic> P;//нагрузки для расчета матрицы X (NxM - N-кол-во образцов, M-кол-во каналов)
		VectorXd Q;//Вектор нагрузок для расчета матрицы Y размерностью A главных компонент
		Matrix <double, Dynamic, Dynamic> W;//Матрица взвешенных нагрузок MxA (M-кол-во каналов, А-кол-во ГК)
		MatrixXd E;//Матрица остатков (ошибок) для X размерностью NxM
		VectorXd F;//Вектор остатков (ошибок) для Y размерностью N
		VectorXd B;//Вектор коэффициентов регрессии размерностью A для предсказания Y
	};

	int Amax;//Максимальное количество Главных компонент
	Matrix<double, Dynamic, Dynamic> Spectra;//Спектральная матрица в формате Eigen 
	VectorXd LE;//Вектор начальных значений LE, инициализированных из файла
	StructPLS LEcoeff;//Искомые параметры разложения матрицы нормирования
	StructPLS XNormcoeff;//Искомые параметры разложения нормированной матрицы спектров

	Matrix<double, Dynamic, CRM_ElementCount> mY;//Матрица концентраций в формате Eigen

	//Методы
	Calibr(const std::string &InitFileName);//Полный путь к файлу инициализации
	~Calibr();

	void LoadSpectra(const std::string &SpectraFileName = "Spectra.dat", const std::string &YFileName="Y.txt", const std::string &LEfilename="LE.dat");
	void SetSumXtoLE();//В качестве параметра инициализации устанавливает все значения вектора LE сумма в каналах
	void SetMeanXtoLE();//В качестве параметра для инициализации устанавливается среднее значение в спектре
	void SetMaxXtoLE();//В качестве параметра для инициализации устанавливается среднее значение в спектре

	/*Калибровка по PLS. Входное значение Y0 - нецентрированный вектор концентраций одного элемента
											A - количество Главных Компонент для разложения
	*/
	void MainCalibrationPLS();//Главный метод по калибровке с нормированием на внутренний стандарт
	void DecomposePLS(const MatrixXd &X0, const VectorXd &Y0, StructPLS &Result);//Возвращаем результат по значению, чтобы не потерять при повторном вызове
	
	//Функция разложения с нормированием на предсказанный внутренний стандарт, будет использоваться для расчета коэффициентов
	void NormDecomposePLS
		(const VectorXd &B_LE, MatrixXd &T, double Ymean, const MatrixXd &X0, const VectorXd &Y0, StructPLS &sX);

	void  ScorePredictPLS(const VectorXd &B, const MatrixXd &T, double Ymean, VectorXd &Y);//PLS прогноз через счета и коэффициенты регрессии. Y-результаты прогноза
	double RMSE(VectorXd const &Y0, VectorXd const &Ycalc);// Расчет параметра градуировки. Минимум RMSE - показатель сходимости
	void SaveResultsPLS(const std::string FileName);
	void MatrixToClipboard(MatrixXd X);//Копирование int матрицы в буфер обмена
	void VectorToClipboard(VectorXd X);//Копирование int вектора в буфер обмена
};

