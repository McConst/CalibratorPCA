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
	std::string CalibrationDataPath;//Путь к файлам с расчетами калибровки
	char SpectaLoadingMethod;//Способ загрузки спектров из файла
	int SpectraCount;//Количество  градуировочных спектров с аттестованными значениями концентрации

	

	//Методы
	void LoadY(const std::string &YFileName, std::vector<std::vector<double> > &Y0);//Загрузка аттестованных концентраций из файла в двумерный вектор
	void LoadNames(std::vector<std::string> &Names, const std::string &NamesFileName = "Names.txt");//Загрузка массива имен стандартных образцов для правильного соотнесения имя/концентрация

	int LoadMatrixLong(const std::string &FileName, MatrixXd &X);//Функция заполнения двумерного массива данными из файла типа long (спектральными интенсивностями)
	int LoadMatrixDouble(const std::string &FileName);
	int LoadVectorLong(const std::string &FileName, VectorXd &X);//Функция инициализация вектора данными из файла

public:
	//Свойства класса

	struct StructPLS//Структура, в которой хранятся результаты PLS-разложения
	{
		Index N;//Количество спектров/образцов
		Index M;//Количество каналов
		Index A;//Количество ГК
		RowVectorXd Xmean;//Средний спектр
		double Ymean;//Средний отклик
		Matrix <double, Dynamic, Dynamic> T;//матрица счетов размерностью NxA (N-кол-во образцов, A - кол-во ГК)
		Matrix <double, Dynamic, Dynamic> P;//нагрузки для расчета матрицы X (MxA - N-кол-во образцов, M-кол-во каналов)
		VectorXd Q;//Вектор нагрузок для расчета матрицы Y размерностью A главных компонент
		Matrix <double, Dynamic, Dynamic> W;//Матрица взвешенных нагрузок MxA (M-кол-во каналов, А-кол-во ГК)
		MatrixXd E;//Матрица остатков (ошибок) для X размерностью NxM
		VectorXd F;//Вектор остатков (ошибок) для Y размерностью N
		VectorXd B;//Вектор коэффициентов регрессии размерностью A для предсказания Y
	};

	std::string WorkingPath;//Путь к каталогу с файлами спектров для анализа
	int Amax;//Максимальное количество Главных компонент
	Matrix<double, Dynamic, Dynamic> Spectra;//Спектральная матрица для калибровки (или для расчетов) в формате Eigen 
	VectorXd LE;//Вектор начальных значений LE, инициализированных из файла
	StructPLS LEcoeff;//Искомые параметры разложения матрицы нормирования
	StructPLS XNormcoeff;//Искомые параметры разложения нормированной матрицы спектров
	Matrix<double, Dynamic, CRM_ElementCount> mY;//Матрица концентраций в формате Eigen
	std::string LEInitialType;//Тут хранится способ инициализации LE для расчетов. Информация будет выводиться в имя файла
	char CalibrMethod[CalibrMethodStringLength];//Метод калибровки: PLS, PCR. 3 символа под метод
	double RMSEC;//Качество калибровки
	int TotalPLSRebuildIterat;//Общее количество итераций PLS
	std::string CalcParametersFile;//Имя файла с результатами вычислений (или промежуточных вычислений)
	char FinalPLS;//Флаг окончания расчетов методом PLS
	char FinalPCR;//Флаг окончания расчетов методом PCR


	//Методы
	Calibr(const std::string &InitFileName);//Полный путь к файлу инициализации
	~Calibr();

	void LoadInitDataForCalibrat(const std::string &SpectraFileName = "Spectra.dat", 
		const std::string &YFileName="Y.txt", 
		const std::string &NamesFile="Names.txt");

	void InitLE_SetSumX();//В качестве параметра инициализации устанавливает все значения вектора LE сумма в каналах
	void InitLE_SetMeanX();//В качестве параметра для инициализации устанавливается среднее значение в спектре
	void InitLE_SetMaxX();//В качестве параметра для инициализации устанавливается среднее значение в спектре
	void InitLE_NonCoherentBackScatter(const std::string LEFileName="LE.dat");//Инициализация по площади некогерентного обратного рассеяния

	/*Калибровка по PLS. Входное значение Y0 - нецентрированный вектор концентраций одного элемента
											A - количество Главных Компонент для разложения
	*/
	void MainCalibrationPLS(int elmnt=0);//Главный метод по калибровке с нормированием на внутренний стандарт
	void DecomposePLS(const MatrixXd &X0, const VectorXd &Y0, StructPLS &Result);//Возвращаем результат по значению, чтобы не потерять при повторном вызове
	
	//Функция разложения с нормированием на предсказанный внутренний стандарт, будет использоваться для расчета коэффициентов
	void NormDecomposePLS
		(const VectorXd &B_LE, MatrixXd &T, double Ymean, const MatrixXd &X0, const VectorXd &Y0, StructPLS &sX);
	
	void  ScorePredictPLS(const VectorXd &B, const MatrixXd &T, double Ymean, VectorXd &Y);//PLS прогноз через счета и коэффициенты регрессии. Y-результаты прогноза
	void PredictPLS(const MatrixXd &X, const StructPLS &Coeff, VectorXd &Ycalc);//Простой расчет откликов Ycalc для неизвестного спектра X по PLS параметрам Coeff
	VectorXd SpectraPredictPLS(const MatrixXd &X, const StructPLS &LEcoeff, const StructPLS &NormXcoeff);//Предсказание концентрации неизвестн. спектра через коэфф. PLS до и после нормирования
	double RMSE(VectorXd const &Y0, VectorXd const &Ycalc);// Расчет параметра градуировки. Минимум RMSE - показатель сходимости
	
	
	void SaveResultsPLS();//Сохранение результатов расчет методом PLS
	ProcessError LoadResultsPLS(const std::string FileName);//Чтение параметров PLS из файла в объект класса
	std::vector<std::string> LoadElvaXSpectra(const std::string &Path, MatrixXd &X);//Загрузка всех спектров каталога в матрицу X
	void LoadElvaXSpectrum(const std::string FullFileName, RowVectorXd &X);//Загрузка спектральных интенсивностей из файла ElvaX в вектор-строку


	void MatrixToClipboard(MatrixXd X);//Копирование int матрицы в буфер обмена
	void VectorToClipboard(VectorXd X);//Копирование int вектора в буфер обмена
};

