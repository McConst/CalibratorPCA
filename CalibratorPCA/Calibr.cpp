
#define _CRT_SECURE_NO_WARNINGS//отключаем предупреждения библиотеки CRT в компиляторе

#include "pch.h"
#include "Calibr.h"
#include "Split.h"
#include "MyFuncs.h"
#include "FileClass.h"//Класс работы с файлами
#include "const.h"//перечень констант
#include <string>
#include <iostream>
#include <stdlib.h>//библиотека для преобразования string в double
#include <fstream>//Библиотека работы с файлами
#include <vector>
#include <Eigen/Dense>//Библиотека Eigen
#include <Windows.h>//WinAPI
#include <omp.h>//Библиотека параллельного программирования OpenMP

using namespace Eigen;

Calibr::Calibr(const std::string &InitFileName)
//В конструкторе открываем файл настроек и присваиваем их значения свойствам класса
{
	std::ifstream file;
	file.open(InitFileName);//Создаем объект file для чтения и открываем по имени
	if (file.is_open())
	{
		std::string str;
		while (!file.eof())
		{
			getline(file, str);
			if (str.substr(0, 2) == "//")//Пропускаем комментарии
			{
				continue;
			}
			else if (str.length() == 0)//Исключаем из обработки пустые строки
			{
				continue;
			}
			//Инициализируем значения
			size_t pos = str.find("\t");
			std::string operat = str.substr(0,pos);
			pos = str.rfind("\t")+1;
			std::string value = str.substr(pos, str.length()-pos);
			if (operat == "Amax")
			{
				Amax = std::stoi(value);
			}
			else if (operat == "SpectraPath")
			{
				SpectraPath = CorrectPath(value);
			}
			else if (operat == "WorkingPath")
			{
				WorkingPath= CorrectPath(value);
			}
			else if (operat == "CalibrationDataPath")
			{
				CalibrationDataPath = CorrectPath(value);
			}
			else if (operat == "CalibrMethod")
			{
				strncpy(CalibrMethod, value.c_str(), CalibrMethodStringLength);
				CalibrMethod[CalibrMethodStringLength];
			}
			else if (operat == "CalcParametersFile")
			{
				CalcParametersFile = value;
			}
			else if (operat == "LEInitialType")
			{
				LEInitialType = value;
			}
		}
	}
	else
	{
		std::cout << "Файл инициализации не найден или не открылся" << std::endl;
	}
	file.close();//Закрываем файл
	FinalPLS=0;
	FinalPCR = 0;
	TotalPLSRebuildIterat = 0;//Нулевое начальное количество итераций
}



Calibr::~Calibr()
{
	//std::cout << "Деструктор Calibr"<<std::endl;
}

void Calibr::LoadInitDataForCalibrat(const std::string &SpectraFileName, const std::string &YFileName)
/*	Функция выполняет инициализацию массива спектров и концентраций по информации из файлов
SpectraFileName - файл массива спектров
YFileName - файл с аттестованными значениями концентраций, соответствующих каждому спектру из массива
LE инициализируется самостоятельными публичными методами на выбор
*/
{
	//Инициализируем временный динамический векторный массив для чтения концентраций из файлов
	std::vector<std::vector<double> > Y0; // Временный векторный массив аттестованных концентраций. Первая размерн. - СО, втор. размерн. - элемент


	std::string FileName = SpectraPath;
	FileName.append(YFileName);
	std::ifstream file;
	file.open(FileName);//Создаем объект file для чтения и открываем по имени
	if (file.is_open())
	{
		std::string str; Split Splitter;
		SpectraCount = 0;
		std::vector<double> Conc(CRM_ElementCount, 0);//Вектор Conc для инициализации двумерного массива концентрациями СО из файла
		while (!file.eof())
		{
			getline(file, str);
			if (str != "\0")
			{
				++SpectraCount;
				Splitter(str, "\t");//Делим строку по символам табуляции
				for (int i = 0; i <= Splitter.Ubound(); ++i)
				{
					Conc[i] = std::stod(Splitter[i], 0);
				}
				Y0.push_back(Conc);
			}
		}
	}
	else
	{
		std::cout << "Файл инициализации не найден или не открылся" << std::endl;
	}
	file.close();//Закрываем файл

	//Перебрасываем концентрации из векторного динамического массива в матрицу Eigen
	mY.resize(SpectraCount, CRM_ElementCount);
	for(int i=0;i<SpectraCount;++i)
		for (int j = 0; j < CRM_ElementCount; ++j)
		{
			mY(i, j) = Y0[i][j];
		}
	
	//Инициализируем спектральный массив из спектра
	FileName = SpectraPath;
	LoadMatrixLong(FileName.append(SpectraFileName), Spectra);//Инициализация матрицы Spectra данными из файла
}



int Calibr::LoadMatrixLong(const std::string &FileName, MatrixXd &X)
//Инициализация двумерной матрицы X данными из файла
{
	int* Arr = new int [SpectraCount*ChannelCount];//Создаем массив для чтения целочисленных спектров

	// Преобразуем строку типа String к w_string
	LPWSTR widestrFileName = StringToW_Char(FileName);

	HANDLE hFile;//Результат выполнения API
	hFile = CreateFile(widestrFileName, GENERIC_READ, 0, 0, OPEN_EXISTING, FILE_ATTRIBUTE_NORMAL, 0);
	delete[] widestrFileName;

	int DataSize = SpectraCount * ChannelCount * sizeof(int);//Размер байт запрашиваемый на чтение из файла
	DWORD nBytesRead;//Размер, полученный чтением из файла
	int res = ReadFile(hFile, Arr, DataSize, &nBytesRead, 0);//В файле массив Spectra расположен поколоночно!!!!
	CloseHandle(hFile);
	X.resize(SpectraCount,ChannelCount);
	for (int i=0;i<SpectraCount*ChannelCount;++i)
	{
		X(i % SpectraCount, i / SpectraCount) = Arr[i];//Инициализация массива Spectra для расчетов
	}
	//std::cout << Spectra.block<10, 10>(0, 0)<<std::endl;//Отображаем блок 10х10 элементов начиная с элемента 0,0
	delete[] Arr;
	return SpectraCount;
}

int Calibr::LoadVectorLong(const std::string &FileName, VectorXd &X)
//Инициализация двумерной матрицы X данными из файла
{
	int* Arr = new int[SpectraCount];//Создаем массив для чтения целочисленных значений


	LPWSTR widestrFileName = StringToW_Char(FileName);

	HANDLE hFile;//Результат выполнения API
	hFile = CreateFile(widestrFileName, GENERIC_READ, 0, 0, OPEN_EXISTING, FILE_ATTRIBUTE_NORMAL, 0);
	delete[] widestrFileName;

	int DataSize = SpectraCount * sizeof(int);//Размер байт запрашиваемый на чтение из файла
	DWORD nBytesRead;//Размер, полученный чтением из файла
	int res = ReadFile(hFile, Arr, DataSize, &nBytesRead, 0);//В файле массив Spectra расположен поколоночно!!!!
	CloseHandle(hFile);
	X.resize(SpectraCount);
	for (int i = 0; i < SpectraCount; ++i)
	{
		X(i) = Arr[i];//Инициализация вектора данными, загруженными из файла
	}
	delete[] Arr;
	return SpectraCount;
}


int Calibr::LoadMatrixDouble(const std::string &FileName)
//Создание и инициализация массива данными из файла
{
	double* Arr = new double[SpectraCount*ChannelCount];//Создаем массив для чтения спектров

// Преобразуем строку типа String к w_string
	LPWSTR widestrFileName = StringToW_Char(FileName);

	HANDLE hFile;//Результат выполнения API
	hFile = CreateFile(widestrFileName, GENERIC_READ, 0, 0, OPEN_EXISTING, FILE_ATTRIBUTE_NORMAL, 0);
	delete[] widestrFileName;

	int DataSize = SpectraCount * ChannelCount * sizeof(double);//Размер байт запрашиваемый на чтение из файла
	DWORD nBytesRead;//Размер, полученный чтением из файла
	int res = ReadFile(hFile, Arr, DataSize, &nBytesRead, 0);//В файле массив Spectra расположен поколоночно!!!!
	CloseHandle(hFile);
	Spectra.resize(SpectraCount, ChannelCount);
	for (int i = 0; i < SpectraCount*ChannelCount; ++i)
	{
		Spectra(i % SpectraCount, i / SpectraCount) = Arr[i];//Инициализация массива Spectra для расчетов
	}
	//std::cout << Spectra.block<10, 10>(0, 0)<<std::endl;//Отображаем блок 10х10 элементов начиная с элемента 0,0
	delete[] Arr;
	return SpectraCount;
}

void Calibr::DecomposePLS(const MatrixXd &X0, const VectorXd &Y0, StructPLS &Result)
//Разложение матриц X и Y по PLS1 NIPALS
// A-количество ГК для разложения
//В переменной Result типа StructPLS возвращаем результаты разложения

{

	//Блок объявления матричных/векторных классов
	//Result.Xmean.resize(ChannelCount);
	MatrixXd Xt, XXt;
	double YtXXtY, c, TtT;
	
	Result.N = Y0.rows();//Количество строк
	Result.M = X0.cols();//Количество каналов
	Result.A = Amax;//Количество ГК
	Result.W.resize(Result.M, Result.A);// Матрица взвешенных нагрузок w
	Result.T.resize(Result.N, Result.A);//Матрица счетов
	Result.P.resize(Result.M, Result.A);//Матрица P-нагрузок
	Result.E.resize(Result.N, ChannelCount);//Матрица спектральных остатков остатков Е
	Result.F.resize(Result.N);//Вектор химических остатков F
	Result.Q.resize(Result.A);//Вектор Q-нагрузок

	Result.Xmean = X0.colwise().mean();//Вектор-строка для центров
	MatrixXd X = X0.rowwise()-Result.Xmean;//Центрированная матрица спектров
	Result.Ymean = Y0.mean();//центр концентраций
	VectorXd Y = Y0.array() - Result.Ymean;//Центрированный вектор концентраций
	VectorXd Y1 = Y;//Запоминаем Y перед циклом итераций. Пригодится при расчете регрессионных коэффициентов

	for (int A = 0; A < Amax; ++A)
	//В цикле последовательно расчитываем Главные Компоненты по PLS1 (NIPALS)
	{
		Xt = X.transpose();
		YtXXtY = Y.transpose()*X*Xt*Y;
		c = 1 / std::sqrt(YtXXtY);
		Result.W.col(A) = c * Xt*Y;
		Result.T.col(A) = X * Result.W.col(A);//Расчитываем счета компоненты А
		TtT = Result.T.col(A).transpose()*Result.T.col(A);
		Result.P.col(A)=Xt*Result.T.col(A)/TtT;
		Result.Q.row(A) = (Y.transpose()*Result.T.col(A))/TtT;
		Result.E = X - Result.T.col(A)*Result.P.col(A).transpose();
		Result.F = Y - Result.T.col(A)*Result.Q(A);
		//Копируем переменные для следующей итерации
		Y = Result.F;
		X = Result.E;
	}
	//Вычисляем коэффициенты регрессии
	MatrixXd Tt=Result.T.transpose();
	Result.B = (Tt*Result.T).inverse()*Tt*Y1;
	return;
}

inline void  Calibr::ScorePredictPLS(const VectorXd &B, const MatrixXd &T, double Ymean, VectorXd &Y)
//Функция расчета прогнозного значения Y через счета T и коэффициенты регрессии B
{
	Y = (T * B).array() + Ymean;//Расчет прогнозных значений
}


void Calibr::MatrixToClipboard(MatrixXd X)
{
	Index Cols = X.cols();
	Index Rows = X.rows();
	std::string str;
	for (int i = 0; i < Rows; ++i)

	{
		for (int j = 0; j < Cols; ++j)
		{
			if (Cols - j > 1)
			{
				str.append(std::to_string(int(X(i, j)))).append("\t");
			}
			else
			{
				if (Rows - i > 1)
				{
					str.append(std::to_string(int(X(i, j)))).append("\r\n");
				}
				else
				{
					str.append(std::to_string(int(X(i, j))));
				}
			}
		}
	}
	StringToClipboard(str);
}

void Calibr::VectorToClipboard(VectorXd X)
{
	Index Cols = X.cols();
	Index Rows = X.rows();
	std::string str;
	for (int i = 0; i < Rows; ++i)

	{
		for (int j = 0; j < Cols; ++j)
		{
			if (Cols - j > 1)
			{
				str.append(std::to_string(int(X(i, j)))).append("\t");
			}
			else
			{
				if (Rows - i > 1)
				{
					str.append(std::to_string(int(X(i, j)))).append("\r\n");
				}
				else
				{
					str.append(std::to_string(int(X(i, j))));
				}
			}
		}
	}
	StringToClipboard(str);
}

double Calibr::RMSE(const VectorXd &Y0, const VectorXd &Ycalc)
//Метод расчета RMSE
//Минимизация RMSE - показатель качества градуировки
//Y0 - вектор реальных значений
//Ycalc - вектор расчетных значений
{
	return std::sqrt((Y0 - Ycalc).array().square().sum());
}

void Calibr::NormDecomposePLS
	(const VectorXd &B_LE, MatrixXd &T, double Ymean, const MatrixXd &X0, const VectorXd &Y0, StructPLS &sX)
/*Метод разложения спектральной матрицы с учетом нормализации
Входные параметры:
B_LE - коэффициенты для нахождения LE
T - счета при разложении для нахождения LE
Ymean - средний LE
X0 - исходная спектральная матрица
Y0 - вектор концентраций, соответствующих образцам в матрице спектров
Выходные данные:
sX - структура с результатами нормализованного разложения
*/
{
	VectorXd LEp;// Вектор предсказанных значений
	ScorePredictPLS(B_LE, T, Ymean, LEp);//Получаем предсказанные вектора нормирования
	//Получаем диагональную матрицу нормирования размерностью NxN
	DiagonalMatrix<double, Dynamic> LEn(LEp.array().inverse().matrix());
	MatrixXd Xn = LEn * X0;//Нормированная матрица спектров
	DecomposePLS(Xn, Y0, sX);//Выполняем PLS разложение после нормирования
}

void Calibr::MainCalibrationPLS(int elmnt)
//Главный метод для PLS калибровки с нормированием на внутренний стандарт
{
	Index N = Spectra.rows();//Количество спектров
	VectorXd Conc0;//Вектор расчетных концентраций для начального значения
	MatrixXd Jac(N, Amax);//Якобиан
	VectorXd NewLE(N);//Вектор новых значений LE
	StructPLS NewXCoeff;//Результаты нового разложения нормированных спектров
	VectorXd NewConc;//Новые концентрации при выполнении численного диффиренцирования
	VectorXd dB;//Приращение коэффициентов в итерации Левенб.-Маркв.
	MatrixXd Jt;//Транспонированный якобиан
	MatrixXd JtJ;//Квадратная матрица AxA (Гессиан). Одна матрица нужна для хранения и восстановления гессиана в цикле подбора величины лямбда
	VectorXd F(N);//Вектор-функции невязки. В каждой строке разность между реальным и расчетным значением (минимизируется при решении)
	double iteratErr;//Ошибка нахождения коэффициентов. Условие выхода из алгоритма Левенберга-Марквардта
	StructPLS OldPLS;//Старые параметры разложения по PLS;
	double F_errLimit;//абсолютная допустимая ошибка вычисления невязки функций
	double F_err;//Значение ошибки невязки функции

	
	
	//Выполняем градуировку для элемента с порядковым номером elmnt
	
	TotalPLSRebuildIterat=0;//Счетчик итераций пересчета PLS scores для поиска LE
	do
		//Цикл перерасчета параметров PLS для поиска коэффициентов нормирования через матрицу X
	{
		double lambda{ StartLambda };//Инициализируем лямбда для коррекции скорости решения регрессии методом Левенберга-Марквардта
		int TotalIterat{ 0 };
		iteratErr = LEcoeff.B.norm()*err_relatCompareResult;
		do
		{
			//Полное разложение для текущих значений коэффициентов нормирования. Коэфф. нормирования сохраняем в свойство класса XNormcoeff
			NormDecomposePLS(LEcoeff.B, LEcoeff.T, LEcoeff.Ymean, Spectra, mY.col(elmnt), XNormcoeff);
			//По концентрационным счетам предсказываем расчетные концентрации в Сonc0
			ScorePredictPLS(XNormcoeff.B, XNormcoeff.T, XNormcoeff.Ymean, Conc0);

			omp_set_nested(1);//Разрешаем вложенный параллелизм
#pragma omp parallel for ordered private(NewXCoeff, NewConc)
			for (int i = 0; i < Amax; ++i)
				//В цикле последовательно рассчитываем колонки якобиана
			{
#pragma omp ordered
				{
					VectorXd B_LE2 = LEcoeff.B;//Восстанавливаем B_LE2 для новой итерации цикла
					double DiffStep = B_LE2(i)*err_relatDiffStep;//Абсолютный шаг численного дифференцирования для колонки якобиана

					//Здесь допущение, что минимальном изменении коэффициентов при счетах
					//почти не влияет на сами счета
					B_LE2(i) += DiffStep;//Меняем очередной коэффициент по которому проводится частное дифференцирование

					//Выполняем полное нормированное разложение при новом коэффициенте B_LE2
					NormDecomposePLS(B_LE2, LEcoeff.T, LEcoeff.Ymean, Spectra, mY.col(elmnt), NewXCoeff);

					//Выполняем расчет новых прогнозных значений согласно новых счетов полученных в полном разложении
					ScorePredictPLS(NewXCoeff.B, NewXCoeff.T, NewXCoeff.Ymean, NewConc);
					VectorXd dYdX = (NewConc - Conc0) / DiffStep;//Численное дифференцирование. i-я колонка якобиана
					Jac.col(i) = dYdX;
				}
			}

			F = mY.col(elmnt) - Conc0;//Задаем вектор-функцию невязки
			F_errLimit = F.norm() / std::sqrt(N)*err_relatCompareResult;//абсолютная допустимая ошибка вычисления невязки функций

			Jt = Jac.transpose();
			JtJ = Jt * Jac;
			int iteratCount = 0;
			bool DoIterat;//Флаг продолжения итерации при большем лямбда, если при текущем расходится
			VectorXd NewB; //Коэффициенты регрессии, которые будут действовать на выходе из цикла
			do
				//Цикл подбора величины Лямбда для регулирования шага сходимости
			{
				//Обсчитываем в параллельном цикле результат сразу для нескольких лямбда.
				lambda /= 2;//Предварительно уменьшаем лямбда на случай, если функция будет сходиться при большем шаге
				int MaxThreads = omp_get_max_threads();//Получаем максимальное количество потоков для данной архитектуры
				VectorXd minF2 = F;//Минимальный вектор F2 среди расчитанных


				//Вход в параллельную секцию
#pragma omp parallel 
				{
					//Во всех потоках инициализируем локальные переменные
					VectorXd NewB_local;//Вектор коэффициентов, рассчитываемый параллельно в цикле for
					VectorXd db_local;//Вектор приращения коэффициентов в параллельном цикле for
					VectorXd F2(N);
					double lambda_local;//Коэффициент лямбда в параллельном потоке

#pragma omp	for ordered private(NewXCoeff, NewConc)
					for (int i = 0; i < MaxThreads; ++i)
#pragma omp ordered
					{
						lambda_local = lambda * pow(2, i);//В потоках лямбда отличается друг от друга кратно 2
						//Cоздаем диагональную матрицу  NxN со значениями Лямбда в диагонали
						MatrixXd lambdI = MatrixXd::Zero(Amax, Amax);
						lambdI.diagonal().array() = lambda_local;//Главной диагонали присвоили значение лямбда
						db_local = (JtJ + lambdI).inverse()*Jt*F; //Вычисляем приращения для вектора коэффициентов
						NewB_local = LEcoeff.B + db_local;//Корректируем коэффициенты B с учетом приращения
						//Снова проводим нормированное PLS-разложение для проверки сходимости невязки к минимуму при новых коэфф. B
						NormDecomposePLS(NewB_local, LEcoeff.T, LEcoeff.Ymean, Spectra, mY.col(elmnt), NewXCoeff);
						//Выполняем расчет новых прогнозных значений при новых коэффициентах
						ScorePredictPLS(NewXCoeff.B, NewXCoeff.T, NewXCoeff.Ymean, NewConc);
						F2 = mY.col(elmnt) - NewConc;//Считаем невязку при новых коэффициентах	
					}

#pragma omp critical(bestLambdaSearch)
					//Из локальных переменных F2 каждого потока выбираем с лучшей сходимостью
					//Если при данном наборе функция расходится, код в блок if не попадет,
					//тогда на выходе из параллельной секции F==minF2
					{
						if (F2.norm() < minF2.norm())
						{
							//Сохраняем лучшие данные
							minF2 = F2;
							NewB = NewB_local;
							lambda = lambda_local;
							dB = db_local;
						}
					}
				}

				//Вышли из параллельной секции

				F_err = std::abs((F.norm() - minF2.norm()) / std::sqrt(N));
				if (minF2 == F)
					//Функция расходится. Повторяем итерацию в цикле Do..While при большем лямбда.
					//С увеличением лямбда приближение идет более мелкими шагами
				{
					lambda *= pow(2, MaxThreads);
					DoIterat = true;
				}
				else
				{
					DoIterat = false;//При данном лямбда регрессия сходится. Можно выходить из цикла
				}

				/*if ((TotalPLSRebuildIterat == 18) && (TotalIterat == 20))
				{
					std::cout << "stop\n";
				}*/
				++iteratCount;
				/*Условия выхода из цикла:
				1. Превышено количество максимальных итераций
				2. Лямбда слишком большое, не удалось подобрать сходимость функции
				3. Невязка при новой лямбда отличается от невязки при старой лямбда меньше чем ошибка вычисления
				*/
			} while ((iteratCount < MaxLambdaIterat) && DoIterat);
			if (iteratCount == MaxLambdaIterat)
				//Функция расходится. Заканчиваем работу
			{
				std::cout << "Функция расходится. Расчет прерван" << std::endl << std::endl;
				return;
			}
			//Присваиваем системе новые коэффициенты
			LEcoeff.B = NewB;
			++TotalIterat;
			RMSEC = F.norm() / std::sqrt(N);
			std::cout << "\r" << "Элемент " << elmnt + 1 << ", Перестроено PLS: " << TotalPLSRebuildIterat << ", Шаг: " << TotalIterat;
			std::cout << ", ошибка dB в LevMaq: " << dB.norm() / std::sqrt(Amax) << ", RMSEC= " << F.norm() / std::sqrt(N) << "   ";
			/*Условия выхода из цикла:
			1. Превышено количество максимальных итераций
			2. Ошибка по коэффициентам не превышает ошибку вычислений
			3. Изменение невязки функций не превышает ошибку вычислений
			*/
		} while ((TotalIterat < MaxIterationsLevMaq) && (dB.norm() / std::sqrt(Amax) > iteratErr) && (F_err > F_errLimit));

		//Перерасчитываем параметры PLS для поиска LE
		++TotalPLSRebuildIterat;
		ScorePredictPLS(LEcoeff.B, LEcoeff.T, LEcoeff.Ymean, LE);//Получили новые значения LE для переразложения
		//Сохраняем предыдущие значения коэффициентов, чтобы выйти из цикла когда итерации сойдутся
		OldPLS = LEcoeff;
		DecomposePLS(Spectra, LE, LEcoeff);//Выполнили переразложение
		std::cout << std::endl << "Ошибка коэффициентов B_LE между итерациями: " << (OldPLS.B - LEcoeff.B).norm() / std::sqrt(Amax) << std::endl;
		
		if (TotalPLSRebuildIterat % 10 == 0)
		//Каждые 10 шагов итерации сохраняемся в файл 
		{
			FinalPLS = 0;//Выход из цикла не выполнен, продолжаем расчет
			SaveResultsPLS();
		}
	} while ((TotalPLSRebuildIterat < MaxPLSRebuildIterat) && ((OldPLS.B - LEcoeff.B).norm() / std::sqrt(Amax) > iteratErr));

	std::cout << std::endl << std::endl << "Коэффициенты LE для " << elmnt + 1 << "-го элемента:" << std::endl << LEcoeff.B << std::endl << std::endl;
	std::cout << "RMSEC= " << RMSEC << std::endl << std::endl;
	FinalPLS = 1;//Расчеты выполнены, сохраняем результаты вычислений в файл с признаком окончания расчетов
	SaveResultsPLS();//Сохраняем информацию о параметрах разложения LE
}


void Calibr::InitLE_SetSumX()
//В качестве начальных значений для LE используется нормирование на сумму всех интенсивностей в спектре
{
	LE = Spectra.rowwise().sum();
	LEInitialType = "SummNorm";//Нормирование на сумму интенсивностей
}

void Calibr::InitLE_SetMeanX()
//В качестве начальных значений для LE используется нормирование на среднюю интенсивность в спектре
{
	LE = Spectra.rowwise().mean();
	LEInitialType = "MeanNorm";//Нормирование на среднюю интенсивность в спектре
}

void Calibr::InitLE_SetMaxX()
//В качестве начальных значений для LE используется нормирование на интенсивность максимальной линии в спектре
{
	LE = Spectra.rowwise().maxCoeff();
	LEInitialType = "MaxNorm";//Нормирование на максимальную интенсивность в спектре
}

void Calibr::InitLE_NonCoherentBackScatter(const std::string LEFileName)
//нормирование на интенсивность некогерентного обратного рассеяния
{
	std::string FileName = SpectraPath;
	LoadVectorLong(FileName.append(LEFileName), LE);//Инициализация вектора LE данными из файла
	LEInitialType = "BckScattrng";//инициализация вектора LE значениями интенсивностей обратного рассеяния
}

void Calibr::SaveResultsPLS()
/*Сохранение результатов вычисления на диск в каталог градуировки
	
*/
{
	//Создаем имя файла
	std::string FileName = "A";
	FileName.append(std::to_string(this->Amax));
	FileName.append(LEInitialType);
	FileName.append("_Iter_");
	FileName.append(std::to_string(TotalPLSRebuildIterat));
	FileName.append("RMSEC");
	FileName.append(std::to_string(RMSEC));

	//Добавляем расширение для файла c данными нормирования по PLS
	FileName.append(".nPLS");
	bool res;
	std::string FullFileName{ CalibrationDataPath };
	FullFileName.append(FileName);
	FileClass ResultPLS;
	res=ResultPLS.OpenForSave(FullFileName);//Открываем файл для записи
	if (!res)
	{
		std::cout << "Не удалось открыть файл для записи параметров градуировки\n";
		std::system("Pause");
		exit(FILE_SAVING_ERROR);
	}
	ResultPLS.SaveObject(CalibrMethod, CalibrMethodStringLength);//Сохраняем в файл тип используемой калибровки
	ResultPLS.SaveObject(LEcoeff);
	ResultPLS.SaveObject(XNormcoeff);
	ResultPLS.SaveObject(TotalPLSRebuildIterat);
	ResultPLS.SaveObject(FinalPLS);
	ResultPLS.Close();
}

ProcessError Calibr::LoadResultsPLS(std::string FileName)
//Подгружаем результаты расчетов из файла
{
	bool res;
	std::string FullFileName{ CalibrationDataPath };
	FullFileName.append(FileName);
	FileClass ResultPLS;
	res = ResultPLS.OpenForRead(FullFileName);//Открываем файл для чтения
	if (!res)
	{
		std::cout << "Не удалось открыть файл\n" << FullFileName<<"\nдля чтения данных PLS-градуировки";
		std::system("Pause");
		return FILE_READING_ERROR;
	}
	ResultPLS.LoadObject(CalibrMethod, CalibrMethodStringLength);//читаем из файла тип используемой калибровки
	if (_stricmp(CalibrMethod, "PLS"))
	//строки не равны. Калибровка не PLS выходим по ошибке
	{
		std::cout << "Файл:\n"<< FullFileName<<"\nпредназначен для калибровки методом "<< CalibrMethod<<"\n"
			<<"Выберите файл с калибровкой по методу PLS\n\n";
		std::system("Pause");
	}
	ResultPLS.LoadObject(LEcoeff);
	ResultPLS.LoadObject(XNormcoeff);
	ResultPLS.LoadObject(TotalPLSRebuildIterat);
	ResultPLS.LoadObject(FinalPLS);
	ResultPLS.Close();
}

void Calibr::LoadElvaXSpectrum(const std::string FullFileName, RowVectorXd &X)
//Загрузка интенсивностей из файла ElvaX в вектор-строку
{
	ifstream spectrum;
	spectrum.open(FullFileName, ios::binary);
	spectrum.seekg(0x0100);
	DWORD ChannelCount;
	spectrum.read((char*) &ChannelCount, sizeof(DWORD));
	//Переопределили вектор-строку под размер спектра в файле. Размер спектра уменьшили на 1, так как последний канал не используется
	X.resize(--ChannelCount);
	DWORD* Arr = new DWORD[ChannelCount];//Определили массив под чтение спектра
	spectrum.read((char*)Arr, ChannelCount * sizeof(DWORD));//Считываем интенсивности в каналах
	spectrum.close();
	for (DWORD i = 0; i < ChannelCount; ++i)
		X(i) = Arr[i];
	delete[] Arr;
	
}

void Calibr::LoadElvaXSpectra(const std::string &Path, MatrixXd &X)
//Загрузка всех спектров каталога в матрицу X
{
	//определяем количество файлов в каталоге
	std::vector<std::string> FileList;//Вектор со списком файлов
	FileClass f;
	std::string mask = Path;
	mask=CorrectPath(mask);
	mask.append("*.evt");
	f.FindFileList(mask, FileList);

	//В цикле создаем матрицу интенсивностей для всех имеющихся файлов
	for (int i=0; i < FileList.size(); ++i)
	{
		RowVectorXd Xrow;
		std::string FullFileName = WorkingPath;
		FullFileName.append(FileList[i]);
		LoadElvaXSpectrum(FullFileName, Xrow);
		if (i == 0)
		{
			X = Xrow;
		}
		else
		{
			X.conservativeResize(X.rows() + 1, X.cols());
			X.row(X.rows() - 1) = Xrow;
		}
	}
}