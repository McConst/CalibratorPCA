
#define _CRT_SECURE_NO_WARNINGS//��������� �������������� ���������� CRT � �����������

#include "pch.h"
#include "Calibr.h"
#include "Split.h"
#include "MyFuncs.h"
#include "FileClass.h"//����� ������ � �������
#include "const.h"//�������� ��������
#include <string>
#include <iostream>
#include <stdlib.h>//���������� ��� �������������� string � double
#include <fstream>//���������� ������ � �������
#include <vector>
#include <Eigen/Dense>//���������� Eigen
#include <Windows.h>//WinAPI
#include <map>//���������� map ��� ������������� �������� �������������� �������������� �� �� �����
#include <omp.h>//���������� ������������� ���������������� OpenMP

using namespace Eigen;

Calibr::Calibr(const std::string &InitFileName)
//� ������������ ��������� ���� �������� � ����������� �� �������� ��������� ������
{
	std::ifstream file;
	file.open(InitFileName);//������� ������ file ��� ������ � ��������� �� �����
	if (file.is_open())
	{
		std::string str;
		while (!file.eof())
		{
			getline(file, str);
			if (str.substr(0, 2) == "//")//���������� �����������
			{
				continue;
			}
			else if (str.length() == 0)//��������� �� ��������� ������ ������
			{
				continue;
			}
			//�������������� ��������
			size_t pos = str.find("\t");
			std::string operat = str.substr(0, pos);
			pos = str.rfind("\t") + 1;
			std::string value = str.substr(pos, str.length() - pos);
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
				WorkingPath = CorrectPath(value);
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
			else if (operat=="SpectaLoadingMethod")
			{
				SpectaLoadingMethod = value[0]-0x30;
			}
		}
	}
	else
	{
		std::cout << "���� ������������� �� ������ ��� �� ��������" << std::endl;
	}
	file.close();//��������� ����
	FinalPLS = 0;
	FinalPCR = 0;
	TotalPLSRebuildIterat = 0;//������� ��������� ���������� ��������
}



Calibr::~Calibr()
{
	//std::cout << "���������� Calibr"<<std::endl;
}


void Calibr::LoadInitDataForCalibrat(const std::string &SpectraFileName, const std::string &YFileName, const std::string &NamesFile)
/*	������� ��������� ������������� ������� �������� � ������������ �� ���������� �� ������
SpectraFileName - ���� ������� �������� 
YFileName - ���� � �������������� ���������� ������������, ��������������� ������� ������� �� �������
LE ���������������� ���������������� ���������� �������� �� �����
*/
{
	if (SpectaLoadingMethod==1)
	//������������� ������������ ������ �� ������ ������������ *.evt � ������ Names.txt � Y.txt
	{
		std::vector<std::vector<double> > Y0;//��������� �������� ������ � �������������� ����������
		LoadY(YFileName, Y0);//���������� � ��������� ������ ���� Vector ������������� ������������ �� �����
		std::vector<std::string> Names;
		LoadNames(Names);
		map <std::string, int> ConcIndex;//��������� ��� ������ ������������� ������������ �� ����� �����
		for (int i = 0; i < Names.size(); ++i)
			ConcIndex[Names[i]] = i;//����������� ����� ������� ������, �� �������� ����� ������ ������������

		//��������� �� ������ *.evt ������������ ������
		//������ �������� ���� �������� � LoadSpectraNames, �� �������� ����� ���������������� ������ ������������ mY
		std::vector < std::string > LoadedSpectraNames= LoadElvaXSpectra(SpectraPath, Spectra);
		mY.resize(SpectraCount, CRM_ElementCount);
		for (int i = 0; i < SpectraCount; ++i)
			for (int j = 0; j < CRM_ElementCount; ++j)
			{
				std::string str = LoadedSpectraNames[i].substr(0, LoadedSpectraNames[i].length() - 4);
				mY(i, j) = Y0[ConcIndex[str]][j];
			}
	}
	else
	//������ ������������ ������ �� ����� SpectraFileName � ������������� ������������ �� ����� YFileName
	{
		//�������������� ��������� ������������ ��������� ������ ��� ������ ������������ �� ������
		std::vector<std::vector<double> > Y0; // ��������� ��������� ������ ������������� ������������. ������ �������. - ��, ����. �������. - �������
		LoadY(YFileName, Y0);//���������� � ��������� ������ ���� Vector ������������� ������������ �� �����

		//������������� ������������ �� ���������� ������������� ������� � ������� Eigen
		mY.resize(SpectraCount, CRM_ElementCount);
		for (int i = 0; i < SpectraCount; ++i)
			for (int j = 0; j < CRM_ElementCount; ++j)
			{
				mY(i, j) = Y0[i][j];
			}

		//�������������� ������������ ������ �� ����� Spectra.dat
		std::string FileName = SpectraPath;
		LoadMatrixLong(FileName.append(SpectraFileName), Spectra);//������������� ������� Spectra ������� �� �����
	}
}

void Calibr::LoadNames(std::vector<std::string> &Names, const std::string &NamesFileName)
//� ������ Names �� ����� NamesFileName ���������� ����� ����������� ��������
//����� ���������� ��� ������ ������������� �������� ��� �������
{
	std::string FileName = SpectraPath;
	FileName.append(NamesFileName);
	std::ifstream file;
	file.open(FileName);//������� ������ file ��� ������ � ��������� �� �����
	if (file.is_open())
	{
		Names.clear();//������� ������ ����� �����������
		std::string str;
		while (!file.eof())
		{
			getline(file, str);
			if (str != "\0")
			{
				Names.push_back(str);
			}
		}
	}
	else
	{
		std::cout << "���� ������������� �������� ��� ����������� �� ������ ��� �� ��������" << std::endl;
	}
	file.close();//��������� ����

}


void Calibr::LoadY(const std::string &YFileName, std::vector<std::vector<double> > &Y0)
//��������� ������������ �� ����� YFileName � ������ ��� ���������� ���������

{
	std::string FileName = SpectraPath;
	FileName.append(YFileName);
	std::ifstream file;
	file.open(FileName);//������� ������ file ��� ������ � ��������� �� �����
	if (file.is_open())
	{
		std::string str; Split Splitter;
		SpectraCount = 0;
		std::vector<double> Conc(CRM_ElementCount, 0);//������ Conc ��� ������������� ���������� ������� �������������� �� �� �����
		while (!file.eof())
		{
			getline(file, str);
			if (str != "\0")
			{
				++SpectraCount;
				Splitter(str, "\t");//����� ������ �� �������� ���������
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
		std::cout << "���� ������������� �������� ��� ����������� �� ������ ��� �� ��������" << std::endl;
	}
	file.close();//��������� ����
}

int Calibr::LoadMatrixLong(const std::string &FileName, MatrixXd &X)
//������������� ��������� ������� X ������� �� �����
{
	int* Arr = new int[SpectraCount*ChannelCount];//������� ������ ��� ������ ������������� ��������

	// ����������� ������ ���� String � w_string
	LPWSTR widestrFileName = StringToW_Char(FileName);

	HANDLE hFile;//��������� ���������� API
	hFile = CreateFile(widestrFileName, GENERIC_READ, 0, 0, OPEN_EXISTING, FILE_ATTRIBUTE_NORMAL, 0);
	delete[] widestrFileName;

	int DataSize = SpectraCount * ChannelCount * sizeof(int);//������ ���� ������������� �� ������ �� �����
	DWORD nBytesRead;//������, ���������� ������� �� �����
	int res = ReadFile(hFile, Arr, DataSize, &nBytesRead, 0);//� ����� ������ Spectra ���������� �����������!!!!
	CloseHandle(hFile);
	X.resize(SpectraCount, ChannelCount);
	for (int i = 0; i < SpectraCount*ChannelCount; ++i)
	{
		X(i % SpectraCount, i / SpectraCount) = Arr[i];//������������� ������� Spectra ��� ��������
	}
	//std::cout << Spectra.block<10, 10>(0, 0)<<std::endl;//���������� ���� 10�10 ��������� ������� � �������� 0,0
	delete[] Arr;
	return SpectraCount;
}

int Calibr::LoadVectorLong(const std::string &FileName, VectorXd &X)
//������������� ��������� ������� X ������� �� �����
{
	int* Arr = new int[SpectraCount];//������� ������ ��� ������ ������������� ��������


	LPWSTR widestrFileName = StringToW_Char(FileName);

	HANDLE hFile;//��������� ���������� API
	hFile = CreateFile(widestrFileName, GENERIC_READ, 0, 0, OPEN_EXISTING, FILE_ATTRIBUTE_NORMAL, 0);
	delete[] widestrFileName;

	int DataSize = SpectraCount * sizeof(int);//������ ���� ������������� �� ������ �� �����
	DWORD nBytesRead;//������, ���������� ������� �� �����
	int res = ReadFile(hFile, Arr, DataSize, &nBytesRead, 0);//� ����� ������ Spectra ���������� �����������!!!!
	CloseHandle(hFile);
	X.resize(SpectraCount);
	for (int i = 0; i < SpectraCount; ++i)
	{
		X(i) = Arr[i];//������������� ������� �������, ������������ �� �����
	}
	delete[] Arr;
	return SpectraCount;
}


int Calibr::LoadMatrixDouble(const std::string &FileName)
//�������� � ������������� ������� ������� �� �����
{
	double* Arr = new double[SpectraCount*ChannelCount];//������� ������ ��� ������ ��������

// ����������� ������ ���� String � w_string
	LPWSTR widestrFileName = StringToW_Char(FileName);

	HANDLE hFile;//��������� ���������� API
	hFile = CreateFile(widestrFileName, GENERIC_READ, 0, 0, OPEN_EXISTING, FILE_ATTRIBUTE_NORMAL, 0);
	delete[] widestrFileName;

	int DataSize = SpectraCount * ChannelCount * sizeof(double);//������ ���� ������������� �� ������ �� �����
	DWORD nBytesRead;//������, ���������� ������� �� �����
	int res = ReadFile(hFile, Arr, DataSize, &nBytesRead, 0);//� ����� ������ Spectra ���������� �����������!!!!
	CloseHandle(hFile);
	Spectra.resize(SpectraCount, ChannelCount);
	for (int i = 0; i < SpectraCount*ChannelCount; ++i)
	{
		Spectra(i % SpectraCount, i / SpectraCount) = Arr[i];//������������� ������� Spectra ��� ��������
	}
	//std::cout << Spectra.block<10, 10>(0, 0)<<std::endl;//���������� ���� 10�10 ��������� ������� � �������� 0,0
	delete[] Arr;
	return SpectraCount;
}

void Calibr::DecomposePLS(const MatrixXd &X0, const VectorXd &Y0, StructPLS &Result)
//���������� ������ X � Y �� PLS1 NIPALS
// A-���������� �� ��� ����������
//� ���������� Result ���� StructPLS ���������� ���������� ����������

{

	//���� ���������� ���������/��������� �������
	//Result.Xmean.resize(ChannelCount);
	MatrixXd Xt, XXt;
	double YtXXtY, c, TtT;

	Result.N = Y0.rows();//���������� �����
	Result.M = X0.cols();//���������� �������
	Result.A = Amax;//���������� ��
	Result.W.resize(Result.M, Result.A);// ������� ���������� �������� w
	Result.T.resize(Result.N, Result.A);//������� ������
	Result.P.resize(Result.M, Result.A);//������� P-��������
	Result.E.resize(Result.N, ChannelCount);//������� ������������ �������� �������� �
	Result.F.resize(Result.N);//������ ���������� �������� F
	Result.Q.resize(Result.A);//������ Q-��������

	Result.Xmean = X0.colwise().mean();//������-������ ��� �������
	MatrixXd X = X0.rowwise() - Result.Xmean;//�������������� ������� ��������
	Result.Ymean = Y0.mean();//����� ������������
	VectorXd Y = Y0.array() - Result.Ymean;//�������������� ������ ������������
	VectorXd Y1 = Y;//���������� Y ����� ������ ��������. ���������� ��� ������� ������������� �������������

	for (int A = 0; A < Amax; ++A)
		//� ����� ��������������� ����������� ������� ���������� �� PLS1 (NIPALS)
	{
		Xt = X.transpose();
		YtXXtY = Y.transpose()*X*Xt*Y;
		c = 1 / std::sqrt(YtXXtY);
		Result.W.col(A) = c * Xt*Y;
		Result.T.col(A) = X * Result.W.col(A);//����������� ����� ���������� �
		TtT = Result.T.col(A).transpose()*Result.T.col(A);
		Result.P.col(A) = Xt * Result.T.col(A) / TtT;
		Result.Q.row(A) = (Y.transpose()*Result.T.col(A)) / TtT;
		Result.E = X - Result.T.col(A)*Result.P.col(A).transpose();
		Result.F = Y - Result.T.col(A)*Result.Q(A);
		//�������� ���������� ��� ��������� ��������
		Y = Result.F;
		X = Result.E;
	}
	return;
}

inline void  Calibr::ScorePredictPLS(const VectorXd &Q, const MatrixXd &T, double Ymean, VectorXd &Y)
//������� ������� ����������� �������� Y ����� ����� T � ������������ ��������� B
{
	Y = (T * Q).array() + Ymean;//������ ���������� ��������
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
				str.append(std::to_string(X(i, j))).append("\t");
			}
			else
			{
				if (Rows - i > 1)
				{
					str.append(std::to_string(X(i, j))).append("\r\n");
				}
				else
				{
					str.append(std::to_string(X(i, j)));
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
				str.append(std::to_string(X(i, j))).append("\t");
			}
			else
			{
				if (Rows - i > 1)
				{
					str.append(std::to_string(X(i, j))).append("\r\n");
				}
				else
				{
					str.append(std::to_string(X(i, j)));
				}
			}
		}
	}
	StringToClipboard(str);
}

double Calibr::RMSE(const VectorXd &Y0, const VectorXd &Ycalc)
//����� ������� RMSE
//����������� RMSE - ���������� �������� �����������
//Y0 - ������ �������� ��������
//Ycalc - ������ ��������� ��������
{
	return std::sqrt((Y0 - Ycalc).array().square().sum());
}

void Calibr::NormDecomposePLS
(const VectorXd &Q_LE, MatrixXd &T, double Ymean, const MatrixXd &X0, const VectorXd &Y0, StructPLS &sX)
/*����� ���������� ������������ ������� � ������ ������������
������� ���������:
B_LE - ������������ ��� ���������� LE
T - ����� ��� ���������� ��� ���������� LE
Ymean - ������� LE
X0 - �������� ������������ �������
Y0 - ������ ������������, ��������������� �������� � ������� ��������
�������� ������:
sX - ��������� � ������������ ���������������� ����������
*/
{
	VectorXd LEp;// ������ ������������� ��������
	ScorePredictPLS(Q_LE, T, Ymean, LEp);//�������� ������������� ������� ������������
	//�������� ������������ ������� ������������ ������������ NxN
	DiagonalMatrix<double, Dynamic> LEn(LEp.array().inverse().matrix());
	MatrixXd Xn = LEn * X0;//������������� ������� ��������
	DecomposePLS(Xn, Y0, sX);//��������� PLS ���������� ����� ������������
}

void Calibr::MainCalibrationPLS(int elmnt)
//������� ����� ��� PLS ���������� � ������������� �� ���������� ��������
{
	Index N = Spectra.rows();//���������� ��������
	VectorXd Conc0;//������ ��������� ������������ ��� ���������� ��������
	MatrixXd Jac(N, Amax);//�������
	VectorXd NewLE(N);//������ ����� �������� LE
	StructPLS NewXCoeff;//���������� ������ ���������� ������������� ��������
	VectorXd NewConc;//����� ������������ ��� ���������� ���������� �����������������
	VectorXd dQ;//���������� ���. �������� � �������� ������.-�����.
	MatrixXd Jt;//����������������� �������
	MatrixXd JtJ;//���������� ������� AxA (�������). ���� ������� ����� ��� �������� � �������������� �������� � ����� ������� �������� ������
	VectorXd F(N);//������-������� �������. � ������ ������ �������� ����� �������� � ��������� ��������� (�������������� ��� �������)
	double iteratErr;//������ ���������� �������������. ������� ������ �� ��������� ����������-����������
	StructPLS OldPLS;//������ ��������� ���������� �� PLS;
	double F_errLimit;//���������� ���������� ������ ���������� ������� �������
	double PLS_errIterat;//���������� ��������� ������������� ����� ���������� PLS. ��� ������ �� ���������� �������� ����� ��������
	double F_err;//�������� ������ ������� �������



	//��������� ����������� ��� �������� � ���������� ������� elmnt

	TotalPLSRebuildIterat = 0;//������� �������� ��������� PLS scores ��� ������ LE
	do
		//���� ����������� ���������� PLS ��� ������ ������������� ������������ ����� ������� X
	{
		double lambda{ StartLambda };//�������������� ������ ��� ��������� �������� ������� ��������� ������� ����������-����������
		int TotalIterat{ 0 };
		iteratErr = LEcoeff.Q.norm()*err_relatCompareResult;
		do
		{
			//������ ���������� ��� ������� �������� ������������� ������������. �����. ������������ ��������� � �������� ������ XNormcoeff
			NormDecomposePLS(LEcoeff.Q, LEcoeff.T, LEcoeff.Ymean, Spectra, mY.col(elmnt), XNormcoeff);
			//�� ���������������� ������ ������������� ��������� ������������ � �onc0
			ScorePredictPLS(XNormcoeff.Q, XNormcoeff.T, XNormcoeff.Ymean, Conc0);

			omp_set_nested(1);//��������� ��������� �����������
#pragma omp parallel for ordered private(NewXCoeff, NewConc)
			for (int i = 0; i < Amax; ++i)
				//� ����� ��������������� ������������ ������� ��������
			{
#pragma omp ordered
				{
					VectorXd Q_LE2 = LEcoeff.Q;//��������������� Q_LE2 ��� ����� �������� �����
					double DiffStep = Q_LE2(i)*err_relatDiffStep;//���������� ��� ���������� ����������������� ��� ������� ��������

					//����� ���������, ��� ����������� ��������� ������������� ��� ������
					//����� �� ������ �� ���� �����
					Q_LE2(i) += DiffStep;//������ ��������� ����������� �� �������� ���������� ������� �����������������

					//��������� ������ ������������� ���������� ��� ����� ������������ B_LE2
					NormDecomposePLS(Q_LE2, LEcoeff.T, LEcoeff.Ymean, Spectra, mY.col(elmnt), NewXCoeff);

					//��������� ������ ����� ���������� �������� �������� ����� ������ ���������� � ������ ����������
					ScorePredictPLS(NewXCoeff.Q, NewXCoeff.T, NewXCoeff.Ymean, NewConc);
					VectorXd dYdX = (NewConc - Conc0) / DiffStep;//��������� �����������������. i-� ������� ��������
					Jac.col(i) = dYdX;
				}
			}

			F = mY.col(elmnt) - Conc0;//������ ������-������� �������
			F_errLimit = F.norm() / std::sqrt(N)*err_relatCompareResult;//���������� ���������� ������ ���������� ������� �������

			Jt = Jac.transpose();
			JtJ = Jac.transpose() * Jac;
			int iteratCount = 0;
			bool DoIterat;//���� ����������� �������� ��� ������� ������, ���� ��� ������� ����������
			VectorXd NewQ; //�������� ������� ��������, ������� ����� ����������� �� ������ �� �����
			do
				//���� ������� �������� ������ ��� ������������� ���� ����������
			{
				//����������� � ������������ ����� ��������� ����� ��� ���������� ������.
				lambda /= 2;//�������������� ��������� ������ �� ������, ���� ������� ����� ��������� ��� ������� ����
				int MaxThreads = omp_get_max_threads();//�������� ������������ ���������� ������� ��� ������ �����������
				VectorXd minF2 = F;//����������� ������ F2 ����� �����������


				//���� � ������������ ������
#pragma omp parallel 
				{
					//�� ���� ������� �������������� ��������� ����������
					VectorXd NewQ_local;//������ �������������, �������������� ����������� � ����� for
					VectorXd NewConc_local;//������ ��������� ������������
					VectorXd dQ_local;//������ ���������� ������������� � ������������ ����� for
					VectorXd F2(N);
					double lambda_local;//����������� ������ � ������������ ������

#pragma omp	for ordered private(NewXCoeff)
					for (int i = 0; i < MaxThreads; ++i)
#pragma omp ordered
					{
						lambda_local = lambda * pow(2, i);//� ������� ������ ���������� ���� �� ����� ������ 2
						//C������ ������������ �������  NxN �� ���������� ������ � ���������
						MatrixXd lambdI = MatrixXd::Zero(Amax, Amax);
						lambdI.diagonal().array() = lambda_local;//������� ��������� ��������� �������� ������
						dQ_local = (JtJ + lambdI).inverse()*Jt*F; //��������� ���������� ��� ������� �������������
						NewQ_local = LEcoeff.Q + dQ_local;//������������ ������ Q � ������ ����������
						//����� �������� ������������� PLS-���������� ��� �������� ���������� ������� � �������� ��� ����� Q
						NormDecomposePLS(NewQ_local, LEcoeff.T, LEcoeff.Ymean, Spectra, mY.col(elmnt), NewXCoeff);
						//��������� ������ ����� ���������� �������� ��� ����� �������������
						ScorePredictPLS(NewXCoeff.Q, NewXCoeff.T, NewXCoeff.Ymean, NewConc_local);
						F2 = mY.col(elmnt) - NewConc_local;//������� ������� ��� ����� �������������	
					}

#pragma omp critical(bestLambdaSearch)
					//�� ��������� ���������� F2 ������� ������ �������� � ������ �����������
					//���� ��� ������ ������ ������� ����������, ��� � ���� if �� �������,
					//����� �� ������ �� ������������ ������ F==minF2
					{
						if (F2.norm() < minF2.norm())
						{
							//��������� ������ ������
							minF2 = F2;
							NewConc = NewConc_local;
							NewQ = NewQ_local;
							lambda = lambda_local;
							dQ = dQ_local;
						}
					}
				}

				//����� �� ������������ ������

				F_err = std::abs((F - minF2).norm() / std::sqrt(N));
				if (minF2 == F)
					//������� ����������. ��������� �������� � ����� Do..While ��� ������� ������.
					//� ����������� ������ ����������� ���� ����� ������� ������
				{
					lambda *= pow(2, MaxThreads);
					DoIterat = true;
				}
				else
				{
					DoIterat = false;//��� ������ ������ ��������� ��������. ����� �������� �� �����
				}

				++iteratCount;
				/*������� ������ �� �����:
				1. ��������� ���������� ������������ ��������
				2. ������ ������� �������, �� ������� ��������� ���������� �������
				3. ������� ��� ����� ������ ���������� �� ������� ��� ������ ������ ������ ��� ������ ����������
				*/
			} while ((iteratCount < MaxLambdaIterat) && DoIterat);
			if (iteratCount == MaxLambdaIterat)
				//������� ����������. ����������� ������
			{
				std::cout << "������� ����������. ������ �������" << std::endl << std::endl;
				return;
			}
			//����������� ������� ����� ������������
			LEcoeff.Q = NewQ;
			++TotalIterat;
			RMSEC = F.norm() / std::sqrt(N);
			std::cout << "\r" << "������� " << elmnt + 1 << ", ����������� PLS: " << TotalPLSRebuildIterat << ", ���: " << TotalIterat;
			std::cout << ", ������ dQ � LevMaq: " << dQ.norm() / std::sqrt(Amax) << ", RMSEC= " << F.norm() / std::sqrt(N) << "   ";

			/*������� ������ �� �����:
			1. ��������� ���������� ������������ ��������
			2. ������ �� ������������� �� ��������� ������ ����������
			3. ��������� ������� ������� �� ��������� ������ ����������
			*/
		} while ((TotalIterat < MaxIterationsLevMaq) && (dQ.norm() / std::sqrt(Amax) > iteratErr) && (F_err > F_errLimit));

		//��������������� ��������� PLS ��� ������ LE
		++TotalPLSRebuildIterat;
		ScorePredictPLS(LEcoeff.Q, LEcoeff.T, LEcoeff.Ymean, LE);//�������� ����� �������� LE ��� ��������������
		//��������� ���������� �������� �������������, ����� ����� �� ����� ����� �������� ��������
		OldPLS = LEcoeff;
		DecomposePLS(Spectra, LE, LEcoeff);//��������� ��������������
		PLS_errIterat = (OldPLS.Q - LEcoeff.Q).norm() / std::sqrt(Amax);//��������� ������������� ����� ���������� PLS ����������
		std::cout << std::endl << "��������� ������� Q_LE ����� ����������: " << PLS_errIterat << std::endl;
		if (TotalPLSRebuildIterat % 10 == 0)
			//������ 10 ����� �������� ����������� � ���� 
		{
			FinalPLS = 0;//����� �� ����� �� ��������, ���������� ������
			SaveResultsPLS();
			cout << "��������� ������������ ������ 25 ��������:\n" << NewConc.block<25, 1>(0, 0) << "\n\n";
		}
	/*������� ������ �� �����
	1. �������� ��������� ������ ��� ������ � ������������
	2. ������ ����� �������������� �������������� � ��������� PLS ������ �������� ������
	*/
	} while ((TotalPLSRebuildIterat < MaxPLSRebuildIterat) &&
		PLS_errIterat>LEcoeff.Q.norm()*err_relatPLSiterat);

	std::cout << std::endl << std::endl << "������ Q_LE ��� " << elmnt + 1 << "-�� ��������:" << std::endl << LEcoeff.Q << std::endl << std::endl;
	std::cout << "RMSEC= " << RMSEC << std::endl << std::endl;
	FinalPLS = 1;//������� ���������, ��������� ���������� ���������� � ���� � ��������� ��������� ��������
	SaveResultsPLS();//��������� ���������� � ���������� ���������� LE
}


void Calibr::InitLE_SetSumX()
//� �������� ��������� �������� ��� LE ������������ ������������ �� ����� ���� �������������� � �������
{
	LE = Spectra.rowwise().sum();
	LEInitialType = "SummNorm";//������������ �� ����� ��������������
}

void Calibr::InitLE_SetMeanX()
//� �������� ��������� �������� ��� LE ������������ ������������ �� ������� ������������� � �������
{
	LE = Spectra.rowwise().mean();
	LEInitialType = "MeanNorm";//������������ �� ������� ������������� � �������
}

void Calibr::InitLE_SetMaxX()
//� �������� ��������� �������� ��� LE ������������ ������������ �� ������������� ������������ ����� � �������
{
	LE = Spectra.rowwise().maxCoeff();
	LEInitialType = "MaxNorm";//������������ �� ������������ ������������� � �������
}

void Calibr::InitLE_NonCoherentBackScatter(const std::string LEFileName)
//������������ �� ������������� �������������� ��������� ���������
{
	std::string FileName = SpectraPath;
	LoadVectorLong(FileName.append(LEFileName), LE);//������������� ������� LE ������� �� �����
	LEInitialType = "BckScattrng";//������������� ������� LE ���������� �������������� ��������� ���������
}

void Calibr::SaveResultsPLS()
/*���������� ����������� ���������� �� ���� � ������� �����������

*/
{
	//������� ��� �����
	std::string FileName = "A";
	FileName.append(std::to_string(this->Amax));
	FileName.append(LEInitialType);
	FileName.append("_Iter_");
	FileName.append(std::to_string(TotalPLSRebuildIterat));
	FileName.append("RMSEC");
	FileName.append(std::to_string(RMSEC));

	//��������� ���������� ��� ����� c ������� ������������ �� PLS
	FileName.append(".nPLS");
	bool res;
	std::string FullFileName{ CalibrationDataPath };
	FullFileName.append(FileName);
	FileClass ResultPLS;
	res = ResultPLS.OpenForSave(FullFileName);//��������� ���� ��� ������
	if (!res)
	{
		std::cout << "�� ������� ������� ���� ��� ������ ���������� �����������\n";
		std::system("Pause");
		exit(FILE_SAVING_ERROR);
	}
	ResultPLS.SaveObject(CalibrMethod, CalibrMethodStringLength);//��������� � ���� ��� ������������ ����������
	ResultPLS.SaveObject(LEcoeff);
	ResultPLS.SaveObject(XNormcoeff);
	ResultPLS.SaveObject(TotalPLSRebuildIterat);
	ResultPLS.SaveObject(FinalPLS);
	ResultPLS.Close();
}

ProcessError Calibr::LoadResultsPLS(std::string FileName)
//���������� ���������� �������� �� �����
{
	bool res;
	std::string FullFileName{ CalibrationDataPath };
	FullFileName.append(FileName);
	FileClass ResultPLS;
	res = ResultPLS.OpenForRead(FullFileName);//��������� ���� ��� ������
	if (!res)
	{
		std::cout << "�� ������� ������� ����\n" << FullFileName << "\n��� ������ ������ PLS-�����������";
		std::system("Pause");
		return FILE_READING_ERROR;
	}
	ResultPLS.LoadObject(CalibrMethod, CalibrMethodStringLength);//������ �� ����� ��� ������������ ����������
	if (_stricmp(CalibrMethod, "PLS"))
		//������ �� �����. ���������� �� PLS ������� �� ������
	{
		std::cout << "����:\n" << FullFileName << "\n������������ ��� ���������� ������� " << CalibrMethod << "\n"
			<< "�������� ���� � ����������� �� ������ PLS\n\n";
		std::system("Pause");
	}
	ResultPLS.LoadObject(LEcoeff);
	ResultPLS.LoadObject(XNormcoeff);
	ResultPLS.LoadObject(TotalPLSRebuildIterat);
	ResultPLS.LoadObject(FinalPLS);
	ResultPLS.Close();
	return ProcessError::OK;
}

void Calibr::LoadElvaXSpectrum(const std::string FullFileName, RowVectorXd &X)
//�������� �������������� �� ����� ElvaX � ������-������
{
	ifstream spectrum;
	spectrum.open(FullFileName, ios::binary);
	if (spectrum.is_open())
	{
		spectrum.seekg(0x0100);
		DWORD ChannelCount;
		spectrum.read((char*)&ChannelCount, sizeof(DWORD));
		//�������������� ������-������ ��� ������ ������� � �����. ������ ������� ��������� �� 1, ��� ��� ��������� ����� �� ������������
		X.resize(--ChannelCount);
		DWORD* Arr = new DWORD[ChannelCount];//���������� ������ ��� ������ �������
		spectrum.read((char*)Arr, ChannelCount * sizeof(DWORD));//��������� ������������� � �������
		spectrum.close();
		for (DWORD i = 0; i < ChannelCount; ++i)
			X(i) = Arr[i];
		delete[] Arr;
	}
	else
	{
		std::cout << "������\n" << FullFileName << "\n�� ������� �������\n";
	}
}

std::vector<std::string> Calibr::LoadElvaXSpectra(const std::string &Path, MatrixXd &X)
//�������� ���� �������� �������� � ������� X
//����� �������� ���������� � ������� �� �������� ��� ������������� ������������ ���� ��������, ���� ����������
{
	//���������� ���������� ������ � ��������
	std::vector<std::string> FileList;//������ �� ������� ������
	FileClass f;
	std::string mask = Path;
	mask = CorrectPath(mask);
	mask.append("*.evt");
	f.FindFileList(mask, FileList);

	//� ����� ������� ������� �������������� ��� ���� ��������� ������
	for (int i = 0; i < FileList.size(); ++i)
	{
		RowVectorXd Xrow;
		std::string FullFileName = Path;
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
	return FileList;
}

VectorXd Calibr::SpectraPredictPLS(const MatrixXd &X, const StructPLS &LEcoeff, const StructPLS &NormXcoeff)
//������������ ������������ ������� ����������� �������� ����� LECoeff - ������������ PLS ��� ������ ������������� ����� ������������
// � NormXcoeff - ������������ PLS ��� �������������� �������
{
	VectorXd LE(X.rows());//��������� ������ LE ��� ��������� ������������� ��������������
	PredictPLS(X, LEcoeff, LE);//���������� ������������� �������������� �� PLS
	//�������� ������������ ������� ������������ ������������ NxN
	DiagonalMatrix<double, Dynamic> LEn(LE.array().inverse().matrix());
	MatrixXd Xnorm = LEn * X;//������������� ������� ��������
	//����� �� �������������� ������, �������� ������ ������������� ������������ � ������ LE
	PredictPLS(Xnorm, NormXcoeff, LE);
	return LE;//���������� �� ������� ��������� ���������� ������������� ������������
}

void Calibr::PredictPLS(const MatrixXd &X0, const StructPLS &Coeff, VectorXd &Ycalc)
//������� ������ �������� Y ��� ������� ����������� �������� X �� PLS ���������� Coeff
{
	MatrixXd T(X0.rows(), Coeff.A);//������� ������ NxA (���-�� �������� x ���-�� ��)
	MatrixXd X = X0.rowwise() - Coeff.Xmean;//������������ ������� ������������� ��������, ��� ������� ����������� �������
	for (int i = 0; i < Coeff.A; ++i)
		//� ����� ��������������� ��������� ������� ������
	{
		T.col(i) = X * Coeff.W.col(i);//����� ����� i+1 ����������
		X = X - T.col(i)*Coeff.P.col(i).transpose();//�������� ����� ������� X
	}

	Ycalc = (T * Coeff.Q).array() + Coeff.Ymean;//�������� ������������ ������� ������� ��� ������� X0 �� PLS
}