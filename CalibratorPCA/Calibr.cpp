
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
		}
	}
	else
	{
		std::cout << "���� ������������� �� ������ ��� �� ��������" << std::endl;
	}
	file.close();//��������� ����
}



Calibr::~Calibr()
{
	//std::cout << "���������� Calibr"<<std::endl;
}

void Calibr::LoadSpectra(const std::string &SpectraFileName, const std::string &YFileName, const std::string &LEfilename)
/*	������� ��������� ������������� ������� �������� � ������������ �� ���������� �� ������
SpectraFileNmae - ���� ������� ��������
YFileName - ���� � �������������� ���������� ������������, ��������������� ������� ������� �� �������
*/
{
	//�������������� ��������� ������������ ��������� ������ ��� ������ ������������ �� ������
	std::vector<std::vector<double> > Y0; // ��������� ��������� ������ ������������� ������������. ������ �������. - ��, ����. �������. - �������


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
		std::cout << "���� ������������� �� ������ ��� �� ��������" << std::endl;
	}
	file.close();//��������� ����

	//������������� ������������ �� ���������� ������������� ������� � ������� Eigen
	mY.resize(SpectraCount, CRM_ElementCount);
	for(int i=0;i<SpectraCount;++i)
		for (int j = 0; j < CRM_ElementCount; ++j)
		{
			mY(i, j) = Y0[i][j];
		}
	
	//�������������� ������������ ������ �� �������
	FileName = SpectraPath;
	LoadMatrixLong(FileName.append(SpectraFileName), Spectra);//������������� ������� Spectra ������� �� �����
	FileName = SpectraPath;
	LoadVectorLong(FileName.append(LEfilename), LE);//������������� ������� LE ������� �� �����
	//LoadMatrixDouble(FileName.append(SpectraFileName));//�������� ������� � ��������������� ���� double � �������

}



int Calibr::LoadMatrixLong(const std::string &FileName, MatrixXd &X)
//������������� ��������� ������� X ������� �� �����
{
	int* Arr = new int [SpectraCount*ChannelCount];//������� ������ ��� ������ ������������� ��������

	// ����������� ������ ���� String � w_string
	LPWSTR widestrFileName = StringToW_Char(FileName);

	HANDLE hFile;//��������� ���������� API
	hFile = CreateFile(widestrFileName, GENERIC_READ, 0, 0, OPEN_EXISTING, FILE_ATTRIBUTE_NORMAL, 0);
	delete[] widestrFileName;

	int DataSize = SpectraCount * ChannelCount * sizeof(int);//������ ���� ������������� �� ������ �� �����
	DWORD nBytesRead;//������, ���������� ������� �� �����
	int res = ReadFile(hFile, Arr, DataSize, &nBytesRead, 0);//� ����� ������ Spectra ���������� �����������!!!!
	CloseHandle(hFile);
	X.resize(SpectraCount,ChannelCount);
	for (int i=0;i<SpectraCount*ChannelCount;++i)
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
	MatrixXd X = X0.rowwise()-Result.Xmean;//�������������� ������� ��������
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
		Result.P.col(A)=Xt*Result.T.col(A)/TtT;
		Result.Q.row(A) = (Y.transpose()*Result.T.col(A))/TtT;
		Result.E = X - Result.T.col(A)*Result.P.col(A).transpose();
		Result.F = Y - Result.T.col(A)*Result.Q(A);
		//�������� ���������� ��� ��������� ��������
		Y = Result.F;
		X = Result.E;
	}
	//��������� ������������ ���������
	MatrixXd Tt=Result.T.transpose();
	Result.B = (Tt*Result.T).inverse()*Tt*Y1;
	return;
}

inline void  Calibr::ScorePredictPLS(const VectorXd &B, const MatrixXd &T, double Ymean, VectorXd &Y)
//������� ������� ����������� �������� Y ����� ����� T � ������������ ��������� B
{
	Y = (T * B).array() + Ymean;//������ ���������� ��������
}


void Calibr::MatrixToClipboard(MatrixXd X)
{
	int Cols = X.cols();
	int Rows = X.rows();
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
	int Cols = X.cols();
	int Rows = X.rows();
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
//����� ������� RMSE
//����������� RMSE - ���������� �������� �����������
//Y0 - ������ �������� ��������
//Ycalc - ������ ��������� ��������
{
	return std::sqrt((Y0 - Ycalc).array().square().sum());
}

void Calibr::NormDecomposePLS
	(const VectorXd &B_LE, MatrixXd &T, double Ymean, const MatrixXd &X0, const VectorXd &Y0, StructPLS &sX)
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
	ScorePredictPLS(B_LE, T, Ymean, LEp);//�������� ������������� ������� ������������
	//�������� ������������ ������� ������������ ������������ NxN
	DiagonalMatrix<double, Dynamic> LEn(LEp.array().inverse().matrix());
	MatrixXd Xn = LEn * X0;//������������� ������� ��������
	DecomposePLS(Xn, Y0, sX);//��������� PLS ���������� ����� ������������
}

void Calibr::MainCalibrationPLS()
//������� ����� ��� PLS ���������� � ������������� �� ���������� ��������
{
	int N = Spectra.rows();//���������� ��������
	VectorXd Conc0;//������ ��������� ������������ ��� ���������� ��������
	MatrixXd Jac(N, Amax);//�������
	VectorXd NewLE(N);//������ ����� �������� LE
	StructPLS NewXCoeff;//���������� ������ ���������� ������������� ��������
	VectorXd NewConc;//����� ������������ ��� ���������� ���������� �����������������
	VectorXd dB;//���������� ������������� � �������� ������.-�����.
	MatrixXd Jt;//����������������� �������
	MatrixXd JtJ;//���������� ������� AxA (�������). ���� ������� ����� ��� �������� � �������������� �������� � ����� ������� �������� ������
	VectorXd F(N);//������-������� �������. � ������ ������ �������� ����� �������� � ��������� ��������� (�������������� ��� �������)
	double iteratErr;//������ ���������� �������������. ������� ������ �� ��������� ����������-����������
	StructPLS OldPLS;//������ ��������� ���������� �� PLS;
	double F_errLimit;//���������� ���������� ������ ���������� ������� �������
	double F_err;//�������� ������ ������� �������

	//for (
	int elmnt = 0;
		//; elmnt < CRM_ElementCount; ++elmnt)
	//��������� ���������������� ����������� ��� ���� ������������� ���������
	//{
		int TotalPLSRebuildIterat{ 0 };//������� �������� ��������� PLS scores ��� ������ LE
		do
			//���� ����������� ���������� PLS ��� ������ ������������� ������������ ����� ������� X
		{
			double lambda{ StartLambda };//�������������� ������ ��� ��������� �������� ������� ��������� ������� ����������-����������
			int TotalIterat{ 0 };
			iteratErr = LEcoeff.B.norm()*err_relatCompareResult;
			do
			{
				//������ ���������� ��� ������� �������� ������������� ������������. �����. ������������ ��������� � �������� ������ XNormcoeff
				NormDecomposePLS(LEcoeff.B, LEcoeff.T, LEcoeff.Ymean, Spectra, mY.col(elmnt), XNormcoeff);
				//�� ���������������� ������ ������������� ��������� ������������ � �onc0
				ScorePredictPLS(XNormcoeff.B, XNormcoeff.T, XNormcoeff.Ymean, Conc0);
				
				omp_set_nested(1);//��������� ��������� �����������
#pragma omp parallel for ordered private(NewXCoeff, NewConc)
				for (int i = 0; i < Amax; ++i)
					//� ����� ��������������� ������������ ������� ��������
				{
#pragma omp ordered
					{
						VectorXd B_LE2 = LEcoeff.B;//��������������� B_LE2 ��� ����� �������� �����
						double DiffStep = B_LE2(i)*err_relatDiffStep;//���������� ��� ���������� ����������������� ��� ������� ��������

						//����� ���������, ��� ����������� ��������� ������������� ��� ������
						//����� �� ������ �� ���� �����
						B_LE2(i) += DiffStep;//������ ��������� ����������� �� �������� ���������� ������� �����������������

						//��������� ������ ������������� ���������� ��� ����� ������������ B_LE2
						NormDecomposePLS(B_LE2, LEcoeff.T, LEcoeff.Ymean, Spectra, mY.col(elmnt), NewXCoeff);

						//��������� ������ ����� ���������� �������� �������� ����� ������ ���������� � ������ ����������
						ScorePredictPLS(NewXCoeff.B, NewXCoeff.T, NewXCoeff.Ymean, NewConc);
						VectorXd dYdX = (NewConc - Conc0) / DiffStep;//��������� �����������������. i-� ������� ��������
						Jac.col(i) = dYdX;
					}
				}

				F = mY.col(elmnt)- Conc0;//������ ������-������� �������
				F_errLimit = F.norm()/std::sqrt(N)*err_relatCompareResult;//���������� ���������� ������ ���������� ������� �������

				Jt = Jac.transpose();
				JtJ = Jt * Jac;
				int iteratCount = 0;
				bool DoIterat;//���� ����������� �������� ��� ������� ������, ���� ��� ������� ����������
				VectorXd NewB; //������������ ���������, ������� ����� ����������� �� ������ �� �����
				do
					//���� ������� �������� ������ ��� ������������� ���� ����������
				{
					//����������� � ������������ ����� ��������� ����� ��� ���������� ������.
					lambda /= 2;//�������������� ��������� ������ �� ������, ���� ������� ����� ��������� ��� ������� ����
					int MaxThreads = omp_get_max_threads();//�������� ������������ ���������� ������� ��� ������ �����������
					VectorXd minF2=F;//����������� ������ F2 ����� �����������

					
					//���� � ������������ ������
					#pragma omp parallel 
					{
						//�� ���� ������� �������������� ��������� ����������
						VectorXd NewB_local;//������ �������������, �������������� ����������� � ����� for
						VectorXd db_local;//������ ���������� ������������� � ������������ ����� for
						VectorXd F2(N);
						double lambda_local;//����������� ������ � ������������ ������

						#pragma omp	for ordered private(NewXCoeff, NewConc)
						for (int i = 0; i < MaxThreads; ++i)
						#pragma omp ordered
						{
							lambda_local = lambda * pow(2, i);//� ������� ������ ���������� ���� �� ����� ������ 2
							//C������ ������������ �������  NxN �� ���������� ������ � ���������
							MatrixXd lambdI = MatrixXd::Zero(Amax, Amax);
							lambdI.diagonal().array() = lambda_local;//������� ��������� ��������� �������� ������
							db_local = (JtJ + lambdI).inverse()*Jt*F; //��������� ���������� ��� ������� �������������
							NewB_local = LEcoeff.B + db_local;//������������ ������������ B � ������ ����������
							//����� �������� ������������� PLS-���������� ��� �������� ���������� ������� � �������� ��� ����� �����. B
							NormDecomposePLS(NewB_local, LEcoeff.T, LEcoeff.Ymean, Spectra, mY.col(elmnt), NewXCoeff);
							//��������� ������ ����� ���������� �������� ��� ����� �������������
							ScorePredictPLS(NewXCoeff.B, NewXCoeff.T, NewXCoeff.Ymean, NewConc);
							F2 = mY.col(elmnt) - NewConc;//������� ������� ��� ����� �������������	
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
								NewB = NewB_local;
								lambda = lambda_local;
								dB = db_local;
							}
						}
					}

					//����� �� ������������ ������

					F_err = std::abs((F.norm() - minF2.norm()) / std::sqrt(N));
					if (minF2 == F)//����� ����� - ������ �� ����� ��������� ��������� �������
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

					/*if ((TotalPLSRebuildIterat == 18) && (TotalIterat == 20))
					{
						std::cout << "stop\n";
					}*/
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
				LEcoeff.B = NewB;
				++TotalIterat;
				std::cout << "\r" << "������� " << elmnt+1 << ", ����������� PLS: " << TotalPLSRebuildIterat << ", ���: " << TotalIterat;
				std::cout << ", ������ dB � LevMaq: " << dB.norm() / std::sqrt(Amax) << ", RMSEC= " << F.norm() / std::sqrt(N)<<"   ";
			/*������� ������ �� �����:
			1. ��������� ���������� ������������ ��������
			2. ������ �� ������������� �� ��������� ������ ����������
			3. ��������� ������� ������� �� ��������� ������ ����������
			*/
			} while ((TotalIterat < MaxIterationsLevMaq) && (dB.norm()/std::sqrt(Amax) > iteratErr) && (F_err > F_errLimit));
			
			//��������������� ��������� PLS ��� ������ LE
			++TotalPLSRebuildIterat;
			ScorePredictPLS(LEcoeff.B, LEcoeff.T, LEcoeff.Ymean, LE);//�������� ����� �������� LE ��� ��������������
			//��������� ���������� �������� �������������, ����� ����� �� ����� ����� �������� ��������
			OldPLS = LEcoeff;
			DecomposePLS(Spectra, LE, LEcoeff);//��������� ��������������
			std::cout << std::endl << "������ ������������� B_LE ����� ����������: " << (OldPLS.B - LEcoeff.B).norm()/std::sqrt(Amax) << std::endl;
		} while ((TotalPLSRebuildIterat < MaxPLSRebuildIterat) && ((OldPLS.B - LEcoeff.B).norm() /std::sqrt(Amax) > iteratErr));

		std::cout << std::endl << std::endl << "������������ LE ��� "<< elmnt+1 << "-�� ��������:" << std::endl << LEcoeff.B << std::endl << std::endl;
		std::cout << "RMSEC= " << F.norm()/std::sqrt(N) << std::endl << std::endl;
	//}
	std::cout << std::endl << std::endl << "������ ��������" << std::endl << LEcoeff.B << std::endl << std::endl;
}


void Calibr::SetSumXtoLE()
{
	LE = Spectra.rowwise().sum();
	//std::cout << LE;
}

void Calibr::SetMeanXtoLE()
{
	LE = Spectra.rowwise().mean();
	//std::cout << LE;
}

void Calibr::SetMaxXtoLE()
{
	LE = Spectra.rowwise().maxCoeff();
	//std::cout << LE;
}

void Calibr::SaveResultsPLS(const std::string FileName)
/*���������� ����������� ���������� �� ���� � ������� �����������
	
*/
{
	bool res;
	std::string FullFileName{ SpectraPath };
	FullFileName.append(FileName);
	FileClass ResultPLS;
	res=ResultPLS.OpenForSave(FullFileName);//��������� ���� ��� ������
	if (!res)
	{
		std::cout << "�� ������� ������� ���� ��� ������ ���������� �����������\n";
	}
	ResultPLS.SaveObject(LEcoeff);
	ResultPLS.SaveObject(XNormcoeff);
	ResultPLS.Close();
}