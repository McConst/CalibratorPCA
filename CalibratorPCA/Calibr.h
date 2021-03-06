#pragma once

#include <string>
#include <vector>
#include "const.h"
#include "Eigen/Dense"

using namespace Eigen;




class Calibr
{
	//�������� ������
	std::string InitFile;// ������ ���� � ����� �������������
	std::string SpectraPath;//���� � ������ �������� � �����������;
	std::string CalibrationDataPath;//���� � ������ � ��������� ����������
	char SpectaLoadingMethod;//������ �������� �������� �� �����
	int SpectraCount;//����������  �������������� �������� � �������������� ���������� ������������

	

	//������
	void LoadY(const std::string &YFileName, std::vector<std::vector<double> > &Y0);//�������� ������������� ������������ �� ����� � ��������� ������
	void LoadNames(std::vector<std::string> &Names, const std::string &NamesFileName = "Names.txt");//�������� ������� ���� ����������� �������� ��� ����������� ����������� ���/������������

	int LoadMatrixLong(const std::string &FileName, MatrixXd &X);//������� ���������� ���������� ������� ������� �� ����� ���� long (������������� ���������������)
	int LoadMatrixDouble(const std::string &FileName);
	int LoadVectorLong(const std::string &FileName, VectorXd &X);//������� ������������� ������� ������� �� �����

public:
	//�������� ������

	struct StructPLS//���������, � ������� �������� ���������� PLS-����������
	{
		Index N;//���������� ��������/��������
		Index M;//���������� �������
		Index A;//���������� ��
		RowVectorXd Xmean;//������� ������
		double Ymean;//������� ������
		Matrix <double, Dynamic, Dynamic> T;//������� ������ ������������ NxA (N-���-�� ��������, A - ���-�� ��)
		Matrix <double, Dynamic, Dynamic> P;//�������� ��� ������� ������� X (MxA - N-���-�� ��������, M-���-�� �������)
		VectorXd Q;//������ �������� ��� ������� ������� Y ������������ A ������� ���������
		Matrix <double, Dynamic, Dynamic> W;//������� ���������� �������� MxA (M-���-�� �������, �-���-�� ��)
		MatrixXd E;//������� �������� (������) ��� X ������������ NxM
		VectorXd F;//������ �������� (������) ��� Y ������������ N
	};

	std::string WorkingPath;//���� � �������� � ������� �������� ��� �������
	int Amax;//������������ ���������� ������� ���������
	Matrix<double, Dynamic, Dynamic> Spectra;//������������ ������� ��� ���������� (��� ��� ��������) � ������� Eigen 
	VectorXd LE;//������ ��������� �������� LE, ������������������ �� �����
	StructPLS LEcoeff;//������� ��������� ���������� ������� ������������
	StructPLS XNormcoeff;//������� ��������� ���������� ������������� ������� ��������
	Matrix<double, Dynamic, CRM_ElementCount> mY;//������� ������������ � ������� Eigen
	std::string LEInitialType;//��� �������� ������ ������������� LE ��� ��������. ���������� ����� ���������� � ��� �����
	char CalibrMethod[CalibrMethodStringLength];//����� ����������: PLS, PCR. 3 ������� ��� �����
	double RMSEC;//�������� ����������
	int TotalPLSRebuildIterat;//����� ���������� �������� PLS
	std::string CalcParametersFile;//��� ����� � ������������ ���������� (��� ������������� ����������)
	char FinalPLS;//���� ��������� �������� ������� PLS
	char FinalPCR;//���� ��������� �������� ������� PCR


	//������
	Calibr(const std::string &InitFileName);//������ ���� � ����� �������������
	~Calibr();

	void LoadInitDataForCalibrat(const std::string &SpectraFileName = "Spectra.dat", 
		const std::string &YFileName="Y.txt", 
		const std::string &NamesFile="Names.txt");

	void InitLE_SetSumX();//� �������� ��������� ������������� ������������� ��� �������� ������� LE ����� � �������
	void InitLE_SetMeanX();//� �������� ��������� ��� ������������� ��������������� ������� �������� � �������
	void InitLE_SetMaxX();//� �������� ��������� ��� ������������� ��������������� ������� �������� � �������
	void InitLE_NonCoherentBackScatter(const std::string LEFileName="LE.dat");//������������� �� ������� �������������� ��������� ���������

	/*���������� �� PLS. ������� �������� Y0 - ���������������� ������ ������������ ������ ��������
											A - ���������� ������� ��������� ��� ����������
	*/
	void MainCalibrationPLS(int elmnt=0);//������� ����� �� ���������� � ������������� �� ���������� ��������
	void DecomposePLS(const MatrixXd &X0, const VectorXd &Y0, StructPLS &Result);//���������� ��������� �� ��������, ����� �� �������� ��� ��������� ������
	
	//������� ���������� � ������������� �� ������������� ���������� ��������, ����� �������������� ��� ������� �������������
	void NormDecomposePLS
		(const VectorXd &B_LE, MatrixXd &T, double Ymean, const MatrixXd &X0, const VectorXd &Y0, StructPLS &sX);
	
	void  ScorePredictPLS(const VectorXd &B, const MatrixXd &T, double Ymean, VectorXd &Y);//PLS ������� ����� ����� � ������������ ���������. Y-���������� ��������
	void PredictPLS(const MatrixXd &X, const StructPLS &Coeff, VectorXd &Ycalc);//������� ������ �������� Ycalc ��� ������������ ������� X �� PLS ���������� Coeff
	VectorXd SpectraPredictPLS(const MatrixXd &X, const StructPLS &LEcoeff, const StructPLS &NormXcoeff);//������������ ������������ ���������. ������� ����� �����. PLS �� � ����� ������������
	double RMSE(VectorXd const &Y0, VectorXd const &Ycalc);// ������ ��������� �����������. ������� RMSE - ���������� ����������
	
	
	void SaveResultsPLS();//���������� ����������� ������ ������� PLS
	ProcessError LoadResultsPLS(const std::string FileName);//������ ���������� PLS �� ����� � ������ ������
	std::vector<std::string> LoadElvaXSpectra(const std::string &Path, MatrixXd &X);//�������� ���� �������� �������� � ������� X
	void LoadElvaXSpectrum(const std::string FullFileName, RowVectorXd &X);//�������� ������������ �������������� �� ����� ElvaX � ������-������


	void MatrixToClipboard(MatrixXd X);//����������� int ������� � ����� ������
	void VectorToClipboard(VectorXd X);//����������� int ������� � ����� ������
};

