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
	int SpectraCount;//���������� �������� � �������������� ���������� ������������
	

	//������
	int LoadMatrixLong(const std::string &FileName, MatrixXd &X);//������� ���������� ���������� ������� ������� �� ����� ���� long (������������� ���������������)
	int LoadMatrixDouble(const std::string &FileName);
	int LoadVectorLong(const std::string &FileName, VectorXd &X);//������� ������������� ������� ������� �� �����

public:
	//�������� ������

	struct StructPLS//���������, � ������� �������� ���������� PLS-����������
	{
		int N;//���������� ��������/��������
		int M;//���������� �������
		int A;//���������� ��
		RowVectorXd Xmean;//������� ������
		double Ymean;//������� ������
		Matrix <double, Dynamic, Dynamic> T;//������� ������ ������������ NxA (N-���-�� ��������, A - ���-�� ��)
		Matrix <double, Dynamic, Dynamic> P;//�������� ��� ������� ������� X (NxM - N-���-�� ��������, M-���-�� �������)
		VectorXd Q;//������ �������� ��� ������� ������� Y ������������ A ������� ���������
		Matrix <double, Dynamic, Dynamic> W;//������� ���������� �������� MxA (M-���-�� �������, �-���-�� ��)
		MatrixXd E;//������� �������� (������) ��� X ������������ NxM
		VectorXd F;//������ �������� (������) ��� Y ������������ N
		VectorXd B;//������ ������������� ��������� ������������ A ��� ������������ Y
	};

	int Amax;//������������ ���������� ������� ���������
	Matrix<double, Dynamic, Dynamic> Spectra;//������������ ������� � ������� Eigen 
	VectorXd LE;//������ ��������� �������� LE, ������������������ �� �����
	StructPLS LEcoeff;//������� ��������� ���������� ������� ������������
	StructPLS XNormcoeff;//������� ��������� ���������� ������������� ������� ��������

	Matrix<double, Dynamic, CRM_ElementCount> mY;//������� ������������ � ������� Eigen

	//������
	Calibr(const std::string &InitFileName);//������ ���� � ����� �������������
	~Calibr();

	void LoadSpectra(const std::string &SpectraFileName = "Spectra.dat", const std::string &YFileName="Y.txt", const std::string &LEfilename="LE.dat");
	void SetSumXtoLE();//� �������� ��������� ������������� ������������� ��� �������� ������� LE ����� � �������
	void SetMeanXtoLE();//� �������� ��������� ��� ������������� ��������������� ������� �������� � �������
	void SetMaxXtoLE();//� �������� ��������� ��� ������������� ��������������� ������� �������� � �������

	/*���������� �� PLS. ������� �������� Y0 - ���������������� ������ ������������ ������ ��������
											A - ���������� ������� ��������� ��� ����������
	*/
	void MainCalibrationPLS();//������� ����� �� ���������� � ������������� �� ���������� ��������
	void DecomposePLS(const MatrixXd &X0, const VectorXd &Y0, StructPLS &Result);//���������� ��������� �� ��������, ����� �� �������� ��� ��������� ������
	
	//������� ���������� � ������������� �� ������������� ���������� ��������, ����� �������������� ��� ������� �������������
	void NormDecomposePLS
		(const VectorXd &B_LE, MatrixXd &T, double Ymean, const MatrixXd &X0, const VectorXd &Y0, StructPLS &sX);

	void  ScorePredictPLS(const VectorXd &B, const MatrixXd &T, double Ymean, VectorXd &Y);//PLS ������� ����� ����� � ������������ ���������. Y-���������� ��������
	double RMSE(VectorXd const &Y0, VectorXd const &Ycalc);// ������ ��������� �����������. ������� RMSE - ���������� ����������
	void SaveResultsPLS(const std::string FileName);
	void MatrixToClipboard(MatrixXd X);//����������� int ������� � ����� ������
	void VectorToClipboard(VectorXd X);//����������� int ������� � ����� ������
};

