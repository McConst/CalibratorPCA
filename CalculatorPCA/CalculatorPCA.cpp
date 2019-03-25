#include "pch.h"
#include <iostream>
#include <string>
#include <locale.h> /* ��� �������� ����� */
#include "const.h"
#include "Calibr.h"
#include "Split.h"
#include "Eigen/Dense"

using namespace std;
using namespace Eigen;

int main(int argc, char* argv[])
{
	//� argv ��������� ����, ������� ����� �������������� ��� ����������� ����� ������������� ��� ������ Calibr
	
	setlocale(LC_ALL, "Rus");//��������� ������� ���������
	//setlocale(LC_NUMERIC, "C");//������������ ������� ����� ������ �����

	string Path = argv[0];
	string::size_type pos = Path.rfind("\\") + 1;
	Path = Path.substr(0, pos); //�������� �������, � ������� ��������� ���� �������������
	Calibr Calc (Path.append(CalculatorInitFile));//�������������� ����� ������ Calc ������� �� �����
	Calc.LoadInitDataForCalibrat();//������ ������������ ������ � ������������� �������� �� ���

}