#pragma once
//����� ��� ������� ������ string �� ��������� �� ����������� delimiter

#include <string>
#include <vector>

class Split
{
	std::string Str;//�������� ������, ������� ���������� ���������
	std::string Dlmtr;//������-�����������;
	int ubound;//������� ������� �������, ���������� ����� ���������� ������
	std::vector<std::string> Result;//������, ������ ��������� ������� ������ �� �������
	void Splitting();//������� ������� ������ �� ���������.

public:
	Split();
	Split(const std::string &Strng, const std::string &delimiter);//������� ������ ����� ��������� ������������

	void operator()(const std::string &Strng, const std::string &delimiter);//������� ����� ������ ��� �� ����������� ������

	const std::string& operator[] (const int index);//���������� ���������� ������. ���������� ��������� �������
	const int Ubound() { return ubound; }
	
	~Split();
};

