#include "pch.h"
#include "Split.h"
#include <string>
#include <iostream>


Split::Split()
{
	//std::cout << "����������� Split" << std::endl;
}

Split::Split(const std::string &Strng, const std::string &delimiter) :Str(Strng), Dlmtr(delimiter)
//������� ������ �� ��������� � ������ ����������� ��������� ��������������� � ������������
{
	
	Splitting();
}

void Split::operator()(const std::string &Strng, const std::string &delimiter)
//������� ������� ����� ������ �� ��������� ��� �� ����������� ������
{
	Result.clear();//������� ������ �� ������ ������
	Str = Strng;
	Dlmtr = delimiter;
	Splitting();
}

const std::string& Split::operator[] (const int index)
{
	if (index > ubound)
	{
		return "\0";
	}
	return Result [index];
}

void Split::Splitting()
//������� ������� ������ �� ���������. �������� � ��������� ����������� ������
{
	ubound = -1;//�� ������� ������� �� ��������� ��������� ������� ���
	if (Str.empty() || Dlmtr.empty())//������ ������, �� ��������� ������� �� ���������
	{

		return;
	}
	size_t StartPos = 0;//��������� ������� ���������
	auto pos = Str.find(Dlmtr);//�������� ������� ���������
	while (pos != std::string::npos)
	{
		ubound++;
		Result.push_back(Str.substr(StartPos, pos - StartPos));//� ������ ���������� string
		StartPos = pos + Dlmtr.size();
		if (StartPos >= Str.size()) break;
		pos = Str.find(Dlmtr, StartPos);
	}
	if (StartPos < Str.size())
		//������������� ������, ���������� ����� ���������� �����������
	{
		++ubound;
		Result.push_back(Str.substr(StartPos, Str.size() - StartPos));
	}
}

Split::~Split()
{
	//std::cout << "���������� Split"<<std::endl;
}
