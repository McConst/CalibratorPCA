#include "pch.h"
#include "Split.h"
#include <string>
#include <iostream>


Split::Split()
{
	//std::cout << "Конструктор Split" << std::endl;
}

Split::Split(const std::string &Strng, const std::string &delimiter) :Str(Strng), Dlmtr(delimiter)
//Деление строки на подстроки с учетом разделителя выполняем непосредственно в конструкторе
{
	
	Splitting();
}

void Split::operator()(const std::string &Strng, const std::string &delimiter)
//Функция деления новой строки на подстроки тем же экземпляром класса
{
	Result.clear();//Очищаем вектор от старых данных
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
//Функция деления строки на подстроки. Работает с закрытыми переменными класса
{
	ubound = -1;//До момента деления на подстроки элементов массива нет
	if (Str.empty() || Dlmtr.empty())//Строка пустая, не выполняем деление на подстроки
	{

		return;
	}
	size_t StartPos = 0;//Начальная позиция подстроки
	auto pos = Str.find(Dlmtr);//Конечная позиция подстроки
	while (pos != std::string::npos)
	{
		ubound++;
		Result.push_back(Str.substr(StartPos, pos - StartPos));//В вектор загоняется string
		StartPos = pos + Dlmtr.size();
		if (StartPos >= Str.size()) break;
		pos = Str.find(Dlmtr, StartPos);
	}
	if (StartPos < Str.size())
		//Прорисовываем строку, оставшуюся после последнего разделителя
	{
		++ubound;
		Result.push_back(Str.substr(StartPos, Str.size() - StartPos));
	}
}

Split::~Split()
{
	//std::cout << "Деструктор Split"<<std::endl;
}
