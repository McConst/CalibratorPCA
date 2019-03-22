#pragma once
//Класс для деления строки string на подстроки по разделителю delimiter

#include <string>
#include <vector>

class Split
{
	std::string Str;//Исходная строка, которую необходимо разделить
	std::string Dlmtr;//Строка-разделитель;
	int ubound;//Верхняя граница массива, получаемая после разделения строки
	std::vector<std::string> Result;//Вектор, хранит результат деления строки на участки
	void Splitting();//Функция деления строки на подстроки.

public:
	Split();
	Split(const std::string &Strng, const std::string &delimiter);//Деление строки через параметры конструктора

	void operator()(const std::string &Strng, const std::string &delimiter);//Деление новой строки тем же экземпляром класса

	const std::string& operator[] (const int index);//Перегрузка квадратных скобок. Возвращает результат деления
	const int Ubound() { return ubound; }
	
	~Split();
};

