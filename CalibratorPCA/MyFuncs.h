#include <string>
#include <Windows.h>

#ifndef MYFUNC_H
#define MYFUNC_H

std::string &CorrectPath(std::string &Path);//������� ��������� ���� � ���� ������, ���� ����� �����������
LPWSTR StringToW_Char(const std::string &Str);//������� �������������� ���������� �������� � ���� w_char
std::string W_charToString(const wchar_t *s, char dfault = '?', const std::locale &loc = std::locale());//�������������� ������ wchar_t � string
void StringToClipboard(const std::string &str);//������ ���� String ������� � ����� ������
#endif

