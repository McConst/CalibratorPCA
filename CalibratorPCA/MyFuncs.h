#include <string>
#include <Windows.h>

#ifndef MYFUNC_H
#define MYFUNC_H

std::string &CorrectPath(std::string &Path);//������� ��������� ���� � ���� ������, ���� ����� �����������
LPWSTR StringToW_Char(const std::string &Str);//������� �������������� ���������� �������� � ���� w_char
void StringToClipboard(const std::string &str);//������ ���� String ������� � ����� ������
#endif

