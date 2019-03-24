#include <string>
#include <Windows.h>

#ifndef MYFUNC_H
#define MYFUNC_H

std::string &CorrectPath(std::string &Path);//Функция добавляет слеш к пути справа, если такой отсутствует
LPWSTR StringToW_Char(const std::string &Str);//Функция преобразования строкового значения к типу w_char
std::string W_charToString(const wchar_t *s, char dfault = '?', const std::locale &loc = std::locale());//Преобразование строки wchar_t к string
void StringToClipboard(const std::string &str);//Строку типа String помщаем в буфер обмена
#endif

