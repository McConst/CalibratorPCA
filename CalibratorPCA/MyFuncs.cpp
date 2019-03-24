#include "pch.h"
#include <string>
#include <locale>
#include <sstream>
#include <Windows.h>

std::string &CorrectPath(std::string &Path)
//Функция добавляет слеш в конец пути, если он отсутствует
{
	if (Path.substr(Path.length() - 1, 1) != "\\") Path.append("\\");
	return Path;
}

LPWSTR StringToW_Char(const std::string &Str)
//Функция конвертирует строковое значение к типу w_char
{
	//Преобразуем строку FileName от типа string к типу widechar через WinAPI
	int bufferlen = ::MultiByteToWideChar(CP_ACP, 0, Str.c_str(), Str.size(), NULL, 0);//Получаем длину буфера
	if (bufferlen == 0)
	{
		// Где-то ошибка, проверьте GetLastError() and log.
		return 0;//Ошибка. Система не выделила буфер под строку. Возвращаем нулевой указатель
	}

	// Создаем буферный массив, который необходимо будет удалить через delete вне функции
	LPWSTR widestr = new WCHAR[bufferlen + 1];
	::MultiByteToWideChar(CP_ACP, 0, Str.c_str(), Str.size(), widestr, bufferlen);//Заполняем буферный массив строкой типа w_char

	// Добавляем нуль в конце, чтобы получить нуль-строку
	widestr[bufferlen] = 0;
	return widestr;
}

std::string W_charToString(const wchar_t *s, char dfault = '?', const std::locale &loc = std::locale())
//Преобразование wchar_t к string
//символы, не существующие в ASCII будут заменены на ?
//для преобрвзования в string используются настройки локали loc (текущей по умолчанию)
{
	std::ostringstream stm;

	while (*s != L'\0') 
	{
		stm << std::use_facet< std::ctype<wchar_t> >(loc).narrow(*s++, dfault);
	}
	return stm.str();
}

void StringToClipboard(const std::string &str)
{
	const char* output = str.c_str();
	const size_t len = strlen(output) + 1;
	HGLOBAL hMem = GlobalAlloc(GMEM_MOVEABLE, len);
	memcpy(GlobalLock(hMem), output, len);
	GlobalUnlock(hMem);
	OpenClipboard(0);
	EmptyClipboard();
	SetClipboardData(CF_TEXT, hMem);
	CloseClipboard();
}