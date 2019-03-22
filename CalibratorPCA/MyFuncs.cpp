#include "pch.h"
#include <string>
#include <Windows.h>

std::string &CorrectPath(std::string &Path)
//������� ��������� ���� � ����� ����, ���� �� �����������
{
	if (Path.substr(Path.length() - 1, 1) != "\\") Path.append("\\");
	return Path;
}

LPWSTR StringToW_Char(const std::string &Str)
//������� ������������ ��������� �������� � ���� w_char
{
	//����������� ������ FileName �� ���� string � ���� widechar ����� WinAPI
	int bufferlen = ::MultiByteToWideChar(CP_ACP, 0, Str.c_str(), Str.size(), NULL, 0);//�������� ����� ������
	if (bufferlen == 0)
	{
		// ���-�� ������, ��������� GetLastError() and log.
		return 0;//������. ������� �� �������� ����� ��� ������. ���������� ������� ���������
	}

	// ������� �������� ������, ������� ���������� ����� ������� ����� delete ��� �������
	LPWSTR widestr = new WCHAR[bufferlen + 1];
	::MultiByteToWideChar(CP_ACP, 0, Str.c_str(), Str.size(), widestr, bufferlen);//��������� �������� ������ ������� ���� w_char

	// ��������� ���� � �����, ����� �������� ����-������
	widestr[bufferlen] = 0;
	return widestr;
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