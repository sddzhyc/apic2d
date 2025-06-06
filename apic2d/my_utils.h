// 需要包含头文件
#include <string>
#ifdef _WIN32
#include <commdlg.h>
#include <windows.h>
#endif
#ifdef _WIN32
#include <shlobj.h>
// 检查目录是否存在，不存在则创建
inline void EnsureDirectoryExists(const std::string& path) {
#ifdef _WIN32
  DWORD ftyp = GetFileAttributesA(path.c_str());
  if (ftyp == INVALID_FILE_ATTRIBUTES || !(ftyp & FILE_ATTRIBUTE_DIRECTORY)) {
    CreateDirectoryA(path.c_str(), NULL);
  }
#else
  struct stat st;
  if (stat(path.c_str(), &st) != 0) {
    mkdir(path.c_str(), 0755);
  }
#endif
}


// 弹出选择文件夹窗口，返回用户选择的绝对路径，取消则返回空字符串
inline std::string SelectFolderDialog(const char* title = "Please select folder") {
  char szPath[MAX_PATH] = {0};
  BROWSEINFOA bi = {0};
  bi.lpszTitle = title;
  bi.ulFlags = BIF_RETURNONLYFSDIRS | BIF_NEWDIALOGSTYLE;
  LPITEMIDLIST pidl = SHBrowseForFolderA(&bi);
  if (pidl != nullptr) {
    SHGetPathFromIDListA(pidl, szPath);
    CoTaskMemFree(pidl);
    return std::string(szPath);
  }
  return std::string();
}
//inline std::string SelectFolderDialog(const wchar_t* title = L"请选择文件夹") {
//  wchar_t szPath[MAX_PATH] = {0};
//  BROWSEINFOW bi = {0};
//  bi.lpszTitle = title;
//  bi.ulFlags = BIF_RETURNONLYFSDIRS | BIF_NEWDIALOGSTYLE;
//  LPITEMIDLIST pidl = SHBrowseForFolderW(&bi);
//  if (pidl != nullptr) {
//    SHGetPathFromIDListW(pidl, szPath);
//    CoTaskMemFree(pidl);
//    // 转为UTF-8
//    int len = WideCharToMultiByte(CP_UTF8, 0, szPath, -1, nullptr, 0, nullptr, nullptr);
//    std::string result(len, 0);
//    WideCharToMultiByte(CP_UTF8, 0, szPath, -1, &result[0], len, nullptr, nullptr);
//    // 去除末尾的\0
//    if (!result.empty() && result.back() == '\0') result.pop_back();
//    return result;
//  }
//  return std::string();
//}
#endif