/***************************************************************
 *
 * DLL export/import definitions for Windows
 *
 ***************************************************************/
#ifndef PSURFACE_WIN_DLL_API_H
#define PSURFACE_WIN_DLL_API_H

#ifdef _WIN32
#  ifdef PSURFACE_EXPORTS
#     define PSURFACE_API __declspec(dllexport)
#  else
#     define PSURFACE_API __declspec(dllimport)
#  endif
#else
#   define PSURFACE_API 
#endif

#endif
