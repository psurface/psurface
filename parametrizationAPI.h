/***************************************************************
 *
 * DLL export/import definitions for Windows
 *
 ***************************************************************/
#ifndef PARAMETRIZATION_WIN_DLL_API_H
#define PARAMETRIZATION_WIN_DLL_API_H

#ifdef _WIN32
#  ifdef PARAMETRIZATION_EXPORTS
#     define PARAMETRIZATION_API __declspec(dllexport)
#  else
#     define PARAMETRIZATION_API __declspec(dllimport)
#  endif
#else
#   define PARAMETRIZATION_API 
#endif

#endif
