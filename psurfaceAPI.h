/***************************************************************
 *
 * DLL export/import definitions for Windows
 *
 ***************************************************************/
#ifndef PSURFACE_WIN_DLL_API_H
#define PSURFACE_WIN_DLL_API_H

#ifdef HX_OS_WIN
    #ifdef PSURFACE_EXPORTS
        #define PSURFACE_API __declspec(dllexport)
    #else
        #define PSURFACE_API __declspec(dllimport)
    #endif

    #define PSURFACE_EXPORT __declspec(dllexport)
#else
    #define PSURFACE_API
    #define PSURFACE_EXPORT
#endif

#endif
