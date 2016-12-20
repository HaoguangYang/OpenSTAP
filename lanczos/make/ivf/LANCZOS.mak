# Microsoft Developer Studio Generated NMAKE File, Based on LANCZOS.dsp
!IF "$(CFG)" == ""
CFG=LANCZOS - Win32 Debug
!MESSAGE No configuration specified. Defaulting to LANCZOS - Win32 Debug.
!ENDIF 

!IF "$(CFG)" != "LANCZOS - Win32 Release" && "$(CFG)" != "LANCZOS - Win32 Debug"
!MESSAGE Invalid configuration "$(CFG)" specified.
!MESSAGE You can specify a configuration when running NMAKE
!MESSAGE by defining the macro CFG on the command line. For example:
!MESSAGE 
!MESSAGE NMAKE /f "LANCZOS.mak" CFG="LANCZOS - Win32 Debug"
!MESSAGE 
!MESSAGE Possible choices for configuration are:
!MESSAGE 
!MESSAGE "LANCZOS - Win32 Release" (based on "Win32 (x86) Console Application")
!MESSAGE "LANCZOS - Win32 Debug" (based on "Win32 (x86) Console Application")
!MESSAGE 
!ERROR An invalid configuration is specified.
!ENDIF 

!IF "$(OS)" == "Windows_NT"
NULL=
!ELSE 
NULL=nul
!ENDIF 

CPP=cl.exe
F90=df.exe
RSC=rc.exe

!IF  "$(CFG)" == "LANCZOS - Win32 Release"

OUTDIR=.\Release
INTDIR=.\Release
# Begin Custom Macros
OutDir=.\Release
# End Custom Macros

ALL : "$(OUTDIR)\LANCZOS.exe"


CLEAN :
	-@erase "$(INTDIR)\LANCZOS.obj"
	-@erase "$(INTDIR)\ORTHMGS.obj"
	-@erase "$(INTDIR)\SSPACE90.obj"
	-@erase "$(INTDIR)\subsp.obj"
	-@erase "$(OUTDIR)\LANCZOS.exe"

"$(OUTDIR)" :
    if not exist "$(OUTDIR)/$(NULL)" mkdir "$(OUTDIR)"

F90_PROJ=/compile_only /include:"$(INTDIR)\\" /nologo /warn:nofileopt /module:"Release/" /object:"Release/" 
F90_OBJS=.\Release/
CPP_PROJ=/nologo /ML /W3 /GX /O2 /D "WIN32" /D "NDEBUG" /D "_CONSOLE" /D "_MBCS" /Fp"$(INTDIR)\LANCZOS.pch" /YX /Fo"$(INTDIR)\\" /Fd"$(INTDIR)\\" /FD /c 
BSC32=bscmake.exe
BSC32_FLAGS=/nologo /o"$(OUTDIR)\LANCZOS.bsc" 
BSC32_SBRS= \
	
LINK32=link.exe
LINK32_FLAGS=kernel32.lib user32.lib gdi32.lib winspool.lib comdlg32.lib advapi32.lib shell32.lib ole32.lib oleaut32.lib uuid.lib odbc32.lib odbccp32.lib /nologo /subsystem:console /incremental:no /pdb:"$(OUTDIR)\LANCZOS.pdb" /machine:I386 /out:"$(OUTDIR)\LANCZOS.exe" 
LINK32_OBJS= \
	"$(INTDIR)\LANCZOS.obj" \
	"$(INTDIR)\ORTHMGS.obj" \
	"$(INTDIR)\subsp.obj" \
	"$(INTDIR)\SSPACE90.obj"

"$(OUTDIR)\LANCZOS.exe" : "$(OUTDIR)" $(DEF_FILE) $(LINK32_OBJS)
    $(LINK32) @<<
  $(LINK32_FLAGS) $(LINK32_OBJS)
<<

!ELSEIF  "$(CFG)" == "LANCZOS - Win32 Debug"

OUTDIR=.\Debug
INTDIR=.\Debug
# Begin Custom Macros
OutDir=.\Debug
# End Custom Macros

ALL : "$(OUTDIR)\LANCZOS.exe"


CLEAN :
	-@erase "$(INTDIR)\LANCZOS.obj"
	-@erase "$(INTDIR)\ORTHMGS.obj"
	-@erase "$(INTDIR)\SSPACE90.obj"
	-@erase "$(INTDIR)\subsp.obj"
	-@erase "$(OUTDIR)\LANCZOS.exe"
	-@erase "$(OUTDIR)\LANCZOS.ilk"
	-@erase "$(OUTDIR)\LANCZOS.pdb"

"$(OUTDIR)" :
    if not exist "$(OUTDIR)/$(NULL)" mkdir "$(OUTDIR)"

F90_PROJ=/check:bounds /compile_only /debug:full /include:"$(INTDIR)\\" /nologo /warn:argument_checking /warn:nofileopt /module:"Debug/" /object:"Debug/" /pdbfile:"Debug/DF60.PDB" 
F90_OBJS=.\Debug/
CPP_PROJ=/nologo /MLd /W3 /Gm /GX /ZI /Od /D "WIN32" /D "_DEBUG" /D "_CONSOLE" /D "_MBCS" /Fp"$(INTDIR)\LANCZOS.pch" /YX /Fo"$(INTDIR)\\" /Fd"$(INTDIR)\\" /FD /GZ /c 
BSC32=bscmake.exe
BSC32_FLAGS=/nologo /o"$(OUTDIR)\LANCZOS.bsc" 
BSC32_SBRS= \
	
LINK32=link.exe
LINK32_FLAGS=kernel32.lib user32.lib gdi32.lib winspool.lib comdlg32.lib advapi32.lib shell32.lib ole32.lib oleaut32.lib uuid.lib odbc32.lib odbccp32.lib /nologo /subsystem:console /incremental:yes /pdb:"$(OUTDIR)\LANCZOS.pdb" /debug /machine:I386 /out:"$(OUTDIR)\LANCZOS.exe" /pdbtype:sept 
LINK32_OBJS= \
	"$(INTDIR)\LANCZOS.obj" \
	"$(INTDIR)\ORTHMGS.obj" \
	"$(INTDIR)\subsp.obj" \
	"$(INTDIR)\SSPACE90.obj"

"$(OUTDIR)\LANCZOS.exe" : "$(OUTDIR)" $(DEF_FILE) $(LINK32_OBJS)
    $(LINK32) @<<
  $(LINK32_FLAGS) $(LINK32_OBJS)
<<

!ENDIF 

.c{$(INTDIR)}.obj::
   $(CPP) @<<
   $(CPP_PROJ) $< 
<<

.cpp{$(INTDIR)}.obj::
   $(CPP) @<<
   $(CPP_PROJ) $< 
<<

.cxx{$(INTDIR)}.obj::
   $(CPP) @<<
   $(CPP_PROJ) $< 
<<

.c{$(INTDIR)}.sbr::
   $(CPP) @<<
   $(CPP_PROJ) $< 
<<

.cpp{$(INTDIR)}.sbr::
   $(CPP) @<<
   $(CPP_PROJ) $< 
<<

.cxx{$(INTDIR)}.sbr::
   $(CPP) @<<
   $(CPP_PROJ) $< 
<<

.SUFFIXES: .fpp

.for{$(F90_OBJS)}.obj:
   $(F90) $(F90_PROJ) $<  

.f{$(F90_OBJS)}.obj:
   $(F90) $(F90_PROJ) $<  

.f90{$(F90_OBJS)}.obj:
   $(F90) $(F90_PROJ) $<  

.fpp{$(F90_OBJS)}.obj:
   $(F90) $(F90_PROJ) $<  


!IF "$(NO_EXTERNAL_DEPS)" != "1"
!IF EXISTS("LANCZOS.dep")
!INCLUDE "LANCZOS.dep"
!ELSE 
!MESSAGE Warning: cannot find "LANCZOS.dep"
!ENDIF 
!ENDIF 


!IF "$(CFG)" == "LANCZOS - Win32 Release" || "$(CFG)" == "LANCZOS - Win32 Debug"
SOURCE=.\LANCZOS.f90

"$(INTDIR)\LANCZOS.obj" : $(SOURCE) "$(INTDIR)"


SOURCE=.\ORTHMGS.f90

"$(INTDIR)\ORTHMGS.obj" : $(SOURCE) "$(INTDIR)"


SOURCE=.\SSPACE90.f90

"$(INTDIR)\SSPACE90.obj" : $(SOURCE) "$(INTDIR)"


SOURCE=.\subsp.f90

"$(INTDIR)\subsp.obj" : $(SOURCE) "$(INTDIR)"



!ENDIF 

