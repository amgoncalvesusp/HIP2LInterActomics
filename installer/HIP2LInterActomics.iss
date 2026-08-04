; Inno Setup script for HIP2LInterActomics (native Windows installer).
;
; Packages the PyInstaller onedir bundle into a single Setup.exe that installs
; per-user (no admin / no UAC), creates Start Menu and Desktop
; shortcuts, and registers an uninstaller. The GUI is self-contained; running
; LUNA analyses still requires the separate "luna-env" conda environment
; (installable from the app's "1. Inicio" tab).
;
; Build:
;   "%LOCALAPPDATA%\Programs\Inno Setup 6\ISCC.exe" installer\HIP2LInterActomics.iss
; Override the bundle/output dirs from the command line if needed, e.g.:
;   ISCC.exe /DBundleDir="C:\path\to\HIP2LInterActomics" /DOutputDir="C:\out" installer\HIP2LInterActomics.iss

#define MyAppName "HIP2LInterActomics"
#define MyAppVersion "1.4.0"
#define MyAppPublisher "Daniel Andres Grajales Ruiz e Adriano Marques Goncalves"
#define MyAppExe "HIP2LInterActomics.exe"

; Defaults can be overridden with ISCC /D switches.
#ifndef BundleDir
  #define BundleDir "C:\luna_build\dist\HIP2LInterActomics"
#endif
#ifndef OutputDir
  #define OutputDir "C:\luna_build\installer_out"
#endif
#ifndef IconFile
  #define IconFile "..\luna_gui\assets\hip2l_interactomics_icon.ico"
#endif

[Setup]
AppId={{8F3A1C42-9B7D-4E16-AE58-2C1D9F4B7A60}
AppName={#MyAppName}
AppVersion={#MyAppVersion}
AppPublisher={#MyAppPublisher}
DefaultDirName={localappdata}\Programs\{#MyAppName}
DefaultGroupName={#MyAppName}
DisableProgramGroupPage=yes
PrivilegesRequired=lowest
OutputDir={#OutputDir}
OutputBaseFilename=HIP2LInterActomics-Setup
SetupIconFile={#IconFile}
UninstallDisplayIcon={app}\{#MyAppExe}
Compression=lzma2
SolidCompression=yes
WizardStyle=modern
ArchitecturesAllowed=x64compatible
ArchitecturesInstallIn64BitMode=x64compatible
ChangesEnvironment=yes

[Languages]
Name: "brazilianportuguese"; MessagesFile: "compiler:Languages\BrazilianPortuguese.isl"
Name: "english"; MessagesFile: "compiler:Default.isl"
Name: "spanish"; MessagesFile: "compiler:Languages\Spanish.isl"

[Files]
; PyInstaller places environment.yml at the bundle root; this recursive rule
; installs it physically beside the application on Windows.
Source: "{#BundleDir}\*"; DestDir: "{app}"; Flags: recursesubdirs createallsubdirs ignoreversion
Source: "hipplinteractomics-terminal.cmd"; DestDir: "{app}\bin"; Flags: ignoreversion
Source: "hipplinteractomics-multiple-run.cmd"; DestDir: "{app}\bin"; Flags: ignoreversion

[Registry]
Root: HKCU; Subkey: "Environment"; ValueType: expandsz; ValueName: "Path"; ValueData: "{olddata};{app}\bin"; Check: NeedsAddPath(ExpandConstant('{app}\bin'))

[Icons]
Name: "{group}\{#MyAppName}"; Filename: "{app}\{#MyAppExe}"
Name: "{group}\Uninstall {#MyAppName}"; Filename: "{uninstallexe}"
Name: "{autodesktop}\{#MyAppName}"; Filename: "{app}\{#MyAppExe}"

[Run]
Filename: "{app}\{#MyAppExe}"; Description: "{cm:LaunchProgram,{#MyAppName}}"; Flags: nowait postinstall skipifsilent

[Code]
function NeedsAddPath(Param: string): Boolean;
var
  CurrentPath: string;
begin
  if not RegQueryStringValue(HKCU, 'Environment', 'Path', CurrentPath) then
    CurrentPath := '';
  Result := Pos(';' + Uppercase(Param) + ';', ';' + Uppercase(CurrentPath) + ';') = 0;
end;
