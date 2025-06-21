REM @echo off

set reflector="\graphQL\src\mysqli\App\Reflector.exe"
set R_src="../CellRender\MySql\cad_lab"

%reflector% --reflects /sql ./cad_lab.sql -o %R_src% /namespace biocad_labModel --language visualbasic /split /auto_increment.disable

REM pause