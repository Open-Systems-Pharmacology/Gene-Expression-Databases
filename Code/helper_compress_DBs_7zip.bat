@echo off
REM Change to PK-Sim DBs directory
cd "../PK-Sim DBs"

REM Find all .expressionDB files and compress each to .7z
for /R %%F in (*ADME_ONLY*.expressionDB) do (
    REM Compress using 7z (overwrite if exists)
    7z a "%%F.7z" "%%F"
)