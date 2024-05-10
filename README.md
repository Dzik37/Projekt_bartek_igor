Projekt 1 - Transformacje

1.Program służy do przeprowadzania transformacji współrzędnych między różnymi układami na kilku dostępnych elipsoidach.

2.Wymagania techniczne:
-system operacyjny Windows 10 lub Windows 11
-python 3.11.5
-zainstalowana biblioteka numpy

3.Transformacje zawarte w programie:
- XYZ -> BLH
- BLH -> XYZ
- XYZ -> NEUp
- BL -> xy2000
- BL -> xy1992

4.Obsługiwane elipsoidy:
- WGS84
- GRS80
- Krasowski (ta elipsoida nie działa prawidłowo dla funkcji: BL -> xy2000 oraz BL -> xy1992)

5.Obsługa programu:
W folderze z programem należy uruchomić konsole.
Aby program zadziałał poprawnie w konsole należy wpisać następującą komende:
python transformacje.py [nazwa transformacji] [nazwa pliku wejściowego] [nazwa elipsoidy] [wsp. X0 Y0 Z0 (tylko w przypadku transformacji XYZ -> NEUp)]

W przypadku kiedy program zadziała poprawnie w konsoli powinien wyświetlić się komunikat, który będzie zawierał 
komendy, które zostały użyte podczas wywoływania programu i nie powinien pojawić się żadnej błąd.

nazwy transformacji:
--xyz2plh 
--plh2xyz
--xyz2neu
--uklad2000
--uklad1992

nazwy elipsoid:
--wgs84
--grs80
--krasowski

Przykładowe wywołanie programu:
python transformacje.py --xyz2plh wsp_xyz_inp.txt --wgs84

Przykładowe wywołanie programu dla transformacji XYZ -> NEUp:
python transformacje.py --xyz2neu wsp_xyz_inp.txt --wgs84 3664940.500 1409153.590 5009571.170

6.Pliki wejściowe:
Plik z danymi powinien mieć format txt oraz znajdować się w tym samym folderze na dysku co program. 
Plik nie powinien mieć żadnego nagłówka. Dane w pliku powinny być współrzędnymi X,Y,Z lub B,L,H oddzielonymi przecinkiem 
a separatorem dziesiętnym powinna być kropka. W każdym wierszu powinny znajdować się współrzędne tylko jednego punktu.
Przykład pliku wejściowego z danymi X,Y,Z:

3664940.500,1409153.590,5009571.170
3664940.510,1409153.580,5009571.167
3664940.520,1409153.570,5009571.167
3664940.530,1409153.560,5009571.168

lub B,L,H:

52.09727,21.03153,141.399 
52.09727,21.03153,141.400 
52.09727,21.03153,141.403 
52.09727,21.03153,141.408 

Przykładowy plik z danymi (wsp_xyz_inp.txt) jest dołączony razem z programem.

7.Pliki wyjściowe:
Plik wyjściowy jest w formacie txt, ma taką samą struktórę jak plik wejściowy i zapisuje się w tym samym folderze,
w którym znajduje się program oraz plik wejściowy.


