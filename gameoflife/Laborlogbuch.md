# Labor 01 Aufgabe 3
* Wenn pX > pY, ist das besser (pX=4, pY=1 ist am besten)
 - Gilt für alle gemessenen / durchgeführten Spielfeldgrößen
 - Aufgrund des Aufwandes für den Speicherzugriffs
* 4 x Größeres Spielfeld braucht 4 x mehr Zeit

## Arbeitsspeicher Berechnet
-------- 2 Spielfelder * Breite * Höhe * Double-Größe --------
2 * 1024 * 1024 * 8 Byte =  16.777.216 B =  16.777 KB =  17 MB
2 * 2048 * 2048 * 8 Byte =  67.108.864 B =  67.109 KB =  67 MB
2 * 4096 * 4096 * 8 Byte = 268.435.456 B = 268.435 KB = 268 MB

## Festplattenplatz Berechnet (1 Timestep)
-------- 1 Spielfeld * Breite * Höhe * Float-Größe --------
1 * 1024 * 1024 * 4 Byte =  4.194.304 B =  4.194 KB =  4 MB
1 * 2048 * 2048 * 4 Byte = 16.777.216 B = 16.777 KB = 17 MB
1 * 4096 * 4096 * 4 Byte = 67.108.864 B = 67.109 KB = 67 MB


## Festplattenplatz gemessen (1 Timestep)
------- Größer der VTI- + PVTI-Dateien -------
1024 x 1024 =  4.195.120 B =  4.195 KB =  4 MB
2048 x 2048 = 16.778.032 B = 16.778 KB = 17 MB
4096 x 4096 = 67.109.680 B = 67.110 KB = 67 MB