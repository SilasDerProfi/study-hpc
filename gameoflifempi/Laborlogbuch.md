# Labor 03
* Welche Probleme ergeben sich im Zusammenhang mit einer 2D Gebietszerlegung?
  * Daten sind nicht aneinanderliegend im Speicher -> Können nicht ohne weiteres an anderen Worker gesendet werden
  * Komplexe MPI-Datentypen benötigt
* Welche Vorraussetzungen müssen an den Aufrufer des Programms gestellt werden?
  * Das Feld sollte in jeder Dimension durch die Anzahl an Threads teilbar sein.
