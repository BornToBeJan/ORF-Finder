# ORF-Finder

ORF Finder by Jan Renziehausen

Funktion: ermöglicht die Identifizierung von offenen Leserastern (Open Reading Frames) in Sequenzdateien im 
	  FASTA-Format. Das Programm identifziert Start- und Stopcodons innerhalb der Sequenz innerhalb aller 3 möglichen 
	  Leseraster pro Strang und matcht diese so zueinander, dass kein Stop-Codon überlesen wird. 
	  (d.h. Read-Through wird nicht berücksichtigt)


Umgang mit Introns: Das Programm ist so ausgelegt, alle Kleinbuchstaben aus der Sequenz zu entfernen.
		    ORF-Längen und Positionen sind also immer ohne Introns angegeben.


Aufruf des Programms mit folgendem Eingabeformat über die Shell:

	python ORF_FinderV8.py Q -z
	
	Q - Pfad der Quelldatei (FASTA-Datei) (notwendiges Argument)
	-z - Pfad der Zieldatei (Bitte Pfad/Dateiname.csv angeben) (optionales Argument)

	Wenn -z nicht spezifiziert wird, wird das Ergebnis in Results.csv im gleichen Ordner wie
	das Python-Programm gespeichert.


Eingabe von Suchparametern für ORFs:

	1) minimale Länge des ORFs: Ganzzahl > 0 und <= 2.5 Mio oder 'a' für default Einstellung
				 default: 30
	2) maximale Länge des ORFs: Ganzzahl >= Minimum und <= 2.5 Mio oder 'a' für default Einstellung
				 default: 2.5 Mio
	
	3) Angabe der Spezies, aus dem die Sequenz stammt:
	'v': Virus
	'p': Prokaryot
	'e': Eukaryot
	
	folgende Start- und Stopcodons werden nach Wahl der Spezies hinterlegt:
	
	Start-Codons:
	Viren: 		["ATG", "CTG", "GTG"]
	Prokaryoten: 	["ATG", "GTG", "TTG"]
	Eukaryoten:	["ATG"]

	Stop-Codons (unabhängig von Spezies):
	["TAG", "TGA", "TTA"]

	4) Beantworten der Frage, ob nach regulatorischen Elementen gesucht werden soll:
	   ('y' für ja, beliebige Eingabe für nein)
	
	Es wird je nach Wahl der Spezies nach folgenden regulatorischen Elementen in einem definierten 
	Abstand vom ORF gesucht.
	(Distanz gemessen von der ersten Base der Konsensus-Sequenz des regulatorischen Elements zur ersten Base 
	des Startcodons des ORFs):

	Viren: 		keine regulatorischen Elemente
	Prokaryoten: 	Pribnow-Box ["TATAAT"] 			   Distanz vom ORF: 10 bis 100 Basenpaare
			Shine-Dalgarno-Sequenz ["AGGAGGT"] 	   Distanz vom ORF: 10 bis 16 Basenpaare 
	Eukaryoten: 	TATA-Boxen 
			[ "TATAAAA","TATAAAT","TATATAA","TATATAT"] Distanz vom ORF: 25 bis 36 Basenpaare
	
	Es wird nur nach den regulatorischen Elementen der angegebenen Spezies gesucht.		 
	
Ausgabe des Ergebnisses im .csv - Format:
	
	Auflistung aller gefundenen ORFs im Tabellen-Format mit folgenden Spalten:
	
	Strang: 	Angabe des Strangs auf dem sich ORF befindet (Forward (Sense-Strang) oder Reverse (Anti-Sense-Strang))
	Index: 		Laufnummer des ORFs (beginnt für jeden Strang neu bei 1)
	Startindex: 	Angabe der Position der ersten Base des Startcodons des ORFs
	Stopindex:	Angabe der Position der letzten Base des Stopcodons des ORFs
	Laenge:		Anzahl der Basenpaare des ORFs (von Startindex bis Stopindex)
	TATA-Box:	Angabe, ob zum ORF eine TATA-Box im oben angegebenen Abstand gefunden wurde
	Shine-Dalgarno:	Angabe, ob zum ORF eine Shine-Dalgarno-Sequenz im oben angegebenen Abstand gefunden wurde
	Pribnow:	Angabe, ob zum ORF eine Pribnow-Box im oben angegebenen Abstand gefunden wurde
	Sequenz:	Ausgabe der Sequenz des gefundenen ORFs
	
