# -*- coding: utf-8 -*-
"""
Created on Mon Mar 14 09:36:44 2022

@author: Jan Renziehausen und Florian Carlo Fischer
about: ORF Finder
"""

import argparse

StartCodonsEuk = ["ATG"]
StartCodonsProk = ["ATG", "GTG", "TTG"]
StartCodonsVir = ["ATG", "CTG", "GTG"]

StopCodons = ["TAG", "TGA", "TAA"]

TATABOX= [ "TATAAAA","TATAAAT","TATATAA","TATATAT"]
ShineDalgano = ["AGGAGGT"]
Pribnow = ["TATAAT"]


##############################################
# Funktion zum Prozessieren des Eingabefiles #
##############################################

def readFile(Pfad):
    
    FastaFile = open(Pfad, 'r')
    Text = FastaFile.readlines()[1:]
    FastaFile.close()
    
    #Zusammenfügen aller Zeilen bis auf die erste und Löschen der NewLine-Zeichen
    Sequence ="".join(Text)
    Sequence = Sequence.replace('\n', "")

    return Sequence



###############################################################
# Funktion zum Erzeugen der Stränge und Entfernen der Introns #
###############################################################

def CutIntronsAndCreateReverseComplement(Sequence):
    
    Sense_strand = ""
    ComplSeq = ''
    
    for letter in Sequence: 
        
        #Komplementärstrang wird parallel gebaut
        if letter == 'A':
             ComplSeq += 'T'
        elif letter == 'C':
             ComplSeq += 'G'
        elif letter == 'T' or letter == 'U':
             ComplSeq += 'A'
        elif letter == 'G':
             ComplSeq += 'C'
        elif ord(letter) > 64 and ord(letter) < 91:
             ComplSeq += letter
        else:
            continue
        
        #Zusammenbauen des Sense-Strangs
        Sense_strand += letter
    
    # Umdrehen des Komplementärstrangs zum Erzeugen des Antisense-Strangs
    Antisense_strand = ComplSeq[::-1]
        
    return [Sense_strand, Antisense_strand]


####################################################
# Funktion zum Handeln und Speichern von Userinput #
####################################################

def askUserforSpecs ():
    
    UserInput = {
        "Min": 30,
        "Max": 2500000,
        "Organism" : StartCodonsEuk,
        "regElements" : False
    }
       
    # Eingabe minimale ORF-Länge
    while True:      
        Eingabe = input("\nMinimale ORF-Laenge \n(Zahl > 0 und =< 2.5 Mio oder 'a' für Standardeinstellung(30)): ")
        if Eingabe == 'a' or Eingabe == 'A':
            break      
        try:
            Eingabe = int(Eingabe)
        except:
            print("\nkeine valide Eingabe!")
            continue
        if Eingabe > 0 and Eingabe <= UserInput["Max"]:
            UserInput["Min"] = Eingabe
            break
        else:
            print("\nkeine valide Eingabe!")
    
     # Eingabe maximale ORF-Länge
    while True:      
        Eingabe = input("\nMaximale ORF-Laenge \n(Zahl >= Minimum oder 'a' für Standardeinstellung(2.5 Mio)): ")
        if Eingabe == 'a' or Eingabe == 'A':
            break      
        try:
            Eingabe = int(Eingabe)
        except:
             print("\nkeine valide Eingabe!")
             continue
        if Eingabe > 0 and Eingabe >= UserInput["Min"]:
            UserInput["Max"] = Eingabe
            break
        else:
            print("\nkeine valide Eingabe!")
            
    # Eingabe Organismus
    while True:      
        Eingabe = input("\nUm welchen Organismus handelt es sich?\n(noetig, damit Start- und Stopcodons adequat gewählt werden.\n'e' fuer Eukaryoten, 'p' fuer Prokaryoten, 'v' fuer Viren): ")
        if Eingabe != 'E' and Eingabe != 'e' and Eingabe != 'p' and Eingabe != 'P' and Eingabe != 'V' and Eingabe != 'v':
            print("\nkeine valide Eingabe!")
        elif Eingabe == 'E' or Eingabe == 'e':
            break
        elif Eingabe == 'p' or Eingabe == 'P':
            UserInput["Organism"] = StartCodonsProk
            break
        elif Eingabe == 'v' or Eingabe == 'V':
            UserInput["Organism"] = StartCodonsVir
            break
    
     # Frage, ob regulatorische Sequenzen gescannt werden sollen 
    Eingabe = input("\nSoll nach regulatorischen Sequenzen gescannt werden?\n(Erhoeht die Laufzeit.)\n('y' für Ja, beliebig fuer Nein): ")
    
    if Eingabe == 'y' or Eingabe == 'Y':
        UserInput["regElements"] = True       
        
    return UserInput


#################################################
# Funktion zum Finden aller Open-Reading-frames #
#################################################

def findOpenReadingFrame(Sequence, StopCodons, Startposition, UserInput):
      
    AllORFs =[]
    ORFStart = []
    ORFStop = []
    
    #find StartCodon
    for i in range(Startposition, len(Sequence)-2, 3):

        #Kreiert ein Codon aus dem nächsten 3 Buchstaben des Strings
        Codon = Sequence[i] +  Sequence[i+1] + Sequence[i+2]
        
        #Vergleicht das Codon mit den StartCodons und speichert die Start-Indizes im String in einer Liste
        for StartCodon in UserInput["Organism"]:
            if StartCodon == Codon:
                ORFStart.append(i)
        
        #Vergleicht das Codon mit den StopCodons und speichert die Stop-Indizes im String in einer Liste
        for StopCodon in StopCodons:
            if StopCodon == Codon:
                ORFStop.append(i)
                
        Codon = ''
                   
    # Wenn kein Startcodon bzw. StopCodon gefunden wird, 
    #verlasse die Funktion und gib Fehlermeldung aus
    if len(ORFStart) < 1 or len(ORFStop) < 1:
        return None

    # fügt die Start- und Stopcodons so zusammen, dass kein Stopcodon überlesen wird   
    j = 0
    for i in range(0, len(ORFStop)):
        while j < len(ORFStart):
          
            if ORFStart[j] < ORFStop[i]:
                
                # Berechnet die ORF-Länge und fügt ORF zu ALLORFs hinzu, wenn innerhalb von Min und Max
                orfLength = ORFStop[i] + 3 - ORFStart[j]
                
                if orfLength >= UserInput['Min'] and orfLength <= UserInput['Max']:
                   
                    # Erstellen eines OpenReadingFrames aus Start und Stopcodon
                    ORF =  {
                        'x': ORFStart[j],
                        'y': ORFStop[i] + 2,
                        'Length': ORFStop[i] + 3 - ORFStart[j],
                        'TATA': None,
                        'Pribnow': None,
                        'Shine': None
                        }
                
                    #Scannen nach regulatorischen Elementen upstream des gefundenen ORFs
                    if UserInput["regElements"] == True:
                        
                        # für Eukaryoten
                        if UserInput["Organism"] == StartCodonsEuk:
                            find_regulatory_Elements(ORF, Sequence, TATABOX)
                            
                        # für Prokaryoten
                        if UserInput["Organism"] == StartCodonsProk:   
                            find_regulatory_Elements(ORF, Sequence, Pribnow)
                            find_regulatory_Elements(ORF, Sequence, ShineDalgano)
            
                    AllORFs.append(ORF)
                j += 1
            else:
                break
    
    # Wenn keine ORFs mit passender Länge gefunden wurden, wird None zurückgegeben
    if len(AllORFs) < 1:
        return None
    
    # Ansonsten gib alle OpenReadingFrames innerhalb eines Shifts zurück
    return AllORFs


#####################################################
# Funktion zum Finden von regulatorischen Elementen #
#####################################################

def find_regulatory_Elements(ORF, Sequence, regElement):
    
    
    # Lage der TATA-BOX upstream relativ zum ORF
    if regElement == TATABOX:
        upperBorder = 35
        lowerBorder = 25
        regElementName = 'TATA'   
    
    # Lage der ShineDalgano-Sequenz upstream relativ zum ORF
    elif regElement == ShineDalgano:
        upperBorder = 16
        lowerBorder = 10
        regElementName = 'Shine'
    
    # Lage der Pribnow-Sequenz upstream relativ zum ORF
    elif regElement == Pribnow:
        upperBorder = 100
        lowerBorder = 10
        regElementName = 'Pribnow'
    
    else:
        print("Unbekanntes regulatorisches Element!")
        return None
          
    Codon = ""
     
    # Suche nach regulatorischen Elementen upstream des betrachteten ORFs
    for i in range(lowerBorder, upperBorder + 1):        
               
        DistanceToORF = ORF['x'] - i

        # Zusammenbauen eines Vergleichcodons vom ORF-Startpunkt aus im entsprechenden Abstand 
        # zum ORF welches die Länge des regulatorischen Elements besitzt
        
        if DistanceToORF >= 0:
            for t in range(0, len(regElement[0])):
                Codon +=  Sequence[DistanceToORF + t]
  
        # Vergleich des Codons mit den Sequenzen des jeweiligen regulatorischen Elements
        # und bei Erfolg Setzen des entsprechenden Eintrags im betrachteten ORF auf True
        for j in regElement:
            if Codon == j:
                ORF[regElementName] = True
                break
        
        # Verlässt die Funktion, sobald regulatorisches Element gefunden wurde, 
        # ansonsten wird regElement auf false gesetzt und weitergesucht
        if ORF[regElementName] == True:
            break
        else:
            ORF[regElementName] = False
        
        Codon = ""
    

######################################################################
# Funktion zum Ausgeben aller Open-Reading-Frames auf dem Bildschirm #
######################################################################

def printAllOpenReadingFrames(AllOpenReadingFrames):
    
    Ausgabestring = ''
    
    for i in range(0,len(AllOpenReadingFrames)):
        if i == 0:
            Ausgabestring = "\nForward Strang: (" + str(len(AllOpenReadingFrames[i])) + " ORFs found):\n"
        if i > 0:
            Ausgabestring = "Reverse Strang: (" + str(len(AllOpenReadingFrames[i])) + " ORFs found):\n"
        if AllOpenReadingFrames[i] != None:
            Ausgabestring += str(AllOpenReadingFrames[i])
        else:
            Ausgabestring += "kein Open Reading Frame gefunden!"
        print(Ausgabestring)
        


#########################################################
# Funktion zum Schreiben der Ergebnisse in Ausgabedatei #
#########################################################

def writeResultsToFile(Pfad, AllOpenReadingFrames, Sequence, revSeq):

    AusgabeFile =open(Pfad, 'w')
        
    # Schreiben des Headers
    AusgabeFile.write("Strang; ORF Nr.;Startindex;Stopindex;Laenge; TATA-BOX; Shine-Dalgano; Pribnow; Sequenz\n")
     
    for i in range(0, len(AllOpenReadingFrames)):
        for j in range (0, len(AllOpenReadingFrames[i])):
            
            if i == 0: # Schreiben des aktuell betrachteten Strangs
                AusgabeFile.write("Forward;")
            else:
                AusgabeFile.write("Reverse;")
                    
            AusgabeFile.write(str(j+1)+ ";") # Schreiben des ORF-Index 
            AusgabeFile.write(str(AllOpenReadingFrames[i][j]['x'] + 1) + ";") #Schreiben von Startindex
            AusgabeFile.write(str(AllOpenReadingFrames[i][j]['y'] + 1) + ";") # Schreiben von Stopindex
            AusgabeFile.write(str(AllOpenReadingFrames[i][j]['Length']) + ";") #Schreiben der ORF-Länge
            
            # Eintragen, ob regulatorische Elemente gefunden wurden 
            #Schreiben, ob TATA ja / nein / n/a
            if AllOpenReadingFrames[i][j]['TATA'] != None:
                if AllOpenReadingFrames[i][j]['TATA'] == True: 
                    AusgabeFile.write("ja;")
                else:
                    AusgabeFile.write("nein;")
            else:
                AusgabeFile.write("n/a;")
            
            #Schreiben, ob Shine-Dalgano ja / nein / n/a
            if AllOpenReadingFrames[i][j]['Shine'] != None:
                if AllOpenReadingFrames[i][j]['Shine'] == True:
                    AusgabeFile.write("ja;")
                else:
                    AusgabeFile.write("nein;")
            else:
                AusgabeFile.write("n/a;")

            #Schreiben, ob Pribnow ja / nein / n/a
            if AllOpenReadingFrames[i][j]['Pribnow'] != None:
                if AllOpenReadingFrames[i][j]['Pribnow'] == True:
                    AusgabeFile.write("ja;")
                else:
                    AusgabeFile.write("nein;")
            else:
                AusgabeFile.write("n/a;")
        
        
            # Schreiben der jeweiligen Sequenz bzw. rev. komplementierten Sequenz
            if i == 0: 
                AusgabeFile.write(Sequence[AllOpenReadingFrames[i][j]['x'] : AllOpenReadingFrames[i][j]['y'] + 1] + "\n")
            else:
                AusgabeFile.write(revSeq[AllOpenReadingFrames[i][j]['x'] : AllOpenReadingFrames[i][j]['y'] + 1] + "\n")
       
    AusgabeFile.close()


###################
# Main - Programm #
###################

#Liste des Sense-Strangs
AllOpenReadingFrames_forward = []
#Liste des Antisense-Strangs
AllOpenReadingFrames_reverse = []

# Liste zum Speichern beider Stränge
AllOpenReadingFrames = []

# Argument-Parsing
parser = argparse.ArgumentParser(description='ORF-Finder')
parser.add_argument('Quellpfad', metavar='Q', type = str, help = 'Pfad der Quelldatei')
parser.add_argument('-z','--Zielpfad', type = str, help = 'Pfad der Zieldatei', default = 'Results.csv')
args = parser.parse_args()

# Abfragen der User-Specs
UserInput = askUserforSpecs()
print("\nProcessing input file...")

# Einlesen der Sequenz aus dem angegebenen Pfad
print("Reading sequence...")
Sequence = readFile(args.Quellpfad)


# Intron Splicing und Erzeugung des Antisense-Strangs
print("Cutting introns and creating both strands...")
strands = CutIntronsAndCreateReverseComplement(Sequence)


# Kreieren aller Open Reading Frames (jeder Strang mit 3 Shifts (0, 1, 2))
print("Finding open reading frames...")
for i in range(0, 3):
    # liefert alle Open-Reading-Frames innerhalb eines Shifts als eine Liste
    Temp_Liste = findOpenReadingFrame(strands[0], StopCodons, i, UserInput)
    
    # kopiert jeden Eintrag aus der Shift-Liste als einzelne Einträge in eine Liste für den gesamten Forward-Strang
    if Temp_Liste != None:
        for j in Temp_Liste:
            AllOpenReadingFrames_forward.append(j)   
    
    # das gleiche wird für den Reverse-Strang gemacht
    Temp_Liste = findOpenReadingFrame(strands[1], StopCodons, i, UserInput)
   
    if Temp_Liste != None:
        for j in Temp_Liste:
            AllOpenReadingFrames_reverse.append(j)

# Sortieren der Listen entsprechend des Startcodon Index
AllOpenReadingFrames_forward = sorted(AllOpenReadingFrames_forward, key = lambda l:l['x'])
AllOpenReadingFrames_reverse = sorted(AllOpenReadingFrames_reverse, key = lambda l:l['x'])

# Hinzufügen von Forward- und Reverseliste zu einer großen Gesamtliste
AllOpenReadingFrames.append(AllOpenReadingFrames_forward)
AllOpenReadingFrames.append(AllOpenReadingFrames_reverse) 

#Ausgabe des Ergebnisses in Ausgabedatei
writeResultsToFile(args.Zielpfad, AllOpenReadingFrames, strands[0], strands[1])
print("\nErgebnisse wurden in " + args.Zielpfad + " exportiert!")
