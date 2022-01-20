# ModelingAndSimulation

## Motivation
Die Strömung inkompressibler Fluide, d.h. Gase und Flüssigkeiten, bei denen Druck keine Volumenänderung bewirkt, wird durch die Impulsgleichung ("Navier-Stokes-Gleichung") und die Kontinuitätsgleichung beschrieben. Diese Gleichungen bilden häufig die Grundlage für Strömungssimulationen, wobei verschiedene numerische Lösungsmethoden eingesetzt werden. Die numerische, näherungsweise Lösung hat aber den Nachteil, dass Rundungsfehler oder numerische Instabilitäten auftreten können. Dadurch kann das Modell wesentliche Eigenschaften des entsprechenden physikalischen Systems verlieren, z.B. Masse- und Impulserhaltung. Ein Ansatz, um dieses Problem zu umgehen, ist, stattdessen ein diskretes Modell zu verwenden, das nicht näherungsweise berechnet werden muss (siehe Abbildung 1). Das HPP- und das FHP-Modell sind solche Modelle. Sie wurden in der zweiten Hälfte des 20. Jahrhunderts entwickelt und gehören zur Klasse der zellulären Automaten.

## Das HPP-Modell
Die Grundidee dieses Modells ist, Bewegung und Interaktion (Kollision) von Molekülen auf einem endlichen Gitter durch einfache Regeln nachzubilden, sodass Impuls und Masse erhalten bleiben. Dazu wird ein zweidimensionales, regelmäßiges, rechteckiges Gitter verwendet, das mit einer Von-Neumann- Nachbarschaft versehen ist. Es werden 16 Zellenzustände unterschieden: eine Zelle ist entweder leer oder enthält bis zu vier Partikel, die durch ihre Bewegungsrichtung identifiziert werden (Rechts, Hinauf, Links, Hinunter). Jede Richtung kann pro Zelle höchstens einmal vorkommen. Die Partikel sind dabei ununterscheidbar. Es gibt nur eine Kollisionsregel (siehe auch Abbildung 2): Befinden sich genau zwei Partikel in einer Zelle und haben diese entgegengesetzte Bewegungsrichtungen, so werden die Partikel jeweils im rechten Winkel abgelenkt und verlassen die Zelle in entgegengesetzte Richtungen. Ansonsten wandern die Partikel - anschaulich gesprochen - in diejenige Nachbarzelle, die der aktuellen Bewegungsrichtung entspricht. Die Kollisionsregel sorgt dafür, dass beim Aufeinandertreffen zweier Partikel mit entgegengesetztem Impuls (Bewegungsrichtung) diese danach wieder entgegengesetzten Impuls haben und außerdem keine Masse (kein Partikel) verloren geht. Findet keine Kollision statt, bleiben Impuls und Masse ebenso erhalten.

## Das FHP-Modell
Dieses Modell baut auf dem HPP-Modell auf und es gibt drei verschiedene Varianten (FHP-I, FHP-II, FHP-III), die sich durch ihre Kollisionsregeln unterscheiden. Alle Varianten verwenden ein zweidimensionales, regelmäßiges, dreieckiges Gitter. Dadurch hat jede Zelle eine direkte Verbindung zu sechs Nachbarzellen. Entsprechend ist auch die Anzahl zulässiger Partikel pro Zelle größer als beim HPP-Modell. FHP-I erlaubt sechs Partikel (eines pro Richtung), FHP-II- und FHP-III sogar sieben, weil zusätzlich ein unbewegtes Partikel erlaubt ist, womit es insgesamt 64 bzw. 128 mögliche Zellenzustände gibt. Je nach Variante gibt es Kollisionsregeln für zwei, drei und vier (in FHP-II und FHP-III) Partikel. Diese Kollisionsregeln sind zusammen mit der HPP-Kollisionsregel in Abbildung 2 zusammengefasst.  


## ToDo's
### Task 1
Implementieren Sie einen zellulären Automaten, der die Gegebenheiten in einem Strömungskanal repräsentiert (vgl. Abbildung 3). Die Ausmaße des Gitters (Standardszenario: 200 x 500-Matrix) sowie die Form (HPP oder FHP) und Platzierung von Hindernissen sollten frei wählbar sein.

### Task 2
Erstellen Sie eine Plotvariante, um die mittlere Bewegungsrichtung in den Zellen als Vektoren anzuzeigen. Verwenden Sie dazu die MATLAB-Funktion quiver und vergleichen Sie die verschiedenen Modelle miteinander.

### Task 3
Erstellen Sie nun eine Plotvariante, um die Verteilung der Partikeldichte auf dem Gitter anzuzeigen. Verwenden Sie dazu Pseudoisolinien, die mit der MATLAB-Funktion contourf erzeugt werden.

### Task 4
Implementieren Sie für HPP und FHP-II eine Plotvariante, um die Bewegung eines Partikels über das Gitter zu verfolgen und den gesamten bisher zurückgelegten Weg anzuzeigen. Vergleichen Sie die beiden Modelle in Bezug auf die Partikelbewegung.

### Task 5
Versuchen Sie mit FHP-III eine Kármánsche Wirbelstraße zu erzeugen.

## Programmiersprache
MATLAB
