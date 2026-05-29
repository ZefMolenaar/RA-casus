# RA-casus transcriptomics analyse veranderingen in genregulatie reumatoïde artritis.
## Project transcripnomics J2P4 / Zef Molenaar / BT2A

## Inhoud

- `Assets/` - Overige documenten voor opmaak
- `Bronnen/` - Gebruikte bronnen
- `Data/Processed` - Verwerkte datasets gegenereerd met scripts 
- `Data/Raw` – De gebruikte ruwe data
- `Uitleg competentie beheren/` - Uitleg competentie beheren 
- `Resultaten/` - Verkregen grafieken en tabellen
- `Script/` – Gebruikte script om de transcriptomics-analyse mee uit te voeren


---


## Inleiding: 
Wereldwijd wordt ongeveer 1% van de bevolking getroffen door de auto-immuunziekte reumatoïde artritis (RA). Het ontstaan van deze aandoening wordt beïnvloed door erfelijke aanleg en omgevingsinvloeden. Op genetisch vlak vormen specifieke HLA-DRB1-allelen de belangrijkste risicofactor. Wat betreft omgevingsvariabelen is er een wisselend effect te zien: zo hebben vrouwelijke hormonen een beschermende werking en verlagen ze de kans op RA, terwijl roken de kans op het ontwikkelen van de ziekte juist aanzienlijk vergroot [(Silman & Pearson, 2002)](Bronnen/SilmanPearson2002.pdf)

Pijnlijke en gezwollen gewrichten zijn de typische kenmerken van reumatoïde artritis (RA). Naast deze directe gewrichtsklachten lopen patiënten een groter risico op ernstige bijkomende aandoeningen, waaronder infecties, bloedarmoede (anemie), hart- en vaatziekten en specifieke typen kanker [(Wilson et al., 2004)](Bronnen/Wilson2004.pdf), [(Sparks, 2019)](Bronnen/Sparks2019.pdf). Wanneer de ontsteking chronisch wordt, kan er onherstelbare schade ontstaan aan het bot en het kraakbeen, wat op den duur kan leiden tot functionele beperkingen of een handicap [(Smolen et al., 2016)](Bronnen/Smolen2016.pdf),[(Karlson et al., 2008)](Bronnen/Karlson2008.pdf).

Omdat er nog geen geneesmiddel is voor RA worden patienten die klachten ervaren persoonlijk behandeld om specifieke symptomen te onderdrukken [(Bullock et al., 2019)](Bronnen/Bullock2019.pdf). Hierdoor kan 90% van de patiënten de gewrichtsschade sterk worden afgeremd of zelfs helemaal worden voorkomen. Dit is van essentieel belang om blijvende schade te voorkomen [(Aletaha & Smolen, 2018)](Bronnen/AletahaSmolen2018.pdf).

Aangezien de exacte oorzaak van RA nog niet volledig is opgehelderd, ligt de focus van dit onderzoek op het in kaart brengen van genen die een verhoogde of verlaagde expressie vertonen bij RA-patiënten. Om dit te realiseren, wordt er een transcriptomics-analyse uitgevoerd op synoviumweefsel. De verkregen data worden vervolgens met behulp van R geanalyseerd om de significante verschillen in genexpressie tussen personen met en zonder RA te identificeren.




## Methode:
## Materiaal en Methode

Er is gebruikgemaakt van acht synoviumbiopten: vier van ACPA-negatieve controlepersonen zonder reumatoïde artritis (RA) en vier van ACPA-positieve (anti-CCP) patiënten met established RA (> 12 maanden na diagnose). De ruwe gegevens zijn te downloaden als zip-bestand onder: [Data/bam files], [Data].

---

## RNA-seq Analyse

De analyse van de sequencingdata is uitgevoerd in R. Met behulp van het pakket **Rsubread** (versie 2.24.0) [(Liao et al., 2019)](Bronnen/Liao2019.pdf) zijn de reads uitgelijnd op het humane referentiegenoom (GCF_000001405.40_GRCh38.p14_genomic). Vervolgens zijn de verkregen BAM-bestanden gesorteerd en geïndexeerd met het pakket **Rsamtools** (versie 2.26.0). Hierbij is gebruikgemaakt van het bijbehorende RefSeq GTF-annotatiebestand (versie GCF_000001405.40) van NCBI.

---

## Genexpressieanalyse

Er is op basis van de BAM-bestanden met **Rsubread** een count matrix gegenereerd. Er is gebruikgemaakt van het GTF-annotatiebestand dat hoort bij het gebruikte humane referentiegenoom. Daarnaast is de differentiële genexpressieanalyse is het pakket **DESeq2** (versie 1.50.2) [(Love et al., 2014)](Bronnen/Love2014.pdf) gebruikt. Visualisatie van de resultaten is gedaan met behulp van **EnhancedVolcano** (versie 1.28.2) waarmee een volcano plot is gemaakt.

---

## Pathway / Functieanalyse
Genoverrepresentatie-analyse is uitgevoerd met het **goseq** pakket (versie 1.62.0) (https://cloud.wikis.utexas.edu/wiki/spaces/bioiteam/pages/47732482/GO+Enrichment+using+goseq). Deze pathview werd gevisualiseerd met de **pathview** (versie 1.50.0)





## Resultaten
+- 200 woorden, inclusief correcte verwijzingen.

## Conclusie:
+- 200 woorden, inclusief aanbevelingen en onderzoek in context plaatsen.

## Uitleg competentie beheren 
(zie voor hulpvragen het voorbeeld):
o File (bijvoorbeeld een md file) met uitleg over Data Stewardship
o File met uitleg over toepassing beheren met GitHub
