# Mars: hoge-orde zwaartekracht profileren

Standaard: dertig simulatiedagen, graad/order 100, alleen dynamica.
Bij de lokale kalibratie kostten tien dagen circa 24,1 seconden integratie.
Korte controles met graad 20, 50 en 100 en met variatievergelijkingen bij graad
100 slaagden. De twee modi leverden bij dezelfde korte graad-100-run dezelfde
gerapporteerde eindstaat. De volledige standaardrun is nog niet uitgevoerd.

`profile_mars_gravity.py` isoleert de zwaartekracht-workload uit het
onderzoeksvoorstel. Het gebruikt het lokaal beschikbare JGMRO120D-model,
expliciet geselecteerd, met instelbare maximale graad **en** order (standaard
100). De satelliet begint in een bijna cirkelvormige baan op circa 300 km boven
de zwaartekrachtreferentiestraal, met inclinatie 93°. Mars roteert volgens de
standaard Tudat/SPICE-instellingen.

Alleen Mars-zwaartekracht werkt op de satelliet. Er zijn geen atmosfeer,
panelen, observaties, schatting of Python-callbacks in de integratie. Dit is een
gerichte zwaartekrachtbenchmark, geen volledig missiemodel. De SPICE-rotatie
en de normale integrator- en opslagkosten blijven onderdeel van de run.

## Korte controle en vergelijking van graden

Start vanuit de repositoryroot met de Python-omgeving van de gewenste
geoptimaliseerde Tudat-build met debug-symbolen:

```bash
conda activate tudatpy-with-debug-symbols
export OPENBLAS_NUM_THREADS=1
export OMP_NUM_THREADS=1

python benchmarks/profile_mars_gravity.py --days 0.1
python benchmarks/profile_mars_gravity.py --days 0.1 --mode variational

for degree in 20 50 100; do
    python benchmarks/profile_mars_gravity.py --degree "$degree" --days 10
done
```

Iedere aanroep is een vers proces. Gebruik dezelfde duur, stap en beginbaan
voor de vergelijking. De baan zelf zal door de verschillende zwaartekracht-
truncaties enigszins verschillen. De standaardintegrator is RKDP87 met vaste
stap van 30 seconden. De duur is instelbaar met `--days`, de stap met `--step`.
Voor 60–180 seconden integratietijd kunt u de duur ongeveer lineair schalen
op basis van de gerapporteerde `wall_seconds`. Verleng de run in plaats van
de stap kunstmatig te verkleinen om hem zwaarder te maken.

De optionele `--mode variational` propageert ook de 6×6 toestandovergangsmatrix
ten opzichte van de beginstaat. Er worden geen duizenden afzonderlijke
zwaartekrachtcoëfficiënten als parameters toegevoegd. Hiermee worden ook de
zwaartekrachtafgeleiden/Hessiaan belast, zonder een sterk vergrote parameterruimte.

## Opnemen en rapporteren

```bash
perf record -o mars-g100-dyn.perf.data -F 199 --call-graph dwarf,32768 -- \
    python benchmarks/profile_mars_gravity.py --degree 100 --days 30

perf report -i mars-g100-dyn.perf.data --stdio --no-children \
    --call-graph none --sort dso,symbol --percent-limit 0.5 \
    > mars-g100-dyn-self.txt

perf report -i mars-g100-dyn.perf.data --stdio --children --percent-limit 1 \
    > mars-g100-dyn-callers.txt

# Apart profiel voor de variatievergelijkingen
perf record -o mars-g100-var.perf.data -F 199 --call-graph dwarf,32768 -- \
    python benchmarks/profile_mars_gravity.py --degree 100 --days 30 --mode variational

perf report -i mars-g100-var.perf.data --stdio --no-children \
    --call-graph none --sort dso,symbol --percent-limit 0.5 \
    > mars-g100-var-self.txt
```

Kijk eerst naar de eigen kosten van Legendre-berekeningen, sferisch-harmonische
versnellingssommatie en, in de variatiemodus, Hessiaanberekeningen. Een dominant
aandeel is een te toetsen hypothese, geen eigenschap die het script afdwingt.
Inclusieve percentages overlappen. De eerdere LEO-opname had onbetrouwbare
callstacks; controleer die kwaliteit hier opnieuw voordat u groepsaandelen
interpreteert. Geen verloren samples garandeert geen correcte stackreconstructie.

## Uitvoer en grenzen

De JSON-uitvoer bevat het pakketpad, model, geladen coëfficiëntmatrixgrootte,
instellingen, opbouwtijd, integratietijd, eindstaat en minimale opgeslagen hoogte
boven de referentiestraal. De variatiemodus controleert daarnaast dat de
eindmatrix 6×6 en eindig is. De eindtijd en succesvolle integratie worden
gecontroleerd. Dit zijn plausibiliteitscontroles; geen bewijs van
missienauwkeurigheid. Controleer desgewenst de eindstaat met gehalveerde stap,
bij dezelfde graad en duur. De referentiestraal is geen topografisch oppervlak.

De integratietimer omvat de native integratie en bijbehorende native opslag.
Python-export en controle van de geschiedenis vinden daarna plaats. `perf`
neemt het hele proces op, inclusief deze export, imports en initialisatie;
maak daarom de integratiefase dominant. De printmarkers schakelen de profiler
niet aan of uit. Geheugengebruik neemt toe met de opgeslagen geschiedenis.
Gebruik dezelfde opslaginstellingen bij prestatievergelijkingen.

Het script rapporteert welk geïnstalleerd pakket wordt gebruikt. Controleer bij
C++-wijzigingen dat dit daadwerkelijk de opnieuw gebouwde versie is. De lokale
broncheckout en de geïnstalleerde debugomgeving hebben niet exact dezelfde API.
