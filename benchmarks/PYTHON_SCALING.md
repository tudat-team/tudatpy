# Strong scaling via de Python-API

`benchmark_python_independent_arcs.py` meet dezelfde variatiepropagatie met
1/2/4/8 native C++-workers. Het gebruikt de publieke
`IndependentArcVariationalSimulator` met geïsoleerde native Mars-modellen.
Het introduceert geen Python-processen, callbacks, nieuwe modelcategorie of
wijziging van het productie-API-/GIL-contract. De standaard blijft een korte
functionele smoke-run. Het bestaande `python_independent_arcs.py` blijft het
kleine gebruiksvoorbeeld.

## Uitvoeren

Gebruik de ontwikkelomgeving en de **nieuw gebouwde** bindings zoals beschreven
in [het API-contract](../docs/python_independent_arcs_api.md). Vanuit de
repositoryroot, op deze machine (de driver gebruikt Python 3.11 of nieuwer):

```bash
conda activate tudatpy-dev
export PYTHONPATH=build-python-arcs/src
export OPENBLAS_NUM_THREADS=1
export OMP_NUM_THREADS=1
python benchmarks/benchmark_python_independent_arcs.py --help
mkdir -p local-profiles/20261112
python benchmarks/benchmark_python_independent_arcs.py \
  --output local-profiles/20261112/python-smoke.jsonl
```

De standaardrun gebruikt acht arcs van 600 s, graad/orde 20, stap 30 s en één
herhaling per workerstand. Hij controleert de werking; overhead kan groter zijn
dan de tijdwinst. Een mededeling over speedup in deze run is geen bewijs van
schaalbaarheid bij zwaar werk.

Een expliciet gekozen middelzware meetreeks voor deze AMD-machine:

```bash
python benchmarks/benchmark_python_independent_arcs.py \
  --degree 50 --days-per-arc 0.25 --repeats 3 \
  --output local-profiles/20261112/python-scaling-quarter-day.jsonl
```

De zwaardere instelling uit het meetplan start **alleen op expliciet verzoek**:

```bash
python benchmarks/benchmark_python_independent_arcs.py --preset representative \
  --output local-profiles/20261112/python-scaling-two-days.jsonl
```

Dat preset gebruikt acht arcs van twee dagen, graad/orde 50 en drie herhalingen.
Dit kan op een andere machine veel langer duren. Er is geen automatische
kalibratie die zelfstandig naar steeds grotere runs opschaalt. CLI-overrides
hebben voorrang op een preset; de uiteindelijke parameters staan in het verslag.

Ook `--arcs`, `--order`, `--step-size`, `--workers` en `--repeats` zijn instelbaar.
De workerlijst moet unieke positieve aantallen bevatten, inclusief 1, en geen
waarde groter dan het arc-aantal. Voor vier arcs is bijvoorbeeld ook
`--workers 1 2 4` nodig. De één-worker-run gebruikt de echte sequentiële default;
alle grotere workerstanden moeten aantoonbaar de native onafhankelijke route
gebruiken. Volgorde wordt per herhaling cyclisch geroteerd.

`--output` opent uitsluitend een **nieuw** bestand en overschrijft niets. Kies
bij herhalen een andere bestandsnaam; een bestaande output levert exitcode 1.
Zonder `--output` schrijft de driver JSONL naar stdout. Voortgang en de leesbare
samenvatting gaan naar stderr. Grote histories worden niet geëxporteerd.
`local-profiles/` is op deze checkout lokaal uitgesloten via `.git/info/exclude`;
die uitsluiting gaat niet mee met een clone. Bewaar uitvoer op een andere machine
bewust lokaal. Kleine meetverslagen kunnen na beoordeling apart in Git worden
opgenomen; deze driver commit of verwijdert niets.

## Model en vergelijkbaarheid

De invoer volgt `benchmarks/cpp/marsArcBenchmark.h:makeSettings`:

- Mars/J2000, Satellite, JGMRO120D en Cowell met vaste RKDP87-stap;
- eerste epoch 778118400 TDB s sinds J2000;
- één uur tussen het einde van een arc en het begin van de volgende;
- Kepler-elementen: halve lange as `reference_radius + 300000` m,
  eccentriciteit 0,001, inclinatie 93 graden, argument van periapsis 0,2 rad,
  rechte klimming van de stijgende knoop 0,3 rad;
- ware anomalie `fmod(0.4 + 0.7 * arc, 2*pi)`;
- zes beginstaatparameters per arc plus één gedeelde Mars-mu.

Conversie naar Cartesische staten gebeurt met Tudats eigen conversieroutine.
De concrete beginstaten, epochs, radius, mu en hash van het gravity-bestand staan
in JSONL. Arcs zijn onafhankelijk en niet-overlappend; het zijn synthetische
schattingsarcs, geen aansluitend missietraject of volledige baanschatting met
observaties/normaalvergelijkingen.

**Vergelijk oude native tijden niet rechtstreeks met deze Python-tijden.**
Naast andere machine-/runcondities zijn er concrete implementatieverschillen:

1. De echte bindings gebruiken `tudat::Time`; de historische native benchmark
   gebruikt `double` voor tijd.
2. De native benchmark slaat ook relatieve posities als afhankelijke variabele
   op en werkt de bijbehorende interpolator bij. De huidige Python-facade biedt
   die extra uitvoerinstelling niet aan.
3. De historische getimede één-worker-run configureert al de geïsoleerde
   factory-route; onze één-worker-run gebruikt de sequentiële default zonder
   workerfactory. De oude tabel vergelijkt dus niet precies dezelfde baseline.
4. Python-opbouw omvat ook de modeldescriptor; initialisatie-/meetgrenzen komen
   niet volledig overeen.

Het JSONL-verslag legt deze verschillen vast. De nieuwe speedup is uitsluitend
binnen één meetreeks berekend, tegen diens eigen één-worker-mediaan. Er wordt
geen verbetering van de C++-kern ten opzichte van de historische tabel geclaimd.

## Tijdgrenzen en numerieke controle

Imports en standaard-SPICE-kernellading gebeuren eenmaal, vóór de metingen.
Elke invocation doet dezelfde korte, ongetimede sequentiële warmup (één arc,
graad/orde 2, 600 s, stap 30 s). Daarna wordt een ongetimede sequentiële
referentie voor de gevraagde volledige workload gemaakt. De reference is voor
correctheid en wordt **niet** als getimede één-worker-baseline hergebruikt.

Iedere getimede run bouwt een verse native modeldescriptor en facade:

| Veld | Meetgrens |
|---|---|
| `setup_seconds` | Modeldescriptor en facade, inclusief configuratie van workers |
| `integration_seconds` | `integrate_full_equations`, inclusief seriële workerfactory en native nabewerking |
| `total_seconds` | Som van bovenstaande fasen |
| `validation_seconds` | Resultaatuitlezing, globale probes en vergelijking; buiten bovenstaande som |

Inputgeneratie, warmup, referentierun, validatie, destructie en JSON-output zitten
niet in de totale tijd. Het is dus geen totale Python-procesduur en geen zuivere
workerkerneltijd. Validatie gebeurt na iedere run met dezelfde werkwijze. Start
geen compilaties of andere zware taken tegelijk met een echte tijdmeting.

Voor elke run worden **alle opgeslagen epochs** van alle arcs vergeleken voor
Cartesische toestanden, lokale 6×6-STM's en 6×1-mu-gevoeligheden. Epochlijsten,
shapes, integratiegrenzen en eindige waarden moeten kloppen. Vergelijking volgt
het native criterium per component:

```text
abs(actual - reference) <= 1e-14 + 1e-11 * abs(reference)
```

Bij drie opgeslagen epochs per arc (begin, midden, einde) controleert de driver
ook of de globale matrix de lokale STM in de juiste zes kolommen zet, alle
inactieve arc-kolommen nul laat en mu in de gedeelde laatste kolom zet. Hij
vergelijkt deze matrices en één extra interpolatiepunt per arc met de
sequentiële referentie. De globale matrix is dus op vier probes per arc
gecontroleerd, niet op alle histories opgeslagen. Voor acht arcs is de shape
6×49. Mu-gevoeligheid op het midden van de arc moet niet identiek nul zijn.

Het verslag bevat maximale absolute componentafwijkingen en maximale ratio's
ten opzichte van de tolerantie, afzonderlijk voor staten, STM's, sensitiviteiten
en globale matrices. Een ratio groter dan 1, een NaN/Inf, een verkeerde route of
een mismatch doet de driver falen. Een gedeeltelijk verslag eindigt dan met een
`error`-record en heeft geen succesvolle `summary`. Deze controles bewijzen
parity met de sequentiële uitvoering, niet de absolute fysieke baannauwkeurigheid.

## Verslag en interpretatie

JSONL bevat achtereenvolgens:

1. `metadata`: werkelijk gebruikte workload, tijdgrenzen, bron-HEAD/branch en
   dirty-status, hash van de driver en tracked diff, systeem/CPU-affiniteit,
   Python/NumPy/pybind11, relevante Conda-/CMakegegevens, compiler,
   geladen modulepad/-grootte/-mtime, threadinstellingen, gravity-hash en invoer;
2. `reference`: controle van de ongetimede sequentiële referentie;
3. `run`: herhaling/positie/workers, gebruikte native route, fase-tijden en
   numerieke controlegegevens;
4. `summary`: mediaan, minimum en maximum per fase/workerstand, plus speedup
   voor integratie en totaal ten opzichte van de mediaan van één worker.

Een SHA-256 en dirty-status maken een ongecommitteerde run herkenbaar, maar
vervangen het bewaren van de betreffende broncode niet. Noteer bij delen welke
build daadwerkelijk geladen was. De driver rapporteert waar mogelijk
bouwgegevens; ontbrekende Conda/CMakegegevens bij een geïnstalleerde package
worden niet verzonnen. De module-mtime is geen bewijs dat iedere bron actueel is.

Er is geen minimale speedup-assertie. Eén herhaling geeft geen nuttige spreiding;
gebruik voor een schaalclaim meerdere herhalingen en voldoende werk. De mediaan
verbergt uitschieters niet: alle losse tijden en de min/max blijven beschikbaar.
De machine kan nog steeds door ander werk of thermische effecten beïnvloed zijn.

## Eerste meting op 12 november 2026

Op dezelfde AMD Ryzen 9 7900X (24 logische CPU's), met de bestaande Release-build,
Python 3.14.6 en GCC 14.4.0, is de hierboven genoemde meetreeks van zes uur per
arc uitgevoerd: acht arcs, graad/orde 50, stap 30 s, drie herhalingen, steeds
1/2/4/8 workers in geroteerde volgorde. `OPENBLAS_NUM_THREADS` en
`OMP_NUM_THREADS` waren beide 1. De tweedaagse preset is niet uitgevoerd.

| Workers | Opbouwmediaan (s) | Integratiemediaan (s) | Totale mediaan (s) | Totale min–max (s) | Totale speedup |
|---:|---:|---:|---:|---:|---:|
| 1 | 0,01334 | 4,07902 | 4,09237 | 4,04569–4,11632 | 1,00× |
| 2 | 0,01359 | 2,09505 | 2,10864 | 2,08911–2,12312 | 1,94× |
| 4 | 0,01330 | 1,10106 | 1,11428 | 1,09283–1,13485 | 3,67× |
| 8 | 0,01339 | 0,74709 | 0,76062 | 0,74646–0,76274 | 5,38× |

De afzonderlijke fase-medianen hoeven niet exact op te tellen tot de mediaan
van het totaal; het totaal wordt per run gemeten. Iedere run controleerde alle
**5.768 arc-tijdstippen** en **32 globale matrixprobes**. De maximale gemeten
afwijkingen voor staten, STM's, sensitiviteiten en globale matrices waren alle
nul in deze meetreeks. Dat is geen belofte van bit-identieke resultaten op
andere platformen of bewijs van absolute baannauwkeurigheid.

Het kleine [JSONL-meetverslag](python-scaling-amd-20261112.jsonl) bevat alle
herhalingen en provenance. De driver was tijdens de meting nog untracked op
basis-HEAD `83c0f289d`; de vastgelegde driverhash hoort bij de broncode in deze
wijziging. Dit is één lokale Python-meetreeks, geen rechtstreeks vergelijk met
de oude native cijfers en geen garantie voor andere modellen of machines.

## Testen

```bash
PYTHONPATH=build-python-arcs/src python -m pytest \
  tests/test_tudatpy/test_python_arc_benchmark.py --no-remote-data -q
```

De tests controleren ongeldige workloads, geforceerde NaN/shape/epoch/numerieke
fouten, alle vergeleken resultaatcategorieën, mediaan/spreiding, bescherming van
bestaande output en een echte korte subprocess-run met beide native routes,
geroteerde workerstanden, provenance en globale-matrixcontroles. De bestaande
22 bindingtests blijven afzonderlijk van toepassing. Voor de passende bredere
offline/serialisatie-OFF-selectie, zie het API-document.

Op 12 november slaagden alle **28 drivertests**, de bestaande 22 bindingtests
en de drie native CTests. De passende bredere Python-selectie slaagde met
**680 passed, 104 skipped, 2 deselected**; de twee serialisatiemodules en de
twee bekende ongemarkeerde externe-downloadtests waren expliciet uitgesloten,
zoals beschreven in het API-document. Dit is geen nieuwe all-features/full-suite
claim. Het korte standaardvoorbeeld is niet gewijzigd.

Na de testselectiefix, later op 12 november, slaagde de **volledige** offline-run
zonder `--ignore` of naamfilters met **694 passed, 329 skipped**. Deze run omvat
ook elf nieuwe tests voor de expliciete serialisatie-bouwvlag en selectie.
Gebruik voortaan:

```bash
PYTHONPATH=build-python-arcs/src python -m pytest tests/test_tudatpy --no-remote-data -q
```

Serialisatie-afhankelijke tests worden nu alleen bij de expliciete OFF-bouwvlag
overgeslagen; de twee downloadtests respecteren de offline-optie. De onafhankelijke
`Time`-pickle-tests blijven actief. De historische handmatige uitsluitingen zijn
dus niet meer nodig. De aparte serialisatie-ON-build is vervolgens eveneens
[gevalideerd](../docs/independent_arc_development.md#validation-evidence-and-limits): 915 passed,
108 skipped, plus vier native serialisatiedoelen. Daarbij zijn verouderde fixtures
en waarde-gelijkheidsbugs gecorrigeerd. Dit betreft correctheidsvalidatie;
de schaalmeting hierboven blijft de eerdere OFF-meting.

De benchmarkdriver, tests en documentatie zijn AI-assisted. Review vooral de
baseline, meetgrenzen, invoerreceptuur en numerieke controles voordat resultaten
breder worden gebruikt. Achtergrond over herhaalde metingen:
[hyperfine](https://github.com/sharkdp/hyperfine) en
[asv](https://github.com/airspeed-velocity/asv). Hier volstaan de bestaande
NumPy-stack en Python-standaardbibliotheek; er is geen extra benchmarkframework
nodig.
