# C++-profilering met een LEO-propagatie

Voor native multi-arc variatiepropagatie en de C++-taakgrens voor
arc-parallellisatie: zie [MULTI_ARC.md](MULTI_ARC.md).
Het [C++-prototype](cpp/README.md) voert onafhankelijke arcs met threads uit.

Voor de hoge-orde Mars-zwaartekrachtbenchmark uit het onderzoeksvoorstel:
zie [MARS_GRAVITY.md](MARS_GRAVITY.md) en `profile_mars_gravity.py`.

`profile_leo.py` draait één baanpropagatie, standaard met variatievergelijkingen.
Het model bevat een satelliet van 500 kg op circa 700 km hoogte, inclinatie 98°,
aardzwaartekracht tot graad/order 8, luchtweerstand, zonnestralingsdruk met
aardverduistering en puntmassazwaartekracht van zon, maan en Mars. De
variatievergelijkingen betreffen de zes beginstaatcomponenten, de
aardzwaartekrachtparameter en de luchtweerstandscoëfficiënt.

Dit is een representatieve workload voor baanpropagatie en sensitiviteiten;
geen universele benchmark voor alle Tudat-toepassingen. Observatiemodellen,
schatting, grote satellietconstellaties en hoge-orde zwaartekracht hebben andere
profielen. Er zijn geen Python-callbacks in de integratie.

## Uitvoeren en duur kiezen

Gebruik de Python-interpreter van de omgeving met de gewenste gecompileerde
Tudat-bibliotheek. Het script rapporteert het geïmporteerde pakketpad. Controleer
dat dit de build is die u wilt optimaliseren; een wijziging in deze checkout
verandert een elders geïnstalleerde `kernel.so` niet automatisch.

```bash
conda activate tudatpy-with-debug-symbols
export OPENBLAS_NUM_THREADS=1
export OMP_NUM_THREADS=1

# Korte controle van beide codepaden
python benchmarks/profile_leo.py --days 0.1 --mode dynamics
python benchmarks/profile_leo.py --days 0.1 --mode variational

# Kalibratie zonder profiler
python benchmarks/profile_leo.py --days 30

# Standaard: 120 dagen, stap 30 s, graad/order 8, variatievergelijkingen
python benchmarks/profile_leo.py
```

Streef naar 60–180 seconden **integratietijd** per profiel. Schaal daarvoor
`--days`: ongeveer `nieuwe_dagen = oude_dagen * gewenste_tijd / gemeten_tijd`.
De standaardduur is een startpunt, geen garantie. Geheugengebruik neemt eveneens
toe met de duur: de native solver bewaart toestanden en variatiematrices.
Laat opslaginstellingen tussen vergelijkingen gelijk.

Bij een lokale controle met de geïnstalleerde `tudatpy-with-debug-symbols`-omgeving
kostten 30 simulatiedagen circa 22,0 s met variatievergelijkingen en 13,1 s zonder.
Beide runs eindigden succesvol met dezelfde gerapporteerde eindstaat en 86.401
opgeslagen toestanden. Dit waren enkele kalibratieruns, geen herhaalde
prestatiemetingen. Een korte `perf`-proefopname slaagde eveneens; deze was te kort
voor conclusies over bottlenecks.

De JSON-uitvoer onderscheidt opbouwtijd, CPU-tijd en verstreken integratietijd.
De markers `INTEGRATION_BEGIN` en `INTEGRATION_END` omringen alleen de native
integratie. Resultaten worden pas daarna naar Python gekopieerd; het script
controleert succesvolle beëindiging, eindtijd en eindige eindstaat en print een
eindstaat voor vergelijkingen. Dat is geen bewijs van numerieke nauwkeurigheid.

Verhoog eerst de duur, niet de zwaartekrachtgraad. `--degree 20` kan later een
afzonderlijke, zwaardere zwaartekracht-workload vertegenwoordigen, maar verandert
de verdeling van het rekenwerk. Controleer bij een echte toepassing de
integratienauwkeurigheid door de stap te verkleinen. De vaste RKDP87-stap maakt
deze benchmark reproduceerbaar en vermijdt wisselende adaptieve stapkeuzes.

## Profiel opnemen en routines rangschikken

Gebruik een geoptimaliseerde build met debug-symbolen, bijvoorbeeld
`RelWithDebInfo`; een ongeoptimaliseerde Debug-build kan andere bottlenecks tonen.
Begin met 199 samples/s om overhead en bestandsgrootte te beperken.

```bash
perf record -o leo-var.perf.data -F 199 --call-graph dwarf,32768 -- \
    python benchmarks/profile_leo.py --days 120 --mode variational

# Platte ranglijst: eigen CPU-aandeel per routine, zonder callgraph
perf report -i leo-var.perf.data --stdio --no-children --call-graph none \
    --sort dso,symbol --percent-limit 0.5 > leo-var-self.txt

# Aanroepketens en inclusief aandeel; niet optellen over parent/child
perf report -i leo-var.perf.data --stdio --children --percent-limit 1 \
    > leo-var-callers.txt

# Interactief rapport; selecteer een hotspot en bekijk annotatie/aanroepers
perf report -i leo-var.perf.data

# Dezelfde workload zonder variatievergelijkingen
perf record -o leo-dyn.perf.data -F 199 --call-graph dwarf,32768 -- \
    python benchmarks/profile_leo.py --days 120 --mode dynamics
perf report -i leo-dyn.perf.data --stdio --no-children --call-graph none \
    --sort dso,symbol --percent-limit 0.5 > leo-dyn-self.txt
```

Deze commando's schrijven nieuwe bestanden en laten het bestaande `perf.data`
ongemoeid. `perf` meet het volledige proces, dus ook imports, opbouw en de
uiteindelijke export. Gebruik voor conclusies over integratie de native
integratie-aanroepketens en maak de integratiefase dominant met een langere run.
Voor exacte fase-isolatie is apart aan-/uitschakelen van de profiler nodig;
de printmarkers filteren de opname niet automatisch.

Een callgraph is een weergave van het profiel, geen vervanging van de platte
ranglijst. Gebruik `perf report` als uitgangspunt: percentages van gprof2dot
kunnen afwijken door filtering en de weging van samples met hun eventperiode.
Ook CPU-aandeel is geen aandeel van verstreken tijd, zeker bij meerdere threads.
Vergelijk voor variatiekosten absolute integratietijden en routineprofielen;
het aftrekken van twee procentkolommen geeft geen absolute kosten.

Controleer verloren samples en stackkwaliteit. Onmogelijke frames zoals
`[unknown] 5` maken de bovenliggende aanroepketen onbetrouwbaar. Een grotere
DWARF-stackbuffer kan truncatie verminderen maar herstelt geen ontbrekende of
onjuiste unwind-informatie. Eigen routinekosten kunnen nog bruikbaar zijn
wanneer de leaf-symbolen goed zijn opgelost. Frame-pointerprofilering is pas
een alternatief als de betrokken bibliotheken daarvoor zijn gebouwd.

Meet optimalisaties ook zonder profiler, meerdere malen achter elkaar op een
rustige machine. Houd model, stap, opslag, compilerinstellingen en threadaantallen
gelijk. Het script gebruikt standaard één BLAS/OpenMP-thread via bovenstaande
omgeving; het verandert die instellingen niet zelf. Kies later bewust andere
waarden wanneer u parallelle implementaties vergelijkt.
