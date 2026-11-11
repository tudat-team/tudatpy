# Overdracht: parallellisatie van Tudat-arcs

**Aanvulling 11 november 2026:** de eerste beperkte Python-API is uitgewerkt en
getest op de AMD-machine. Zie [python_independent_arcs_api.md](python_independent_arcs_api.md)
voor ontwerp, modelbeperkingen, lifetime/GIL/exception-afspraken, reproduceerbare
commando's en testresultaten. Op de branch `parallellization_inquiry` zijn nu
ook `fb07784bf` (bouwfix zonder serialisatie) en `7688a026b` (Python-API) gecommit.
De vijf basiscommits hieronder zijn behouden. De oorspronkelijke tekst
hieronder beschrijft de overdrachtsstand van 9 november en is historische context.

De overdracht en het Intel-reproductieverslag worden nu ook in Git opgenomen.
De oude profileringsbestanden zijn lokaal verzameld onder
`local-profiles/20261111/`, met behoud van de oorspronkelijke bestandsnamen en
de submappen `Old/` en `oud/`. Een lokaal manifest bevat bestandsgroottes en
SHA-256-controlesommen. De profileringsmap en agent-tooling zijn uitsluitend
via `.git/info/exclude` uitgesloten; die bestanden blijven op deze machine
aanwezig en worden niet met een clone of Git-bundle overgedragen.

Opgesteld op 9 november 2026 voor een nieuwe Codex-sessie op een andere machine.
Dit is projectcontext en een reproduceerbaar startpunt; controleer de actuele
repositorystatus voordat u verder werkt.

**Bijgewerkt na de sessie van 9 november 2026 op de Intel-machine. De
eerstvolgende taak is de native C++-parallellisatie toegankelijk maken vanuit
Python (Tudatpy). De gebruiker wil dit in een latere sessie op een andere
machine oppakken. Begin met een concreet Python-API-ontwerp voor één ondersteund
voorbeeld en werk daarna de bindings uit. De native acceptatietests zijn op de
Intel-machine inmiddels succesvol gereproduceerd; zie hieronder.**

## Doel en stand van zaken

Het onderzoeksvoorstel richt zich op zware baananalyses met veel ruimtevaartuigen,
lange tijdvakken en hoge nauwkeurigheid. We hebben prioriteit gegeven aan
parallellisatie **over onafhankelijke arcs**, boven parallellisatie van losse
krachtberekeningen binnen een arc.

Er is inmiddels een opt-in route in de bestaande C++-klassen
`MultiArcDynamicsSimulator` en `MultiArcVariationalEquationsSolver`. Deze gebruikt
threads voor complete arc-integraties en behoudt de native resultaatcontainers
en gezamenlijke resultaatverwerking. Het is meer dan het oorspronkelijke losse
prototype, maar nog geen generieke automatische modelkopie of Python-feature.

Belangrijkste grenzen:

- Standaard blijft uitvoering sequentieel.
- De parallelle route ondersteunt voorlopig Cowell-translatiepropagatie.
- De aanroeper levert een expliciete C++-factory voor geïsoleerde arc-modellen.
- Toestandsoverdracht/NaN-beginstaten houdt de solver conservatief sequentieel,
  ook bij herintegratie nadat de NaNs zijn ingevuld.
- Andere toestandsrepresentaties vallen terug op sequentiële uitvoering.
- Consider parameters worden bij configuratie van de variatiefactory afgewezen.
- Python-bindings voor deze configuratie zijn **nog niet toegevoegd**.
- Er is geen algemene thread-safety-garantie voor willekeurige Tudat-modellen.

## Git en overdracht naar de andere machine

Tak: `parallellization_inquiry` (let op deze exacte spelling).
Gecontroleerde basis-HEAD bij het schrijven: **`9a4c3365d`**.

| Commit | Inhoud |
|---|---|
| `2873e55d4` | Los parallel C++-prototype, controles en eerste meetresultaten |
| `5efa213d2` | Zelfstandige thread-safety-bugfix voor Runge–Kutta-tabellen |
| `a7718ee1d` | Integratie in de native C++-multi-arc solvers |
| `fc9d93882` | Single-arc-profileringsscripts en LEO-analyse |
| `9a4c3365d` | Historische onderzoeksnotitie met verwijzing naar huidige implementatie |

Controleer of deze commits op de andere machine beschikbaar zijn. De
beschikbaarheid op een remote is niet gecontroleerd; er is hier geen push gedaan.
Dit overdrachtsdocument is na bovenstaande HEAD aangemaakt en moet zelf ook
worden meegenomen (los kopiëren of eerst committen).
Hetzelfde geldt voor `docs/parallel_arcs_reproduction_20261109.md`: beide
documenten zijn aan het einde van deze sessie nog untracked. Er is in deze
sessie geen commit of push gedaan en geen wetenschappelijke broncode gewijzigd.

Bijvoorbeeld: maak desgewenst een bundle van de tak nadat de gewenste documenten
zijn gecommit. Dit kan een groot bestand worden, omdat het bereikbare historie
bevat. Gebruik een gekozen overdrachtspad buiten de bronbestanden:

```bash
git bundle create /gekozen/pad/tudat-parallel-arcs.bundle parallellization_inquiry
# Op de andere machine, in een nieuwe directory:
git clone -b parallellization_inquiry /gekozen/pad/tudat-parallel-arcs.bundle tudatpy
```

Een bundle bevat geen untracked bestanden, dependency-installaties of
submodulerepositories. Gebruik waar nodig de normale repository/submodule-setup;
de afwijkende lokale submodulestand hieronder hoort niet bij deze feature.

De oude werkmap had daarnaast bestaande, bewust niet opgenomen wijzigingen:

- `examples/tudatpy`: submodule van `0113d3586…` naar `8227687d5…`; deze wijziging
  is niet door dit werk gemaakt en mag niet stilzwijgend worden teruggedraaid.
- Untracked profielbestanden (`perf.data*`, `*.perf.data`, callgraphs, rapporten),
  `Old/`, `oud/` en `gprof-environment-dev.yaml`. Ze zijn behouden; ze zijn niet
  nodig voor de native acceptatietest en horen niet automatisch in een commit.
- Alle eerdere benchmarkbestanden onder `benchmarks/` zijn inmiddels gecommit.

## Eerst lezen en belangrijkste bestanden

Lees [independent_arc_execution.md](independent_arc_execution.md) voor het
actuele API-contract. [parallellization_inquiry.md](parallellization_inquiry.md)
is uitsluitend historische context, geen actuele implementatiestatus.

| Bestand | Rol |
|---|---|
| `include/tudat/simulation/propagation_setup/multiArcDynamicsSimulator.h` | Factory, aliascontroles, scheduling, fallback en native resultaatverwerking |
| `include/tudat/simulation/estimation_setup/multiArcVariationalEquationsSolver.h` | Variatiefactory, actuele parameterwaarden en layoutcontroles |
| `include/tudat/simulation/propagation_setup/singleArcDynamicsSimulator.h` | `integrateIndependentArc`: schrijven naar caller-owned arc-resultaten |
| `include/tudat/simulation/propagation_setup/independentArcExecutor.h` | Thread-taakverdeler, join en deterministische fouttoerekening |
| `src/tudat/math/integrators/rungeKuttaCoefficients.cpp` | Eenmalige, constante Runge–Kutta-tabellen |
| `tests/test_tudat/src/math/integrators/unitTestRungeKuttaCoefficients.cpp` | Zelfstandige concurrentietest, zonder arc-executor |
| `benchmarks/cpp/nativeMultiArc.cpp` | Native acceptatietest, model-factory en schaalmeting |
| `benchmarks/cpp/marsArcBenchmark.h` | Gedeelde Mars-opbouw en SPICE-wrappers voor de benchmarks |
| `benchmarks/cpp/parallelMultiArc.cpp` | Oorspronkelijk los prototype; bewaren als referentie |
| `benchmarks/profile_multi_arc.py` | Sequentiële native multi-arc Python-benchmark |
| `benchmarks/cpp/native-scaling-forward.jsonl` | Vastgelegde native schaalmeting |

De overige profileringsscripts en uitleg staan in
[benchmarks/README.md](../benchmarks/README.md),
[MARS_GRAVITY.md](../benchmarks/MARS_GRAVITY.md) en
[MULTI_ARC.md](../benchmarks/MULTI_ARC.md).

## Ontwerp dat behouden moet blijven

De C++-ingang is `setIndependentArcExecution(workers, factory)`, op beide
multi-arc klassen. Configureer vóór integratie en schakel integratie bij
constructie uit. `clearIndependentArcExecution()` herstelt de gewone route.

De factory wordt **serieel, op iedere integratie opnieuw** aangeroepen, vóór
workerthreads starten. De variatiefactory krijgt arc-index, oorspronkelijke
arc-instellingen en actuele lokale parameterwaarden. Zo worden gewijzigde
beginstaten en gedeelde parameters bij herintegratie meegenomen.

De factory moet dezelfde fysica en instellingen reconstrueren met eigen
veranderlijke bodies, modellen, caches en parameterobjecten. Een kopie van
`SystemOfBodies` is geen deep copy. Bekende aliases worden afgewezen, maar
controle van een willekeurige volledige pointer-/callbackgraph is niet mogelijk.

De Mars-benchmark beschermt SPICE-ephemeris- en rotatieaanroepen met één gedeelde
mutex. Die wrappers staan in benchmarkcode: Tudat heeft hiermee **geen algemene
globale SPICE-lock** gekregen. Voor een andere factory moet externe gedeelde
toestand opnieuw correct worden behandeld.

Workers schrijven rechtstreeks in de native per-arc resultaten. Na join volgt
de bestaande sequentiële nabewerking: ephemeriden, afhankelijke variabelen en
gezamenlijke STM/sensitiviteitsinterface. Er is geen handmatig samengevoegde,
dichte globale matrixhistorie zoals in het eerste prototype. De native globale
parameterindeling koppelt alle lokale mu-sensitiviteiten aan dezelfde parameter.

Exceptions worden na join met arc-index gemeld; ook niet-succesvolle terminatie
geldt als fout. Bij fouten volgt geen globale nabewerking/succesmarkering.
Oorspronkelijke derivative-model-scratchcaches zijn geen snapshots van workers;
gebruik de resultaatobjecten voor histories en evaluatietellingen.

## Bekende bugfix en bewijsgrenzen

Een incidentele afwijking in de dynamica-only-test leidde tot inspectie van
gedeelde integratortoestand. De DP87-getter controleerde `higherOrder != 7`,
terwijl de initializer orde 8 instelde. Iedere lookup herschreef daardoor een
statische tabel, mogelijk terwijl een andere thread die kopieerde.

De getter gebruikt nu voor alle 21 methoden function-local `static const`
tabellen met C++ thread-safe initialisatie. De numerieke initializers en
coëfficiënten zijn niet gewijzigd. De code-inspectie toont de race; zonder
ThreadSanitizer is niet bewezen dat juist die race de eerste incidentele
afwijking veroorzaakte. Maak dat onderscheid in toekomstige verslaglegging.

## Verse bouw op een andere machine

De onderstaande shellcommando's zijn voor Linux/POSIX. Bouw opnieuw; neem geen
`CMakeCache.txt`, statische libraries of uitvoerbare bestanden van de oude
machine over. De oude map `build/` bevatte libraries van augustus bij headers
van november en bleek niet compatibel.

De repository heeft `environment.yaml` met onder meer CMake, Make, compilers,
Eigen, Boost, Cereal, CSPICE, SOFA, NRLMSISE, Tudat Resources en pybind11. Dit is
geen versielockbestand. Een mogelijke basis is:

```bash
conda env create -f environment.yaml
conda activate tudatpy-dev
```

Gebruik een andere omgevingsnaam als deze al bestaat. Controleer de gekozen
compiler en dependency-paden in de CMake-uitvoer. Tudat Resources moet de
standaard SPICE-kernels en het JGMRO120D Mars-model kunnen vinden. Vul ontbrekende
wetenschappelijke data niet aan met verzonnen/vereenvoudigde coëfficiënten.

```bash
# Vanuit de repositoryroot; gebruik op deze machine een NIEUWE buildmap.
cmake -S . -B build-parallel-arcs -G "Unix Makefiles" \
    -DCMAKE_PREFIX_PATH="$CONDA_PREFIX" \
    -DCMAKE_BUILD_TYPE=Release \
    -DTUDAT_BUILD_TESTS=ON \
    -DTUDAT_BUILD_WITH_SERIALIZATION=OFF \
    -DTUDAT_BUILD_WITH_EXTENDED_PRECISION_PROPAGATION_TOOLS=OFF \
    -DTUDAT_BUILD_EXPLICIT_INSTANTIATIONS=ON

cmake --build build-parallel-arcs --target copy_test_data \
    test_propagators_IndependentArcExecution \
    test_integrators_RungeKuttaCoefficients \
    test_integrators_RungeKutta87DormandPrinceIntegrator -j 4
```

Pas het aantal bouwjobs aan het beschikbare geheugen aan; templates compileren
zwaar. Dit bouwt C++-doelen, geen nieuwe geïnstalleerde Python-module. Het
native acceptatieprogramma gebruikt `double` voor tijd en toestand.

Oude lokale bijzonderheden die niet naar de andere machine mogen worden gekopieerd:

- Compiler/dependencies: `/home/hanno/miniforge3/envs/tudatpy-dev`.
- Cereal ontbrak in die installatie en werd tijdelijk onder
  `/tmp/tudat-cereal-prefix` geïnstalleerd. Een nieuwe correcte omgeving kan
  Cereal rechtstreeks leveren; kopieer die `/tmp`-prefix niet in nieuwe commando's.
- De Python-profielen gebruikten daarnaast `tudatpy-with-debug-symbols` met een
  oudere geïnstalleerde interface. Dat pakket was niet deze nieuwe C++-build.
- `build_prototype.py` leest `flags.make`/`link.txt` uit een Makefiles-build.
  Gebruik bij voorkeur het native CMake-testdoel hierboven; die helper is geen
  universeel Ninja/Windows-bouwsysteem.

## Tests en meting reproduceren

```bash
export OPENBLAS_NUM_THREADS=1
export OMP_NUM_THREADS=1

ctest --test-dir build-parallel-arcs --output-on-failure \
    -R '^(test_propagators_IndependentArcExecution|test_integrators_RungeKuttaCoefficients|test_integrators_RungeKutta87DormandPrinceIntegrator)$'

ctest --test-dir build-parallel-arcs --output-on-failure \
    -R '^test_propagators_IndependentArcExecution$' --repeat until-fail:30

ctest --test-dir build-parallel-arcs --output-on-failure \
    -R '^test_integrators_RungeKuttaCoefficients$' --repeat until-fail:20

# Eerst een korte functionele schaalproef, daarna pas de langere meting.
build-parallel-arcs/tests/test_propagators_IndependentArcExecution \
    --arcs 4 --degree 20 --days-per-arc 0.02
build-parallel-arcs/tests/test_propagators_IndependentArcExecution \
    --arcs 8 --degree 50 --days-per-arc 2
```

De executable voert eerst de korte regressiecontroles uit en meet daarna
hetzelfde werk met 1/2/4/8 workers. Vergelijk met de sequentiële native referentie,
niet alleen met een tweede parallelle run. Stop compilaties tijdens tijdmetingen.
Sla nieuwe machine-uitkomsten onder een nieuwe bestandsnaam op; overschrijf de
historische meetbestanden niet automatisch.

Getest op de oude machine: AMD Ryzen 9 7900X, 12 cores/24 logische CPU's,
x86_64 Linux, conda-forge GCC 14.4.0, Release. Native meting:

| Workers | Totale tijd | Versnelling t.o.v. 1 |
|---:|---:|---:|
| 1 | 31,070 s | 1,00 |
| 2 | 15,615 s | 1,99 |
| 4 | 7,977 s | 3,90 |
| 8 | 4,717 s | 6,59 |

Dit is één lokale native schaalmeting, geen universele garantie. Inclusief
native opbouw, factory-opbouw en nabewerking; exclusief validatie en destructie.
Alle 46.088 opgeslagen arc-tijdstippen zijn numeriek vergeleken. De tests dekken
ook gedeelde parameterupdates, resultaatidentiteit, ephemeriden, interpolatoren,
fouten/herstel, dynamics-only en gekoppelde herintegratie. Er is geen volledige
ThreadSanitizer-run van Tudat gedaan. De afzonderlijke executor is met
AddressSanitizer/UBSan getest, met leak detection uit wegens ptrace-beperkingen.

## Reproductie op de Intel-machine: voltooid op 9 november 2026

De huidige sessie gebruikte een Intel Core i7-8550U (4 cores / 8 logische CPU's),
16 GB RAM en de Conda-omgeving `tudatpy-dev`: GCC 14.4.0, Eigen 5.0.1,
Boost 1.85.0, CSPICE 67 en Tudat Resources 2.4. Er is een verse Release-build
gemaakt in `build-parallel-arcs-local`, met twee bouwjobs en dezelfde opties als
hierboven. Cereal ontbrak; uitsluitend Cereal 1.3.2 is toegevoegd, zonder andere
pakketten bij te werken. De benodigde SPICE-kernels en JGMRO120D waren aanwezig.

- Alle drie geselecteerde native C++-tests slaagden.
- De native multi-arc-acceptatietest slaagde 30 opeenvolgende herhalingen.
- De Runge–Kutta-coëfficiëntentest slaagde 20 opeenvolgende herhalingen.
- De korte proef met vier arcs, graad 20 en 0,02 dag per arc valideerde alle
  236 opgeslagen tijdstippen bij 1/2/4/8 ingestelde workers tegen de native
  sequentiële referentie. Totale tijden: respectievelijk 0,2351 / 0,2816 /
  0,2197 / 0,2230 seconde. Dit is te weinig werk voor een overtuigende
  schaalbaarheidsclaim; met vier arcs zijn maximaal vier workers actief.

De bestaande toleranties zijn ongewijzigd. De multi-arc-vergelijking controleert
per element `abs(actual - expected) <= 1e-14 + 1e-11 * abs(expected)` en wijst
niet-eindige waarden af. Dit toetst overeenkomst met de sequentiële referentie,
niet de absolute fysieke baannauwkeurigheid. Maximale afwijkingen worden bij
succes niet gerapporteerd.

De lange schaalmeting, volledige Tudat-testsuite en ThreadSanitizer zijn hier
niet uitgevoerd. De build gaf een waarschuwing over een mogelijk dangling
reference in `observationOutput.cpp:542`; de benchmark meldde deprecated gebruik
van `getNumericalVariationalEquationsSolution`. De tests slaagden desondanks.

Zie [het reproductieverslag](parallel_arcs_reproduction_20261109.md) voor
commando's en beperkingen. Logs en pakketmanifest staan lokaal onder
`build-parallel-arcs-local/reproduction/` (gitignored, niet automatisch
overgedragen). Neem geen binaries of CMake-cache over naar een andere machine.

De toegevoegde rseng-praktijkregels staan in `AGENTS.md` en `.agents/skills/`;
deze waren lokaal untracked, evenals `rseng-check/` en agentconfiguratie.
Fasedekking staat in `.rseng-agent-skills-coverage.md`. De zelfcontrole meldt
nog ontbrekende `CITATION.cff` en agent-directoryregels in `.gitignore`; dit is
vastgelegd en niet als onderdeel van de testreproductie aangepast.

## Voorgesteld vervolg voor de nieuwe sessie

1. **Eerstvolgende ontwikkeltaak: ontsluit de native C++-parallellisatie vanuit
   Python.** De gebruiker heeft deze vervolgstap gekozen, maar uitgesteld tot
   een latere sessie op een andere machine. Werk eerst het Python-API-contract
   concreet uit aan de hand van één ondersteund voorbeeld: hoe levert de
   gebruiker geïsoleerde arc-modellen aan en wie bezit die objecten? Een
   worker-aantal zonder modelisolatie is geen veilige uitbreiding.
2. Controleer bij de start branch/commits en lokale wijzigingen. Bouw op een
   nieuwe machine zo nodig opnieuw en voer de korte native tests uit als
   omgevingscontrole. De reproductie op de Intel-machine is al voltooid;
   een lange schaalmeting is geen voorwaarde voor het Python-API-ontwerp.
3. Werk vervolgens de Python-bindings uit binnen het afgesproken contract.
   Behandel factory-lifetime, GIL, exceptions en eventuele Python-callbacks in
   modellen expliciet; geef alleen C++-workerwerk vrij waarvan de voorwaarden
   duidelijk zijn. De oude Python-benchmark activeert nog geen parallellisatie.
4. Overweeg meer modellen/platforms en een gerichte ThreadSanitizer-build vóór
   bredere thread-safety-claims. Uitbreiding buiten Cowell/zonder factory en
   productie-integratie van SPICE-synchronisatie zijn afzonderlijke ontwerpstappen.

Gebruikersvoorkeuren uit deze sessie: communiceer in het Nederlands, werk met
concrete metingen, houd prototype en native implementatie duidelijk uit elkaar,
en maak afgebakende commits. Commit/push geen bestaande submodulewijziging,
profieldata of andere lokale bestanden zonder passende opdracht.

### Starttekst om in de nieuwe sessie te plakken

> Lees `docs/codex_handoff_parallel_arcs.md` en daarna
> `docs/independent_arc_execution.md`. De eerstvolgende taak is de native
> C++-parallellisatie toegankelijk maken vanuit Python. Begin met een concreet
> API-ontwerp en één ondersteund voorbeeld voor geïsoleerde arc-modellen; werk
> daarna de bindings uit, met expliciete afspraken over lifetime, GIL,
> exceptions en Python-callbacks. Controleer de branch `parallellization_inquiry`
> en de vijf genoemde commits. De native acceptatietests zijn op 9 november
> 2026 op de Intel-machine succesvol gereproduceerd; lees ook
> `docs/parallel_arcs_reproduction_20261109.md`. Doe op een nieuwe machine de
> benodigde bouw-/omgevingscontrole. Behoud sequentiële defaults,
> isolatievoorwaarden, native resultaatinterfaces en sequentiële fallback voor
> gekoppelde arcs. Voer geen lange benchmark uit als voorwaarde voor het
> API-ontwerp.
