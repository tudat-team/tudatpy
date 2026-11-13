# Overdracht: parallellisatie van Tudat-arcs

Bijgewerkt op **13 november 2026** (Europe/Amsterdam), na de commits voor
serialisatie/testselectie en de Python-schaalbenchmark. De eerste Python-bindings
én de aanbevolen benchmark zijn **afgerond en gecommit**. De gebruiker heeft op
dezelfde AMD-machine verder gewerkt; dit document blijft bruikbaar voor een
volgende sessie op een andere machine.

De actuele benchmark en het meet-/validatiecontract staan in
[PYTHON_SCALING.md](../benchmarks/PYTHON_SCALING.md). Acht arcs van zes uur,
graad 50 en drie herhalingen geven lokaal **5,38× totale versnelling bij acht
workers**, met nul waargenomen verschillen ten opzichte van de sequentiële
referentie. De 28 drivertests slagen, opnieuw gecontroleerd op 13 november.
De tweedaagse preset is niet uitgevoerd. De CPU heeft **12 fysieke kernen en
24 logische processors**; de aanvankelijke aanname van 24 fysieke kernen was
onjuist. Een aanvullende gebruikersrun met 24 arcs schaalde het beste tot acht
workers (5,40×); de uitvoer blijft lokaal in
`local-profiles/20261112/python-scaling-24-arcs.jsonl`. De oorzaak van het afvlakken
is nog niet door profiling vastgesteld.

De volledige offline-suite is op 12 november zonder handmatige uitsluitingen
gevalideerd: **694 passed, 329 skipped** met de behouden serialisatie-OFF-binary
en **915 passed, 108 skipped** met een afzonderlijke serialisatie-ON-build.
Deze aantallen omvatten de 28 benchmarktests. Daarnaast slaagden 235 gerichte
serialisatie-/policytests (twee refactoringskips) en vier native
serialisatiedoelen. Zie [het ON-verslag](serialization_validation_20261112.md)
voor omgeving, opdrachten, gevonden fixtures en gelijkheidsbugs. De oude
uitgesloten testselecties verderop zijn uitsluitend historische informatie.

## Begin hier

1. Controleer branch `parallellization_inquiry`, de commits hieronder en lokale
   wijzigingen. Neem geen bestaande binaries, buildmappen of CMake-cache over.
2. Lees [python_independent_arcs_api.md](python_independent_arcs_api.md) voor het
   Python-contract en [independent_arc_execution.md](independent_arc_execution.md)
   voor de native C++-route. Het [Intel-reproductieverslag](parallel_arcs_reproduction_20261109.md)
   beschrijft de eerdere native controle op 9 november.
3. Doe op de andere machine eerst de bouw-/omgevingscontrole en de korte tests
   hieronder. Verifieer dat Python de nieuw gebouwde module importeert.
4. **Een volgende ontwikkeltaak is nog niet gekozen.** Bespreek die na de korte
   reproductie met de gebruiker. Bouw de bestaande Python-API niet opnieuw en
   start geen lange benchmark als voorwaarde voor dat gesprek.

De gebruiker werkt in het Nederlands, wil concrete metingen en afgebakende
commits, en richt zich op parallellisatie van onderliggende Tudat-routines.
Prioriteit is gegeven aan onafhankelijke arcs boven losse krachtberekeningen
binnen één arc. Prototype, native solver en Python-interface moeten duidelijk
onderscheiden blijven.

## Wat werkt nu

De bestaande C++-klassen `MultiArcDynamicsSimulator` en
`MultiArcVariationalEquationsSolver` hebben een opt-in route met native
**C++-threads** voor onafhankelijke arcs. Dit gebruikt geen Python
`multiprocessing`. Resultaten en globale matrixinterfaces blijven native Tudat.

Python biedt nu een bewust beperkte ingang:

- `simulator.MarsGravityArcModel(degree=..., order=...)`;
- `simulator.IndependentArcVariationalSimulator(...)`;
- één ondersteund model: Mars/JGMRO120D, Cowell, vaste RKDP87-stap, zes
  beginstaatparameters per arc en één gedeelde Mars-mu;
- eigen native modellen per arc, opgebouwd door een C++-factory;
- sequentiële default (`number_of_workers=1`), expliciet parallel opt-in;
- geen willekeurige Python-factory, Python-callbacks of toegang tot de
  onderliggende veranderlijke solver/bodies vanuit deze wrapper.

Het [uitvoerbare voorbeeld](../benchmarks/python_independent_arcs.py) heeft twee
arcs van elk **600 seconden**, stap **30 seconden**, zwaartekrachtgraad **20** en
**twee workers**. De gebruiker heeft het zelf uitgevoerd. Het is terecht zeer
snel klaar: dit is een functioneel voorbeeld, geen zware snelheidsbenchmark.

```bash
# Vanuit de repositoryroot, na een eigen verse build:
conda activate tudatpy-dev
PYTHONPATH=build-python-arcs/src python benchmarks/python_independent_arcs.py
```

Verwachte hoofduitvoer:

```text
Native independent execution: True
Arc 0: 21 saved states
Arc 1: 21 saved states
Global state-transition/sensitivity matrix: (6, 13)
```

## Git en overdracht

Deze twaalf commits staan op `parallellization_inquiry`, plus de actualisering
van dit overdrachtsdocument. Code-/benchmark-HEAD vóór deze documentatiecommit:
`bc682fefc`.

| Commit | Inhoud |
|---|---|
| `2873e55d4` | Los parallel C++-prototype en eerste controles/metingen |
| `5efa213d2` | Thread-safe initialisatie van Runge–Kutta-tabellen |
| `a7718ee1d` | Onafhankelijke arcs in de bestaande native multi-arc solvers |
| `fc9d93882` | Single-arc-profileringsscripts en LEO-analyse |
| `9a4c3365d` | Historische onderzoeksnotitie |
| `fb07784bf` | Bouwfix voor Python met serialisatie uitgeschakeld |
| `7688a026b` | Python-API, native Mars-model, voorbeeld, tests en API-contract |
| `990b7b4fa` | Overdracht en Intel-reproductieverslag in Git |
| `9be245b4c` | Actualisering overdracht voor de volgende machine |
| `83c0f289d` | Aanbeveling voor de Python-schaalbenchmark |
| `eba2dca47` | Serialisatievalidatie, gelijkheidsfixes en serialisatiecontract |
| `bc682fefc` | Python-schaalbenchmark, 28 drivertests en klein AMD-meetverslag |

Er is in deze sessie **niet gepusht**. Beschikbaarheid op een remote is niet
gecontroleerd. Zorg dus dat ook de laatste documentatiecommit wordt meegenomen.
Bijvoorbeeld via een bundle, na het committen van alle gewenste wijzigingen:

```bash
# Bronmachine; kies zelf een overdrachtspad:
git bundle create /gekozen/pad/tudat-parallel-arcs.bundle parallellization_inquiry
# Andere machine, in een nieuwe directory:
git clone -b parallellization_inquiry /gekozen/pad/tudat-parallel-arcs.bundle tudatpy
cd tudatpy
```

Een bundle bevat geen submodulerepositories, omgevingsinstallaties, ignored of
untracked bestanden. Regel submodules waar nodig afzonderlijk. Op 11 november
stond `examples/tudatpy` op `0113d35864519867da66acd0a3ca3a3f2f2e5d63`, zonder
lokale submodulewijziging. Een oudere werkmap had een afwijkende submodulestand;
neem die historische afwijking niet automatisch over en herstel geen bestaand
werk op de nieuwe machine zonder opdracht.

### Wat uitsluitend lokaal blijft

De gebruiker heeft expliciet gekozen de agent-tooling **alleen lokaal** te
houden. Deze paden zijn via `.git/info/exclude` uitgesloten:

```text
/.agents/
/.codex/
/AGENTS.md
/rseng-check/
/.rseng-agent-skills-*
/.rseng-agent-skills.codex.json
/local-profiles/
```

Die uitsluitingen zijn lokaal; `.gitignore` is niet gewijzigd. De nieuwe machine
krijgt de tooling en uitsluitingen niet automatisch. Gebruik daar de aanwezige
agentinstructies; voeg de lokale tooling niet alsnog aan Git toe.

Alle **24** oude profilerings-/omgevingsbestanden (**2,71 GiB**) zijn op verzoek
verplaatst naar `local-profiles/20261111/`. Bestandsnamen en submappen `Old/` en
`oud/` zijn behouden. Een `manifest.json` bevat oorspronkelijke paden, groottes,
mtime en SHA-256; controles vóór en na verplaatsing kwamen overeen. Een lokale
README licht de indeling toe. Niets is verwijderd. Ook `gprof-environment-dev.yaml`
en de drie `*-self.txt`-rapporten staan daar. De YAML beschrijft een historische
Python 3.11-omgeving en is niet de omgeving voor de nieuwe bindings.

Deze data is niet nodig voor de korte reproductie. Kopieer de map afzonderlijk
als de nieuwe sessie de oorspronkelijke profielen moet onderzoeken. Vlak na de
opruimcommit was `git status --short` leeg.

## Bouw-/omgevingscontrole op de andere machine

Gebruik een verse build en een eigen ontwikkelomgeving. `environment.yaml`
bevat de dependencies, maar is geen versielockbestand. Controleer onder andere
CMake, compilers, Eigen, Boost, **Cereal**, CSPICE, SOFA, NRLMSISE, Tudat Resources,
pybind11, NumPy, pytest en de gedeclareerde Python-dependencies.

```bash
# Alleen indien een geschikte omgeving nog niet bestaat:
conda env create -f environment.yaml
conda activate tudatpy-dev
```

De AMD-machine gebruikte op 11 november Python 3.14.6, GCC 14.4.0,
pybind11 3.0.3 en Tudat Resources 2.4 (`h68c723b_5`). Cereal 1.3.2 kwam hier uit
het reeds aanwezige `/tmp/tudat-cereal-prefix`; **neem dit tijdelijke pad niet
over**. Zorg dat Cereal op de nieuwe machine correct is geïnstalleerd.

De publieke package-import vereiste hier ook herstel van de al gedeclareerde
`atdf2ascii`-dependency: 1.0.0, upstream commit
`ae11cb1438d2843923272caee4eda69ea476fa93`, plus `bitstring` 5.0.0,
`termcolor` 3.3.0 en `tibs` 2.0.2. Geen bestaande pakketten zijn opgewaardeerd.
Zie het API-document voor het volledige lokale bouw-/omgevingsverslag.

De standaard SPICE-kernels en `gravity_models/Mars/jgmro120d.txt` moeten
beschikbaar zijn. Gebruik echte resourcebestanden; vervang ze niet door
verzonnen of vereenvoudigde coëfficiënten. Kopieer geen oude geïnstalleerde
`tudatpy-with-debug-symbols`-package als vervanging voor de nieuwe build.

Onderstaande compilerpaden zijn voor Linux met de gebruikte Conda-toolchain;
pas ze op een ander platform aan. Kies een andere nieuwe buildmap als
`build-python-arcs` daar al bestaat. Pas bouwjobs aan het geheugen aan; op de
Intel-machine met 16 GB waren twee jobs gebruikt.

```bash
cmake -S . -B build-python-arcs -G 'Unix Makefiles' \
  -DCMAKE_CXX_COMPILER="$CONDA_PREFIX/bin/x86_64-conda-linux-gnu-c++" \
  -DCMAKE_C_COMPILER="$CONDA_PREFIX/bin/x86_64-conda-linux-gnu-cc" \
  -DCMAKE_PREFIX_PATH="$CONDA_PREFIX" \
  -DCMAKE_BUILD_TYPE=Release -DTUDAT_BUILD_TESTS=ON \
  -DTUDAT_BUILD_WITH_SERIALIZATION=OFF \
  -DTUDAT_BUILD_WITH_EXTENDED_PRECISION_PROPAGATION_TOOLS=OFF \
  -DTUDAT_BUILD_EXPLICIT_INSTANTIATIONS=ON \
  -DCCACHE_PROGRAM=CCACHE_PROGRAM-NOTFOUND
cmake --build build-python-arcs --target kernel -j 2
PYTHONPATH=build-python-arcs/src python -c 'import tudatpy.kernel; print(tudatpy.kernel.__file__)'
PYTHONPATH=build-python-arcs/src python benchmarks/python_independent_arcs.py
PYTHONPATH=build-python-arcs/src python -m pytest \
  tests/test_tudatpy/test_independent_arc_execution.py --no-remote-data -q

cmake --build build-python-arcs --target copy_test_data \
  test_propagators_IndependentArcExecution test_integrators_RungeKuttaCoefficients \
  test_integrators_RungeKutta87DormandPrinceIntegrator -j 2
ctest --test-dir build-python-arcs \
  -R '^test_(propagators_IndependentArcExecution|integrators_RungeKuttaCoefficients|integrators_RungeKutta87DormandPrinceIntegrator)$' \
  --output-on-failure
```

De geprinte module moet uit de eigen build komen. De native acceptatietest
gebruikt `double` voor tijd/toestand; de echte Python-build gebruikt de
bindingsconfiguratie met `tudat::Time`. De native test alleen vervangt daarom
niet de Python-build en de 22 bindingtests.

## Bewijs en grenzen

De actuele resultaten van 12–13 november staan bovenaan; onderstaande
11-novemberresultaten beschrijven de toestand vóór de test-/serialisatiefix.

Op de **AMD-machine, 11 november 2026**:

- Verse Release-kernel, publiek Python-voorbeeld en **22 nieuwe Python-tests**
  geslaagd. Dekking: sequentieel/parallel, mu-updates, gekopieerde invoer,
  lifetime van native resultaten/interfaces, deferred reset, dynamics-only,
  gekoppelde herintegratie, validatie en workerfouten/herstel.
- Alle **drie native CTests** geslaagd.
- Volledige Python-run met `--no-remote-data`: **660 passed, 209 failed,
  114 skipped, 1 xfailed**. Daarvan vereisen 207 fouten de uitgeschakelde
  serialisatie; twee andere tests downloaden ongemarkeerd externe bestanden
  en faalden op sandbox-DNS. Dit is geen volledig groene all-features-suite.
- Passende selectie voor serialisatie-OFF/offline: **652 passed, 104 skipped,
  2 deselected**. Exacte uitsluitingen en commando staan in het API-document.
- Black, clang-format en whitespacecontrole geslaagd bij die eerdere controles.

Op de **Intel-machine, 9 november 2026**: drie native tests geslaagd, native
acceptatie 30 keer en coëfficiëntencontrole 20 keer achtereen geslaagd. De korte
schaalproef met vier arcs was te licht voor een overtuigende snelheidsclaim.
De Python-bindings bestonden toen nog niet; ze zijn daar nog niet getest.
Zie het afzonderlijke reproductieverslag voor details.

De eerdere native AMD-meting (acht arcs, graad 50, twee dagen per arc) gaf
31,070 / 15,615 / 7,977 / 4,717 s bij 1/2/4/8 workers (maximaal 6,59×).
Zie `benchmarks/cpp/native-scaling-forward.jsonl`. Dit is historische native
meetdata, **geen nieuwe meting van de Python-API**. Op 11 november is geen nieuwe
lange benchmark uitgevoerd. Gebruik de korte demonstratie niet als speeduptest.

Numerieke parity gebruikt de bestaande tolerantie
`abs(actual - expected) <= 1e-14 + 1e-11 * abs(expected)` en eindige waarden.
Dit toetst overeenkomst met de sequentiële route, geen absolute fysieke
baannauwkeurigheid. Er is geen volledige ThreadSanitizer-run of onafhankelijke
menselijke/scientifieke goedkeuring van de nieuwe bindings uitgevoerd.

## Ontwerpafspraken die behouden moeten blijven

- Native C++: `setIndependentArcExecution(workers, factory)`; de factory wordt
  **serieel vóór workers, opnieuw bij elke integratie** aangeroepen met actuele
  lokale parameters. Modellen/caches/integrators moeten echt geïsoleerd zijn;
  `SystemOfBodies` kopiëren is geen deep copy. Bekende aliases worden gecheckt,
  maar niet iedere willekeurige pointer-/callbackgraph.
- Sequentiële defaults, Cowell-beperking en native resultaatinterfaces blijven.
  Toestandsoverdracht/NaN-afhankelijkheid blijft voor de solver gelatcht op
  sequentiële uitvoering, ook na het oplossen van de NaNs. Andere native
  toestandsrepresentaties vallen terug; consider parameters worden afgewezen.
- C++-workers vullen bestaande arc-resultaten; na join volgt native globale
  verwerking van ephemeriden en STM/sensitiviteiten. Gedeelde mu blijft één
  globale parameter. Workerfouten worden na join met arc-index doorgegeven.
- Python: invoer/configuratie wordt gekopieerd; workerfactories bevatten geen
  geleende Python-objecten, Python-callbacks of verwijzing naar de facade.
  Native resultaatobjecten hebben live/update-semantiek, geen snapshots.
- **GIL blijft vastgehouden** tijdens Python-integratie. Native workers hebben
  hem niet nodig; Python-threads kunnen hierdoor niet gelijktijdig live
  resultaten lezen/wijzigen. Geen `gil_scoped_release` toevoegen zonder een
  apart synchronisatie-/snapshotontwerp voor alle toegankelijke objecten.
- Free-threaded Python is besproken, **niet geïmplementeerd of getest**. De
  C++-workers zijn al parallel. Zonder GIL is expliciete bescherming nodig,
  ook voor eerder verkregen native resultaten. Alleen een mutex om de facade
  of een modulevlag is onvoldoende; willekeurige callbacks blijven uitgesloten.
- `reintegrate=False` in de nieuwe facade stelt alle integratie uit. De oudere
  native `resetParameterEstimate(values, false)` doet juist dynamics-only;
  verwar die betekenissen niet. De facade synchroniseert parameters expliciet.
- De ondersteunde Mars-ephemeris-/rotatiedelegates gebruiken een gedeelde mutex.
  Dit is **geen algemene globale SPICE-lock** voor alle Tudat-code of modellen.
- De RK-bugfix initialiseert alle 21 tabellen éénmaal als `static const`. De
  oude DP87-getter testte orde 7 terwijl de initializer 8 instelde en herschreef
  daardoor gedeelde tabellen. De race is uit code-inspectie aangetoond; dat hij
  de oorspronkelijke incidentele afwijking veroorzaakte is zonder TSan niet bewezen.

## Codekaart en mogelijk vervolg

| Bestand | Rol |
|---|---|
| `src/tudatpy/dynamics/simulator/independent_arc_models.h` | Native model, eigendom, validatie en facade |
| `src/tudatpy/dynamics/simulator/expose_simulator_independent_arcs.cpp` | Python-bindings; bewust GIL vasthouden |
| `tests/test_tudatpy/test_independent_arc_execution.py` | 22 korte binding-/paritytests |
| `benchmarks/python_independent_arcs.py` | Kort werkend Python-voorbeeld |
| `include/tudat/simulation/propagation_setup/multiArcDynamicsSimulator.h` | Native scheduling, aliascontroles, fallback en resultaatverwerking |
| `include/tudat/simulation/estimation_setup/multiArcVariationalEquationsSolver.h` | Actuele lokale parameters en variatiefactory |
| `include/tudat/simulation/propagation_setup/singleArcDynamicsSimulator.h` | Schrijven naar native per-arc resultaten |
| `include/tudat/simulation/propagation_setup/independentArcExecutor.h` | Taakverdeling, join en fouten |
| `benchmarks/cpp/nativeMultiArc.cpp` | Native acceptatie/factory/schaalmeting |
| `benchmarks/cpp/marsArcBenchmark.h` | Native benchmarkmodel en SPICE-wrappers |
| `benchmarks/cpp/parallelMultiArc.cpp` | Eerder zelfstandig prototype; referentie behouden |
| `src/tudat/math/integrators/rungeKuttaCoefficients.cpp` | Constante RK-tabellen |

### Afgeronde benchmark en mijn aanbeveling voor vervolg

De eerdere aanbeveling om een Python-schaalbenchmark te bouwen is uitgevoerd
in `bc682fefc`. De driver ondersteunt instelbare arcs, graad/orde, duur, stap,
workers en herhalingen, een korte smoke-default en een expliciete zwaardere
preset. Eén worker gebruikt de echte sequentiële default. Iedere gemeten run
krijgt verse modellen en een facade; alle opgeslagen toestanden, STM's en
gevoeligheden worden met een sequentiële referentie vergeleken. De globale
matrixinterface wordt eveneens gecontroleerd. Opbouw, integratie en totale tijd
worden afzonderlijk gerapporteerd, met mediaan en spreiding. Het kleine
[AMD-meetverslag](../benchmarks/python-scaling-amd-20261112.jsonl) staat bewust
in Git; grotere uitvoer blijft lokaal.

Voor gebruik en exacte meetgrenzen geldt het benchmarkdocument. De historische
native tijden zijn niet rechtstreeks vergelijkbaar met deze Python-tijden:
onder meer tijdtype, opgeslagen afhankelijke variabelen, sequentiële route en
meetgrenzen verschillen.

**Mijn aanbeveling voor de volgende ontwikkeltaak is gericht onderzoeken waarom
de schaalbaarheid boven acht workers afvlakt.** Dit is een voorstel, nog geen
opdracht voor een nieuwe meting of codewijziging. Begin met de bestaande
24-arc-uitvoer en CPU-topologie. Profileer vervolgens, indien gekozen, een
beperkte run met acht en twaalf workers om seriële modelopbouw, workerwerk,
SPICE-synchronisatie en native nabewerking te onderscheiden. Deel geen caches
en versoepel geen isolatie- of GIL-afspraken om alleen de speedup te verhogen.
Kies pas een optimalisatie wanneer de meting het knelpunt aantoont.

Het serialisatiecontract staat in het API-document: integratie is synchroon;
na worker-join en native nabewerking mogen resultaten sequentieel worden
opgeslagen. Verkregen resultaatobjecten zijn live objecten, geen snapshots.
Gelijktijdig muteren en serialiseren van dezelfde objectgraph is niet
ondersteund. De workers gebruiken geen Cereal-archives; parallelle integratie
vereist daarom op zichzelf geen parallelle serialisatie. Checkpoints,
achtergrondserialisatie en snapshots vragen een afzonderlijk ontwerp.

Meer modellen, gerichte ThreadSanitizer-controle, andere platformen en generieke
Python-factories zijn afzonderlijke vervolgstappen. De oude sequentiële
`benchmarks/profile_multi_arc.py` activeert de nieuwe parallelle API niet.
[parallellization_inquiry.md](parallellization_inquiry.md) is historische context.

Lokale rseng-tooling meldde nog ontbrekende `CITATION.cff` en agent-ignore-regels
in de gedeelde `.gitignore`. Agent-tooling is inmiddels bewust lokaal uitgesloten;
projectbrede metadata/tooling is geen nieuw gekozen ontwikkeldoel. Code, tests
en documentatie zijn AI-assisted; menselijke review van modelbeperking, GIL,
lifetime en parameter-resetcontract blijft gewenst vóór bredere opname.

## Starttekst voor een volgende sessie

> Lees `docs/codex_handoff_parallel_arcs.md`, daarna
> `docs/python_independent_arcs_api.md`, `docs/independent_arc_execution.md`,
> `docs/serialization_validation_20261112.md` en `benchmarks/PYTHON_SCALING.md`.
> Controleer branch `parallellization_inquiry`, de twaalf genoemde commits en
> de laatste overdrachtscommit. Python-bindings, serialisatie-/testfixes en de
> schaalbenchmark zijn gecommit. Doe op een andere machine eerst een verse
> bouw-/omgevingscontrole en reproduceer het korte Python-voorbeeld, de 22
> bindingtests, de 28 benchmarktests en de drie native arc/RK-CTests.
> Voer de volledige offline-suite uit met `--no-remote-data`, zonder handmatige
> uitsluitingen; gebruik de expliciete serialisatiebouwvlag voor testselectie.
> Valideer ON en OFF in afzonderlijke buildmappen als beide nodig zijn.
> Behoud sequentiële defaults, modelisolatie, native resultaten, gelatchte
> fallback voor gekoppelde arcs en het GIL/callback-/serialisatiecontract.
> Een volgende ontwikkeltaak is nog niet gekozen. Bespreek de aanbeveling om
> het afvlakken boven acht workers gericht te profilen. Start geen lange
> benchmark automatisch. Agent-tooling en grote profieldata blijven lokaal;
> het kleine gekozen benchmarkverslag staat in Git. Communiceer in het
> Nederlands en maak afgebakende commits.
