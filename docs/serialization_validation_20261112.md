# Validatie van de serialisatie-ON-build — 12 november 2026

Op verzoek is op dezelfde AMD-machine een afzonderlijke, verse build gemaakt
met `TUDAT_BUILD_WITH_SERIALIZATION=ON`. Doel is de werkelijk gecompileerde
Cereal-/Python-roundtrips te testen, naast de al geteste selectielogica voor
optionele features. De bestaande OFF-build staat in `build-python-arcs` en wordt
niet vervangen of geïnstalleerd.

## Omgeving en bouw

- AMD Ryzen 9 7900X, 12 fysieke kernen / 24 logische CPU's, circa 64 GB RAM.
- Bestaande Conda-omgeving `tudatpy-dev`: Python 3.14.6, GCC 14.4.0,
  pybind11 3.0.3 en dezelfde bestaande wetenschappelijke dependencies/resources.
- Cereal 1.3.2 via de reeds aanwezige lokale prefix `/tmp/tudat-cereal-prefix`.
  Op een andere machine moet Cereal correct geïnstalleerd zijn; dit tijdelijke
  pad is geen portable installatie-instructie.
- Verse Release-buildmap `build-python-arcs-serialization-on`; expliciete
  instantiaties aan, uitgebreide propagatieprecisie uit, tests aan.
- Geen afhankelijkheden opgewaardeerd, geen installatie over de OFF-package.
  De numerieke tolerantie is niet aangepast.
- Basis-HEAD `83c0f289d`, met de ongecommitteerde benchmark- en testselectiefix
  uit deze sessie. De interne boolean `_build_with_serialization` komt uit
  de werkelijk geladen kernel. Gebruik een verse Python-processaanroep met
  het juiste `PYTHONPATH`; laad de twee kernels niet samen in één interpreter.

De gebruikte configuratie op deze machine:

```bash
cmake -S . -B build-python-arcs-serialization-on -G 'Unix Makefiles' \
  -DCMAKE_CXX_COMPILER="$CONDA_PREFIX/bin/x86_64-conda-linux-gnu-c++" \
  -DCMAKE_C_COMPILER="$CONDA_PREFIX/bin/x86_64-conda-linux-gnu-cc" \
  -DCMAKE_PREFIX_PATH="$CONDA_PREFIX;/tmp/tudat-cereal-prefix" \
  -DCMAKE_BUILD_TYPE=Release -DTUDAT_BUILD_TESTS=ON \
  -DTUDAT_BUILD_WITH_SERIALIZATION=ON \
  -DTUDAT_BUILD_WITH_EXTENDED_PRECISION_PROPAGATION_TOOLS=OFF \
  -DTUDAT_BUILD_EXPLICIT_INSTANTIATIONS=ON \
  -DCCACHE_PROGRAM=CCACHE_PROGRAM-NOTFOUND
cmake --build build-python-arcs-serialization-on --target kernel -j 4
```

Pas compiler-/dependency-paden en het aantal jobs aan het eigen systeem aan.
Alle libraries en de echte Python-module worden voor deze configuratie opnieuw
gebouwd; de OFF-libraries worden niet gekopieerd of in de ON-build gemengd.

## Python-controles

Controleer expliciet de geladen module en de gecompileerde featurevlag:

```bash
PYTHONPATH=build-python-arcs-serialization-on/src python -c \
  'from tudatpy import kernel; assert kernel._build_with_serialization is True; print(kernel.__file__)'
```

De vlag moet `True` zijn en het modulepad moet uit de ON-build komen. De gerichte
roundtrip-/selectietests en de volledige offline-suite draaien vervolgens met:

```bash
PYTHONPATH=build-python-arcs-serialization-on/src python -m pytest \
  tests/test_tudatpy/test_serialization_file_io.py \
  tests/test_tudatpy/test_serialization_pickle.py \
  tests/test_tudatpy/test_serialization_build_policy.py --no-remote-data -q
PYTHONPATH=build-python-arcs-serialization-on/src python -m pytest \
  tests/test_tudatpy --no-remote-data -q
PYTHONPATH=build-python-arcs-serialization-on/src python benchmarks/python_independent_arcs.py
```

Deze tests controleren JSON-, binaire en Python-pickle-roundtrips van beschikbare
instellingen en resultaatobjecten, gelijkheid/typebehoud, polymorfe dispatch,
verkeerde loader-types en verwachte fouten bij niet-serialiseerbare functies.
De policytest met een kunstmatige ON/OFF-vlag is op zichzelf geen bewijs van
echte serialisatie; bovenstaande aanroepen gebruiken de ON-binary.

### Bevindingen tijdens de echte ON-controle

De eerste bestaande testselectie gaf 226 passed, 10 skipped en 1 xpassed.
Inspectie leidde tot de volgende gerichte correcties:

- Acht tests zochten een los Python-package `spice`, hoewel ze Tudats eigen
  SPICE-module gebruiken. Die tests bevatten ook oude body-/variatiesolver-API's.
  Een gedeelde fixture bouwt nu via de actuele publieke API een korte echte
  propagatie, eventueel met variatievergelijkingen of hybride terminatie.
  Deze tests worden daadwerkelijk uitgevoerd en slaan geen ontbrekende externe
  `spice`-module meer als reden aan. Ze slaan ook een relatieve positie op.
- De mode-coupled Love-number-roundtrip slaagt. De achterhaalde `xfail` is
  verwijderd; de oorspronkelijke waarde-/typeasserties blijven staan.
- Polymorfe filetests schrijven nu expliciet via dezelfde base-class-root als
  waarmee ze lezen. De root-type-header controleert dat terecht; schrijven via
  een afgeleid type en lezen via een andere root is geen geldige testaanroep.
- Sommige C++-subtypen hebben geen eigen Python-binding. Een base-loader kan
  daarom een andere Python-fallbackwrapper kiezen dan de originele factory.
  De helpers controleren de virtuele C++-gelijkheid in beide richtingen,
  zodat echte slicing/waardeverlies niet via een oppervlakkige base-vergelijking
  wordt gemist. Expliciet gebonden resultaattypen houden daarnaast strikte
  Python-type- en historycontroles.
- De geactiveerde resultaat-roundtrips vonden een echte fout in
  `SingleArcSimulationResults::equals`: `outputSettings_` en
  `dependentVariableInterface_` werden op pointeridentiteit vergeleken.
  Deserialisatie maakt nieuwe objecten met nieuwe adressen. Beide vergelijkingen
  gebruiken nu de onderliggende waarden, met behoud van null-afhandeling,
  overeenkomstig de bestaande vergelijkingen voor multi-/hybrid-arc objecten.
  Dit verandert geen opgeslagen schema, krachtmodel, integrator of tolerantie.
- De daaropvolgende waarde-vergelijking maakte nog twee configuratieproblemen
  zichtbaar: de standaard print-/opslagintervallen gebruiken NaN voor
  "niet ingesteld", maar hun gelijkheidsoperatoren vergeleken die met gewone
  floating-point `==`. Alleen deze twee sentinelvelden behandelen nu twee
  ontbrekende intervallen als gelijk. De native regressiecase
  `test_DefaultProcessingSettingsEquality` controleert reflexiviteit, twee
  defaults en beide richtingen van default versus een eindig interval van 30 s.
  Numerieke histories en de NaN/Inf-validatie daarvan worden niet versoepeld.

De resultaat-tests blijven objectgelijkheid eisen en vergelijken daarnaast
typen en opgeslagen numerieke histories. Alleen de twee bestaande
observatiecollectie-filetests onder refactoring blijven expliciet overgeslagen.

## Native serialisatiecontroles

```bash
cmake --build build-python-arcs-serialization-on --target copy_test_data \
  test_io_SettingsSerialization test_io_ObservationCollectionSerialization \
  test_io_PropagationResultsSerialization test_io_PodInputOutputTypesSerialization -j 2
ctest --test-dir build-python-arcs-serialization-on \
  -R '^test_io_(SettingsSerialization|ObservationCollectionSerialization|PropagationResultsSerialization|PodInputOutputTypesSerialization)$' \
  --output-on-failure
```

Alle vier native CTests zijn geslaagd. Zij dekken instellingen,
observatiecollecties, propagatieresultaten en schattingsinvoer/-uitvoer. Dit is
geen volledige native Tudat-testsuite.

## Eindresultaten

Na de bovenstaande test-/gelijkheidsfixes en herbouw:

| Controle | Resultaat |
|---|---|
| Gecompileerde ON-vlag en importpad | `True`, module uit de aparte ON-build |
| Gerichte Python-serialisatie-/policytests | **235 passed, 2 skipped** |
| Volledige Python-offline-suite, ON | **915 passed, 108 skipped** |
| Vier native serialisatie-CTest-doelen | **4/4 passed**, inclusief nieuwe default-/sentinelregressie |
| Kort Python-voorbeeld, ON | Geslaagd; twee onafhankelijke arcs, 21 toestanden per arc, matrix 6×13 |
| Volledige Python-offline-suite, behouden OFF-binary | **694 passed, 329 skipped** |
| OFF-binary SHA-256, grootte en mtime | Ongewijzigd |

De volledige offline-runs gebruiken geen `--ignore` of naamfilters. Er zijn geen
onverwachte passes of nieuwe xfail-markeringen. De twee gerichte skips betreffen
de bestaande `SingleObservationSet`-/`ObservationCollection`-filetests onder
refactoring; het C++-observatiecollectiedoel is wel getest. De overige skips in
de volledige suite volgen de bestaande netwerk-/optionele-testvoorwaarden.
Bestaande deprecation-waarschuwingen blijven aanwezig.

De volledige Python-suitetotalen bevatten ook de 28 schaalbenchmarktests die
tijdens de validatie al in de werkboom stonden. Die benchmarkwijzigingen krijgen
een aparte commit; de gerichte 235 tests en vier native doelen horen volledig
bij deze serialisatie-/testfix.

De behouden OFF-binary is bewust niet overschreven of herbouwd. De nieuwe
C++-waarde-gelijkheidsfixes zijn hier in de ON-build gecompileerd en gevalideerd;
een toekomstige OFF-herbouw zal die eveneens opnemen. Geen nieuwe lange
benchmark of ON/OFF-prestatievergelijking is uitgevoerd.

## Vastlegging en bewijsgrenzen

### Contract met de parallelle arc-uitvoering

De ondersteunde volgorde is **succesvol integreren, alle workers afronden en
joinen, daarna de voltooide resultaten sequentieel opslaan**. De workers
gebruiken geen Cereal-archives en modelisolatie gebeurt via de C++-factory,
niet door modellen met Cereal te kopiëren. `TUDAT_BUILD_WITH_SERIALIZATION=ON`
schakelt I/O in; het is geen verklaring dat willekeurige gelijktijdige
serialisatie veilig is.

Een resultaatreferentie blijft een live object, geen snapshot. Tijdens het
doorlopen van een objectgraaf voor opslag mag die niet door integratie,
parameterupdates of andere mutaties veranderen. Na een mislukte integratie
mag een bewaarde referentie niet als een compleet nieuw resultaat worden
beschouwd. De huidige GIL voorkomt zulke gelijktijdige toegang door gewone
Python-threads tijdens de integratie-aanroep. Externe C++-threads en echt
free-threaded Python worden daarmee niet beschermd.

De hier geteste Cereal-configuratie gebruikt de standaard
`CEREAL_THREAD_SAFE=0`; er is geen gelijktijdige archieftoegang nodig. Opslaan
tijdens propagatie, achtergrond-I/O en checkpoints vergen eerst een apart
snapshot-/synchronisatieontwerp. Ook met `CEREAL_THREAD_SAFE=1` mag één archive
niet gelijktijdig door meerdere threads worden gebruikt. Die optie moet voor
parallelle afzonderlijke archives in alle betrokken compilatie-eenheden vóór
de Cereal-headers ingesteld zijn en beschermt alleen Cereals interne toestand,
niet de Tudat-objecten of gedeelde streams. Zie
[Cereals documentatie](https://uscilab.github.io/cereal/thread_safety.html) en
[het API-contract](python_independent_arcs_api.md#serialization-contract-for-independent-arcs).

Deze commit documenteert de bestaande gebruiksvoorwaarde. Hij introduceert
geen gelijktijdige opslag, extra locks of snapshots. Opgeslagen resultaten zijn
geen herstartbaar checkpoint van de integrator en zijn runtime-callbacks.

### Lokale gegevens en validatiescope

Lokale gegevens staan onder
`build-python-arcs-serialization-on/validation/`: build-/broncontext,
OFF-modulefingerprint en testlogs. Deze buildmap is genegeerd en gaat niet mee
met Git. Een checksum en bestandstijd/-grootte van de OFF-module dienen om te
controleren dat de ON-validatie die binary ongemoeid laat.

Dit is een controle op dezelfde machine en binnen dezelfde build-/schemafamilie;
geen bewijs van archiefcompatibiliteit tussen willekeurige versies of platformen.
Alleen vertrouwde, tijdens de tests gemaakte bestanden worden teruggelezen.
De documentatie en validatiehandelingen zijn AI-assisted; automatische
roundtrips vervangen geen onafhankelijke inhoudelijke review.
