# Native C++-acceptatietests op de huidige machine

Uitgevoerd op 9 november 2026 door Codex, op verzoek van de gebruiker.
Branch: `parallellization_inquiry`; HEAD: `9a4c3365d`.
Er zijn geen wetenschappelijke bronbestanden of toleranties aangepast.

## Omgeving en bouw

- Intel Core i7-8550U, 4 cores / 8 logische CPU's, x86_64 Linux, 16 GB RAM.
- Conda-omgeving: `/home/hanno/miniforge3/envs/tudatpy-dev`.
- GCC 14.4.0, Eigen 5.0.1, Boost 1.85.0, CSPICE 67, Tudat Resources 2.4.
- Nieuwe buildmap: `build-parallel-arcs-local`; bestaande builds niet hergebruikt.
- Cereal ontbrak. Via Conda uitsluitend `cereal=1.3.2` toegevoegd
  (`--no-deps`); geen andere pakketten gewijzigd. Een eerste dry-run met
  `--freeze-installed` stelde ook certificaat-/OpenSSL-updates voor, die niet
  zijn uitgevoerd.
- SPICE-kernels en `gravity_models/Mars/jgmro120d.txt` waren aanwezig onder
  `/home/hanno/.tudat/resource`.

Vanuit de repositoryroot, met de genoemde Conda-omgeving geactiveerd:

```bash
cmake -S . -B build-parallel-arcs-local -G "Unix Makefiles" \
    -DCMAKE_PREFIX_PATH="$CONDA_PREFIX" \
    -DCMAKE_BUILD_TYPE=Release \
    -DTUDAT_BUILD_TESTS=ON \
    -DTUDAT_BUILD_WITH_SERIALIZATION=OFF \
    -DTUDAT_BUILD_WITH_EXTENDED_PRECISION_PROPAGATION_TOOLS=OFF \
    -DTUDAT_BUILD_EXPLICIT_INSTANTIATIONS=ON
cmake --build build-parallel-arcs-local --target copy_test_data \
    test_propagators_IndependentArcExecution \
    test_integrators_RungeKuttaCoefficients \
    test_integrators_RungeKutta87DormandPrinceIntegrator -j 2
export OPENBLAS_NUM_THREADS=1
export OMP_NUM_THREADS=1
ctest --test-dir build-parallel-arcs-local --output-on-failure \
    -R '^(test_propagators_IndependentArcExecution|test_integrators_RungeKuttaCoefficients|test_integrators_RungeKutta87DormandPrinceIntegrator)$'
ctest --test-dir build-parallel-arcs-local --output-on-failure \
    -R '^test_propagators_IndependentArcExecution$' --repeat until-fail:30
ctest --test-dir build-parallel-arcs-local --output-on-failure \
    -R '^test_integrators_RungeKuttaCoefficients$' --repeat until-fail:20
build-parallel-arcs-local/tests/test_propagators_IndependentArcExecution \
    --arcs 4 --degree 20 --days-per-arc 0.02
```

De uitvoering gebruikte de absolute paden naar CMake en CTest in bovenstaande
omgeving; de configuratie koos de compilers uit dezelfde omgeving.

## Resultaten

Alle drie geselecteerde tests slaagden. De native acceptatietest slaagde tevens
30 keer achter elkaar; de coëfficiëntentest 20 keer. De eerste native test en
de 30 herhalingen liepen terwijl de resterende integratortests nog compileerden;
de gezamenlijke drie-test-run en de schaalproef liepen na afloop van de build.

De korte schaalproef valideerde bij ieder worker-aantal alle 236 opgeslagen
arc-tijdstippen tegen de native sequentiële referentie:

| Ingestelde workers | Totale tijd (s) | Validatie |
|---:|---:|---|
| 1 | 0,2351 | geslaagd |
| 2 | 0,2816 | geslaagd |
| 4 | 0,2197 | geslaagd |
| 8 | 0,2230 | geslaagd |

Bij vier arcs zijn maximaal vier workers tegelijk actief. Deze kleine workload
en één meting per instelling onderbouwen geen schaalbaarheidsclaim. De lange
8-arc/2-daagse meting is niet uitgevoerd. Ook de volledige Tudat-testsuite,
Python-bindings en ThreadSanitizer vallen buiten deze reproductie.

## Logs en aandachtspunten

Lokale logs en het pakketmanifest staan in de gitignored map
`build-parallel-arcs-local/reproduction/`: `environment.json`, `configure.log`,
`build.log`, `acceptance-first.log`, `acceptance-repeat30.log`,
`tests-three.log`, `coefficients-repeat20.log`, `short-scaling.jsonl` en
`short-scaling.stderr.log`. Deze bestanden worden niet automatisch met Git
overgedragen. De historische meetbestanden zijn niet overschreven.

De compiler gaf een waarschuwing over een mogelijk dangling reference in
`src/tudat/simulation/estimation_setup/observationOutput.cpp:542`.
De acceptatietest gaf deprecated-API-waarschuwingen voor
`getNumericalVariationalEquationsSolution`. Beide verhinderden slagen niet;
ze zijn niet als onderdeel van deze reproductie opgelost.

De verplichte rseng-zelfcontrole is uitgevoerd. Fasedekking is vastgelegd in
`.rseng-agent-skills-coverage.md`; de controle meldt daarnaast bestaande
projectbrede hiaten: ontbrekende `CITATION.cff` en agent-directoryregels in
`.gitignore`. Deze metadata-/beleidswijzigingen zijn uitgesteld buiten de
afgebakende testreproductie. Dit verslag is door de agent opgesteld en nog
niet door de gebruiker gereviewd.
