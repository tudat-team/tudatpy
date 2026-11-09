# Eerste onderzoek naar parallellisatie

> **Historische onderzoeksnotitie — status bijgewerkt op 9 november 2026.**
> Dit document bewaart de eerste verkenning. De profielcijfers hieronder komen
> uit een oudere opname waarvan scenario en bouwconfiguratie onvoldoende zijn
> vastgelegd; ze beschrijven niet de huidige benchmarks. Uitspraken over nog
> ontbrekende workloads of een ongewijzigde rekenkern gelden voor die eerste
> onderzoeksfase. Inmiddels zijn representatieve benchmarks, een parallel
> prototype en een opt-in route in de native C++-multi-arc solver toegevoegd.
> Zie [Independent arcs in the native C++ multi-arc solver](independent_arc_execution.md)
> voor de huidige implementatie, tests, resultaten en beperkingen. De tekst
> hieronder blijft als oorspronkelijke onderzoekscontext behouden.

## Status en beschikbare meting

Onderzoekstak: `parallellization_inquiry`. Nog geen wijzigingen aan de rekenkern.
De bestaande lokale bestanden `graph.dot`, `callgraph.svg`, `out.perf` en
`perf.data.old` zijn behouden. De onderstaande getallen komen uit de knooplabels
in `graph.dot`, niet uit een nieuwe benchmark. Het scenario, de broncommit,
buildopties, hardware en profileringsopdracht van die opname zijn onbekend.
Daarom zijn dit aanwijzingen voor vervolgmetingen, geen algemene Tudat-benchmark.
De grafiek verwijst naar een geïnstalleerde `tudatpy/kernel.so` in `tudatpy-dev`;
die hoeft niet overeen te komen met deze checkout.

| Functie | Inclusief aangeroepen functies | Eigen aandeel |
| --- | ---: | ---: |
| `DynamicsStateDerivativeModel::computeStateDerivative` | 78,70% | 1,23% |
| `EnvironmentUpdater::updateEnvironment` | 50,84% | 0,14% |
| `SpiceEphemeris::getCartesianState` | 31,41% | 0,03% |
| `SpiceRotationalEphemeris::getRotationToBaseFrame` | 7,98% | 0,00% |
| `blas_thread_server` | 7,19% | 7,04% |

Deze aandelen overlappen langs aanroepketens en mogen niet worden opgeteld.
Het zijn aandelen van de profileringssamples, geen gemeten wandkloktijden.
Het BLAS-aandeel vraagt om een vergelijking met één BLAS-thread: de grafiek
alleen vertelt niet hoeveel nuttig rekenwerk of wachtwerk die threads doen.

## Advies op basis van de broncode

1. **Begin met onafhankelijke simulaties**, als het doel Monte Carlo,
   parameterstudies of optimalisatie is. Verdeel complete runs over processen.
   Bouw per proces de eigen bodies, modellen en simulator op; geef eenvoudige
   configuratie en compacte resultaten door. Meet opstart-, geheugen- en
   overdrachtskosten. Stel aanvankelijk BLAS/OpenMP-threadaantallen per worker
   op één in en meet alternatieven om overbezetting te voorkomen.
2. **Onderzoek minder SPICE-werk binnen één simulatie.** De bestaande
   `InterpolatedSpiceEphemerisSettings` in
   `include/tudat/simulation/environment_setup/createEphemeris.h` biedt hiervoor
   een kandidaat. Vergelijk directe en geïnterpoleerde ephemeriden inclusief
   tabelopbouw, geheugengebruik en fout over het gehele relevante tijdsinterval.
   Dit vervangt niet automatisch de rotatiemodellen. Kies interpolatiestap en
   dekking op basis van de vereiste nauwkeurigheid.
3. **Onderzoek interne threads pas na gerichte meting.**
   `DynamicsStateDerivativeModel::computeStateDerivative` wijzigt eigen buffers,
   werkt de omgeving bij en werkt vervolgens modellen bij
   (`include/tudat/astro/propagators/dynamicsStateDerivativeModel.h:137`).
   `Body` bewaart en wijzigt actuele toestanden
   (`include/tudat/simulation/environment_setup/body.h:184`). Een parallelle lus
   over deze bestaande objecten vereist dus eerst een analyse van afhankelijkheden
   en gedeelde toestand. Controleer ook expliciet de threadveiligheid van externe
   bibliotheken en Python-callbacks voordat threads worden gebruikt.
4. **Verdeel opeenvolgende integratiestappen niet zomaar.** Elke stap gebruikt
   de voorgaande toestand; ook arcs kunnen van elkaar afhangen. Selecteer alleen
   aantoonbaar onafhankelijke taken. Kleine krachtberekeningen kunnen minder
   opleveren dan de overhead van scheduling en synchronisatie.

## Reproduceerbaar meetplan

1. Kies het concrete script en de representatieve invoer. Leg commit,
   geïmporteerde modulelocatie, compiler/buildtype, bibliotheken, CPU,
   threadinstellingen, kernels, integrator en toleranties vast.
2. Meet een geoptimaliseerde build met symbolen. Scheid imports, laden van
   kernels, modelopbouw, integratie/variatievergelijkingen, schatting en export.
   Rapporteer zowel totale tijd als rekentijd; houd uitvoerinstellingen gelijk.
3. Doe een opwarmrun en ten minste vijf herhalingen; rapporteer mediaan en
   spreiding van wandkloktijd, CPU-tijd en piekgeheugen. Verzamel daarnaast een
   native samplingprofiel met herkenbare C++-symbolen. Een Python-profiel alleen
   kan de tijd binnen `kernel.so` niet voldoende uitsplitsen.
4. Vergelijk afzonderlijk: bestaande configuratie, één BLAS-thread,
   geïnterpoleerde ephemeriden, en onafhankelijke runs met 1/2/4/8 processen
   voor zover hardware en geheugen dit toelaten. Gebruik steeds hetzelfde werk.
5. Controleer toestand/resultaten op dezelfde evaluatietijden tegen vooraf
   gekozen foutgrenzen. Meet versnelling als T(1)/T(p), naast totale doorvoer en
   geheugen. Accepteer alleen winst die groter is dan de meetspreiding.

## Nog nodig

Het script/scenario dat sneller moet worden en de gewenste uitkomst:
een kortere individuele simulatie of meer simulaties per uur. Ook moet worden
vastgesteld of de bestaande profileringsopname bij dat scenario hoort.
Zonder die koppeling is een onderbouwde keuze voor een parallelle implementatie
nog niet mogelijk.

## Aanvulling: Eigen en parallelle BLAS/LAPACK-backends

Het door een collega aangeleverde Claude-antwoord noemt een bruikbare optie:
Eigen behouden en geselecteerde lineaire algebra via een externe backend laten
uitvoeren. De officiële Eigen-documentatie bevestigt beperkte interne
OpenMP-parallellisatie (onder meer matrixproducten, PartialPivLU en enkele
sparse iteratieve solvers). BLAS/LAPACKE-delegatie vereist daarnaast expliciete
compileermacro's en het linken van de betreffende bibliotheken.

Bij de eerste inspectie van CMakeLists.txt, cmake en configuration zijn geen
expliciete EIGEN_USE_BLAS/EIGEN_USE_LAPACKE/EIGEN_USE_MKL_ALL-instellingen
gevonden. Extern meegegeven compileropties zijn daarmee niet uitgesloten.
OpenBLAS in het bestaande profiel bewijst niet dat Tudats Eigen-code daarvan
gebruikmaakt; ook Python-afhankelijkheden kunnen deze bibliotheek laden.

Een backendexperiment is vooral relevant als grote matrixbewerkingen veel
rekentijd innemen. De aanwezige SPICE-hotspots worden er niet rechtstreeks door
versneld. Kleine vaste 3x3/6x6-bewerkingen zijn geen overtuigende reden voor
een migratie. PARDISO is een kandidaat voor grote sparse stelsels, waarvoor
in deze eerste inspectie nog geen dominante workload is aangetoond.

In src/tudat/math/basic/leastSquaresEstimation.cpp staat een JacobiSVD-pad.
Eigen waarschuwt dat EIGEN_USE_LAPACKE dit door gesvd kan vervangen, met een
andere robuustheid; EIGEN_USE_LAPACKE_STRICT voorkomt die specifieke vervanging.
Een benchmark moet dus ook slecht geconditioneerde problemen en numerieke
uitkomsten controleren. Een andere backend is geen zuivere threadinstelling.

Bronnen (officiële ontwikkelversiedocumentatie; controleer ook de daadwerkelijk
gebruikte Eigen-versie):

- https://libeigen.gitlab.io/eigen/docs-nightly/TopicMultiThreading.html
- https://libeigen.gitlab.io/eigen/docs-nightly/TopicUsingIntelMKL.html

De licentie- en AMD-prestatieclaims uit het aangeleverde antwoord zijn hier niet
afzonderlijk geverifieerd en vormen geen basis voor een backendkeuze.
