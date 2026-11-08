# Multi-arc variatiepropagatie: benchmark en C++-analyse

`profile_multi_arc.py` roept eenmaal Tudats native multi-arc variatiesolver aan.
De Python-lus bouwt alleen instellingen op; hij plant geen parallelle taken en
integreert geen afzonderlijke arcs. Dit is de sequentiële uitgangsmeting voor
latere parallellisatie in C++.

Standaard: acht niet-overlappende arcs van elk twee dagen rond Mars, met een
uur tussenruimte, graad/order 50, vaste RKDP87-stap 30 s. Iedere arc heeft zes
eigen beginstaatparameters. Eén Mars-zwaartekrachtparameter is gemeenschappelijk:
in totaal 6N+1 parameters. De toestanden worden expliciet voorgeschreven en
verschillen per arc; er wordt geen toestand overgenomen uit een vorige arc.
Dit zijn synthetische orbit-determination-arcs, geen aansluitend missietraject.
Er worden nog geen observaties of normaalvergelijkingen gemaakt.

## Uitvoeren

```bash
conda activate tudatpy-with-debug-symbols
export OPENBLAS_NUM_THREADS=1
export OMP_NUM_THREADS=1

# Korte functionele controle
python benchmarks/profile_multi_arc.py --arcs 2 --days-per-arc 0.01 --degree 20

# Sequentiële uitgangsmeting
python benchmarks/profile_multi_arc.py

# Werk per arc gelijk houden; totaal werk neemt toe met aantal arcs
for arcs in 1 2 4 8; do
    python benchmarks/profile_multi_arc.py --arcs "$arcs" --days-per-arc 2 --degree 50
done

perf record -o multi-arc-var.perf.data -F 199 --call-graph dwarf,32768 -- \
    python benchmarks/profile_multi_arc.py --days-per-arc 6
perf report -i multi-arc-var.perf.data --stdio --no-children --call-graph none \
    --sort dso,symbol --percent-limit 0.5 > multi-arc-var-self.txt
```

Kies indien nodig meer dagen per arc voor een opname van 60–180 s. Het verhogen
van `--arcs` verhoogt nu de hoeveelheid sequentieel werk; het activeert geen
threads. Voor latere strong-scalingmetingen moet juist hetzelfde aantal arcs,
dezelfde duur en hetzelfde model bij verschillende worker-aantallen worden
gebruikt. Meet meerdere runs zonder profiler, ook de totale duur en het geheugen.

Lokale controle: de standaardrun met acht arcs van twee dagen slaagde in circa
26,0 s native solverduur; alle acht 6×49-matrices voldeden aan de blokcontroles.
Een korte vergelijking van één versus twee arcs gaf voor arc 0 dezelfde STM en
mu-gevoeligheid (rtol 1e-12, atol 1e-15). Zes dagen per arc is daarom een
startpunt voor ongeveer 78 s profieltijd, geen gegarandeerde duur.

De timer omvat native solverconstructie, alle integraties en native
nabewerking. Sommige geïnstalleerde versies exposen voor multi-arc alleen de
basisklasse; daarom gebruikt het script de factory met integratie bij creatie.
De tijd is dus geen zuivere integratietijd zoals in de single-arc benchmarks.
De opbouw van Python-instellingen en de uitvoer/controles vallen erbuiten.
`perf` meet het volledige proces. JSON bevat configuratie, arcgrenzen,
CPU/walltijd en de uiteindelijke STM en mu-gevoeligheid per arc.

De controle vraagt voor iedere arc op de eindtijd een 6×(6N+1)-matrix op:
de eigen beginstaatkolommen zijn actief, die van andere arcs zijn nul,
de laatste kolom bevat de gedeelde mu-gevoeligheid. Alle waarden moeten eindig
zijn. Dit controleert de blokindeling, niet de volledige schattingsketen of
numerieke nauwkeurigheid. Er wordt geen gedeelde parameter per arc dubbel geteld.

## Concrete ingang in de lokale C++-broncode

1. `include/tudat/simulation/estimation_setup/multiArcVariationalEquationsSolver.h`:
   `integrateVariationalAndDynamicalEquations` roept de multi-arc dynamicasolver
   aan met een variatieresultatencontainer en initial-state-provider.
2. `include/tudat/simulation/propagation_setup/multiArcDynamicsSimulator.h`:
   de template-overload van `integrateEquationsOfMotion` bevat de sequentiële
   lus over `singleArcDynamicsSimulators_`. De aanroep per arc schrijft naar
   `propagationResults->getSingleArcResults().at(i)`.
3. Daarna worden de gezamenlijke status, eventueel bijgewerkte beginstaten en
   geïntegreerde ephemeriden verwerkt via
   `processNumericalEquationsOfMotionSolution`.
4. De variatiesolver maakt per-arc interpolatoren en vervolgens de gezamenlijke
   toestandovergangs-/sensitiviteitsinterface in
   `resetVariationalEquationsInterpolators`.

Dit is de concrete taakgrens: **één volledige arc-integratie per taak, gevolgd
door gezamenlijke resultaatverwerking**. De laatste stap blijft eerst serieel.

## Wat een parallelle implementatie moet scheiden

De single-arc simulators worden nu met dezelfde `SystemOfBodies` gebouwd.
Toestanden, tijden, rotaties en caches veranderen tijdens evaluaties. De
afzonderlijke modellen in het benchmarkscript maken deze omgeving niet
automatisch thread-safe. Ook de arc-lus heeft gedeelde scratchvariabelen en een
groeiende lijst met beginstaten. Een `parallel for` om de bestaande lus volstaat
daarom niet.

Een eerste C++-implementatie zou expliciete onafhankelijke beginstaten moeten
vereisen, per worker een geïsoleerde veranderlijke omgeving en solver gebruiken,
arcresultaten op vaste indices opslaan, fouten/status per arc verzamelen, en
pas na een join de globale ephemeriden en variatie-interface bijwerken.
Identieke onveranderlijke coëfficiënten kunnen waar verantwoord gedeeld blijven.
Controleer ook bibliotheekcalls zoals SPICE; een gekopieerde bodies-container
garandeert geen isolatie van externe bibliotheektoestand of gedeelde pointers.

Voor gedeelde parameters moet iedere worker dezelfde waarde krijgen, maar de
sensitiviteitskolom moet bij samenvoegen naar **dezelfde globale parameterkolom**
wijzen. De zes beginstaatkolommen horen juist bij de betreffende arc. Een
latere gezamenlijke schatting bevat daardoor ook koppelingen tussen arcs via
die gedeelde parameters; die worden niet met afzonderlijke fits vervangen.

Niet iedere multi-arc-configuratie is onafhankelijk: state transfer, ontbrekende
beginstaten die uit een vorige arc worden afgeleid, hybride dynamica en andere
expliciete koppelingen moeten worden herkend en voorlopig sequentieel blijven.

## Acceptatie van de eerste parallelle versie

Vergelijk bij vast model alle arc-eindstaten, STM's, sensitiviteiten,
arcgrenzen en globale parameterindeling met de sequentiële versie; test ook
foutafhandeling. De JSON-matrices zijn bruikbaar als eerste regressiereferentie.
Meet daarna dezelfde workload met 1, 2, 4 en 8 workers, inclusief opbouw en
samenvoegen. Verwacht geen versnelling groter dan het aantal gelijktijdige arcs;
geheugen, taakbalans en sequentiële nabewerking beperken de winst.

De bronanalyse betreft deze checkout. De geïnstalleerde debugbibliotheek heeft
een oudere Python-interface (geen geëxposeerde MultiArcVariationalSimulator).
De factory-route is getest; vóór C++-wijzigingen moet ook de gebruikte
bron/buildversie worden vastgesteld. Het Python-script blijft de sequentiële
referentie. Er is inmiddels een afzonderlijk [C++-threadprototype](cpp/README.md)
met geïsoleerde omgevingen, vergelijking met de native referentie en instelbaar
worker-aantal. De publieke Tudat-solver en Python-API zijn nog niet aangepast.
