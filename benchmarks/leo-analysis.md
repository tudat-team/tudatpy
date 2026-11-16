# Analyse van leo-var.perf.data

De opname bevat 18.133 samples over circa 90,9 seconden, zonder verloren
samples. De platte percentages zijn gewogen met de eventperiode.

## Gevonden samenhangende groep: aerodynamische toestandsafgeleiden

De lokale C++-broncode bevat de volgende keten:

```
DynamicsStateDerivativeModel::computeStateDerivative
  VariationalEquations::updatePartials
    AerodynamicAccelerationPartial::update
      zes toestandcomponenten, telkens positieve en negatieve verstoring
        vluchtcondities/caches resetten
        voertuigtoestand vervangen door verstoorde toestand
        AtmosphericFlightConditions::updateConditions
          aerodynamische hoeken en frametransformaties
        AerodynamicAcceleration::updateMembers
      nominale toestand en modellen herstellen
```

`AerodynamicAccelerationPartial::update` berekent centraal numerieke
afgeleiden: twaalf versnellingsberekeningen en één aanvullende berekening
voor herstel van de nominale toestand per aanroep. De resultaten zijn de zes
kolommen van de 3×6 aerodynamische toestands-Jacobiaan.

Hiermee vormen de afzonderlijk zichtbare rotaties, vluchtcondities,
cache-resets en versnellingsevaluaties één concreet optimalisatiegebied.
`FlightConditions::resetCurrentTime` kopieert bovendien een `vector<bool>`
met cacheflags. Dat is een broncodekandidaat voor de gemeten kopieerkosten,
geen bewezen toerekening van alle samples in `vector<bool>::operator=`.

## Profielbewijs en beperking

Eigen aandelen in het bestaande rapport:

| Routine | Eigen aandeel |
|---|---:|
| aerodynamische frametransformatie | 4,63% |
| `vector<bool>::operator=` | 2,56% |
| `AerodynamicAngleCalculator::update` | 1,34% |
| `AerodynamicAcceleration::scaleAerodynamicAcceleration` | 1,20% |
| `AerodynamicAcceleration::updateMembers` | 0,69% |

Deze routines kunnen ook buiten de afgeleidenberekening worden aangeroepen;
hun aandelen mogen niet volledig aan deze groep worden toegeschreven.

Een inclusief rapport uit dezelfde opname geeft 6,29% voor
`AerodynamicAccelerationPartial::update`, 8,73% voor `updatePartials` en
15,98% voor `computeStateDerivative`. Deze waarden zijn **geen betrouwbare
totale groepskosten**: zelfs de bovenliggende integratieroutine wordt maar
in circa 16% van de gewogen samples teruggevonden. Veel stacks zijn afgebroken
of bevatten onmogelijke/ongekende frames. Geen verloren samples betekent niet
dat stack-unwinding correct is. De opname bewijst dus geen dominante groep
van meer dan 50%.

## Optimalisatie en parallellisatie

De twaalf verstoorde evaluaties zijn wiskundig onafhankelijk bij een vaste
nominale toestand en tijd. De huidige implementatie wijzigt echter dezelfde
voertuigtoestand, vluchtcondities, caches en versnellingsobjecten. Een OpenMP-lus
om de bestaande code zou daarom races veroorzaken.

Een mogelijke implementatie gebruikt per worker eigen evaluatieobjecten en
scratchgeheugen. Iedere taak levert één Jacobiaankolom (twee evaluaties);
de nominale gedeelde omgeving blijft tijdens deze taken onveranderd. Alle
bereikte modelobjecten en eventuele bibliotheekaanroepen moeten daarbij
worden meegenomen; alleen de toestandsvector kopiëren volstaat niet.
De evaluaties zijn klein, dus threadpoolgebruik en voldoende werk per taak
zijn bepalend voor de winst.

Een alternatief is minder werk uitvoeren: analytische of semi-analytische
afgeleiden voor de gebruikte constante-coëfficiëntenconfiguratie, of hergebruik
van grootheden die bij een gegeven verstoring werkelijk onveranderd blijven.
Dat vereist behoud van alle relevante afhankelijkheden en numerieke controles
tegen de bestaande afgeleiden, ook buiten dit ene benchmarkgeval.

## Eerstvolgende meting

Meet in de daadwerkelijk gebruikte geoptimaliseerde build met symbolen de
totale tijd en aantallen van `computeStateDerivative`, omgevingsupdates,
versnellingsupdates en `updatePartials`; splits binnen de laatste de
aerodynamische en zwaartekrachtafgeleiden. Gebruik bijvoorbeeld cumulatieve
C++-timers met slechts één uitvoer na afloop, en meet de timeroverhead.
Geneste timers overlappen; houd inclusief en exclusief werk gescheiden.
Deze meting omzeilt de defecte stackreconstructie en geeft het aandeel waarop
een versnellingsschatting gebaseerd kan worden.

Voor een hersteld samplingprofiel moeten unwind-informatie, debugbestanden
en de daadwerkelijk geladen bibliotheken worden gecontroleerd. Een grotere
stackbuffer alleen garandeert geen herstel. De geïnstalleerde Tudat-build
heeft bovendien een andere interfaceversie dan deze checkout; controleer de
bronversie voordat wijzigingen worden gebouwd en vergeleken.

De SPICE-routines blijven een tweede belangrijk gebied (samen minstens
circa 28% eigen aandeel boven de rapportgrens), maar de bestaande stacks
ondersteunen geen betrouwbare toerekening aan specifieke Tudat-aanroepers.
