# La Chimica del Vaso

Simulatore interattivo dei vasi comunicanti applicato alla capacità tampone acido-base e al potenziale chimico, pensato per lo studio della chimica al liceo.

Ogni vaso rappresenta una coppia acido/base coniugata: la sua forma è la curva di capacità tampone β(pH). I vasi si collegano con valvole come dei veri vasi comunicanti, condividendo lo stesso pH quando sono uniti.

## Come si usa

Apri [`index.html`](index.html) in un browser (anche con un semplice doppio clic, senza bisogno di un server). Il pulsante "Istruzioni per l'uso" nell'app spiega le quattro modalità (Modello libero, Preparazione tampone, Titolazione, Specie predominanti) e tutti i comandi.

## Struttura del progetto

- `index.html` — markup della pagina
- `style.css` — tutti gli stili
- `app.js` — il modello chimico e tutta la logica dell'interfaccia
- `pka-database.js` — sostanze e pK<sub>a</sub> dal database CurTiPot (file generato, vedi sotto)
- `tools/build_pka_database.py` — genera `pka-database.js` dal file Excel di CurTiPot

Nessuna dipendenza da installare: è HTML/CSS/JS puro, nessun passaggio di build.

## Dati dei pKa e citazione

L'elenco di circa 250 sostanze e 9 indicatori con i loro pK<sub>a</sub> proviene dal foglio "Database" di **CurTiPot**, il foglio di calcolo gratuito per Excel di Ivano G. R. Gutz (Instituto de Química, Universidade de São Paulo, Brasile), versione 4.4.0 (gennaio 2024): <https://www.iq.usp.br/gutz>. La maggior parte dei valori è a 25 °C e forza ionica 0, senza garanzia di correttezza (lo dice lo stesso database).

Se usi o riporti questi dati, cita CurTiPot:

> Gutz, I. G. R. *CurTiPot*, versione 4.4.0 (gennaio 2024). Instituto de Química, Universidade de São Paulo. https://www.iq.usp.br/gutz

`pka-database.js` è generato da `tools/build_pka_database.py` (solo libreria standard, legge i dati senza eseguire le macro). Il file originale `curtipot_i.xlsm` non è incluso nel repository (cartella `Source/`, ignorata da git):

```bash
python3 -I tools/build_pka_database.py Source/curtipot_i.xlsm pka-database.js
```
