# La Chimica del Vaso

Simulatore interattivo dei vasi comunicanti applicato alla capacità tampone acido-base e al potenziale chimico, pensato per lo studio della chimica al liceo.

Ogni vaso rappresenta una coppia acido/base coniugata: la sua forma è la curva di capacità tampone β(pH). I vasi si collegano con valvole come dei veri vasi comunicanti, condividendo lo stesso pH quando sono uniti.

## Come si usa

Apri [`index.html`](index.html) in un browser (anche con un semplice doppio clic, senza bisogno di un server). Il pulsante "Istruzioni per l'uso" nell'app spiega le quattro modalità (Modello libero, Preparazione tampone, Titolazione, Specie predominanti) e tutti i comandi.

## Struttura del progetto

- `index.html` — markup della pagina
- `style.css` — tutti gli stili
- `app.js` — il modello chimico e tutta la logica dell'interfaccia

Nessuna dipendenza da installare: è HTML/CSS/JS puro, nessun passaggio di build.
