## 📐 Geometria (Geometry Engine)

- [ ]  **NACA a 5 cifre**: Implementare la formulazione matematica per la generazione di profili NACA della serie a 5 cifre.
- [ ]  **Profili asimmetrici**: Sviluppare il supporto per la definizione di curve indipendenti per dorso e ventre.
- [ ]  **Chiusura bordo d'uscita**: Fissare lo standard geometrico di chiusura del bordo d'uscita (aperto vs chiuso) in funzione dei requisiti per la successiva generazione della mesh.

## 🌪️ Aerodinamica e Modelli Fisici

- [ ]  **Strato limite (Prandtl)**: Integrare il calcolo dello strato limite nel solutore in-house per stimare il Cd. Obiettivo: utilizzare la medesima funzione di costo di XFOIL.
- [ ]  **Correzione di Prandtl-Glauert**: Aggiungere la correzione per effetti di compressibilità subsonica (baseline: il numero di Mach incide per circa -0.02 sul Cl a 10°).
- [ ]  **Analisi Transizione (Spike)**: Analizzare l'implementazione del criterio di transizione in XFOIL (metodo e^N o similare) e definire l'architettura per replicarlo/adattarlo nel codice.
- [ ]  **Gestione Separazione (Spike)**: Strutturare la logica algoritmica per il rilevamento e l'integrazione degli effetti della separazione del flusso.
- [ ]  **Linee di corrente**: implementare l'output con linee di corrente per il profilo vincente

## ⚙️ Architettura I/O e Ottimizzazione

- [ ]  **Inizializzazione Seed**: Definire la policy per il seed del generatore (valutare seed fisso per riproducibilità dei test vs casuale per esplorazione completa dello spazio).
- [ ]  **Tuning Obiettivo 2**: Ricalibrare il limite di Cd di default (la soglia attuale di 0.02 a 4° risulta matematicamente ininfluente sull'ottimizzazione).
- [ ]  **Programma unico**: Unificare i due `run.py` in un unico entry point con selezione del solutore all'avvio (in-house / XFOIL). I due file sono oggi quasi identici: le differenze reali sono il controllo di XFOIL, le domande sull'obiettivo e i risultati finali a video. Da fare prima dello strato limite: quando l'in-house calcolerà il Cd, i due flussi diventeranno identici.
- [ ]  **Cartella risultati unica**: Spostare l'output in una singola cartella `Results/` in radice, con nome `Re<Re>_Alpha<alpha>_<obiettivo>_<solutore>_seed<N>` e suffissi `_inhouse` / `_xfoil`.
- [ ]  **Selezione solutore da CLI**: Valutare se affiancare alla domanda interattiva un'opzione da riga di comando (es. `--solver xfoil`), utile per lanciare più casi in automatico.
- [ ]  **Test di non regressione sul terminale**: Valutare se rendere permanente il confronto dell'output prima/dopo (a seme fisso) come test in `tests/`, oltre alla baseline numerica esistente.

## 📄 Gestione Repository

- [ ]  **Licenza**: Selezionare e applicare una licenza open source (es. MIT o GPL) tramite l'inserimento del file `LICENSE` per definire legalmente l'uso del codice da parte di terzi.
