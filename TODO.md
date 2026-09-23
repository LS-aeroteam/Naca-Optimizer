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

- [ ]  **Export CAD**: Sviluppare la funzione di esportazione della geometria del profilo ottimizzato in formati CAD standard (es. STEP, IGES, o dat).
- [ ]  **Inizializzazione Seed**: Definire la policy per il seed del generatore (valutare seed fisso per riproducibilità dei test vs casuale per esplorazione completa dello spazio).
- [ ]  **Tuning Obiettivo 2**: Ricalibrare il limite di Cd di default (la soglia attuale di 0.02 a 4° risulta matematicamente ininfluente sull'ottimizzazione).

## 📄 Gestione Repository

- [ ]  **Licenza**: Selezionare e applicare una licenza open source (es. MIT o GPL) tramite l'inserimento del file `LICENSE` per definire legalmente l'uso del codice da parte di terzi.
