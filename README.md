# OPTIMIZATION FOR CYBERSECURITY

# DESCRIZIONE:

Progetto realizzato per il corso di Optimization Methods for CyberSecurity con l'obiettivo di:
1) Generare le matrici A e Aeq ed i vettori b, beq, c, lb ed ub del seguente problema MILP:
   > min c^T x
   > 
   > Ax <= b
   > 
   > Aeqx = beq
   > 
   > lb <= x <= ub

2) Popola il problema in tre differenti modi:
   1) Per colonna
   2) Per riga
   3) Per coefficienti diversi a zero

3) Chiama CPLEX per la sua soluzione

# ESTRAZIONE DEI DATI
Per generare le matrici e i vettori del problema MILP a partire dal file MPS letto ho utilizzato i metodi forniti dalla libreria CPLEX in Python:

_**prob.linear_constraints.get_rows():**_
> Restituisce le righe delle matrici A e Aeq in formato di lista di tuple, dove ogni tupla contiene gli indici delle colonne non nulle e i relativi valori.


_**prob.linear_constraints.get_rhs():**_
> Restituisce il vettore b.

_**prob.quadratic_constraints.get_rows():**_
> Restituisce le righe della matrice Aeq in formato di lista di tuple.

_**prob.quadratic_constraints.get_rhs():**_
> Restituisce il vettore beq.

_**prob.objective.get_linear():**_
> Restituisce il vettore c.

_**prob.variables.get_lower_bounds():**_
> Restituisce il vettore lb.

_**prob.variables.get_upper_bounds():**_
> Restituisce il vettore ub.

# POPOLAMENTO PER COLONNA
1. Per ogni colonna, vengono creati tre elenchi vuoti per: variabili coinvolte, coefficienti di ogni variabile nel vincolo e limiti inferiore e superiore delle variabili coinvolte e il tipo di variabile.

2. Per i vincoli di tipo "<=" viene popolato l'elenco delle variabili, dei coefficienti e dei limiti inferiori e superiori.

3. Per i vincoli di tipo "==" viene popolato l'elenco delle variabili, dei coefficienti e dei limiti inferiori e superiori.

4. Viene aggiunta una variabile al problema mediante il metodo "my_prob.variables.add()", specificando il coefficiente dell'obiettivo, i limiti inferiore e superiore della variabile, la tipologia della variabile e il nome della variabile.

5. Viene aggiunto un vincolo lineare di tipo "<=" al problema mediante il metodo "my_prob.linear_constraints.add()", specificando la coppia di valori (variabili e coefficienti) per la parte lineare, il segno di disuguaglianza, il lato destro e il nome del vincolo.

6. Viene aggiunto un vincolo quadratico di tipo "==" al problema mediante il metodo "my_prob.quadratic_constraints.add()", specificando la coppia di valori (variabili e coefficienti) per la parte lineare, la lista vuota per la parte quadratica, il segno di uguaglianza, il lato destro e il nome del vincolo.


# POPOLAMENTO PER RIGA
1. Verifica se le dimensioni delle matrici di input corrispondono ai vincoli del problema.

2. Inizializza tre liste vuote per contenere: indici delle righe, indici delle colonne ed i valori della matrice A del problema di programmazione lineare.

3. Determina il numero di colonne e righe nel problema.

4. Popola i vincoli lineari di tipo "<=" o ">=" in base al valore di my_sense.

5. Aggiunge i vincoli lineari di tipo "==" al problema.

6. Popola i vincoli quadratici di tipo "==".

7. Aggiunge le variabili del problema, distinguendo tra variabili continue e variabili binarie/interne.

# POPOLAMENTO PER COEFFICIENTI DIVERSI DA ZERO
1. Controlla che my_beq non sia vuota.

2. Controlla che la lunghezza di my_beq sia maggiore o uguale a numrows.

3. Per ogni riga, se my_sense[i] è 'E' (cioè l'uguaglianza), aggiunge il vincolo quadratico di tipo '==' al problema lineare, utilizzando my_Aeq[i] come gli indici emy_beq[i] come la soluzione.
Se invece my_sense[i] è '<=' o '>=', aggiunge il vincolo lineare corrispondente usando my_A[i] come gli indici e my_b[i] come la soluzione.

4. Per ogni colonna, se il coefficiente my_c[j] è diverso da 0.0, aggiunge la variabile di decisione, utilizzando my_c[j] come il coefficiente obiettivo, my_colnames[j] come il nome della variabile, e le rispettive limitazioni superiori e inferiori (se definite).

# OUTPUT
Questo script in Python risolve un problema di ottimizzazione lineare CPLEX, e restituisce una serie di informazioni sull'output della soluzione.

In particolare, questa funzione estrae diverse informazioni dal problema di ottimizzazione, come i limiti delle variabili, i vincoli, la funzione obiettivo e le informazioni sul tipo di variabile.

Risolve il problema usando il metodo _**solve()**_ fornito da CPLEX.

Infine, la funzione stampa le informazioni sulla soluzione, inclusi lo stato della soluzione, il valore della funzione obiettivo e i valori delle variabili ottimali.

Inoltre, la funzione stampa le variabili slack e i valori delle variabili di decisione. 

## Linguaggi e Strumenti Utilizzati

◉ Python 3.11.3

## Authors

- [@DanieleCampopiano](https://github.com/DanieleCampopiano)