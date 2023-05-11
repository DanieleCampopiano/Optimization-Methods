# Optimization Methods for Cybersecurity
# Campopiano Daniele - 174624
# Monday 08/05/2023


import cplex
from cplex.exceptions import CplexError
import os


# POPOLAMENTO DEL PROBLEMA PER COLONNA: VINCOLI => RIGHE | VARIABILI => COLONNE
def populateByColumn(my_prob, my_A, my_Aeq, my_c, my_lb, my_ub, my_ctype, my_b, my_sense, my_num_cols):

    for j in range(my_num_cols):
        variables = []                      # variabili coinvolte nel vincolo
        coefficients = []                   # coefficienti di ogni variabile nel vincolo
        lower_bounds = []                   # limite inferiore delle variabili coinvolte
        upper_bounds = []                   # limite superiore delle variabili coinvolte
        types = []                          # tipologia delle variabili coinvolte

        # popolamento delle variabili, dei coefficienti e dei limiti per i vincoli di tipo "<="
        for i in range(len(my_A)):
            variables.append(i)
            coefficients.append(my_A[i][j])
            lower_bounds.append(0.0)
            upper_bounds.append(1.0)
            types.append(my_ctype[j])

        # popolamento delle variabili, dei coefficienti e dei limiti per i vincoli di tipo "=="
        for i in range(len(my_Aeq)):
            variables.append(i + len(my_A))
            coefficients.append(my_Aeq[i][j])
            lower_bounds.append(0.0)
            upper_bounds.append(1.0)
            types.append(my_ctype[j])

        colname = "x_" + str(j)                                 # nome della variabile

        try:
            my_prob.variables.add(
                obj=my_c[j],                                    # coefficiente dell'obiettivo
                lb=my_lb[j],                                    # limite inferiore della variabile
                ub=my_ub[j],                                    # limite superiore della variabile
                types=my_ctype[j],                              # tipologia della variabile
                names=[colname]                                 # nome della variabile
            )
        except TypeError:
            break

        # aggiunta del vincolo lineare di tipo <= al problema
        my_prob.linear_constraints.add(
            lin_expr=[cplex.SparsePair(ind=variables, val=coefficients)],
            senses=my_sense,
            rhs=[my_b],
            names=[colname]
        )

        # aggiunta del vincolo quadratico di tipo == al problema
        my_prob.quadratic_constraints.add(
            lin_expr=[cplex.SparsePair(ind=variables, val=coefficients)],
            quad_expr=[[], [], []],
            sense='E',
            rhs=[0.0],
            names=[colname]
        )


# POPOLAMENTO DEL PROBLEMA PER RIGA: 'L' = "<=" | 'E' = "==" | 'G' = ">="
def populateByRow(my_prob, my_A, my_Aeq, my_c, my_lb, my_ub, my_ctype, my_colnames, my_b, my_beq, my_sense, my_rownames):

    if len(my_A) != len(my_b) or len(my_A[0]) != len(my_c):
        return

    i = 0
    rows = []
    cols = []
    vals = []

    numcols = len(my_c)
    numrows = len(my_b)

    # Popolare i vincoli lineari di tipo "<=" o ">="
    for i in range(numrows):
        row = []
        row_eq = []
        for j in range(numcols):
            val = my_A[i][j]
            if val != 0:
                rows.append(i)
                cols.append(j)
                vals.append(val)
                row.append(j)
        if my_sense[i] == 'E':
            # Aggiunta del vincolo lineare di tipo "=="
            my_prob.linear_constraints.add(lin_expr=[cplex.SparsePair(ind=row, val=vals)],
                                           senses=[my_sense[i]], rhs=[my_beq[i]], names=[my_rownames[i]])
        else:
            # Aggiunta del vincolo lineare di tipo "<=" o ">="
            my_prob.linear_constraints.add(lin_expr=[cplex.SparsePair(ind=row, val=vals)],
                                           senses=[my_sense[i]], rhs=[my_b[i]], names=[my_rownames[i]])
        vals.clear()
        row.clear()

    # Ciclo for per popolare i vincoli quadratici di tipo "=="
    for i in range(len(my_Aeq)):
        row_eq = []
        for j in range(numcols):
            val = my_Aeq[i][j]
            if val != 0:
                row_eq.append(j)
        # Aggiunta del vincolo quadratico di tipo "=="
        my_prob.quadratic_constraints.add(lin_expr=[cplex.SparsePair(ind=row_eq, val=my_Aeq[i])],
                                          rhs=[my_beq[i]], sense='E', name='qc' + str(i + 1))
        row_eq.clear()

    # Ciclo for per popolare le variabili
    for j in range(numcols):
        if my_ctype[j] == 'C':
            # Aggiunta della variabile continua
            my_prob.variables.add(obj=[my_c[j]], lb=[my_lb[j]], ub=[my_ub[j]], types=[my_ctype[j]],
                                  names=[my_colnames[j]])
        else:
            # Aggiunta della variabile binaria o intera
            my_prob.variables.add(obj=[my_c[j]], lb=[my_lb[j]], ub=[my_ub[j]], names=[my_colnames[j]])



# POPOLAMENTO DEL PROBLEMA PER COEFFICIENTI DIVERSI DA ZERO
def populateByNonZero(my_prob, my_A, my_Aeq, my_c, my_colnames, my_b, my_beq, my_sense, my_rownames):

    # Verifica che my_beq non sia vuota
    if not my_beq:
        return

    # Verifica che la lunghezza di my_beq sia maggiore o uguale a numrows
    numrows = len(my_sense)
    if len(my_beq) < numrows:
        return

    # Popola il problema con coefficienti diversi da zero
    for i in range(numrows):
        if my_sense[i] == 'E':
            my_prob.linear_constraints.add(lin_expr=[cplex.SparsePair(ind=my_Aeq[i], val=[1.0] * len(my_Aeq[i]))],
                                           rhs=my_beq[i], names=[my_rownames[i]])
        else:
            my_prob.linear_constraints.add(lin_expr=[cplex.SparsePair(ind=my_A[i], val=[1.0] * len(my_A[i]))],
                                           rhs=my_b[i], names=[my_rownames[i]])
    for j in range(len(my_c)):
        if my_c[j] != 0.0:
            my_prob.objective.set_linear([(my_colnames[j], my_c[j])])


# RISOLUZIONE DEL PROBLEMA
def solveProblem(popMethod, prob):

    try:
        beq = prob.quadratic_constraints.get_rhs()               # Estrae i termini quadrati
        b = prob.linear_constraints.get_rhs()                    # Estrae i termini lineari
        c = prob.objective.get_linear()                          # Estrae la funzione obiettivo

        lb = prob.variables.get_lower_bounds()                   # Estrae i limiti inferiori
        ub = prob.variables.get_upper_bounds()                   # Estrae i limiti superiori

        numcols = prob.variables.get_num()                       # Estrae il numero di colonne
        colnames = prob.variables.get_names()                    # Estrae i nomi delle colonne
        numrows = prob.linear_constraints.get_num()              # Estrae il numero di righe
        rownames = prob.linear_constraints.get_names()           # Estrae i nomi delle righe

        ctype = prob.variables.get_types()                        # Estrae il tipo di variabile
        sense = prob.linear_constraints.get_senses()              # Estrae il tipo di vincolo (<=, >=, =)

        rows = []
        cols = []
        vals = []                                                    # Lista dei valori della matrice
        A = []                                                       # Matrice dei coefficienti di disuguaglianza
        Aeq = []                                                     # Matrice dei coefficienti di uguaglianza

        # Estrae la matrice del vincolo in forma di triple (riga, colonna, valore)
        for i in range(numrows):
            row = prob.linear_constraints.get_rows(i)
            for j in range(len(row.ind)):
                rows.append(i)
                cols.append(row.ind[j])
                vals.append(row.val[j])
            if sense[i] == 'E':
                Aeq.append(vals.copy())
            else:
                A.append(vals.copy())
            vals.clear()

        # Popola il problema di ottimizzazione per riga, colonna o coefficienti diversi da zero
        if popMethod == "r":
            populateByRow(prob, A, Aeq, c, lb, ub, ctype, colnames, b, beq, sense, rownames)
        elif popMethod == "c":
            populateByColumn(prob, A, Aeq, c, lb, ub, ctype, b, sense, numcols)
        elif popMethod == "n":
            populateByNonZero(prob, A, Aeq, c, colnames, b, beq, sense, rownames)

        prob.solve()                                 # Risolve il problema di ottimizzazione

    except CplexError as exc:
        print(exc)
        return

    print("\n")
    print("Solution status = ", prob.solution.get_status(), ":", end=' ')
    print(prob.solution.status[prob.solution.get_status()])
    print("Solution value  = ", prob.solution.get_objective_value())               # Valore della funzione obiettivo
    print("\n")

    slack = prob.solution.get_linear_slacks()
    x = prob.solution.get_values()

    for j in range(numrows):
        print("Row %d:  Slack = %10f" % (j, slack[j]))
    for j in range(numcols):
        print("Column %d:  Value = %10f" % (j, x[j]))


if __name__ == "__main__":

    myProb = cplex.Cplex()

    mps_file = "Mps Files/blend2.mps"

    if os.path.exists(mps_file):
        print("\nAnalizzo il file nella path " + mps_file)
    else:
        print("\nFile non trovato nella path " + mps_file)

    myProb.read(mps_file)

    # COLONNA
    print("\n      c          generate problem by column\n")
    solveProblem("c", myProb)

    # RIGA
    print("\n      r          generate problem by row\n")
    #solveProblem("r", myProb)

    # COEFFICIENTI DIVERSI DA ZERO
    print("\n      n          generate problem by nonzero\n")
    #solveProblem("n", myProb)
